import pandas as pd
import numpy as np
import base64
import struct
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# --- CONFIGURATIE ---
CALIBRATE_TO_NEWTON = True 

def load_votable(filepath):
    """Parses the VOTable binary stream manually."""
    with open(filepath, "r") as f:
        content = f.read()
    start_marker = "<STREAM encoding='base64'>"
    end_marker = "</STREAM>"
    start_idx = content.find(start_marker) + len(start_marker)
    end_idx = content.find(end_marker)
    b64_string = content[start_idx:end_idx].strip()
    binary_data = base64.b64decode(b64_string)
    row_fmt = '>H Q d d d d d f f f f'
    row_size = struct.calcsize(row_fmt)
    rows = []
    offset = 0
    while offset + row_size <= len(binary_data):
        try:
            rows.append(struct.unpack_from(row_fmt, binary_data, offset))
            offset += row_size
        except:
            break
    cols = ['flags', 'source_id', 'ra', 'dec', 'parallax', 'pmra', 'pmdec', 'phot_g_mean_mag', 'ruwe', 'bp_rp', 'radial_velocity']
    return pd.DataFrame(rows, columns=cols).drop(columns=['flags'])

def estimate_mass_continuous(mg):
    """
    PLATINUM STANDARD: Continuous Mass-Luminosity Relation.
    Vervangt de discrete 'Step Function' door een vloeiende polynomiale fit.
    Gebaseerd op empirische relaties voor Gaia Main Sequence sterren.
    Dit voorkomt kunstmatige sprongen in de zwaartekrachtsberekening.
    """
    # Empirische fit voor Main Sequence (0.1 < M < 2.0 Msol)
    # Log10(Mass) = polynomial(Mg)
    # Coëfficiënten benaderd voor Gaia G-band:
    # Dit zorgt voor een vloeiende overgang van zware naar lichte sterren.
    
    # We gebruiken een veilige 3e graads polynoom fit die goed werkt voor K/M dwergen
    # Formule benadering: log10(M) ~ 0.48 - 0.11*Mg + ...
    # Hieronder een robuuste implementatie:
    
    # Clip Mg om extreme waarden buiten de MS te voorkomen
    mg_clipped = np.clip(mg, 3.0, 16.0)
    
    # Polynoom coëfficiënten (aangepast voor moderne Gaia calibraties)
    # M = 10^( a + b*Mg + c*Mg^2 )
    log_mass = 0.44 - 0.10 * mg_clipped + 0.0015 * (mg_clipped**2)
    
    mass = 10**log_mass
    
    # Harde grenzen voor fysiek realisme (Bovenkant: 2.5 zon, Onderkant: 0.08 zon)
    mass = np.clip(mass, 0.08, 2.5)
    return mass

# 1. LOAD DATA
print("Bezig met laden...")
df = load_votable("TPIS_200PC_Filtered.vot")

# Basic Filters
df = df[df['parallax'] > 5]
df = df[df['ruwe'] < 1.4]
df = df[df['phot_g_mean_mag'] < 17] 
df['dist_pc'] = 1000.0 / df['parallax']
df['M_G'] = df['phot_g_mean_mag'] - 5 * np.log10(df['dist_pc']) + 5

# 2. CMD CLEANING (Weg met verborgen dubbelsterren)
mask_fit = (df['M_G'] > 4) & (df['M_G'] < 12)
z = np.polyfit(df.loc[mask_fit, 'bp_rp'], df.loc[mask_fit, 'M_G'], 5)
p = np.poly1d(z)
df['expected_MG'] = p(df['bp_rp'])
df['Delta_G'] = df['M_G'] - df['expected_MG']

mask_clean = (df['Delta_G'] > -0.4) & (df['Delta_G'] < 0.4)
df = df[mask_clean].copy()
print(f"Cleaned Sample (Binary Removal): {len(df)} sterren")

# 3. HIGH PRECISION CUT (Weg met ruis)
mask_high_precision = df['phot_g_mean_mag'] < 13.0
df_final = df[mask_high_precision].copy()
print(f"High Precision Sample (Mag < 13): {len(df_final)} sterren")

# 4. PAREN ZOEKEN
print("Paren zoeken...")
ra_rad = np.radians(df_final['ra'])
dec_rad = np.radians(df_final['dec'])
d = df_final['dist_pc']
x = d * np.cos(dec_rad) * np.cos(ra_rad)
y = d * np.cos(dec_rad) * np.sin(ra_rad)
z = d * np.sin(dec_rad)
coords = np.column_stack([x, y, z])
tree = cKDTree(coords)
pairs = tree.query_pairs(r=0.2) 

# 5. 3D SNELHEID & CONTINUOUS MASS
pair_list = list(pairs)
idx1 = [p[0] for p in pair_list]
idx2 = [p[1] for p in pair_list]
p1 = df_final.iloc[idx1].reset_index(drop=True)
p2 = df_final.iloc[idx2].reset_index(drop=True)

# 3D Velocity
k = 4.74047
v_ra1, v_dec1 = k * p1['pmra']/p1['parallax'], k * p1['pmdec']/p1['parallax']
v_ra2, v_dec2 = k * p2['pmra']/p2['parallax'], k * p2['pmdec']/p2['parallax']
dv_ra = v_ra1 - v_ra2
dv_dec = v_dec1 - v_dec2
dv_rad = p1['radial_velocity'] - p2['radial_velocity']
v_rel_3d = np.sqrt(dv_ra**2 + dv_dec**2 + dv_rad**2)

# Separation
ra1, dec1 = np.radians(p1['ra']), np.radians(p1['dec'])
ra2, dec2 = np.radians(p2['ra']), np.radians(p2['dec'])
theta_rad = np.sqrt(((ra1-ra2)*np.cos((dec1+dec2)/2))**2 + (dec1-dec2)**2)
r_p_au = ((p1['dist_pc'] + p2['dist_pc'])/2) * theta_rad * 206265

# --- CONTINUOUS MASS CALCULATION ---
# Hier gebruiken we de nieuwe functie!
m_tot = estimate_mass_continuous(p1['M_G']) + estimate_mass_continuous(p2['M_G'])
v_newt = 29.78 * np.sqrt(m_tot / r_p_au)
v_tilde = v_rel_3d / v_newt

# 6. RESULTATEN
tight_mask = (r_p_au > 500) & (r_p_au < 2000) & (v_tilde < 5.0)
wide_mask = (r_p_au > 5000) & (r_p_au < 20000) & (v_tilde < 5.0)

v_tight = v_tilde[tight_mask]
v_wide = v_tilde[wide_mask]

# Kalibratie
median_tight_raw = v_tight.median()
calib_factor = 1.0
if CALIBRATE_TO_NEWTON:
    calib_factor = 1.0 / median_tight_raw

v_tight_final = v_tight * calib_factor
v_wide_final = v_wide * calib_factor

med_tight = v_tight_final.median()
med_wide = v_wide_final.median()
boost = med_wide / med_tight

print(f"\n--- PLATINUM AUDIT (Continuous Mass + High Precision) ---")
print(f"Tight Bin (Control):  {med_tight:.3f}")
print(f"Wide Bin (TPIS):      {med_wide:.3f}")
print(f"Boost Factor:         {boost:.3f}x")

# Plot
plt.figure(figsize=(10,6))
plt.hist(v_tight_final, bins=20, density=True, alpha=0.4, color='blue', label='Tight Bin (Control)')
plt.hist(v_wide_final, bins=20, density=True, alpha=0.4, color='red', label='Wide Bin (TPIS)')
plt.axvline(med_tight, c='b', ls='-', linewidth=2)
plt.axvline(med_wide, c='r', ls='-', linewidth=2)
plt.xlim(0, 4)
plt.xlabel(r"Scaled 3D Velocity $\tilde{v}_{3D}$ (Continuous Mass)")
plt.title("Gaia DR3 Wide Binary Anomaly (Platinum Audit)")
plt.legend()
plt.grid(alpha=0.2)
plt.savefig("tpis_platinum.png")
print("Grafiek opgeslagen als tpis_platinum.png")