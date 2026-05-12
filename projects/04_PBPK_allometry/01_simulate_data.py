"""
01_simulate_data.py — Project 04: Minimal PBPK simulation
Generates synthetic adult IV PK observations from the PBPK model and
saves them with realistic lognormal noise for use in 02_fit_model.py.

═══════════════════════════════════════════════════════════════════════════════
TRUE PARAMETER VALUES
═══════════════════════════════════════════════════════════════════════════════
Physiological (ICRP 2002 reference, 70 kg adult):
  QC   = 15 * 70^0.74 ≈ 87 L/hr   (cardiac output)
  Organ flows: QL=25%, QK=19%, Qrest=56% of QC
  Organ volumes: V_liver=1.82 L, V_kidney=0.308 L, V_lung=0.532 L

Drug-specific (assumed in-silico compound):
  CLint   = 50 μL/min/mg prot  (in vitro hepatic intrinsic clearance)
  fu      = 0.10               (free fraction in plasma)
  BP      = 1.20               (blood:plasma concentration ratio)
  Kp_lung = 2.0, Kp_liver = 5.0, Kp_kidney = 3.0, Kp_rest = 1.5

Route: IV bolus 1 mg/kg  (single dose)
Noise: lognormal, CV = 15%
Subjects: n=10 (population-level noise to simulate study variability)
          + one mean-parameter reference simulation
═══════════════════════════════════════════════════════════════════════════════
"""

import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join('..', '..', 'python'))
from pbpk_helpers import MinimalPBPK, allometric_scale, ICRP_ADULT

np.random.seed(2027)

# ── True parameters (adult 70 kg reference) ───────────────────────────────────
BW_ADULT = 70.0      # kg
DOSE_MGKG = 1.0      # mg/kg IV bolus
N_SUBJECTS = 10
CV_OBS = 0.15        # observational noise (lognormal)

# In-vitro drug parameters (compound-specific — not scaled allometrically)
CLINT_INVITRO = 50.0   # μL/min/mg microsomal protein
FU_PLASMA     = 0.10
BP_RATIO      = 1.20
KP_LUNG       = 2.0
KP_LIVER      = 5.0
KP_KIDNEY     = 3.0
KP_REST       = 1.5

# ── Build reference adult parameters ─────────────────────────────────────────
adult_params = MinimalPBPK.build_params(
    bw        = BW_ADULT,
    CLint     = CLINT_INVITRO,
    fu        = FU_PLASMA,
    BP        = BP_RATIO,
    Kp_lung   = KP_LUNG,
    Kp_liver  = KP_LIVER,
    Kp_kidney = KP_KIDNEY,
    Kp_rest   = KP_REST,
)

print(f"[01_simulate_data.py] Adult PBPK parameters:")
for k in ['QC', 'Q_liver', 'Q_kidney', 'V_liver', 'V_kidney', 'CLint']:
    print(f"  {k:12s} = {adult_params[k]:.4f}")

# ── Sampling times ────────────────────────────────────────────────────────────
t_span  = (0.0, 72.0)
t_dense = np.linspace(0.001, 72.0, 500)   # dense for integration accuracy
t_obs   = np.array([0.083, 0.25, 0.5, 1, 2, 4, 8, 12, 24, 36, 48, 72])

# ── Reference simulation (population mean parameters) ────────────────────────
model_ref = MinimalPBPK(adult_params)
sol_ref   = model_ref.solve(t_span, t_dense, dose_mg = DOSE_MGKG * BW_ADULT)
df_ref    = model_ref.to_dataframe(sol_ref)

# Sample at observed time points
idx_obs   = np.searchsorted(df_ref['time'].values, t_obs)
idx_obs   = np.clip(idx_obs, 0, len(df_ref) - 1)
ref_obs   = df_ref.iloc[idx_obs].copy().reset_index(drop=True)
ref_obs['time_obs'] = t_obs

print(f"[01_simulate_data.py] Reference AUC0-72: {np.trapz(df_ref['C_plasma'], df_ref['time']):.4f} mg·h/L")
print(f"[01_simulate_data.py] Reference Cmax:    {df_ref['C_plasma'].max():.4f} mg/L")

# ── Population simulations (IIV on CLint ±20% CV lognormal) ──────────────────
omega_CLint = 0.04   # variance (CV ≈ 20%)
all_sims    = []

for subj_i in range(N_SUBJECTS):
    params_i = adult_params.copy()
    # Subject-level variability on hepatic CLint only (most important driver)
    params_i['CLint'] = adult_params['CLint'] * np.exp(
        np.random.normal(0, np.sqrt(omega_CLint))
    )
    model_i = MinimalPBPK(params_i)
    sol_i   = model_i.solve(t_span, t_dense, dose_mg = DOSE_MGKG * BW_ADULT)
    df_i    = model_i.to_dataframe(sol_i)

    # Sample at obs times + add lognormal noise
    idx_i   = np.searchsorted(df_i['time'].values, t_obs)
    idx_i   = np.clip(idx_i, 0, len(df_i) - 1)
    obs_i   = df_i.iloc[idx_i].copy().reset_index(drop=True)
    obs_i['time_obs']  = t_obs
    obs_i['ID']        = subj_i + 1
    obs_i['BW']        = BW_ADULT
    obs_i['DOSE_mg']   = DOSE_MGKG * BW_ADULT
    obs_i['C_obs']     = obs_i['C_plasma'] * np.exp(
        np.random.normal(0, CV_OBS, len(obs_i))
    )
    obs_i['C_obs']     = np.maximum(obs_i['C_obs'], 0)
    all_sims.append(obs_i)

pop_df = pd.concat(all_sims, ignore_index=True)

# ── Pediatric prediction (6yo, 25 kg) — using allometric scaling ─────────────
BW_PED = 25.0   # kg, 6 years old
ped_params = allometric_scale(adult_params, ref_bw=BW_ADULT, target_bw=BW_PED)

model_ped  = MinimalPBPK(ped_params)
sol_ped    = model_ped.solve(t_span, t_dense, dose_mg = DOSE_MGKG * BW_PED)
df_ped     = model_ped.to_dataframe(sol_ped)

AUC_adult = np.trapz(df_ref['C_plasma'], df_ref['time'])
AUC_ped   = np.trapz(df_ped['C_plasma'],  df_ped['time'])
AUC_ratio = AUC_ped / AUC_adult

print(f"[01_simulate_data.py] Pediatric (25 kg) AUC at same mg/kg dose:")
print(f"  Adult AUC:    {AUC_adult:.4f} mg·h/L")
print(f"  Pediatric AUC:{AUC_ped:.4f} mg·h/L  (ratio = {AUC_ratio:.3f})")

# ── Save data ─────────────────────────────────────────────────────────────────
os.makedirs('data', exist_ok=True)

# Adult population observations (for fitting CLint)
adult_obs = pop_df[['ID','time_obs','C_obs','BW','DOSE_mg']].copy()
adult_obs.columns = ['ID','time','conc','BW','dose_mg']
adult_obs.to_csv('data/adult_obs.csv', index=False)

# Reference PK profile (dense, for plotting)
df_ref['ID']    = 0
df_ref['BW']    = BW_ADULT
df_ref.to_csv('data/adult_ref_pk.csv', index=False)

# Pediatric predicted profile
df_ped['ID']   = 999
df_ped['BW']   = BW_PED
df_ped.to_csv('data/pediatric_pred.csv', index=False)

# Summary of true parameter values for reference
summary_df = pd.DataFrame([
    {'Parameter': 'BW_adult',    'Value': BW_ADULT,      'Units': 'kg'},
    {'Parameter': 'BW_pediatric','Value': BW_PED,        'Units': 'kg'},
    {'Parameter': 'CLint_invitro','Value': CLINT_INVITRO,'Units': 'uL/min/mg'},
    {'Parameter': 'fu',          'Value': FU_PLASMA,     'Units': ''},
    {'Parameter': 'BP',          'Value': BP_RATIO,      'Units': ''},
    {'Parameter': 'Kp_liver',    'Value': KP_LIVER,      'Units': ''},
    {'Parameter': 'QC_adult',    'Value': adult_params['QC'],  'Units': 'L/hr'},
    {'Parameter': 'CLint_invivo','Value': adult_params['CLint'],'Units': 'L/hr'},
    {'Parameter': 'AUC_adult',   'Value': round(AUC_adult, 4), 'Units': 'mg·h/L'},
    {'Parameter': 'AUC_ped_same_mgkg','Value': round(AUC_ped, 4),'Units': 'mg·h/L'},
    {'Parameter': 'AUC_ratio',   'Value': round(AUC_ratio, 4), 'Units': 'ped/adult'},
])
summary_df.to_csv('data/true_parameters.csv', index=False)

print(f"[01_simulate_data.py] Saved: adult_obs.csv ({len(adult_obs)} rows), "
      f"adult_ref_pk.csv, pediatric_pred.csv")
