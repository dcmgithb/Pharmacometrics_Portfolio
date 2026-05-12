"""
02_fit_model.py — Project 04: Fit CLint to adult PBPK observations
Uses scipy.optimize.minimize (Nelder-Mead) to estimate hepatic intrinsic
clearance from the simulated adult IV PK observations.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join('..', '..', 'python'))
from pbpk_helpers import MinimalPBPK, allometric_scale, ICRP_ADULT, plot_pk_profile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(42)
os.makedirs('results', exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
adult_obs = pd.read_csv('data/adult_obs.csv')
ref_pk    = pd.read_csv('data/adult_ref_pk.csv')
ped_pk    = pd.read_csv('data/pediatric_pred.csv')
true_pars = pd.read_csv('data/true_parameters.csv').set_index('Parameter')

BW_ADULT = 70.0
BW_PED   = 25.0
DOSE_MG  = BW_ADULT   # 1 mg/kg × 70 kg

# ── Objective function: weighted SSE on log-scale ────────────────────────────
# Log-scale SSE is appropriate for concentration data spanning orders of magnitude:
# it gives equal weight to relative errors regardless of absolute magnitude.

def objective(log_CLint, obs_df):
    CLint_val = np.exp(log_CLint[0])
    params = MinimalPBPK.build_params(
        bw=BW_ADULT, CLint=CLint_val,
        fu=0.10, BP=1.20,
        Kp_lung=2.0, Kp_liver=5.0, Kp_kidney=3.0, Kp_rest=1.5
    )
    model = MinimalPBPK(params)
    t_obs = obs_df['time'].values
    try:
        sol   = model.solve((0.0, t_obs.max() + 1.0), t_obs, dose_mg=DOSE_MG)
        df    = model.to_dataframe(sol)
        idx   = np.searchsorted(df['time'].values, t_obs)
        idx   = np.clip(idx, 0, len(df) - 1)
        C_pred = df['C_plasma'].values[idx]
        C_obs  = obs_df['conc'].values
        # Ignore near-zero predictions to avoid log(0)
        valid  = (C_pred > 1e-9) & (C_obs > 1e-9)
        sse    = np.sum((np.log(C_obs[valid]) - np.log(C_pred[valid]))**2)
    except Exception:
        sse = 1e9
    return sse

# Use population mean observations (mean across subjects per timepoint)
mean_obs = adult_obs.groupby('time')['conc'].mean().reset_index()

print("[02_fit_model.py] Fitting CLint by Nelder-Mead minimisation...")
CLint_init = 50.0  # starting guess (μL/min/mg) — true value is 50
result = minimize(
    objective,
    x0      = [np.log(CLint_init)],
    args    = (mean_obs,),
    method  = 'Nelder-Mead',
    options = {'xatol': 1e-4, 'fatol': 1e-6, 'maxiter': 5000}
)

CLint_est = np.exp(result.x[0])
print(f"[02_fit_model.py] CLint estimate: {CLint_est:.2f} μL/min/mg  "
      f"(true: {float(true_pars.loc['CLint_invitro','Value']):.1f})")
print(f"[02_fit_model.py] Optimiser converged: {result.success}  SSE={result.fun:.4f}")

# ── Build fitted model and simulate ──────────────────────────────────────────
adult_params_fit = MinimalPBPK.build_params(
    bw=BW_ADULT, CLint=CLint_est,
    fu=0.10, BP=1.20,
    Kp_lung=2.0, Kp_liver=5.0, Kp_kidney=3.0, Kp_rest=1.5
)
model_fit = MinimalPBPK(adult_params_fit)

t_eval    = np.linspace(0.001, 72.0, 500)
sol_fit   = model_fit.solve((0.0, 72.0), t_eval, dose_mg=DOSE_MG)
df_fit    = model_fit.to_dataframe(sol_fit)

AUC_fit   = np.trapz(df_fit['C_plasma'], df_fit['time'])
Cmax_fit  = df_fit['C_plasma'].max()
AUC_true  = float(true_pars.loc['AUC_adult', 'Value'])

print(f"[02_fit_model.py] Fitted AUC0-72: {AUC_fit:.4f} mg·h/L  (true: {AUC_true:.4f})")

# ── Pediatric dose recommendation ────────────────────────────────────────────
ped_params_fit = allometric_scale(adult_params_fit, ref_bw=BW_ADULT, target_bw=BW_PED)
model_ped_fit  = MinimalPBPK(ped_params_fit)

DOSE_PED_MG_KG = 1.0  # mg/kg — same as adult weight-normalised dose
sol_ped_same   = model_ped_fit.solve((0.0, 72.0), t_eval,
                                      dose_mg=DOSE_PED_MG_KG * BW_PED)
df_ped_same    = model_ped_fit.to_dataframe(sol_ped_same)
AUC_ped_same   = np.trapz(df_ped_same['C_plasma'], df_ped_same['time'])
AUC_ratio_same = AUC_ped_same / AUC_fit

# Adjusted dose to match adult AUC
dose_ped_adj_mgkg = DOSE_PED_MG_KG / AUC_ratio_same
sol_ped_adj        = model_ped_fit.solve((0.0, 72.0), t_eval,
                                          dose_mg=dose_ped_adj_mgkg * BW_PED)
df_ped_adj         = model_ped_fit.to_dataframe(sol_ped_adj)
AUC_ped_adj        = np.trapz(df_ped_adj['C_plasma'], df_ped_adj['time'])

print(f"[02_fit_model.py] AUC ratio (ped/adult, same mg/kg): {AUC_ratio_same:.3f}")
print(f"[02_fit_model.py] Adjusted pediatric dose: {dose_ped_adj_mgkg:.3f} mg/kg")
print(f"[02_fit_model.py] Adjusted pediatric AUC:  {AUC_ped_adj:.4f} mg·h/L  "
      f"(adult: {AUC_fit:.4f})")

# ── Adult fit plot ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(df_fit['time'], df_fit['C_plasma'], 'b-', linewidth=2, label='PBPK fit')
for subj_id in adult_obs['ID'].unique():
    sub = adult_obs[adult_obs['ID'] == subj_id]
    ax.scatter(sub['time'], sub['conc'], color='grey', alpha=0.4, s=25, zorder=3)
ax.scatter(mean_obs['time'], mean_obs['conc'],
           color='black', s=60, zorder=5, label='Mean observed')
ax.set_yscale('log')
ax.set_xlabel('Time (h)')
ax.set_ylabel('Plasma concentration (mg/L)')
ax.set_title('Adult PBPK Model Fit — IV Bolus 1 mg/kg')
ax.legend()
plt.tight_layout()
plt.savefig('results/adult_fit.png', dpi=300, bbox_inches='tight')
plt.savefig('results/adult_fit.svg', bbox_inches='tight')
plt.close()

# ── Pediatric comparison plot ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(df_fit['time'],     df_fit['C_plasma'],     'b-',  linewidth=2,
        label=f'Adult 1 mg/kg (AUC={AUC_fit:.2f})')
ax.plot(df_ped_same['time'], df_ped_same['C_plasma'], 'r--', linewidth=2,
        label=f'Ped 1 mg/kg (AUC={AUC_ped_same:.2f}; ratio={AUC_ratio_same:.2f})')
ax.plot(df_ped_adj['time'],  df_ped_adj['C_plasma'],  'g-',  linewidth=2,
        label=f'Ped {dose_ped_adj_mgkg:.2f} mg/kg (AUC≈adult)')
ax.set_yscale('log')
ax.set_xlabel('Time (h)')
ax.set_ylabel('Plasma concentration (mg/L)')
ax.set_title('Adult vs Pediatric (6yo, 25 kg) — PBPK Predicted PK')
ax.legend()
plt.tight_layout()
plt.savefig('results/pediatric_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('results/pediatric_comparison.svg', bbox_inches='tight')
plt.close()

# ── Save results ───────────────────────────────────────────────────────────────
fit_summary = pd.DataFrame([
    {'Parameter': 'CLint_est_uLminmg', 'Value': round(CLint_est, 2),
     'True': float(true_pars.loc['CLint_invitro','Value']),
     'Units': 'uL/min/mg'},
    {'Parameter': 'AUC_adult_fit',     'Value': round(AUC_fit, 4),
     'True': AUC_true, 'Units': 'mg·h/L'},
    {'Parameter': 'AUC_ratio_ped_adult','Value': round(AUC_ratio_same, 4),
     'True': float(true_pars.loc['AUC_ratio','Value']), 'Units': ''},
    {'Parameter': 'Dose_ped_adj_mgkg', 'Value': round(dose_ped_adj_mgkg, 3),
     'True': None, 'Units': 'mg/kg'},
])
fit_summary.to_csv('results/fit_summary.csv', index=False)

df_fit['source']    = 'adult_fitted'
df_ped_same['source'] = 'ped_1mgkg'
df_ped_adj['source']  = f'ped_{dose_ped_adj_mgkg:.3f}mgkg'
combined = pd.concat([df_fit, df_ped_same, df_ped_adj])
combined.to_csv('results/pk_profiles.csv', index=False)

ped_rec = pd.DataFrame([{
    'BW_kg': BW_PED, 'Age_years': 6,
    'Adult_dose_mgkg': 1.0, 'Adult_AUC': round(AUC_fit, 4),
    'Ped_dose_same_mgkg': 1.0, 'Ped_AUC_same': round(AUC_ped_same, 4),
    'AUC_ratio_same_dose': round(AUC_ratio_same, 4),
    'Ped_dose_adj_mgkg': round(dose_ped_adj_mgkg, 3),
    'Ped_AUC_adj': round(AUC_ped_adj, 4),
}])
ped_rec.to_csv('results/pediatric_dose_recommendation.csv', index=False)

# ── Save recommendation text ───────────────────────────────────────────────────
rec_text = (
    f"## Pediatric Dose Recommendation\n\n"
    f"Fitted hepatic CLint: **{CLint_est:.1f} μL/min/mg** "
    f"(true: {float(true_pars.loc['CLint_invitro','Value']):.0f})\n\n"
    f"At the same mg/kg dose as the adult (1 mg/kg), the 6-year-old (25 kg) "
    f"child has a predicted AUC ratio of **{AUC_ratio_same:.2f}** "
    f"(pediatric/adult), driven by lower allometrically-scaled hepatic "
    f"blood flow and CLint.\n\n"
    f"To achieve AUC equivalence with the adult, the recommended pediatric "
    f"dose is **{dose_ped_adj_mgkg:.2f} mg/kg** "
    f"({dose_ped_adj_mgkg * BW_PED:.1f} mg absolute for a 25 kg child).\n\n"
    f"This recommendation should be validated in a dedicated pediatric PK "
    f"study. Ontogeny corrections for CYP enzyme maturation were not applied."
)
with open('results/pediatric_recommendation.txt', 'w') as f:
    f.write(rec_text)

print("[02_fit_model.py] Done. Results saved to results/")
