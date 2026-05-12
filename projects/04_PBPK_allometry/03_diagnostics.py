"""
03_diagnostics.py — Project 04: PBPK model diagnostics
Produces GOF plot, sensitivity analysis (CLint and fu), and sensitivity table.
"""

import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join('..', '..', 'python'))
from pbpk_helpers import MinimalPBPK, allometric_scale, sensitivity_analysis

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

os.makedirs('results', exist_ok=True)

# ── Load fit results ──────────────────────────────────────────────────────────
adult_obs  = pd.read_csv('data/adult_obs.csv')
pk_profiles= pd.read_csv('results/pk_profiles.csv')
fit_summary= pd.read_csv('results/fit_summary.csv').set_index('Parameter')

CLint_est  = float(fit_summary.loc['CLint_est_uLminmg', 'Value'])
BW_ADULT   = 70.0
DOSE_MG    = BW_ADULT

adult_params = MinimalPBPK.build_params(
    bw=BW_ADULT, CLint=CLint_est,
    fu=0.10, BP=1.20,
    Kp_lung=2.0, Kp_liver=5.0, Kp_kidney=3.0, Kp_rest=1.5
)

t_eval  = np.linspace(0.001, 72.0, 500)
model   = MinimalPBPK(adult_params)
sol     = model.solve((0.0, 72.0), t_eval, dose_mg=DOSE_MG)
df_pred = model.to_dataframe(sol)

# ── 1. GOF: predicted vs observed ────────────────────────────────────────────
mean_obs = adult_obs.groupby('time')['conc'].mean().reset_index()
t_obs    = mean_obs['time'].values
idx      = np.clip(np.searchsorted(df_pred['time'].values, t_obs), 0, len(df_pred)-1)
C_pred_at_obs = df_pred['C_plasma'].values[idx]

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Predicted vs observed
axes[0].scatter(C_pred_at_obs, mean_obs['conc'], color='steelblue', s=60, zorder=5)
lim = max(C_pred_at_obs.max(), mean_obs['conc'].max()) * 1.1
axes[0].plot([0, lim], [0, lim], 'r-', linewidth=1)
axes[0].set_xlabel('Predicted C_plasma (mg/L)')
axes[0].set_ylabel('Observed mean concentration (mg/L)')
axes[0].set_title('Predicted vs Observed')
axes[0].set_xscale('log'); axes[0].set_yscale('log')

# Residuals vs time
residuals = np.log(mean_obs['conc'].values) - np.log(C_pred_at_obs)
axes[1].scatter(t_obs, residuals, color='steelblue', s=60)
axes[1].axhline(0, color='tomato', linewidth=1)
axes[1].axhline( 1.96 * mean_obs['conc'].std() / mean_obs['conc'].mean(),
                 linestyle='--', color='grey', linewidth=0.8)
axes[1].axhline(-1.96 * mean_obs['conc'].std() / mean_obs['conc'].mean(),
                 linestyle='--', color='grey', linewidth=0.8)
axes[1].set_xlabel('Time (h)')
axes[1].set_ylabel('log(Obs) – log(Pred)')
axes[1].set_title('Log-scale Residuals vs Time')

plt.suptitle('PBPK Model GOF — Adult IV Data', fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig('results/gof_adult.png', dpi=300, bbox_inches='tight')
plt.savefig('results/gof_adult.svg', bbox_inches='tight')
plt.close()
print("[03_diagnostics.py] GOF plot saved.")

# ── 2. Sensitivity analysis ───────────────────────────────────────────────────
# One-at-a-time (OAT): vary CLint and fu ±50% in steps, compute AUC
factors = [0.50, 0.75, 1.00, 1.25, 1.50]

sa_CLint = sensitivity_analysis(adult_params, 'CLint', factors,
                                 dose_mg=DOSE_MG, t_span=(0.0, 72.0), metric='AUC')
sa_CLint['Parameter'] = 'CLint'
sa_CLint['Factor_pct'] = [f'{f*100:.0f}%' for f in factors]

sa_fu = sensitivity_analysis(adult_params, 'fu', factors,
                              dose_mg=DOSE_MG, t_span=(0.0, 72.0), metric='AUC')
sa_fu['Parameter'] = 'fu'
sa_fu['Factor_pct'] = [f'{f*100:.0f}%' for f in factors]

# Kp_liver sensitivity
sa_Kp = sensitivity_analysis(adult_params, 'Kp_liver', factors,
                              dose_mg=DOSE_MG, t_span=(0.0, 72.0), metric='AUC')
sa_Kp['Parameter'] = 'Kp_liver'
sa_Kp['Factor_pct'] = [f'{f*100:.0f}%' for f in factors]

sa_all = pd.concat([sa_CLint, sa_fu, sa_Kp], ignore_index=True)

# Express as % change in AUC relative to base (factor=1.0)
base_auc = sa_CLint.loc[sa_CLint['factor'] == 1.0, 'AUC'].values[0]
sa_all['AUC_pct_change'] = (sa_all['AUC'] / base_auc - 1.0) * 100

sa_table = sa_all[['Parameter', 'Factor_pct', 'AUC', 'AUC_pct_change']].copy()
sa_table.columns = ['Parameter', 'Factor', 'AUC (mg·h/L)', 'AUC change (%)']
sa_table = sa_table.round({'AUC (mg·h/L)': 3, 'AUC change (%)': 1})
sa_table.to_csv('results/sensitivity_table.csv', index=False)

# Tornado plot
fig, ax = plt.subplots(figsize=(8, 5))
params_list = ['CLint', 'fu', 'Kp_liver']
colours     = ['steelblue', 'darkorange', 'mediumseagreen']
y_pos       = np.arange(len(params_list))

for i, (param, col) in enumerate(zip(params_list, colours)):
    sub = sa_all[sa_all['Parameter'] == param].sort_values('factor')
    lo  = sub.loc[sub['factor'] == 0.50, 'AUC_pct_change'].values[0]
    hi  = sub.loc[sub['factor'] == 1.50, 'AUC_pct_change'].values[0]
    ax.barh(y_pos[i], hi - lo, left=lo, height=0.5, color=col, alpha=0.8,
            label=param)
    ax.text(hi + 0.5, y_pos[i], f'+50%: {hi:+.1f}%', va='center', fontsize=9)
    ax.text(lo - 0.5, y_pos[i], f'−50%: {lo:+.1f}%', va='center', ha='right', fontsize=9)

ax.axvline(0, color='black', linewidth=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(params_list)
ax.set_xlabel('AUC change from base (%)')
ax.set_title('Sensitivity of AUC to ±50% Parameter Variation')
plt.tight_layout()
plt.savefig('results/sensitivity_tornado.png', dpi=300, bbox_inches='tight')
plt.savefig('results/sensitivity_tornado.svg', bbox_inches='tight')
plt.close()
print("[03_diagnostics.py] Sensitivity analysis complete.")

# ── 3. Summary stats ──────────────────────────────────────────────────────────
rmse_log = np.sqrt(np.mean(
    (np.log(mean_obs['conc'].values) - np.log(C_pred_at_obs))**2
))
mpe      = np.mean((C_pred_at_obs - mean_obs['conc'].values) / mean_obs['conc'].values * 100)
print(f"[03_diagnostics.py] RMSE(log): {rmse_log:.4f}")
print(f"[03_diagnostics.py] MPE:       {mpe:+.2f}%")

stats_df = pd.DataFrame([{'RMSE_log': round(rmse_log, 4), 'MPE_pct': round(mpe, 2)}])
stats_df.to_csv('results/gof_stats.csv', index=False)

print("[03_diagnostics.py] Project 04 diagnostics complete.")
