"""
python/pbpk_helpers.py — shared utilities for the minimal PBPK project.

Classes / functions exported:
    MinimalPBPK         — five-compartment PBPK ODE model
    allometric_scale    — scale physiological parameters from ref BW to target BW
    plot_pk_profile     — matplotlib concentration–time plot
    sensitivity_analysis — one-at-a-time sensitivity on a scalar output
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
from typing import Callable

# Use a non-interactive backend so scripts run headless (no display needed)
matplotlib.use("Agg")

# ── Consistent plot style ────────────────────────────────────────────────────

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": True,
        "grid.color": "0.88",
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 10,
        "figure.dpi": 150,
    }
)


# ── Physiological reference values (ICRP 2002, 70 kg adult) ─────────────────

ICRP_ADULT = {
    # Organ volumes as fraction of body weight
    "fV_lung":   0.0076,
    "fV_liver":  0.0260,
    "fV_kidney": 0.0044,
    "fV_adipose":0.0214,
    # Venous blood pool: remaining after above + red blood cells (~8% BW)
    "fV_blood":  0.0790,
    # Cardiac output coefficient: QC = QC_coeff * BW^QC_exp (L/hr)
    "QC_coeff":  15.0,
    "QC_exp":    0.74,
    # Organ flows as fraction of cardiac output
    "fQ_liver":  0.250,
    "fQ_kidney": 0.190,
    "fQ_adipose":0.050,
    # Remainder goes to the 'rest' compartment
}


# ── Allometric scaling ────────────────────────────────────────────────────────


def allometric_scale(ref_params: dict, ref_bw: float, target_bw: float) -> dict:
    """
    Scale physiological parameters from a reference body weight to a target.

    Volumes scale linearly with BW (exponent 1.0).
    Blood flows scale allometrically with BW^0.74 (Lindstedt 1981).
    Drug-specific parameters (Kp, CLint, fu) are unchanged — these are
    physicochemical properties of the compound, not of the organism.

    Parameters
    ----------
    ref_params : dict  Physiological parameters at ref_bw (from MinimalPBPK.build_params).
    ref_bw     : float Reference body weight (kg).
    target_bw  : float Target body weight (kg).

    Returns
    -------
    dict with scaled parameters.
    """
    scaled = ref_params.copy()
    bw_ratio = target_bw / ref_bw

    # Volumes: linear with BW
    for key in ["V_lung", "V_liver", "V_kidney", "V_adipose", "V_blood", "V_rest"]:
        if key in scaled:
            scaled[key] = ref_params[key] * bw_ratio

    # Blood flows: allometric (exponent 0.74 on BW, same as cardiac output)
    for key in ["QC", "Q_liver", "Q_kidney", "Q_adipose", "Q_rest"]:
        if key in scaled:
            scaled[key] = ref_params[key] * (bw_ratio ** ICRP_ADULT["QC_exp"])

    # Hepatic intrinsic clearance: scales with liver mass (BW^0.75 per Mahmood)
    # This reflects that metabolic enzyme activity scales with liver size.
    if "CLint" in scaled:
        scaled["CLint"] = ref_params["CLint"] * (bw_ratio ** 0.75)

    scaled["BW"] = target_bw
    return scaled


# ── Minimal PBPK model ───────────────────────────────────────────────────────


class MinimalPBPK:
    """
    Five-compartment minimal PBPK model.

    Compartments (index → name):
        0: venous blood (sampling compartment)
        1: lung (arterialises blood)
        2: liver (elimination via CLint)
        3: kidney (distribution only — elimination handled separately if needed)
        4: rest  (all remaining tissues lumped)

    Elimination: hepatic only (well-stirred liver model).
    Route: IV bolus into venous blood (dose added to y[0] at t=0).

    Parameters
    ----------
    params : dict
        Must contain: V_blood, V_lung, V_liver, V_kidney, V_rest,
                      QC, Q_liver, Q_kidney, Q_rest,
                      Kp_lung, Kp_liver, Kp_kidney, Kp_rest,
                      CLint, fu, BP.
    """

    COMPARTMENTS = ["venous", "lung", "liver", "kidney", "rest"]

    def __init__(self, params: dict) -> None:
        self.params = params
        self._validate()

    def _validate(self) -> None:
        required = [
            "V_blood", "V_lung", "V_liver", "V_kidney", "V_rest",
            "QC", "Q_liver", "Q_kidney", "Q_rest",
            "Kp_lung", "Kp_liver", "Kp_kidney", "Kp_rest",
            "CLint", "fu", "BP",
        ]
        missing = [k for k in required if k not in self.params]
        if missing:
            raise ValueError(f"MinimalPBPK: missing parameters: {missing}")

    @staticmethod
    def build_params(
        bw: float = 70.0,
        CLint: float = 50.0,   # μL/min/mg protein → scaled to L/hr below
        fu: float = 0.10,
        BP: float = 1.20,
        Kp_lung: float = 2.0,
        Kp_liver: float = 5.0,
        Kp_kidney: float = 3.0,
        Kp_rest: float = 1.5,
        icrp: dict = ICRP_ADULT,
    ) -> dict:
        """
        Construct a parameter dict for a given body weight from ICRP fractions.

        CLint units: μL/min/mg microsomal protein.
        MPPGL (microsomal protein per gram liver) = 40 mg/g (typical human).
        Conversion: CLint_invivo = CLint_invitro [μL/min/mg] × MPPGL [mg/g]
                                   × liver_weight [g] × 60/1e6 [L/hr per μL/min]
        """
        QC      = icrp["QC_coeff"] * bw ** icrp["QC_exp"]  # L/hr
        Q_liver = icrp["fQ_liver"]  * QC
        Q_kidney= icrp["fQ_kidney"] * QC
        Q_adipose = icrp["fQ_adipose"] * QC
        Q_rest  = QC - Q_liver - Q_kidney - Q_adipose

        V_liver_L  = icrp["fV_liver"]  * bw  # L
        V_liver_g  = V_liver_L * 1000         # g (density ≈ 1 g/mL)

        MPPGL      = 40.0  # mg microsomal protein per g liver
        # CLint invitro → invivo via IVIVE: units L/hr
        CLint_invivo = CLint * MPPGL * V_liver_g * 60.0 / 1e6

        return {
            "BW":        bw,
            "QC":        QC,
            "Q_liver":   Q_liver,
            "Q_kidney":  Q_kidney,
            "Q_adipose": Q_adipose,
            "Q_rest":    Q_rest,
            "V_blood":   icrp["fV_blood"]   * bw,
            "V_lung":    icrp["fV_lung"]    * bw,
            "V_liver":   V_liver_L,
            "V_kidney":  icrp["fV_kidney"]  * bw,
            "V_adipose": icrp["fV_adipose"] * bw,
            "V_rest":    bw * (1 - sum(icrp[k] for k in
                              ["fV_lung","fV_liver","fV_kidney","fV_adipose","fV_blood"])),
            "Kp_lung":   Kp_lung,
            "Kp_liver":  Kp_liver,
            "Kp_kidney": Kp_kidney,
            "Kp_rest":   Kp_rest,
            "CLint":     CLint_invivo,
            "fu":        fu,
            "BP":        BP,
        }

    def odes(self, t: float, y: np.ndarray) -> list[float]:
        """
        ODE right-hand side for scipy.integrate.solve_ivp.

        y[0] = C_venous   (blood total concentration, mg/L)
        y[1] = C_lung
        y[2] = C_liver
        y[3] = C_kidney
        y[4] = C_rest
        """
        p = self.params
        C_ven, C_lung, C_liver, C_kidney, C_rest = y

        # Arterial concentration = blood leaving lung
        C_art = C_lung / p["Kp_lung"]

        # Well-stirred liver model: CLh accounts for limited extraction
        # CLh = QH * fuB * CLint / (QH + fuB * CLint)
        # fuB (free fraction in blood) = fu / BP
        fuB  = p["fu"] / p["BP"]
        QH   = p["Q_liver"]
        CLh  = QH * fuB * p["CLint"] / (QH + fuB * p["CLint"])

        dC_ven    = (
            p["Q_liver"]  * C_liver  / p["Kp_liver"]
            + p["Q_kidney"] * C_kidney / p["Kp_kidney"]
            + p["Q_rest"]   * C_rest   / p["Kp_rest"]
            - p["QC"]       * C_ven
        ) / p["V_blood"]

        dC_lung   = (p["QC"] * C_ven - p["QC"] * C_lung / p["Kp_lung"]) / p["V_lung"]

        dC_liver  = (
            p["Q_liver"] * C_art
            - p["Q_liver"] * C_liver / p["Kp_liver"]
            - CLh * C_liver / p["Kp_liver"]
        ) / p["V_liver"]

        dC_kidney = (
            p["Q_kidney"] * C_art - p["Q_kidney"] * C_kidney / p["Kp_kidney"]
        ) / p["V_kidney"]

        dC_rest   = (
            p["Q_rest"] * C_art - p["Q_rest"] * C_rest / p["Kp_rest"]
        ) / p["V_rest"]

        return [dC_ven, dC_lung, dC_liver, dC_kidney, dC_rest]

    def solve(
        self,
        t_span: tuple[float, float],
        t_eval: np.ndarray,
        dose_mg: float,
        route: str = "iv",
    ) -> "scipy.integrate.OdeSolution":
        """
        Solve the PBPK model for a single IV bolus or first-order oral input.

        For IV: dose is added directly to venous blood compartment at t = 0.
        For oral: a depot compartment is prepended; Ka is assumed constant.

        Parameters
        ----------
        t_span   : (t_start, t_end) in hours.
        t_eval   : time points at which to evaluate (hours).
        dose_mg  : total dose in mg.
        route    : 'iv' (default) or 'oral' (uses Ka = 1.0 1/hr).

        Returns
        -------
        scipy OdeSolution object with .t and .y attributes.
        """
        p = self.params

        if route == "iv":
            # Initial concentration in venous blood = dose / V_blood
            C0_ven = dose_mg / p["V_blood"]
            y0 = [C0_ven, 0.0, 0.0, 0.0, 0.0]
            sol = solve_ivp(
                fun    = self.odes,
                t_span = t_span,
                y0     = y0,
                t_eval = t_eval,
                method = "LSODA",   # auto-switches between stiff/non-stiff — appropriate
                                    # for PBPK where fast lung equilibration creates stiffness
                rtol   = 1e-6,
                atol   = 1e-8,
                dense_output=False,
            )
        elif route == "oral":
            Ka = 1.0  # 1/hr (illustrative; could be parameterised)
            F  = 0.85  # assumed oral bioavailability

            def odes_oral(t, y):
                depot = y[0]
                state = y[1:]
                absorption_rate = Ka * depot * F
                dydt_pbpk = self.odes(t, state)
                # Absorbed drug enters venous blood
                dydt_pbpk[0] += absorption_rate / p["V_blood"]
                return [-Ka * depot] + dydt_pbpk

            y0 = [dose_mg, 0.0, 0.0, 0.0, 0.0, 0.0]  # depot + 5 compartments
            sol_oral = solve_ivp(
                fun    = odes_oral,
                t_span = t_span,
                y0     = y0,
                t_eval = t_eval,
                method = "LSODA",
                rtol   = 1e-6,
                atol   = 1e-8,
            )
            # Strip depot row so shape matches IV output
            from types import SimpleNamespace
            sol = SimpleNamespace(
                t      = sol_oral.t,
                y      = sol_oral.y[1:],  # rows 1-5 are the PBPK compartments
                success= sol_oral.success,
                message= sol_oral.message,
            )
        else:
            raise ValueError(f"route must be 'iv' or 'oral', got '{route}'")

        return sol

    def to_dataframe(self, sol) -> pd.DataFrame:
        """Convert scipy solution to a tidy DataFrame."""
        df = pd.DataFrame(
            sol.y.T,
            columns=["C_venous", "C_lung", "C_liver", "C_kidney", "C_rest"],
        )
        df.insert(0, "time", sol.t)
        # Compute total plasma concentration (venous / BP to convert blood → plasma)
        df["C_plasma"] = df["C_venous"] / self.params["BP"]
        return df


# ── Plotting ─────────────────────────────────────────────────────────────────


def plot_pk_profile(
    df: pd.DataFrame,
    tissue: str = "C_plasma",
    obs: pd.DataFrame | None = None,
    label: str = "Predicted",
    title: str = "PK Profile",
    ax: plt.Axes | None = None,
    colour: str = "#2171b5",
    save_path: str | None = None,
) -> plt.Axes:
    """
    Plot a PBPK predicted concentration–time profile.

    Parameters
    ----------
    df        : DataFrame from MinimalPBPK.to_dataframe().
    tissue    : Column name to plot (default 'C_plasma').
    obs       : Optional DataFrame with 'time' and 'conc' columns for observations.
    label     : Legend label for the predicted line.
    title     : Plot title.
    ax        : Existing Axes to draw on; creates new figure if None.
    colour    : Line colour.
    save_path : If given, saves PNG and SVG (by replacing extension).

    Returns
    -------
    matplotlib Axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(df["time"], df[tissue], colour=colour, linewidth=2, label=label)

    if obs is not None and not obs.empty:
        ax.scatter(obs["time"], obs["conc"],
                   color="black", zorder=5, s=40, label="Observed", alpha=0.75)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Concentration (mg/L)")
    ax.set_title(title)
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
    plt.tight_layout()

    if save_path is not None:
        png_path = save_path if save_path.endswith(".png") else save_path + ".png"
        svg_path = png_path.replace(".png", ".svg")
        plt.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.savefig(svg_path, bbox_inches="tight")

    return ax


# ── Sensitivity analysis ──────────────────────────────────────────────────────


def sensitivity_analysis(
    base_params: dict,
    param_name: str,
    factors: list[float],
    dose_mg: float = 70.0,
    t_span: tuple[float, float] = (0.0, 72.0),
    n_eval: int = 200,
    metric: str = "AUC",
) -> pd.DataFrame:
    """
    One-at-a-time sensitivity analysis on a PBPK scalar output.

    Parameters
    ----------
    base_params : Reference parameter dict from MinimalPBPK.build_params().
    param_name  : Name of the parameter to vary.
    factors     : Multiplicative factors to apply (e.g. [0.5, 0.75, 1.0, 1.25, 1.5]).
    dose_mg     : IV bolus dose in mg.
    t_span      : Simulation time window (hours).
    n_eval      : Number of evaluation time points.
    metric      : 'AUC' (trapezoidal) or 'Cmax'.

    Returns
    -------
    pd.DataFrame with columns: factor, param_value, AUC_or_Cmax.
    """
    if param_name not in base_params:
        raise ValueError(f"'{param_name}' not found in base_params")

    t_eval = np.linspace(t_span[0], t_span[1], n_eval)
    rows   = []

    for f in factors:
        test_params                = base_params.copy()
        test_params[param_name]    = base_params[param_name] * f
        model                      = MinimalPBPK(test_params)
        sol                        = model.solve(t_span, t_eval, dose_mg)
        df                         = model.to_dataframe(sol)

        if metric == "AUC":
            val = float(np.trapz(df["C_plasma"], df["time"]))
        elif metric == "Cmax":
            val = float(df["C_plasma"].max())
        else:
            raise ValueError(f"metric must be 'AUC' or 'Cmax', got '{metric}'")

        rows.append(
            {
                "factor":      f,
                "param_value": test_params[param_name],
                metric:        val,
            }
        )

    return pd.DataFrame(rows)
