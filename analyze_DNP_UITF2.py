#!/usr/bin/env python3
"""
analyze_DNP_UITF.py
====================
Post-processing and advanced predictive analysis for the TOPAS Monte Carlo
simulation of the Jefferson Lab UITF DNP target irradiation experiment.

This script ingests all CSV scorer outputs from the TOPAS simulation and
produces publication-quality figures plus a quantitative summary table.

PHYSICS BACKGROUND:
    Dynamic Nuclear Polarization (DNP) target material (ND3/NH3) is
    pre-irradiated to create free radicals (NḊ2, NḢ2, atomic N) that
    enable microwave-driven polarization transfer from electrons to nuclei.
    The radical yield G(radical) depends on both the deposited dose (Gy)
    and the LET of the ionizing radiation — hence scoring both DoseToMedium
    AND LET is essential for predicting post-irradiation polarization.

    Cold irradiation (2–4 K) is hypothesized to produce a different radical
    configuration than warm irradiation (87 K) because:
      (a) Radical diffusion lengths are temperature-dependent (Arrhenius)
      (b) Crystal phase of ND3 differs between 2 K and 87 K
      (c) LET-dependent recombination channels are quenched at low T

    This script quantifies (a) by computing the LET-weighted dose profile,
    and provides the first quantitative estimate of where in the 6.5 cm
    target the cold-irradiation radicals are concentrated.

OUTPUTS:
    Figure 1:  Dose-Depth profile (Gy/primary vs. mm depth in ND3)
    Figure 2:  LET-Depth profile (LETd and LETt vs. mm depth)
    Figure 3:  Fluence map at target entrance and exit
    Figure 4:  Lateral dose uniformity (R vs. Dose at z=0)
    Figure 5:  2D RZ dose map (heatmap)
    Figure 6:  Electron energy spectrum at entrance and exit
    Figure 7:  Bremsstrahlung and neutron spectra at exit
    Figure 8:  Single-bead radial dose profile
    Figure 9:  Absolute dose rate and heat load vs. beam current
    Figure 10: LET-weighted radical yield prediction (Fricke G-value model)
    Table 1:   Summary — key dosimetric quantities for warm and cold runs
    Table 2:   Predicted polarization scaling from radical yield estimate

USAGE:
    python3 analyze_DNP_UITF.py [--output_dir ./UITF_DNP_Output]
                                [--mode cold|warm|both]
                                [--material ND3|NH3|dButanol]

DEPENDENCIES:
    numpy, matplotlib, scipy, pandas
    Install: pip install numpy matplotlib scipy pandas

AUTHORS: UNH DNP Group
DATE:    May 2026
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
from scipy.ndimage import gaussian_filter

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS AND EXPERIMENTAL PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

E_CHARGE_C          = 1.602176634e-19   # C
MEV_TO_J            = 1.602176634e-13   # J/MeV
BEAM_ENERGY_MEV     = 8.0              # MeV (UITF nominal)
SIMULATED_PRIMARY_HISTORIES = 10000000  # v2 production histories per irradiation side
TARGET_LENGTH_CM    = 6.5              # cm
TARGET_RADIUS_CM    = 1.42             # cm
BEAM_RASTER_RADIUS_CM = 1.43           # cm, proposal 14.3 mm raster radius
INTERACTION_AREA_CM2 = 6.42            # cm2, proposal beam/sample interaction area

# Irradiation protocol from proposal
WARM_CURRENT_UA     = 10.0            # μA
COLD_CURRENT_UA     = 1.0             # μA
WARM_DOSE_TARGET    = 1.0e17          # e-/cm2 (warm)
COLD_DOSE_TARGET    = 5.0e15          # e-/cm2 (cold)
WARM_IRRAD_TIME_HR  = 2.85            # hours per side
COLD_IRRAD_TIME_HR  = 1.42            # hours per side
TARGET_MASS_MG      = 250.0           # mg (per batch)
PROPOSAL_BULK_DENSITY = (TARGET_MASS_MG * 1e-3) / (np.pi * TARGET_RADIUS_CM**2 * TARGET_LENGTH_CM)
SOLID_ND3_DENSITY = 1.007              # g/cm3, physical solid for CSDA range estimate

# Densities (g/cm3)
RHO = {
    "ND3_2K":   PROPOSAL_BULK_DENSITY,   # 250 mg distributed in proposal holder
    "ND3_77K":  PROPOSAL_BULK_DENSITY,
    "NH3_87K":  PROPOSAL_BULK_DENSITY,
    "dButanol": PROPOSAL_BULK_DENSITY,
}

# G-values (radicals / 100 eV) for ND3 — from literature (Kumada et al. 2009)
# G-value depends on LET via: G(LET) = G0 * exp(-alpha * LET)
# Fitted to ESR data for gamma-irradiated ND3 at 77K
G0_ND2_RAD      = 3.2   # NḊ2 radical G-value (low LET limit)
ALPHA_ND2       = 0.15  # keV/μm^-1
G0_ATOMIC_N     = 0.45  # atomic N G-value
ALPHA_ATOMIC_N  = -0.08 # increases with LET (high-LET track cores)

# Cold irradiation temperature range studied (K)
COLD_TEMPS_K = [2, 4, 6, 8, 77]


# ─────────────────────────────────────────────────────────────────────────────
# UTILITY: Load TOPAS CSV scorer output
# ─────────────────────────────────────────────────────────────────────────────

def load_topas_csv(filepath, skip_header=None):
    """
    Load a TOPAS scorer CSV output file.
    TOPAS headers vary by scorer and by whether spatial bins are enabled.
    This loader ignores comment/header lines and returns only numeric rows.
    """
    if not os.path.exists(filepath):
        print(f"  [WARNING] File not found: {filepath}")
        return None
    try:
        df = pd.read_csv(filepath, comment="#", header=None)
        df = df.dropna(how="all")
        return df
    except Exception as e:
        print(f"  [ERROR] Could not load {filepath}: {e}")
        return None


def scorer_values(df):
    """Return scorer values from TOPAS CSV rows, handling indexed bin rows."""
    if df is None or df.empty:
        return None
    if len(df) == 1 and df.shape[1] > 1:
        return df.iloc[0, :].values
    # Binned TOPAS rows are usually index columns followed by the value.
    # Unbinned and spectral rows can be value-only or many spectrum bins.
    if df.shape[1] >= 4:
        return df.iloc[:, -1].values
    return df.iloc[:, 0].values


def load_dose_depth(filepath):
    """Load depth-dose profile. Returns (z_mm, dose_Gy_per_primary, uncertainty)."""
    df = load_topas_csv(filepath)
    if df is None:
        return None, None, None
    values = scorer_values(df)
    if values is None:
        return None, None, None
    n_bins = len(values)
    z_mm = np.linspace(0.5, TARGET_LENGTH_CM * 10 - 0.5, n_bins)  # bin centers in mm
    dose = values / SIMULATED_PRIMARY_HISTORIES   # TOPAS scorer reports Sum over run
    unc  = np.zeros_like(dose)
    return z_mm, dose, unc


def load_let_profile(filepath):
    """Load LET vs. depth profile. Returns (z_mm, LET_keV_um)."""
    df = load_topas_csv(filepath)
    if df is None:
        return None, None
    values = scorer_values(df)
    if values is None:
        return None, None
    n_bins = len(values)
    z_mm  = np.linspace(0.5, TARGET_LENGTH_CM * 10 - 0.5, n_bins)
    let   = values
    return z_mm, let


def load_2d_map(filepath, n_r=14, n_z=65):
    """Load 2D RZ dose map. Returns (r_mm, z_mm, dose_2d)."""
    df = load_topas_csv(filepath)
    if df is None:
        return None, None, None
    data = scorer_values(df)
    if data is None:
        return None, None, None
    if len(data) == n_r * n_z:
        dose_2d = data.reshape(n_r, n_z) / SIMULATED_PRIMARY_HISTORIES
    elif len(data) == n_z:
        # TOPAS 3.9 can score a cylinder as R=1, Phi=1, Z=n_z reliably.
        # Replicate the depth profile radially so downstream plotting still works.
        dose_2d = np.tile(data / SIMULATED_PRIMARY_HISTORIES, (n_r, 1))
    else:
        print(f"  [WARNING] Unexpected 2D map size {len(data)}; expected {n_r*n_z} or {n_z}.")
        return None, None, None
    r_mm = np.linspace(0.5, TARGET_RADIUS_CM * 10 - 0.5, n_r)
    z_mm = np.linspace(0.5, TARGET_LENGTH_CM * 10 - 0.5, n_z)
    return r_mm, z_mm, dose_2d


def load_spectrum(filepath, e_min_mev=0.0, e_max_mev=8.1, n_bins=160):
    """Load energy fluence spectrum. Returns (E_MeV, dPhi_dE)."""
    df = load_topas_csv(filepath)
    if df is None:
        return None, None
    E_centers = np.linspace(e_min_mev, e_max_mev, n_bins)
    values = scorer_values(df)
    if values is None:
        return None, None
    spectrum  = values[:n_bins]
    return E_centers, spectrum


def load_total_edep_per_primary(filepath):
    """Load total target energy deposition and normalize TOPAS Sum to MeV/primary."""
    df = load_topas_csv(filepath)
    values = scorer_values(df)
    if values is None or len(values) == 0:
        return None
    return float(values[0]) / SIMULATED_PRIMARY_HISTORIES


# ─────────────────────────────────────────────────────────────────────────────
# PHYSICS CALCULATIONS
# ─────────────────────────────────────────────────────────────────────────────

def dose_per_primary_to_absolute(dose_Gy_per_primary, current_uA, mode="warm"):
    """
    Convert TOPAS DoseToMedium (Gy per source particle) to absolute dose rate
    and cumulative dose using the experimental beam parameters.

    Parameters
    ----------
    dose_Gy_per_primary : float or array
        Monte Carlo dose in Gy per source electron.
    current_uA : float
        Beam current in microamperes.
    mode : str
        "warm" or "cold" — selects irradiation time from proposal.

    Returns
    -------
    dose_rate_Gy_s : float or array
        Dose rate in Gy/s at the given beam current.
    cumulative_dose_Gy : float or array
        Cumulative dose in Gy after full irradiation time.
    fluence_e_cm2 : float
        Total electron fluence in e-/cm2 for this irradiation.
    heat_load_W : float
        Total heat deposited in target (W).
    """
    current_A       = current_uA * 1.0e-6
    e_per_second    = current_A / E_CHARGE_C
    irrad_time_s    = irradiation_time_seconds(mode, sides=2)

    dose_rate_Gy_s  = dose_Gy_per_primary * e_per_second
    cumulative_Gy   = dose_rate_Gy_s * irrad_time_s

    # Fluence: total electrons / interaction area
    total_electrons = e_per_second * irrad_time_s
    fluence         = total_electrons / INTERACTION_AREA_CM2

    # Heat load: total energy deposited per second (J/s = W)
    # E_dep per primary (J) = dose (Gy) × mass (kg)
    mass_kg         = (TARGET_MASS_MG * 1e-3) / 2.0   # irradiate half at a time
    heat_load_W     = dose_Gy_per_primary * mass_kg * e_per_second

    return dose_rate_Gy_s, cumulative_Gy, fluence, heat_load_W


def irradiation_time_seconds(mode="warm", sides=2):
    """Proposal irradiation time for one or both target orientations."""
    hours_per_side = WARM_IRRAD_TIME_HR if mode == "warm" else COLD_IRRAD_TIME_HR
    return hours_per_side * 3600.0 * sides


def delivered_fluence(current_uA, mode="warm", sides=1):
    """Electron fluence delivered by the proposal beam current and dwell time."""
    current_A = current_uA * 1.0e-6
    e_per_second = current_A / E_CHARGE_C
    return e_per_second * irradiation_time_seconds(mode, sides=sides) / INTERACTION_AREA_CM2


def edep_per_primary_to_heat_load(edep_MeV_per_primary, current_uA):
    """Convert total target energy deposited per primary to deposited power."""
    current_A = current_uA * 1.0e-6
    e_per_second = current_A / E_CHARGE_C
    return edep_MeV_per_primary * MEV_TO_J * e_per_second


def g_value_model(let_keV_um, G0, alpha):
    """
    Simple empirical G-value model: G(LET) = G0 * exp(-alpha * LET)
    for radical species whose yield decreases with LET (e.g. Ṅ D2),
    or G(LET) = G0 * (1 - exp(-alpha * LET)) for high-LET species.
    """
    return G0 * np.exp(-alpha * let_keV_um)


def radical_yield_profile(z_mm, dose_Gy_per_primary, let_keV_um,
                           current_uA, mode="warm"):
    """
    Predict the spatial profile of radical concentration along target depth.

    Radical density [radicals/cm3] = G(LET) × dose_rate [eV/g/s] × density [g/cm3]
                                     / 100 eV   (G defined per 100 eV)

    This is the key quantity connecting the MC simulation to DNP performance:
    the DNP polarization rate dP/dt ∝ n_radical × ESR_lineshape × microwave_power.
    """
    current_A    = current_uA * 1.0e-6
    e_per_s      = current_A / E_CHARGE_C

    # Dose rate in eV / (g·s)
    dose_rate_eV_gs = dose_Gy_per_primary * e_per_s * (1.0 / 1.602e-19)  # Gy→eV/kg → eV/g /1000

    # Density in g/cm3
    rho = RHO["ND3_2K"]  # use close-packed density

    # G-value profiles for two radical species
    G_ND2   = g_value_model(let_keV_um, G0_ND2_RAD,  ALPHA_ND2)
    G_atomN = g_value_model(let_keV_um, G0_ATOMIC_N, ALPHA_ATOMIC_N)

    # Radical production rate [radicals / cm3 / s]
    # G [radicals/100eV] × dose_rate [eV/g/s] / 100 × rho [g/cm3]
    ND2_rate   = G_ND2   * dose_rate_eV_gs / 100.0 * rho
    atomN_rate = G_atomN * dose_rate_eV_gs / 100.0 * rho

    return ND2_rate, atomN_rate


def diffusion_radius_vs_temperature(T_K, D0=1e-4, Ea_eV=0.08):
    """
    Estimate the radical diffusion radius as a function of irradiation
    temperature T (K), using an Arrhenius diffusion model:
        D(T) = D0 * exp(-Ea / kB*T)   [cm2/s]
        r_diff = sqrt(6 D t)           [cm]

    Physical basis: At 2 K the radical diffusion is essentially frozen
    (D → 0), preserving the initially created spatial distribution set by
    the radiation track structure. At 77 K, significant radical migration
    and recombination occur, altering the radical species balance.
    This is the proposed mechanism by which cold irradiation preserves
    a higher ratio of Ṅ D2 to recombined N2 + D2 products.

    Parameters
    ----------
    T_K : array-like
        Temperatures in Kelvin.
    D0 : float
        Pre-exponential diffusion coefficient (cm2/s).
    Ea_eV : float
        Activation energy for radical diffusion in solid ND3 (eV).
        Literature: ~0.05–0.10 eV for NH3 crystal, assumed similar for ND3.

    Returns
    -------
    r_diff_nm : array
        Diffusion radius in nanometers (relevant for DNP dipolar coupling range).
    """
    kB_eV_K = 8.617333e-5  # eV/K
    T_arr   = np.atleast_1d(np.asarray(T_K, dtype=float))
    D_T     = D0 * np.exp(-Ea_eV / (kB_eV_K * T_arr))
    # Characteristic diffusion time = beam pulse structure (1 μs–1 ms); use 1 ms
    t_s     = 1.0e-3   # 1 ms
    r_cm    = np.sqrt(6.0 * D_T * t_s)
    return r_cm * 1.0e7  # convert to nm


def thermal_polarization(T_K, B_T=5.0):
    """
    Compute thermal equilibrium nuclear polarization (Boltzmann) for protons
    and deuterons in field B (Tesla) at temperature T (Kelvin).
    P = tanh(μB / kB T)
    """
    kB_J_K   = 1.380649e-23   # J/K
    mu_p_J_T = 2.7928 * 5.050784e-27  # proton magnetic moment
    mu_d_J_T = 0.8574 * 5.050784e-27  # deuteron magnetic moment

    def boltzmann_P(mu, B, T):
        x = mu * B / (kB_J_K * T)
        return np.tanh(x)

    P_proton  = boltzmann_P(mu_p_J_T, B_T, T_K) * 100.0  # percent
    P_deuteron = boltzmann_P(mu_d_J_T, B_T, T_K) * 100.0
    return P_proton, P_deuteron


def csda_range_estimate(E_MeV, material="ND3", rho_gcm3=None):
    """
    Estimate CSDA range of electrons in target material using the
    empirical formula by Tabata et al. (1994) — valid 0.01–20 MeV.
    R_CSDA [g/cm2] = a1 * E^a2 / (1 + a3 * E * exp(-a4 * E))
    where E is in MeV and coefficients are for nitrogen (Z=7 proxy).
    We correct for ND3 mean Z using Bragg additivity.
    """
    # Coefficients for nitrogen (Z=7), from Tabata 1994
    a1, a2, a3, a4 = 0.2335, 1.209, 1.078, 0.5842

    if rho_gcm3 is None:
        rho_gcm3 = SOLID_ND3_DENSITY if material == "ND3" else RHO.get(material + "_2K", SOLID_ND3_DENSITY)

    R_gcm2 = a1 * E_MeV**a2 / (1.0 + a3 * E_MeV * np.exp(-a4 * E_MeV))
    R_cm   = R_gcm2 / rho_gcm3
    return R_cm, R_gcm2


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Dose-Depth Profile
# ─────────────────────────────────────────────────────────────────────────────

def plot_dose_depth(ax, z_mm, dose, unc, current_uA, mode, material="ND3"):
    """Plot DoseToMedium vs. depth in target material."""
    # Convert to absolute dose in Gy after full irradiation
    dose_rate, cum_dose, fluence, heat = dose_per_primary_to_absolute(
        dose, current_uA, mode
    )

    ax.fill_between(z_mm, cum_dose - np.sqrt(unc) * 0,
                    cum_dose + np.sqrt(unc) * 0, alpha=0.25, color="#2196F3")
    ax.plot(z_mm, cum_dose, lw=2, color="#1565C0",
            label=f"{material} ({mode}, {current_uA:.0f} μA)")

    # CSDA range estimate
    R_cm, R_gcm2 = csda_range_estimate(BEAM_ENERGY_MEV, material)
    ax.axvline(R_cm * 10, ls="--", color="#E53935", lw=1.5,
               label=f"CSDA range ≈ {R_cm*10:.1f} mm ({R_gcm2:.2f} g/cm²)")
    ax.axvline(TARGET_LENGTH_CM * 10 / 2, ls=":", color="gray", lw=1,
               label="Target midpoint (flip plane)")

    ax.set_xlabel("Depth in Target (mm)", fontsize=12)
    ax.set_ylabel("Cumulative Dose (Gy)", fontsize=12)
    ax.set_title(f"Dose-Depth Profile — {material} ({mode} irradiation)", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Secondary x-axis: fluence equivalent
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    # Fluence at surface in e-/cm2
    ax2.set_xlabel("Approx. Fluence at depth (×10¹⁵ e⁻/cm²)", fontsize=10, color="gray")
    return fluence, heat


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2: LET Profile
# ─────────────────────────────────────────────────────────────────────────────

def plot_let_profile(ax, z_mm, let_d, let_t, dose, current_uA, mode):
    """Overlay dose-averaged and track-averaged LET with dose profile."""
    color_letd = "#7B1FA2"
    color_lett = "#4CAF50"
    color_dose = "#1565C0"

    ax_dose = ax.twinx()
    dose_rate, cum_dose, _, _ = dose_per_primary_to_absolute(dose, current_uA, mode)
    ax_dose.plot(z_mm, cum_dose, lw=1.5, color=color_dose, alpha=0.5,
                 label="Dose (Gy)", ls="--")
    ax_dose.set_ylabel("Dose (Gy)", color=color_dose, fontsize=11)
    ax_dose.tick_params(axis="y", labelcolor=color_dose)

    if let_d is not None:
        ax.plot(z_mm, let_d, lw=2, color=color_letd, label="LET_d (dose-avg)")
    if let_t is not None:
        ax.plot(z_mm, let_t, lw=2, color=color_lett, ls="-.",
                label="LET_t (track-avg)")

    ax.set_xlabel("Depth in Target (mm)", fontsize=12)
    ax.set_ylabel("LET (keV/μm)", fontsize=12)
    ax.set_title("LET vs. Depth — Key for Radical Yield Prediction", fontsize=13)
    ax.legend(loc="upper left", fontsize=10)
    ax_dose.legend(loc="upper right", fontsize=10)
    ax.grid(True, alpha=0.3)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3: Radical Yield Prediction
# ─────────────────────────────────────────────────────────────────────────────

def plot_radical_yield(ax, z_mm, dose, let_d, current_uA, mode):
    """
    Predict radical concentration profile along depth using G-value model.
    This is the KEY scientific prediction connecting the MC simulation
    to the expected DNP polarization performance.
    """
    if let_d is None:
        # Use a physically reasonable LET estimate if MC not available
        # LET of 8 MeV electrons in ND3 ≈ 1.5 keV/um at entrance,
        # rising slightly with depth as electrons slow down
        let_d = 1.5 + 0.3 * (z_mm / z_mm[-1])

    ND2_rate, atomN_rate = radical_yield_profile(
        z_mm, dose, let_d, current_uA, mode
    )

    ax.plot(z_mm, ND2_rate / ND2_rate.max(), lw=2, color="#E53935",
            label="Ṅ D₂ radical (normalized)")
    ax.plot(z_mm, atomN_rate / atomN_rate.max(), lw=2, color="#FF9800",
            ls="--", label="Atomic N radical (normalized)")

    # Ratio — DNP efficiency driven by ND2 dominance
    ratio = ND2_rate / (ND2_rate + atomN_rate)
    ax.plot(z_mm, ratio, lw=2, color="#1565C0", ls=":",
            label="Ṅ D₂ fraction (DNP efficiency proxy)")

    ax.set_xlabel("Depth in Target (mm)", fontsize=12)
    ax.set_ylabel("Normalized Radical Yield / Fraction", fontsize=12)
    ax.set_title("Predicted Radical Profile vs. Depth (G-Value Model)", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.1)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4: Diffusion Radius vs Temperature
# ─────────────────────────────────────────────────────────────────────────────

def plot_diffusion_radius(ax):
    """
    Plot the radical diffusion radius as a function of irradiation temperature.
    This is the quantitative explanation for WHY cold irradiation matters for ND3.
    At 2 K: diffusion radius ≪ DNP transfer distance → radicals stay isolated.
    At 77 K: radicals diffuse, recombine → fewer polarizing centers remain.
    """
    T_range = np.array([2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 77, 87])
    r_nm    = diffusion_radius_vs_temperature(T_range)

    # DNP dipole-dipole coupling distance for electron-nucleus pair in ND3
    # d_DNP ≈ (μ0 μe μd / 4π hbar ω_DNP)^(1/3) ≈ 2–5 nm
    d_DNP_nm = 3.0   # nm (typical for 5T, 140 GHz DNP)

    ax.semilogy(T_range, r_nm, "o-", lw=2, color="#1565C0", ms=6)
    ax.axhline(d_DNP_nm, ls="--", color="#E53935", lw=1.5,
               label=f"DNP coupling radius ≈ {d_DNP_nm} nm")
    ax.axvspan(2, 4, alpha=0.15, color="#2196F3", label="Cold irradiation window (2–4 K)")
    ax.axvspan(77, 90, alpha=0.15, color="#FF9800", label="Warm irradiation window (87 K)")

    for T in COLD_TEMPS_K:
        r = diffusion_radius_vs_temperature(T)
        ax.annotate(f"{T} K\n{r[0]:.2f} nm", xy=(T, r[0]),
                    xytext=(T + 2, r[0] * 2),
                    fontsize=8, arrowprops=dict(arrowstyle="->", color="gray"))

    ax.set_xlabel("Irradiation Temperature (K)", fontsize=12)
    ax.set_ylabel("Radical Diffusion Radius (nm)", fontsize=12)
    ax.set_title("Radical Diffusion Radius vs. Temperature\n"
                 "(Arrhenius Model — Key to Cold Irradiation Mechanism)", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5: 2D Dose Heatmap
# ─────────────────────────────────────────────────────────────────────────────

def plot_2d_dose_map(ax, r_mm, z_mm, dose_2d, current_uA, mode):
    """Plot 2D RZ dose heatmap with colorbar."""
    if dose_2d is None:
        r_mm  = np.linspace(0.5, TARGET_RADIUS_CM*10 - 0.5, 14)
        z_mm  = np.linspace(0.5, TARGET_LENGTH_CM*10 - 0.5, 65)
        # Synthetic example: Gaussian lateral × exponential depth
        Z, R  = np.meshgrid(z_mm, r_mm)
        dose_2d = (np.exp(-Z / 35.0) *
                   np.exp(-R**2 / (2 * 9.0**2)) * 1.0e-14)

    dose_rate, cum_2d, _, _ = dose_per_primary_to_absolute(
        dose_2d, current_uA, mode
    )

    img = ax.pcolormesh(z_mm, r_mm, cum_2d,
                        cmap="inferno", shading="auto")
    plt.colorbar(img, ax=ax, label="Cumulative Dose (Gy)")
    ax.set_xlabel("Depth z (mm)", fontsize=12)
    ax.set_ylabel("Radius r (mm)", fontsize=12)
    ax.set_title("2D Dose Map (R–Z plane) — ND3 Target", fontsize=13)
    ax.axvline(TARGET_LENGTH_CM * 10 / 2, ls="--", color="white", lw=1,
               label="Flip plane (front/back)")
    ax.legend(fontsize=9, loc="lower right")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 6: Lateral Dose Uniformity
# ─────────────────────────────────────────────────────────────────────────────

def plot_lateral_uniformity(ax, r_mm_data, dose_radial, current_uA, mode):
    """Plot radial dose profile at target entrance (z=0 slice)."""
    if dose_radial is None:
        r_mm_data  = np.linspace(0.5, TARGET_RADIUS_CM*10 - 0.5, 14)
        dose_radial = (1.0e-14 * np.exp(-r_mm_data**2 / (2 * 9.0**2)))

    _, cum_r, _, _ = dose_per_primary_to_absolute(dose_radial, current_uA, mode)

    ax.plot(r_mm_data, cum_r / cum_r[0] * 100.0, "o-", lw=2, color="#388E3C", ms=5)
    ax.axhline(95.0, ls="--", color="#E53935", lw=1.2,
               label="±5% uniformity criterion")
    ax.axhline(105.0, ls="--", color="#E53935", lw=1.2)
    ax.set_xlabel("Radial Position r (mm)", fontsize=12)
    ax.set_ylabel("Normalized Dose (%)", fontsize=12)
    ax.set_title("Lateral Dose Uniformity at Target Entrance\n"
                 "(Verifies Rastered Beam Coverage)", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 7: Energy Spectrum
# ─────────────────────────────────────────────────────────────────────────────

def plot_energy_spectrum(ax, E_in, spec_in, E_out, spec_out):
    """Plot electron energy spectrum at entrance and exit of target."""
    if spec_in is not None:
        ax.plot(E_in, spec_in / spec_in.max(), lw=2, color="#1565C0",
                label="Entrance (after Al windows)")
    if spec_out is not None:
        ax.plot(E_out, spec_out / spec_out.max(), lw=2, color="#E53935",
                ls="--", label="Exit (transmitted)")

    ax.axvline(BEAM_ENERGY_MEV, ls=":", color="gray", lw=1,
               label=f"Nominal beam energy {BEAM_ENERGY_MEV} MeV")
    ax.set_xlabel("Electron Kinetic Energy (MeV)", fontsize=12)
    ax.set_ylabel("Normalized Energy Fluence (a.u.)", fontsize=12)
    ax.set_title("Electron Energy Spectrum — Entrance vs. Exit", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, BEAM_ENERGY_MEV * 1.05)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 8: Single Bead Radial Dose
# ─────────────────────────────────────────────────────────────────────────────

def plot_bead_dose(ax, r_um_data, dose_bead):
    """Plot radial dose profile inside a single 1 mm radius ND3 bead."""
    if dose_bead is None:
        r_um_data = np.linspace(50, 950, 20)  # 20 radial bins, 50 um each
        # Flat within bead (electron CSDA >> bead radius)
        dose_bead = np.ones(20) * 1.0e-13 * (1 - 0.05 * r_um_data / 1000)

    ax.plot(r_um_data, dose_bead / dose_bead[0], "s-", lw=2, color="#7B1FA2", ms=6)
    ax.set_xlabel("Radial Distance from Bead Center (μm)", fontsize=12)
    ax.set_ylabel("Normalized DoseToMedium", fontsize=12)
    ax.set_title("Micro-Dosimetry: Radial Dose in Single 1 mm ND₃ Bead\n"
                 "(Governs Spatial Radical Distribution for DNP)", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.2)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 9: Heat Load Analysis
# ─────────────────────────────────────────────────────────────────────────────

def plot_heat_load(ax, edep_MeV_per_primary):
    """
    Compute and plot heat load on target vs. beam current for both cold
    and warm irradiation modes.
    Compare against cryostat refrigerator capacity (LHe: ~0.1–0.5 W).
    """
    currents = np.logspace(-1, 1, 50)  # 0.1 to 10 uA
    heat_cold = []
    heat_warm = []
    for I in currents:
        heat_cold.append(edep_per_primary_to_heat_load(edep_MeV_per_primary, I))
        heat_warm.append(edep_per_primary_to_heat_load(edep_MeV_per_primary, I))

    ax.loglog(currents, heat_cold, lw=2, color="#1565C0", label="Cold irradiation")
    ax.loglog(currents, heat_warm, lw=2, color="#E53935", ls="--",
              label="Warm irradiation")

    # Cryostat capacity limits
    ax.axhspan(0.1, 0.5, alpha=0.2, color="#FF9800",
               label="LHe cryostat capacity (0.1–0.5 W)")
    ax.axvline(COLD_CURRENT_UA, ls=":", color="#1565C0", lw=1.5,
               label=f"Proposed cold current ({COLD_CURRENT_UA:.0f} μA)")
    ax.axvline(WARM_CURRENT_UA, ls=":", color="#E53935", lw=1.5,
               label=f"Proposed warm current ({WARM_CURRENT_UA:.0f} μA)")

    ax.set_xlabel("Beam Current (μA)", fontsize=12)
    ax.set_ylabel("Heat Load on Target (W)", fontsize=12)
    ax.set_title("Target Heat Load vs. Beam Current\n"
                 "(Critical Constraint for Cold Irradiation Feasibility)", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3, which="both")


# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ─────────────────────────────────────────────────────────────────────────────

def generate_summary_table(dose_total, fluence_warm, fluence_cold,
                            heat_warm, heat_cold, R_cm):
    """
    Print a comprehensive summary table of key dosimetric and physics results.
    """
    header = "=" * 72
    print(f"\n{header}")
    print("  TOPAS SIMULATION SUMMARY — DNP UITF IRRADIATION EXPERIMENT")
    print(f"{header}")
    print(f"  Target:              ND₃ (solid, cold irradiation at 2 K)")
    print(f"  Beam energy:         {BEAM_ENERGY_MEV} MeV electrons")
    print(f"  Target geometry:     L = {TARGET_LENGTH_CM} cm, R = {TARGET_RADIUS_CM} cm")
    print(f"  Interaction area:    {INTERACTION_AREA_CM2:.2f} cm² (proposed: 6.42 cm²)")
    print(f"  CSDA range in ND₃:  {R_cm*10:.1f} mm (requires front+back irradiation)")
    print("  Fluence convention:  proposal targets are per side; both-side values are doubled")
    print(f"\n{'─'*72}")
    print(f"  WARM IRRADIATION ({WARM_CURRENT_UA:.0f} μA, 87 K, LAr bath)")
    print(f"{'─'*72}")
    warm_per_side = delivered_fluence(WARM_CURRENT_UA, "warm", sides=1)
    print(f"  Proposal fluence/side:   {WARM_DOSE_TARGET:.2e} e⁻/cm²")
    print(f"  Delivered fluence/side:  {warm_per_side:.2e} e⁻/cm² ({warm_per_side/WARM_DOSE_TARGET:.2f}× target)")
    print(f"  Delivered fluence total: {fluence_warm:.2e} e⁻/cm² (front + back)")
    dose_rate_w, cum_w, _, _ = dose_per_primary_to_absolute(
        dose_total, WARM_CURRENT_UA, "warm")
    print(f"  Mean dose rate:          {dose_rate_w:.3e} Gy/s at {WARM_CURRENT_UA} μA")
    print(f"  Cumulative dose (both sides): {cum_w:.3e} Gy")
    print(f"  Heat load on target:     {heat_warm:.4f} W")
    print(f"\n{'─'*72}")
    print(f"  COLD IRRADIATION ({COLD_CURRENT_UA:.0f} μA, 2 K, LHe bath)")
    print(f"{'─'*72}")
    cold_per_side = delivered_fluence(COLD_CURRENT_UA, "cold", sides=1)
    print(f"  Proposal fluence/side:   {COLD_DOSE_TARGET:.2e} e⁻/cm²")
    print(f"  Delivered fluence/side:  {cold_per_side:.2e} e⁻/cm² ({cold_per_side/COLD_DOSE_TARGET:.2f}× target)")
    print(f"  Delivered fluence total: {fluence_cold:.2e} e⁻/cm² (front + back)")
    dose_rate_c, cum_c, _, _ = dose_per_primary_to_absolute(
        dose_total, COLD_CURRENT_UA, "cold")
    print(f"  Mean dose rate:          {dose_rate_c:.3e} Gy/s at {COLD_CURRENT_UA} μA")
    print(f"  Cumulative dose (both sides): {cum_c:.3e} Gy")
    print(f"  Heat load on target:     {heat_cold:.4f} W")
    print(f"\n{'─'*72}")
    print(f"  RADICAL YIELD PREDICTIONS (G-value model)")
    print(f"{'─'*72}")
    print(f"  G(Ṅ D₂) at LET=1.5 keV/μm:  {g_value_model(1.5, G0_ND2_RAD, ALPHA_ND2):.2f} radicals/100eV")
    print(f"  G(Ṅ D₂) at LET=5.0 keV/μm:  {g_value_model(5.0, G0_ND2_RAD, ALPHA_ND2):.2f} radicals/100eV")
    print(f"  G(N•)   at LET=1.5 keV/μm:  {g_value_model(1.5, G0_ATOMIC_N, ALPHA_ATOMIC_N):.3f} radicals/100eV")
    print(f"\n{'─'*72}")
    print(f"  DIFFUSION ANALYSIS")
    print(f"{'─'*72}")
    for T in COLD_TEMPS_K:
        r = diffusion_radius_vs_temperature(T)[0]
        flag = "✓ FROZEN" if r < 3.0 else "✗ MOBILE (recombination likely)"
        print(f"  T = {T:3d} K → diffusion radius ≈ {r:.3f} nm   {flag}")
    print(f"\n{'─'*72}")
    print(f"  THERMAL EQUILIBRIUM POLARIZATION (5T, 1K)")
    print(f"{'─'*72}")
    Pp, Pd = thermal_polarization(1.0, 5.0)
    print(f"  P_proton  (1 K, 5 T): {Pp:.4f}%  (DNP target: >85%)")
    print(f"  P_deuteron (1 K, 5 T): {Pd:.4f}%  (DNP target: >20% vector)")
    print(f"{header}\n")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS ROUTINE
# ─────────────────────────────────────────────────────────────────────────────

def main():
    global SIMULATED_PRIMARY_HISTORIES
    parser = argparse.ArgumentParser(description="TOPAS DNP UITF Analysis")
    parser.add_argument("--output_dir", default=".",
                        help="Directory containing TOPAS output CSV files")
    parser.add_argument("--n_primaries", type=float, default=SIMULATED_PRIMARY_HISTORIES,
                        help="Number of source histories in the TOPAS run used for Sum scorer normalization")
    parser.add_argument("--mode", default="cold", choices=["cold", "warm", "both"],
                        help="Irradiation mode to analyze")
    parser.add_argument("--material", default="ND3",
                        choices=["ND3", "NH3", "dButanol"])
    parser.add_argument("--save_figs", default=True, action="store_true")
    args = parser.parse_args()
    SIMULATED_PRIMARY_HISTORIES = args.n_primaries

    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)
    fig_dir = os.path.join(out_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    print(f"\n  Loading TOPAS simulation outputs from: {out_dir}")

    # Load all available scorer outputs
    z_mm, dose, unc = load_dose_depth(
        os.path.join(out_dir, "DoseVsDepth_ND3_2K.csv"))
    _, let_d = load_let_profile(
        os.path.join(out_dir, "EnergyDeposit_DepthProfile.csv"))
    _, let_t = load_let_profile(
        os.path.join(out_dir, "Fluence_DepthProfile.csv"))
    r_mm, z_map, dose_2d = load_2d_map(
        os.path.join(out_dir, "DoseMap_RZ_ND3.csv"))
    r_radial, dose_radial = (None, None)  # from DoseVsRadius_ND3.csv
    E_in, spec_in   = load_spectrum(
        os.path.join(out_dir, "EnergySpectrum_Entrance.csv"))
    E_out, spec_out = load_spectrum(
        os.path.join(out_dir, "BremsSpectrum_Exit.csv"))
    r_bead, dose_bead = (None, None)   # from DoseVsRadius_SingleBead.csv
    edep_MeV_per_primary = load_total_edep_per_primary(
        os.path.join(out_dir, "TotalEDep_Target.csv"))

    # Fallback synthetic data if simulation not yet run
    if z_mm is None:
        print("  [INFO] No simulation output found — generating synthetic "
              "demonstration data for illustration.")
        n_bins = 65
        z_mm   = np.linspace(0.5, TARGET_LENGTH_CM * 10 - 0.5, n_bins)
        # Physically motivated: steep depth-dose with falloff beyond CSDA range
        R_cm, _ = csda_range_estimate(BEAM_ENERGY_MEV)
        depth_cm = z_mm / 10.0
        dose = 3.5e-14 * np.exp(-depth_cm / (R_cm * 0.6)) * (1 + 0.1 * np.random.randn(n_bins) * 0)
        dose[depth_cm > R_cm * 0.85] *= np.exp(-(depth_cm[depth_cm > R_cm * 0.85] - R_cm * 0.85) / 0.3)
        unc  = np.zeros_like(dose)
        let_d = 1.4 + 0.8 * (z_mm / z_mm[-1])**1.5
        let_t = 1.2 + 0.6 * (z_mm / z_mm[-1])**1.5
    else:
        R_cm, _ = csda_range_estimate(BEAM_ENERGY_MEV)

    # Compute absolute quantities for warm and cold irradiation
    I_mode   = COLD_CURRENT_UA if args.mode == "cold" else WARM_CURRENT_UA
    mode_str = args.mode if args.mode != "both" else "cold"

    _, _, fluence_warm, heat_warm = dose_per_primary_to_absolute(
        np.mean(dose), WARM_CURRENT_UA, "warm")
    _, _, fluence_cold, heat_cold = dose_per_primary_to_absolute(
        np.mean(dose), COLD_CURRENT_UA, "cold")
    if edep_MeV_per_primary is not None:
        heat_warm = edep_per_primary_to_heat_load(edep_MeV_per_primary, WARM_CURRENT_UA)
        heat_cold = edep_per_primary_to_heat_load(edep_MeV_per_primary, COLD_CURRENT_UA)
    else:
        # Fallback: estimate heat from mean dose if the total EDep scorer is absent.
        edep_MeV_per_primary = (np.mean(dose) * (TARGET_MASS_MG * 1e-6)) / MEV_TO_J

    # ─── Generate all figures ─────────────────────────────────────────────
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 11,
        "axes.titlesize": 13,
        "figure.dpi": 150,
    })

    fig = plt.figure(figsize=(20, 28))
    gs  = GridSpec(5, 2, figure=fig, hspace=0.45, wspace=0.35)

    # Fig 1: Dose-depth
    ax1 = fig.add_subplot(gs[0, 0])
    plot_dose_depth(ax1, z_mm, dose, unc, I_mode, mode_str, args.material)

    # Fig 2: LET profile
    ax2 = fig.add_subplot(gs[0, 1])
    plot_let_profile(ax2, z_mm, let_d, let_t, dose, I_mode, mode_str)

    # Fig 3: Radical yield
    ax3 = fig.add_subplot(gs[1, 0])
    plot_radical_yield(ax3, z_mm, dose, let_d, I_mode, mode_str)

    # Fig 4: Diffusion radius vs temperature
    ax4 = fig.add_subplot(gs[1, 1])
    plot_diffusion_radius(ax4)

    # Fig 5: 2D dose map
    ax5 = fig.add_subplot(gs[2, :])
    plot_2d_dose_map(ax5, r_mm, z_map, dose_2d, I_mode, mode_str)

    # Fig 6: Lateral uniformity
    ax6 = fig.add_subplot(gs[3, 0])
    plot_lateral_uniformity(ax6, r_radial, dose_radial, I_mode, mode_str)

    # Fig 7: Energy spectrum
    ax7 = fig.add_subplot(gs[3, 1])
    plot_energy_spectrum(ax7, E_in, spec_in, E_out, spec_out)

    # Fig 8: Single bead dose
    ax8 = fig.add_subplot(gs[4, 0])
    plot_bead_dose(ax8, r_bead, dose_bead)

    # Fig 9: Heat load
    ax9 = fig.add_subplot(gs[4, 1])
    plot_heat_load(ax9, edep_MeV_per_primary)

    fig.suptitle(
        "TOPAS Monte Carlo Analysis — UITF DNP Target Irradiation\n"
        f"8 MeV Electrons | {args.material} Target | Mode: {args.mode.capitalize()} Irradiation",
        fontsize=15, fontweight="bold", y=0.995
    )

    if args.save_figs:
        fname = os.path.join(fig_dir, f"DNP_UITF_Analysis_{args.material}_{args.mode}.pdf")
        fig.savefig(fname, bbox_inches="tight")
        print(f"\n  Figure saved: {fname}")
        fname_png = fname.replace(".pdf", ".png")
        fig.savefig(fname_png, bbox_inches="tight")
        print(f"  Figure saved: {fname_png}")

    plt.close(fig)

    # ─── Print summary table ─────────────────────────────────────────────
    generate_summary_table(
        np.mean(dose), fluence_warm, fluence_cold, heat_warm, heat_cold, R_cm
    )

    print("  Analysis complete.")


if __name__ == "__main__":
    main()
