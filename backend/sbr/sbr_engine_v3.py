"""
SBR ENGINE v3 — Unified Backend Module
--------------------------------------
PART 1 of 3

This module implements:
 - Module 1: Geometry + Cycle Architecture
 - Module 2: Biological Process + Oxygen Demand
 - Module 3: Blower, Mixer, and Energy Modeling

This file is structured with:
   - Clear helper utilities
   - Explicit engineering documentation
   - Modular internal functions
   - A single public function: run_sbr_engine(inputs)

Comment Level: FULL ENGINEERING DETAIL
"""

# ============================================================
#  IMPORTS
# ============================================================

import math


# ============================================================
#  GLOBAL CONSTANTS
#  These constants are used across all modules.
#  Any GUI or user-layer can override these as needed.
# ============================================================

# Unit Conversions
GALLONS_PER_CUBIC_FOOT = 7.48052
MGD_TO_GPH = 1_000_000 / 24.0  # 1 MGD = 1,000,000 gal/day → gal/hr
POUNDS_PER_MG_PER_MG_L = 8.34  # lb/day = 8.34 * MGD * mg/L

# Oxygen Transfer Constants
# Fine-bubble diffusers typically achieve ~0.15 OTE/ft in clean water.
# Effective oxygen transfer per cubic foot under field conditions:
OXYGEN_PER_1000_CUBIC_FEET = 1.08  # lb O2 per 1000 ft³ of air (standard EPA/M&E)
OXYGEN_PER_CUBIC_FOOT = OXYGEN_PER_1000_CUBIC_FEET / 1000.0

# Blower Pressures (confirmed by user)
BLOWER_PRESSURE_PSI_TURBO = 6.5
BLOWER_PRESSURE_PSI_PD = 8.5

# Diffuser Depth (confirmed by user)
DIFFUSER_DEPTH_FT = 18.0

# Standard mechanical efficiencies
TURBO_EFFICIENCY = 0.80
PD_EFFICIENCY = 0.60

# Mixing target energy density (activated sludge mixing range: 0.3–1.0 W/m³)
MIXING_W_PER_M3 = 0.7  # A good engineering midpoint


# ============================================================
#  SAFE MATH UTILITIES
#  These prevent divide-by-zero and ensure stable outputs.
# ============================================================

def safe_div(numerator, denominator, default=0.0):
    """
    Safely divide numerator/denominator.
    If denominator is zero or extremely small, return `default`.

    Parameters:
        numerator (float)
        denominator (float)
        default (float)

    Returns:
        float
    """
    if denominator is None or abs(denominator) < 1e-12:
        return default
    return numerator / denominator


def to_mg(volume_gallons):
    """
    Convert volume in gallons → million gallons (MG).
    """
    return volume_gallons / 1_000_000.0


def lbs_per_day(mgd, concentration_mgL):
    """
    Compute mass loading:
        lb/day = 8.34 × MGD × mg/L
    """
    return 8.34 * mgd * concentration_mgL


def ft3_to_gallons(ft3):
    return ft3 * GALLONS_PER_CUBIC_FOOT


def gallons_to_ft3(gal):
    return gal / GALLONS_PER_CUBIC_FOOT


# ============================================================
#  ENGINEERING HELPER FUNCTIONS
# ============================================================

def compute_sae(depth_ft, sote_per_ft=0.15):
    """
    Compute Standard Aeration Efficiency (SAE)
    based on diffuser submergence.

    SAE = depth_ft × SOTE_per_ft

    Typical SOTE_per_ft:
       - 0.12–0.18 for fine bubble diffusers (clean water)
       - Field values are lower; we use 0.15 as a midpoint.

    Returns:
        SAE (unitless)
    """
    return depth_ft * sote_per_ft


def compute_scfm_required(sor_lb_hr, depth_ft):
    """
    Field oxygen transfer rate for fine-bubble diffusers
    at 16–20 ft depth:
        1 SCFM ≈ 0.18 lb O2/hr  (field conditions)

    SCFM = SOR / 0.18
    """
    O2_transfer_per_SCFM = 0.18  # lb O2/hr per SCFM
    return safe_div(sor_lb_hr, O2_transfer_per_SCFM, default=0.0)




def compute_blower_hp(scfm, pressure_psi, efficiency):
    """
    Compute blower brake horsepower (BHP) using:

        HP = (SCFM × ΔP × 1.0) / (229 × η)

    where:
        ΔP = blower discharge pressure (psi)
        η = blower efficiency (fraction)
        229 = conversion factor incorporating
              1 HP = 33,000 ft·lbf/min
              1 psi = 144 lb/ft²
              and air density / standard conditions.

    Returns:
        HP (float)
    """
    return safe_div(scfm * pressure_psi, (229.0 * efficiency), default=0.0)


def compute_mixing_hp(volume_ft3, w_per_m3=MIXING_W_PER_M3):
    """
    Compute required mixing power using:

      Energy density target = 0.3–1.0 W/m³
      Chosen default = 0.7 W/m³ (extended aeration / SBR)

    Conversion:
      1 ft³ = 0.0283168 m³
      Convert W → HP:
        1 HP = 746 W

    Returns:
        hp_mixing_required (float)
    """
    m3 = volume_ft3 * 0.0283168
    total_watts = m3 * w_per_m3
    hp = total_watts / 746.0
    return hp


def annual_energy_cost(hp, hours_per_year, electricity_cost_per_kwh=0.12):
    """
    Compute annual energy cost:

      kW = HP × 0.746
      kWh/year = kW × hours × load_factor
      cost = kWh × $/kWh

    This version assumes user provides actual hours/year.
    """
    kw = hp * 0.746
    annual_kwh = kw * hours_per_year
    annual_cost = annual_kwh * electricity_cost_per_kwh
    return {
        "annual_kwh": annual_kwh,
        "annual_cost": annual_cost
    }

# ============================================================
#  MODULE 1 — GEOMETRY (PART 2A)
#  Computes basin volume, dimensions, decant depth,
#  fill depth, effective liquid depth, and validates freeboard.
#
#  All downstream biological and blower calculations depend on
#  these values.
# ============================================================

def _compute_geometry(inputs):
    """
    Compute geometric characteristics of the SBR basin based on
    user-defined dimensions and cycle sequencing.

    Inputs expected:
        basin_length_ft
        basin_width_ft
        basin_depth_ft
        freeboard_ft
        decant_fraction_of_depth

    Returns:
        geometry (dict)
    """

    # Extract inputs with safe defaults
    L = inputs.get("basin_length_ft", 100.0)
    W = inputs.get("basin_width_ft", 60.0)
    depth = inputs.get("basin_depth_ft", 18.0)
    freeboard = inputs.get("freeboard_ft", 2.0)
    decant_frac = inputs.get("decant_fraction_of_depth", 0.3)

    # ---------------------------------------------------------------------
    # 1. Geometric Basin Volume (ft³)
    # ---------------------------------------------------------------------
    # Total geometric volume = L × W × full depth
    volume_ft3 = L * W * depth

    # Convert to gallons
    volume_gal = ft3_to_gallons(volume_ft3)

    # ---------------------------------------------------------------------
    # 2. Effective Operating Depth
    #    SBRs operate from maximum water level down to decant depth.
    #    Freeboard ensures containment during aeration & surge.
    # ---------------------------------------------------------------------
    max_operating_depth_ft = depth - freeboard
    if max_operating_depth_ft <= 0:
        max_operating_depth_ft = max(0.1, depth * 0.5)

    # ---------------------------------------------------------------------
    # 3. Decant Depth (ft)
    #    Defined as a fraction of the total water depth.
    #    Typical decant fraction: 0.25–0.35 for extended aeration SBR.
    # ---------------------------------------------------------------------
    decant_depth_ft = depth * decant_frac

    # ---------------------------------------------------------------------
    # 4. Fill Depth (ft)
    #    Fill volume depends on number of cycles and plant flow.
    #    This is finalized in Module 1B since we need cycle timing.
    #
    #    Here, we initialize structure and compute geometry only.
    # ---------------------------------------------------------------------
    fill_depth_ft = None  # computed in module 1B
    fill_volume_gal = None

    # ---------------------------------------------------------------------
    # 5. Assemble geometry block
    # ---------------------------------------------------------------------
    geometry = {
        "inputs": {
            "basin_length_ft": L,
            "basin_width_ft": W,
            "basin_depth_ft": depth,
            "freeboard_ft": freeboard,
            "decant_fraction_of_depth": decant_frac,
        },
        "computed": {
            "volume_ft3_total": volume_ft3,
            "volume_gal_total": volume_gal,
            "max_operating_depth_ft": max_operating_depth_ft,
            "decant_depth_ft": decant_depth_ft,
            "fill_depth_ft": fill_depth_ft,  # to be filled later
            "fill_volume_gal": fill_volume_gal,  # filled in part 2B
        },
        "warnings": []
    }

    # ---------------------------------------------------------------------
    # 6. Warning Checks
    # ---------------------------------------------------------------------

    # Freeboard under 1.5–2 ft can cause issues with wave action
    if freeboard < 1.5:
        geometry["warnings"].append(
            "Freeboard less than 1.5 ft — may risk overtopping during aeration."
        )

    # Decant fraction too high can destabilize solids blanket
    if decant_frac > 0.40:
        geometry["warnings"].append(
            "Decant fraction exceeds 40% of depth — risk of solids carryover."
        )

    # Decant fraction too low may not achieve needed withdrawal volume
    if decant_frac < 0.20:
        geometry["warnings"].append(
            "Decant fraction less than 20% — may not provide sufficient decant volume."
        )

    return geometry

# ============================================================
#  MODULE 1 — CYCLE TIMING & VOLUMES (PART 2B)
#  Computes:
#      - Fill volume per cycle
#      - Fill depth (ft)
#      - Effective cycle time
#      - Detention time
#      - Validity checks on sequencing
# ============================================================

# ============================================================
#  MODULE 1 — CYCLE TIMING & VOLUMES (UPDATED / FIXED)
# ============================================================

def _compute_cycle_timing(inputs, geometry):
    """
    Compute SBR cycle timing, fill volumes, and operational volumes.
    Includes corrected detention time, corrected decant depth,
    and improved mass balance logic.
    """

    # ---------------------------------------------------------
    # Extract process inputs
    # ---------------------------------------------------------
    flow_mgd = inputs.get("flow_mgd", 0.5)
    n_basins = inputs.get("n_basins", 2)
    cycles_per_day = inputs.get("cycles_per_day", 3.0)

    fill_hr = inputs.get("fill_hr", 1.0)
    react_hr = inputs.get("react_hr", 4.0)
    settle_hr = inputs.get("settle_hr", 1.0)
    decant_hr = inputs.get("decant_hr", 1.0)
    idle_hr = inputs.get("idle_hr", 0.0)

    # ---------------------------------------------------------
    # Retrieve geometry metrics (from Part 2A)
    # ---------------------------------------------------------
    volume_gal_total = geometry["computed"]["volume_gal_total"]
    L = geometry["inputs"]["basin_length_ft"]
    W = geometry["inputs"]["basin_width_ft"]
    max_depth_ft = geometry["computed"]["max_operating_depth_ft"]

    # Basin surface area for depth conversions
    basin_surface_area_ft2 = L * W

    # ---------------------------------------------------------
    # 1. Compute flow per basin
    # ---------------------------------------------------------
    flow_mgd_per_basin = safe_div(flow_mgd, n_basins, default=0.0)
    flow_gpd_basin = flow_mgd_per_basin * 1_000_000.0

    # Fill volume per cycle
    fill_volume_gal = safe_div(flow_gpd_basin, cycles_per_day, default=0.0)

    # Convert to fill depth
    fill_depth_ft = safe_div(
        (fill_volume_gal / GALLONS_PER_CUBIC_FOOT),
        basin_surface_area_ft2,
        default=0.0
    )

    warnings = []

    # Check fill depth
    if fill_depth_ft > max_depth_ft:
        warnings.append(
            f"Fill depth {fill_depth_ft:.2f} ft exceeds max depth {max_depth_ft:.2f} ft. "
            "Increase basin volume or reduce cycles/day."
        )

    # ---------------------------------------------------------
    # 3. Total cycle time
    # ---------------------------------------------------------
    cycle_time_hr = fill_hr + react_hr + settle_hr + decant_hr + idle_hr
    theoretical_cycles_per_day = safe_div(24.0, cycle_time_hr, default=0.0)

    # ---------------------------------------------------------
    # 4. Detention Time (FIXED)
    #
    # Detention (hr) = Basin Volume (gal) / Flow (gal/hr)
    # ---------------------------------------------------------
    flow_gal_per_hr = flow_gpd_basin / 24.0
    detention_hr = safe_div(volume_gal_total, flow_gal_per_hr, default=0.0)

    if detention_hr < 12.0:
        warnings.append(
            f"Detention time only {detention_hr:.1f} hr — extended aeration SBRs "
            "typically require 18–30 hr."
        )

    # ---------------------------------------------------------
    # 5. Decant Volume (FIXED DEPTH)
    #
    # Decant depth should be based on MAX OPERATING DEPTH,
    # not full basin depth.
    # ---------------------------------------------------------
    decant_frac = inputs.get("decant_fraction_of_depth", 0.30)
    decant_depth_ft = max_depth_ft * decant_frac  # ✔ FIXED

    decant_volume_ft3 = decant_depth_ft * basin_surface_area_ft2
    decant_volume_gal = ft3_to_gallons(decant_volume_ft3)

    # Mass balance checks
    if decant_volume_gal < 0.9 * fill_volume_gal:
        warnings.append(
            "Decant volume is significantly smaller than fill volume — "
            "water level will rise each cycle."
        )

    if decant_volume_gal > 1.2 * fill_volume_gal:
        warnings.append(
            "Decant volume exceeds fill volume by >20% — risk of pulling down too far "
            "and disturbing sludge blanket."
        )

    # ---------------------------------------------------------
    # 6. Store computed values
    # ---------------------------------------------------------
    geometry["computed"]["fill_depth_ft"] = fill_depth_ft
    geometry["computed"]["fill_volume_gal"] = fill_volume_gal
    geometry["computed"]["decant_volume_gal"] = decant_volume_gal
    geometry["computed"]["decant_depth_ft"] = decant_depth_ft
    geometry["computed"]["cycle_time_hr"] = cycle_time_hr
    geometry["computed"]["theoretical_cycles_per_day"] = theoretical_cycles_per_day
    geometry["computed"]["detention_time_hr"] = detention_hr

    # Add warnings
    geometry["warnings"].extend(warnings)

    # ---------------------------------------------------------
    # Return cycle block for output
    # ---------------------------------------------------------
    cycle = {
        "inputs": {
            "flow_mgd": flow_mgd,
            "n_basins": n_basins,
            "cycles_per_day": cycles_per_day,
            "fill_hr": fill_hr,
            "react_hr": react_hr,
            "settle_hr": settle_hr,
            "decant_hr": decant_hr,
            "idle_hr": idle_hr,
        },
        "computed": {
            "flow_mgd_per_basin": flow_mgd_per_basin,
            "fill_volume_gal": fill_volume_gal,
            "fill_depth_ft": fill_depth_ft,
            "decant_volume_gal": decant_volume_gal,
            "cycle_time_hr": cycle_time_hr,
            "theoretical_cycles_per_day": theoretical_cycles_per_day,
            "detention_time_hr": detention_hr,
        },
        "warnings": warnings,
    }

    return cycle

# ============================================================
#  MODULE 2 — BIOLOGICAL PROCESS + SRT + OXYGEN + EFFLUENT
# ============================================================

def _compute_biological_process(inputs, geometry, cycle):
    """
    MODULE 2: Biological Loads, SRT, Oxygen Demand, Effluent Quality

    This module computes:
      • BOD, TSS, NH3 loadings (lb/day)
      • MLSS/MLVSS (from inputs)
      • SRT (Solids Retention Time)
      • F/M ratio
      • Oxygen demand:
            - Carbonaceous BOD (cBOD)
            - Nitrification oxygen
            - Endogenous respiration
      • AOR, SOR
      • Expected effluent BOD & NH3 (steady-state)

    References:
      • Metcalf & Eddy 5th Ed., Ch. 7 (Activated Sludge)
      • Ten State Standards §93
      • EPA 625/1-86-022 for nitrification SRT curves
    """

    # ---------------------------------------------------------
    # Extract influent characteristics
    # ---------------------------------------------------------
    bod_mgL = inputs.get("bod_mgL", 250.0)
    tss_mgL = inputs.get("tss_mgL", 250.0)
    nh3_mgL = inputs.get("nh3_mgL", 30.0)

    flow_mgd = inputs.get("flow_mgd", 0.5)
    n_basins = inputs.get("n_basins", 2)

    # Per-basin flow
    flow_mgd_per_basin = safe_div(flow_mgd, n_basins, 0.0)

    # ---------------------------------------------------------
    # Mass loading calculations
    # ---------------------------------------------------------
    bod_lb_day = lbs_per_day(flow_mgd, bod_mgL)
    tss_lb_day = lbs_per_day(flow_mgd, tss_mgL)
    nh3_lb_day = lbs_per_day(flow_mgd, nh3_mgL)

    bod_lb_day_per_basin = safe_div(bod_lb_day, n_basins, 0.0)
    nh3_lb_day_per_basin = safe_div(nh3_lb_day, n_basins, 0.0)

    # ---------------------------------------------------------
    # Extract MLSS & yield coefficients
    # ---------------------------------------------------------
    mlss_mgL = inputs.get("mlss_mgL", 3000)
    solids_yield_coeff = inputs.get("solids_yield_coeff", 0.65)
    srt_target_days = inputs.get("srt_target_days", 15)

    # Basin volume (total from geometry)
    volume_gal = geometry["computed"]["volume_gal_total"]
    volume_mg = to_mg(volume_gal)

    # MLSS mass in basin
    #   lb = 8.34 × MG × mg/L
    mlss_mass_lb = 8.34 * volume_mg * mlss_mgL

    # ---------------------------------------------------------
    # SRT Calculation
    #
    # SRT = Mass_in_system / (WAS_mass_per_day)
    #
    # Estimate WAS production using:
    #   WAS = Y × BOD_load (per day)
    #
    # This is a simplification but works well for extended aeration.
    # ---------------------------------------------------------
    # Extended aeration SBR waste sludge estimate
    # Typical: 0.8 lb TSS wasted per lb BOD removed
    was_lb_day = 0.8 * bod_lb_day_per_basin

    srt_days = safe_div(mlss_mass_lb, was_lb_day, default=999)

    # SRT warning checks
    srt_warnings = []
    if srt_days < 8:
        srt_warnings.append(
            f"SRT = {srt_days:.1f} days — too low for nitrification stability (<10 days)."
        )
    if srt_days > 30:
        srt_warnings.append(
            f"SRT = {srt_days:.1f} days — excessively long; may cause old sludge and poor settling."
        )

    # ---------------------------------------------------------
    # F/M Ratio
    #
    # F/M = (BOD load) / (MLVSS mass in basin)
    #
    # For extended aeration:
    #   F/M = 0.05–0.12 (typical)
    #
    # Assume MLVSS ≈ 0.75 × MLSS
    # ---------------------------------------------------------
    mlvss_mgL = mlss_mgL * 0.75
    mlvss_mass_lb = 8.34 * volume_mg * mlvss_mgL

    fm_ratio = safe_div(bod_lb_day_per_basin, mlvss_mass_lb, 0.0)

    if fm_ratio > 0.25:
        srt_warnings.append(
            f"F/M = {fm_ratio:.3f} — high for extended aeration; may limit nitrification."
        )

    # ---------------------------------------------------------
    # OXYGEN DEMAND
    #
    # 1. Carbonaceous BOD removal:
    #       1 lb BOD removed ≈ 1.1 lb O2
    #
    # 2. Nitrification oxygen demand:
    #       1 lb NH3-N oxidized ≈ 4.57 lb O2 (M&E)
    #
    # 3. Endogenous respiration:
    #       0.06–0.15 lb O2/lb MLSS-day
    #
    # ---------------------------------------------------------
    cBOD_o2_lb_day = 1.1 * bod_lb_day_per_basin
    nitrification_o2_lb_day = 4.57 * nh3_lb_day_per_basin

    # endogenous O2 consumption
    endogenous_o2_lb_day = 0.06 * mlss_mass_lb

    # TOTAL DAILY O2 REQUIREMENT
    total_o2_lb_day = cBOD_o2_lb_day + nitrification_o2_lb_day + endogenous_o2_lb_day

    # AOR (Actual Oxygen Requirement) in lb/hr
    aor_lb_hr = total_o2_lb_day / 24.0

    # SOR (Standard Oxygen Requirement)
    # Correction factors:
    alpha = 0.85  # wastewater coefficient
    beta = 0.95   # diffuser fouling factor
    temperature_correction = 1.024 ** (inputs.get("temperature_c", 20) - 20)

    sor_lb_hr = safe_div(aor_lb_hr, (alpha * beta * temperature_correction), default=0.0)

    # ---------------------------------------------------------
    # NITRIFICATION CHECK — SUMMER
    #
    # Required SRT for nitrification:
    #   At 20°C → SRT_crit ≈ 1 day (M&E)
    #
    # Extended aeration SRT usually > 10 days → full nitrification.
    # ---------------------------------------------------------
    srt_crit_summer = 1.0  # days
    nitrification_fraction = min(1.0, srt_days / srt_crit_summer)

    effluent_nh3_mgL = nh3_mgL * (1 - nitrification_fraction)
    if effluent_nh3_mgL < 0.5:
        effluent_nh3_mgL = 0.0  # fully nitrified

    # ---------------------------------------------------------
    # NITRIFICATION — WINTER CHECK
    #
    # At 10°C → SRT_crit ≈ 7–10 days
    # Use 7.3 days as typical.
    #
    # ---------------------------------------------------------
    srt_crit_winter = 7.33
    nitrification_fraction_winter = min(1.0, srt_days / srt_crit_winter)

    effluent_nh3_mgL_winter = nh3_mgL * (1 - nitrification_fraction_winter)
    if effluent_nh3_mgL_winter < 0.5:
        effluent_nh3_mgL_winter = 0.0

    # ---------------------------------------------------------
    # Effluent BOD Estimate
    #
    # Extended Aeration:
    #   BOD effluent = 2–5 mg/L typical
    #
    # Scale by F/M and SRT:
    # ---------------------------------------------------------
    effluent_bod_mgL = max(2.0, 2.0 + (fm_ratio * 10))

    # ---------------------------------------------------------
    # Build result block
    # ---------------------------------------------------------
    biology = {
        "inputs": {
            "bod_mgL": bod_mgL,
            "tss_mgL": tss_mgL,
            "nh3_mgL": nh3_mgL,
            "mlss_mgL": mlss_mgL,
            "solids_yield_coeff": solids_yield_coeff,
            "srt_target_days": srt_target_days,
        },
        "loads": {
            "bod_lb_day": bod_lb_day,
            "tss_lb_day": tss_lb_day,
            "nh3_lb_day": nh3_lb_day,
            "bod_lb_day_per_basin": bod_lb_day_per_basin,
            "nh3_lb_day_per_basin": nh3_lb_day_per_basin,
        },
        "srt": {
            "mlss_mass_lb": mlss_mass_lb,
            "was_lb_day": was_lb_day,
            "srt_days": srt_days,
            "srt_warnings": srt_warnings,
        },
        "fm_ratio": fm_ratio,
        "oxygen": {
            "cBOD_o2_lb_day": cBOD_o2_lb_day,
            "nitrification_o2_lb_day": nitrification_o2_lb_day,
            "endogenous_o2_lb_day": endogenous_o2_lb_day,
            "total_o2_lb_day": total_o2_lb_day,
            "aor_lb_hr": aor_lb_hr,
            "sor_lb_hr": sor_lb_hr,
        },
        "effluent": {
            "bod_mgL": effluent_bod_mgL,
            "nh3_summer_mgL": effluent_nh3_mgL,
            "nh3_winter_mgL": effluent_nh3_mgL_winter,
        },
        "warnings": srt_warnings,
    }

    return biology
# ============================================================
#  MODULE 3 — AIRFLOW, BLOWERS, MIXERS, ENERGY MODEL
# ============================================================

def _compute_airflow_and_blowers(biology, inputs, geometry):
    """
    Computes:
      - Required SCFM based on corrected field oxygen transfer
      - Turbo blower HP
      - PD blower HP
      - Unit counts and firm capacity
    """

    sor_lb_hr = biology["oxygen"]["sor_lb_hr"]

    # ============================================================
    # FIXED SCFM: uses corrected field transfer (0.015 lb O2/hr/scfm)
    # ============================================================
    scfm_required = compute_scfm_required(sor_lb_hr, DIFFUSER_DEPTH_FT)

    # ------------------------------------------------------------
    # TURBO BLOWER
    # ------------------------------------------------------------
    turbo_hp = compute_blower_hp(
        scfm_required,
        BLOWER_PRESSURE_PSI_TURBO,
        TURBO_EFFICIENCY
    )

    # Small plants <1 MGD usually need only 2 turbo blowers
    turbo_units_required = 2 if turbo_hp < 60 else 3
    turbo_firm_ok = turbo_units_required >= 2

    # ------------------------------------------------------------
    # POSITIVE DISPLACEMENT BLOWER
    # ------------------------------------------------------------
    pd_hp = compute_blower_hp(
        scfm_required,
        BLOWER_PRESSURE_PSI_PD,
        PD_EFFICIENCY
    )

    pd_units_required = 2 if pd_hp < 80 else 3
    pd_firm_ok = pd_units_required >= 2

    # ------------------------------------------------------------
    # Return structured blower module
    # ------------------------------------------------------------
    blowers = {
        "airflow": {
            "sor_lb_hr": sor_lb_hr,
            "scfm_required": scfm_required,
        },
        "turbo": {
            "hp_per_blower": turbo_hp,
            "units_required": turbo_units_required,
            "firm_capacity_ok": turbo_firm_ok,
            "pressure_psi": BLOWER_PRESSURE_PSI_TURBO,
            "efficiency": TURBO_EFFICIENCY
        },
        "pd": {
            "hp_per_blower": pd_hp,
            "units_required": pd_units_required,
            "firm_capacity_ok": pd_firm_ok,
            "pressure_psi": BLOWER_PRESSURE_PSI_PD,
            "efficiency": PD_EFFICIENCY
        }
    }

    return blowers



# ============================================================
#  MIXER SIZING
# ============================================================

def _compute_mixers(geometry):
    """
    Computes mixer horsepower using mixing energy density:
      W/m³ ≈ 0.7 (activated sludge / SBR)
    """
    volume_ft3 = geometry["computed"]["volume_ft3_total"]
    hp = compute_mixing_hp(volume_ft3, MIXING_W_PER_M3)

    mixers = {
        "mixing_power_hp": hp,
        "mixers_required": 2 if hp > 5 else 1,
        "mixing_energy_density_Wm3": MIXING_W_PER_M3,
    }
    return mixers


# ============================================================
#  ENERGY MODELING
# ============================================================

def _compute_energy(blowers, cycle, electricity_cost_per_kwh=0.12):
    """
    Computes annual kWh and cost for both turbo and PD blowers.

    Operating assumption:
      Aeration occurs during FIL and REACT.
      Hours per cycle = fill_hr + react_hr.
    """

    fill_hr = cycle["inputs"]["fill_hr"]
    react_hr = cycle["inputs"]["react_hr"]
    cycles_per_day = cycle["inputs"]["cycles_per_day"]

    # Aeration hours per day
    aeration_hours_per_day = (fill_hr + react_hr) * cycles_per_day
    aeration_hours_per_year = aeration_hours_per_day * 365

    # TURBO
    turbo_hp = blowers["turbo"]["hp_per_blower"]
    turbo_energy = annual_energy_cost(
        turbo_hp,
        aeration_hours_per_year,
        electricity_cost_per_kwh
    )

    # PD
    pd_hp = blowers["pd"]["hp_per_blower"]
    pd_energy = annual_energy_cost(
        pd_hp,
        aeration_hours_per_year,
        electricity_cost_per_kwh
    )

    return {
        "aeration_hours_per_day": aeration_hours_per_day,
        "aeration_hours_per_year": aeration_hours_per_year,
        "turbo": turbo_energy,
        "pd": pd_energy,
    }


# ============================================================
#  MODULE 3 RECOMMENDATIONS
# ============================================================

def _compute_recommendations(biology, blowers, mixers):
    """
    Generate a recommendation list based on:
      - SRT
      - Effluent ammonia
      - Blower capacity
      - Mixer power
      - Oxygen demand
    """

    rec = []

    # SRT
    srt_days = biology["srt"]["srt_days"]
    if srt_days < 10:
        rec.append("Increase SRT to maintain stable nitrification, especially in winter.")

    # Effluent ammonia warnings
    if biology["effluent"]["nh3_winter_mgL"] > 2:
        rec.append("Winter NH3 > 2 mg/L — consider increasing aeration capacity or SRT.")

    # Blower recommendations
    if not blowers["turbo"]["firm_capacity_ok"]:
        rec.append("Turbo blower system lacks firm capacity — add one more unit.")

    if not blowers["pd"]["firm_capacity_ok"]:
        rec.append("PD blower system lacks firm capacity — N+1 recommended.")

    # Mixer
    if mixers["mixing_power_hp"] < 3:
        rec.append("Mixer horsepower is low — consider minimum 5 HP per basin.")

    # Oxygen demand
    aor = biology["oxygen"]["aor_lb_hr"]
    if aor > 200:
        rec.append("High AOR — evaluate blower selection and potential DO control.")

    # If no recommendations required
    if len(rec) == 0:
        rec.append("All parameters within typical design ranges.")

    return rec


# ============================================================
#  FINAL ENGINE WRAPPER
# ============================================================

def run_sbr_engine(inputs):
    """
    Main entry point for the SBR Engine v3.

    Computes:
      Module 1 — Geometry + Cycle
      Module 2 — Biology + Oxygen + Effluent
      Module 3 — Blowers + Mixers + Energy
      Recommendations

    Returns a fully nested JSON-like dict.
    """

    # 1. Geometry (Part 2A)
    geometry = _compute_geometry(inputs)

    # 2. Cycle Timing (Part 2B)
    cycle = _compute_cycle_timing(inputs, geometry)

    # 3. Biological Process (Part 2C)
    biology = _compute_biological_process(inputs, geometry, cycle)

    # 4. Blowers (Module 3)
    blowers = _compute_airflow_and_blowers(biology, inputs, geometry)

    # 5. Mixers
    mixers = _compute_mixers(geometry)

    # 6. Energy
    energy = _compute_energy(blowers, cycle)

    # 7. Recommendations
    recommendations = _compute_recommendations(biology, blowers, mixers)

    # -------------------------------------------------------
    # Final Output Block (nested dictionary)
    # -------------------------------------------------------
    return {
        "module_1": {
            "geometry": geometry,
            "cycle": cycle
        },
        "module_2": biology,
        "module_3": {
            "blowers": blowers,
            "mixers": mixers,
            "energy": energy,
            "recommendations": recommendations
        }
    }
