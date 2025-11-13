# backend/sbr/design.py

from __future__ import annotations

from typing import Any, Dict, Optional

# --------------------------------------------------------------------
# 1. CONSTANTS (override with backend.core.constants if available)
# --------------------------------------------------------------------

# Reasonable defaults, will be overwritten if backend.core.constants exists.
MGD_TO_GPD = 1_000_000.0                # gal/day per MGD
FT3_PER_MG = 133_680.0                  # ft³ per MG (1 MG ≈ 133,680 ft³)
LB_PER_MG_PER_MG_L = 8.34               # lb per (MGD · mg/L)
HOURS_PER_DAY = 24.0

# Oxygen requirements (typical design values; will be overridden if defined)
O2_REQUIRED_PER_LB_BOD = 1.1           # lb O₂ / lb BOD removed (includes endogenous)
O2_REQUIRED_PER_LB_NH3 = 4.6           # lb O₂ / lb NH₃-N oxidized

# Default influent characteristics (medium strength domestic)
DEFAULT_BOD_MG_L = 250.0
DEFAULT_TSS_MG_L = 250.0
DEFAULT_NH3_MG_L = 30.0

try:
    # Try to use project-wide constants if they exist
    from backend.core import constants as C  # type: ignore[attr-defined]

    MGD_TO_GPD = getattr(C, "MGD_TO_GPD", MGD_TO_GPD)
    FT3_PER_MG = getattr(C, "FT3_PER_MG", FT3_PER_MG)
    LB_PER_MG_PER_MG_L = getattr(C, "LB_PER_MG_PER_MG_L", LB_PER_MG_PER_MG_L)
    HOURS_PER_DAY = getattr(C, "HOURS_PER_DAY", HOURS_PER_DAY)

    O2_REQUIRED_PER_LB_BOD = getattr(
        C, "O2_REQUIRED_PER_LB_BOD", O2_REQUIRED_PER_LB_BOD
    )
    O2_REQUIRED_PER_LB_NH3 = getattr(
        C, "O2_REQUIRED_PER_LB_NH3", O2_REQUIRED_PER_LB_NH3
    )

    DEFAULT_BOD_MG_L = getattr(C, "DEFAULT_BOD_MG_L", DEFAULT_BOD_MG_L)
    DEFAULT_TSS_MG_L = getattr(C, "DEFAULT_TSS_MG_L", DEFAULT_TSS_MG_L)
    DEFAULT_NH3_MG_L = getattr(C, "DEFAULT_NH3_MG_L", DEFAULT_NH3_MG_L)
except Exception:
    # If backend.core.constants is missing, just use the defaults above.
    pass

# --------------------------------------------------------------------
#  Biological process defaults (Metcalf & Eddy typical values)
# --------------------------------------------------------------------
DEFAULT_MLSS = 3000        # mg/L typical SBR MLSS
DEFAULT_VSS_FRACTION = 0.75  # 75% of MLSS is volatile
DEFAULT_EFF_TSS = 10       # mg/L effluent suspended solids
DEFAULT_TARGET_SRT = 10    # days
DEFAULT_WAS_SOLIDS = 0.8   # WAS solids concentration = 0.8 * MLSS

# --------------------------------------------------------------------
# 2. CORE DESIGN FUNCTION
# --------------------------------------------------------------------


def design_sbr(
    flow_mgd: float,
    bod_mgL: float = DEFAULT_BOD_MG_L,
    tss_mgL: float = DEFAULT_TSS_MG_L,
    nh3_mgL: float = DEFAULT_NH3_MG_L,
    n_basins: int = 2,
    cycles_per_day: Optional[float] = None,
    # Ten States default cycle:
    fill_hr: float = 1.0,
    react_hr: float = 4.0,
    settle_hr: float = 1.5,
    decant_hr: float = 0.75,
    idle_hr: float = 0.5,
    # Volume / geometry assumptions
    design_detention_hr: float = 24.0,    # ≥ 24 hr required by spec
    basin_depth_ft: float = 18.0,         # typical SBR sidewater depth
    decant_fraction_of_depth: float = 0.30,
    solids_storage_fraction: float = 0.15,
) -> Dict[str, Any]:
    """
    RealWW SBR design engine.

    Returns a full design dictionary with:
      - inputs
      - loads
      - cycle
      - basin
      - decanter
      - oxygen
      - ten_states (compliance checks + advisories)

    This function NEVER returns None.
    """

    # Guard against pathological input, but still return a valid dict
    if flow_mgd <= 0.0:
        return {
            "inputs": {
                "flow_mgd": flow_mgd,
                "bod_mgL": bod_mgL,
                "tss_mgL": tss_mgL,
                "nh3_mgL": nh3_mgL,
                "n_basins": n_basins,
            },
            "loads": {
                "bod_lb_day": 0.0,
                "tss_lb_day": 0.0,
                "nh3_lb_day": 0.0,
            },
            "cycle": {
                "fill_hr": fill_hr,
                "react_hr": react_hr,
                "settle_hr": settle_hr,
                "decant_hr": decant_hr,
                "idle_hr": idle_hr,
                "total_cycle_time_hr": fill_hr + react_hr + settle_hr + decant_hr + idle_hr,
                "cycles_per_day": 0.0,
                "notes": "Flow ≤ 0. Design not meaningful; all capacities set to zero.",
            },
            "basin": {},
            "decanter": {},
            "oxygen": {},
            "ten_states": {
                "overall_pass": False,
                "checks": {
                    "min_two_basins": False,
                    "independent_operation_assumed": False,
                    "cycle_time_adequate": False,
                    "min_24hr_detention": False,
                    "max_decant_rate_4x_avg": False,
                    "max_drawdown_2ft_hr": False,
                    "min_react_three_hr": react_hr >= 3.0,
                },
                "required_decant_hr_for_2ft_hr_limit": None,
                "required_decant_hr_for_4x_avg_limit": None,
                "recommended_decant_hr": None,
                "advisories": [
                    "Flow is non-positive; Ten States compliance cannot be demonstrated."
                ],
                "notes": "Flow is non-positive; Ten States compliance cannot be demonstrated.",
            },
        }

    # ----------------------------------------------------------------
    # 2A. Influent Loads
    # ----------------------------------------------------------------
    Q_mgd = flow_mgd
    Q_MG_per_day = Q_mgd  # since units are MG/day already

    bod_lb_day = Q_MG_per_day * bod_mgL * LB_PER_MG_PER_MG_L
    tss_lb_day = Q_MG_per_day * tss_mgL * LB_PER_MG_PER_MG_L
    nh3_lb_day = Q_MG_per_day * nh3_mgL * LB_PER_MG_PER_MG_L

    loads = {
        "bod_lb_day": bod_lb_day,
        "tss_lb_day": tss_lb_day,
        "nh3_lb_day": nh3_lb_day,
    }

    # ----------------------------------------------------------------
    # 2B. SBR Cycle Design
    # ----------------------------------------------------------------
    total_cycle_time_hr = fill_hr + react_hr + settle_hr + decant_hr + idle_hr

    if total_cycle_time_hr <= 0.0:
        # Defensive; but still keep a meaningful cycle
        total_cycle_time_hr = 8.0

    if cycles_per_day is None or cycles_per_day <= 0.0:
        cycles_per_day = HOURS_PER_DAY / total_cycle_time_hr

    cycle_time_from_cycles_hr = HOURS_PER_DAY / cycles_per_day

    if abs(cycle_time_from_cycles_hr - total_cycle_time_hr) > 1e-3:
        cycle_note = (
            "Cycle components and cycles_per_day are not exactly consistent; "
            "cycles_per_day derived from component times."
        )
        cycles_per_day = HOURS_PER_DAY / total_cycle_time_hr
    else:
        cycle_note = "Cycle components and cycles_per_day are consistent."

    cycle_info = {
        "fill_hr": fill_hr,
        "react_hr": react_hr,
        "settle_hr": settle_hr,
        "decant_hr": decant_hr,
        "idle_hr": idle_hr,
        "total_cycle_time_hr": total_cycle_time_hr,
        "cycles_per_day": cycles_per_day,
        "cycle_time_from_cycles_hr": cycle_time_from_cycles_hr,
        "notes": cycle_note,
    }

    # ----------------------------------------------------------------
    # 2C. Basin Volume Requirements
    # ----------------------------------------------------------------
    # Required total working volume for ≥ 24 hr detention at average flow
    design_detention_days = design_detention_hr / HOURS_PER_DAY
    required_total_working_volume_MG = Q_mgd * design_detention_days

    # Include additional storage for mixed liquor / solids in react + settle
    total_volume_including_solids_MG = required_total_working_volume_MG * (
        1.0 + solids_storage_fraction
    )

    # Enforce minimum 2 basins by design
    n_basins = max(2, int(n_basins))

    volume_per_basin_MG = total_volume_including_solids_MG / n_basins
    volume_per_basin_ft3 = volume_per_basin_MG * FT3_PER_MG

    # Sidewater depth, surface area, and decant drawdown
    depth_ft_max = basin_depth_ft
    if depth_ft_max <= 0.0:
        depth_ft_max = 18.0

    surface_area_ft2 = volume_per_basin_ft3 / depth_ft_max

    # Assume a certain fraction of the depth is drawn down during decant
    decant_drawdown_depth_ft = depth_ft_max * decant_fraction_of_depth
    depth_ft_min = depth_ft_max - decant_drawdown_depth_ft

    # Actual detention time provided at this volume
    actual_total_working_volume_MG = volume_per_basin_MG * n_basins
    actual_detention_hr = (actual_total_working_volume_MG / Q_mgd) * HOURS_PER_DAY

    basin_info = {
        "n_basins": n_basins,
        "required_total_working_volume_MG": required_total_working_volume_MG,
        "total_volume_including_solids_MG": total_volume_including_solids_MG,
        "volume_per_basin_MG": volume_per_basin_MG,
        "volume_per_basin_ft3": volume_per_basin_ft3,
        "basin_depth_ft_max": depth_ft_max,
        "basin_depth_ft_min": depth_ft_min,
        "basin_decant_drawdown_depth_ft": decant_drawdown_depth_ft,
        "surface_area_ft2": surface_area_ft2,
        "design_detention_hr_target": design_detention_hr,
        "actual_detention_hr": actual_detention_hr,
    }

    # ----------------------------------------------------------------
    # 2D. Air / Oxygen Requirements
    # ----------------------------------------------------------------
    # BOD + nitrification demand
    o2_bod_lb_day = bod_lb_day * O2_REQUIRED_PER_LB_BOD
    o2_nh3_lb_day = nh3_lb_day * O2_REQUIRED_PER_LB_NH3
    o2_total_lb_day = o2_bod_lb_day + o2_nh3_lb_day

    # Aeration only runs during react-period
    total_react_hr_per_day = react_hr * cycles_per_day
    if total_react_hr_per_day <= 0.0:
        total_react_hr_per_day = 1.0  # defensive

    o2_lb_per_hr_during_react = o2_total_lb_day / total_react_hr_per_day

    oxygen_info = {
        "o2_bod_lb_day": o2_bod_lb_day,
        "o2_nh3_lb_day": o2_nh3_lb_day,
        "o2_total_lb_day": o2_total_lb_day,
        "total_react_hr_per_day": total_react_hr_per_day,
        "o2_lb_per_hr_during_react": o2_lb_per_hr_during_react,
    }

    # ----------------------------------------------------------------
    # 2E. Decanter Sizing
    # ----------------------------------------------------------------
    # Volume decanted per event is based on average daily flow
    decant_events_per_day = n_basins * cycles_per_day
    if decant_events_per_day <= 0.0:
        decant_events_per_day = 1.0

    volume_per_decant_MG = Q_mgd / decant_events_per_day  # MG/event

    # Decant rate as MGD and gpm
    decant_time_days = decant_hr / HOURS_PER_DAY
    if decant_time_days <= 0.0:
        decant_time_days = 1.0 / HOURS_PER_DAY  # 1 hr; defensive

    decant_rate_MGD = volume_per_decant_MG / decant_time_days
    decant_rate_gpm = decant_rate_MGD * (MGD_TO_GPD / 1440.0)

    # Drawdown velocity (ft/hr)
    drawdown_velocity_ft_hr = (
        decant_drawdown_depth_ft / decant_hr if decant_hr > 0.0 else 0.0
    )

    decanter_info = {
        "decant_events_per_day": decant_events_per_day,
        "volume_per_decant_MG": volume_per_decant_MG,
        "decant_rate_MGD": decant_rate_MGD,
        "decant_rate_gpm": decant_rate_gpm,
        "drawdown_velocity_ft_hr": drawdown_velocity_ft_hr,
        "limits": {
            "max_drawdown_ft_hr": 2.0,
            "max_decant_rate_multiplier_of_avg": 4.0,
        },
    }

    # ----------------------------------------------------------------
    # 2F. Biological Process Calculations (MLSS, SRT, F:M)
    # ----------------------------------------------------------------

    # MLSS and MLVSS
    mlss = DEFAULT_MLSS  # user input could be added in future
    mlvss = mlss * DEFAULT_VSS_FRACTION

    # Total active basin volume (reactor volume) = total working volume
    reactor_volume_MG = actual_total_working_volume_MG
    reactor_volume_ft3 = reactor_volume_MG * FT3_PER_MG

    # F/M ratio
    # BOD load (lb/day) / (MLVSS mass in reactor, lb)
    # MLVSS mass = volume * mg/L * 8.34
    mlvss_mass_lb = reactor_volume_MG * mlvss * LB_PER_MG_PER_MG_L
    f_m_ratio = bod_lb_day / mlvss_mass_lb

    # SRT calculation
    # WAS solids concentration (mg/L)
    Xw = mlss * DEFAULT_WAS_SOLIDS
    Xe = DEFAULT_EFF_TSS

    # SRT target → compute required WAS flow
    # SRT = (V * MLSS) / (Qw*Xw + Qe*Xe)
    # Solve for Qw:
    # Qw = (V*MLSS / SRT - Qe*Xe) / Xw
    V_MG = reactor_volume_MG
    V_MG_per_L = V_MG * 1_000_000  # convert MG → gallons → liters? NO.
    # Actually: SRT formula uses:
    # lb solids = MLSS (mg/L) * 8.34 lbs/MG/L * MG volume

    total_solids_lb = V_MG * mlss * LB_PER_MG_PER_MG_L
    eff_flow_mgd = Q_mgd  # effluent flow ~ influent flow

    required_was_flow_mgd = (
        (total_solids_lb / DEFAULT_TARGET_SRT) - (eff_flow_mgd * Xe * LB_PER_MG_PER_MG_L)
    ) / (Xw * LB_PER_MG_PER_MG_L)

    if required_was_flow_mgd < 0:
        required_was_flow_mgd = 0

    biology_info = {
        "mlss_mgL": mlss,
        "mlvss_mgL": mlvss,
        "reactor_volume_MG": reactor_volume_MG,
        "reactor_volume_ft3": reactor_volume_ft3,
        "f_m_ratio": f_m_ratio,
        "srt_days": DEFAULT_TARGET_SRT,
        "was_flow_mgd": required_was_flow_mgd,
        "was_flow_gpd": required_was_flow_mgd * 1_000_000,
    }


    # ----------------------------------------------------------------
    # 2F. Ten States Standards Checks + Advisory Calculations
    # ----------------------------------------------------------------
    checks: Dict[str, bool] = {}

    # 1) ≥ 2 basins required
    checks["min_two_basins"] = n_basins >= 2

    # 2) Basins must be capable of independent operation
    checks["independent_operation_assumed"] = n_basins >= 2

    # 3) Cycle time adequacy — practical range
    checks["cycle_time_adequate"] = (
        4.0 <= total_cycle_time_hr <= 12.0
        and 3.0 <= cycles_per_day <= 12.0
    )

    # 4) ≥ 24 hr total detention time at average flow
    checks["min_24hr_detention"] = actual_detention_hr >= 24.0

    # 5) Decant rate ≤ 4× average flow
    max_decant_rate_allowed_MGD = 4.0 * Q_mgd
    checks["max_decant_rate_4x_avg"] = decant_rate_MGD <= max_decant_rate_allowed_MGD

    # 6) Drawdown velocity ≤ 2 ft/hr
    max_drawdown_allowed_ft_hr = 2.0
    checks["max_drawdown_2ft_hr"] = drawdown_velocity_ft_hr <= max_drawdown_allowed_ft_hr

    # 7) Minimum aeration period ≥ 3 hr react
    checks["min_react_three_hr"] = react_hr >= 3.0

    # ---- Advisory calculations ----

    # For 2 ft/hr drawdown limit:
    required_decant_hr_for_2ft_hr_limit = (
        decant_drawdown_depth_ft / max_drawdown_allowed_ft_hr
        if max_drawdown_allowed_ft_hr > 0
        else decant_hr
    )

    # For 4× average flow limit:
    # decant_rate_MGD = volume_per_decant_MG / (decant_hr / 24)
    # Solve for decant_hr:
    required_decant_hr_for_4x_avg_limit = (
        (volume_per_decant_MG / max_decant_rate_allowed_MGD) * HOURS_PER_DAY
        if max_decant_rate_allowed_MGD > 0
        else decant_hr
    )

    # Controlling requirement
    recommended_decant_hr = max(
        required_decant_hr_for_2ft_hr_limit,
        required_decant_hr_for_4x_avg_limit,
    )

    advisories = []

    if not checks["max_drawdown_2ft_hr"]:
        advisories.append(
            f"Increase decant time to ≥ {required_decant_hr_for_2ft_hr_limit:.2f} hr "
            "to satisfy drawdown velocity (≤ 2 ft/hr)."
        )

    if not checks["max_decant_rate_4x_avg"]:
        advisories.append(
            f"Increase decant time to ≥ {required_decant_hr_for_4x_avg_limit:.2f} hr "
            "to satisfy decant rate (≤ 4× average flow)."
        )

    advisories.append(
        f"Recommended decant time: {recommended_decant_hr:.2f} hr."
    )

    overall_pass = all(checks.values())

    ten_states_info = {
        "overall_pass": overall_pass,
        "checks": checks,
        "required_decant_hr_for_2ft_hr_limit": required_decant_hr_for_2ft_hr_limit,
        "required_decant_hr_for_4x_avg_limit": required_decant_hr_for_4x_avg_limit,
        "recommended_decant_hr": recommended_decant_hr,
        "advisories": advisories,
        "notes": (
            "Ten States checks include advisory values indicating the minimum "
            "design changes required to achieve compliance."
        ),
    }

    # ----------------------------------------------------------------
    # 2G. Assemble full dictionary and return (never None)
    # ----------------------------------------------------------------
    results: Dict[str, Any] = {
        "inputs": {
            "flow_mgd": flow_mgd,
            "bod_mgL": bod_mgL,
            "tss_mgL": tss_mgL,
            "nh3_mgL": nh3_mgL,
            "n_basins": n_basins,
            "cycles_per_day": cycles_per_day,
            "fill_hr": fill_hr,
            "react_hr": react_hr,
            "settle_hr": settle_hr,
            "decant_hr": decant_hr,
            "idle_hr": idle_hr,
            "design_detention_hr": design_detention_hr,
            "basin_depth_ft": basin_depth_ft,
            "decant_fraction_of_depth": decant_fraction_of_depth,
            "solids_storage_fraction": solids_storage_fraction,
        },
        "loads": loads,
        "cycle": cycle_info,
        "basin": basin_info,
        "decanter": decanter_info,
        "oxygen": oxygen_info,
        "ten_states": ten_states_info,
    }

    return results
# --------------------------------------------------------------------
# AUTO-TUNE FUNCTION
# --------------------------------------------------------------------

def design_sbr_autotune(**kwargs) -> Dict[str, Any]:
    """
    Runs design_sbr(), reads recommended decant time from Ten States,
    re-runs design with corrected decant time, and returns both models.
    """

    original = design_sbr(**kwargs)

    ten_states = original.get("ten_states", {})
    recommended = ten_states.get("recommended_decant_hr", None)
    original_decant_hr = original["inputs"]["decant_hr"]


    changes = []

    if recommended is not None and original_decant_hr is not None:
        if recommended > original_decant_hr + 1e-6:
            changes.append(
                f"Decant time increased from {original_decant_hr:.2f} hr to {recommended:.2f} hr."
            )
            kwargs["decant_hr"] = recommended

    autotuned = design_sbr(**kwargs)

    return {
        "original_design": original,
        "autotuned_design": autotuned,
        "changes_made": changes,
        "notes": "Auto-tune adjusted parameters to satisfy Ten States constraints."
    }
# --------------------------------------------------------------------
# FULL CYCLE AUTO-TUNE (Option 1)
# --------------------------------------------------------------------

def design_sbr_autotune_full_cycle(**kwargs) -> Dict[str, Any]:
    """
    Automatically adjusts:
      - Decant time (first)
      - Idle time (down to 0)
      - Fill time (down to 0.2 hr)
      - Settle time (down to 1.0 hr)
    until:
      - cycles/day ≥ 3
      - all Ten States requirements pass (or no further adjustments possible)

    Returns:
      {
        "original_design": ...,
        "autotuned_design": ...,
        "changes_made": [...],
        "notes": "..."
      }
    """

    MIN_FILL = 0.2       # hr
    MIN_SETTLE = 1.0     # hr
    MIN_IDLE = 0.0       # hr
    MIN_REACT = 3.0      # hr  (cannot reduce below this)

    changes = []

    # ---- PART 1: initial run ----
    original = design_sbr(**kwargs)

    # Extract original cycle values
    cycle = original["cycle"]
    fill_hr = cycle["fill_hr"]
    react_hr = cycle["react_hr"]
    settle_hr = cycle["settle_hr"]
    decant_hr = cycle["decant_hr"]
    idle_hr = cycle["idle_hr"]

    # ---- PART 2: decant auto-tune first ----
    reco_decant = original["ten_states"]["recommended_decant_hr"]
    if reco_decant > decant_hr + 1e-6:
        changes.append(
            f"Decant time increased from {decant_hr:.2f} hr to {reco_decant:.2f} hr."
        )
        decant_hr = reco_decant
        kwargs["decant_hr"] = decant_hr

    # Run again with updated decant
    tuned = design_sbr(**kwargs)

    # ---- PART 3: Ensure cycles/day ≥ 3 ----
    # Try reducing idle → fill → settle
    def recompute():
        kw = dict(kwargs)
        kw["fill_hr"] = fill_hr
        kw["react_hr"] = react_hr
        kw["settle_hr"] = settle_hr
        kw["decant_hr"] = decant_hr
        kw["idle_hr"] = idle_hr
        return design_sbr(**kw)

    tuned = recompute()
    cycles_per_day = tuned["cycle"]["cycles_per_day"]

    # Adjust cycle elements until cycles/day ≥ 3
    while cycles_per_day < 3.0 - 1e-6:

        # 1. Reduce idle first
        if idle_hr > MIN_IDLE:
            old = idle_hr
            idle_hr = max(MIN_IDLE, idle_hr - 0.25)
            changes.append(f"Idle reduced from {old:.2f} hr to {idle_hr:.2f} hr.")
        # 2. Then reduce fill
        elif fill_hr > MIN_FILL:
            old = fill_hr
            fill_hr = max(MIN_FILL, fill_hr - 0.25)
            changes.append(f"Fill reduced from {old:.2f} hr to {fill_hr:.2f} hr.")
        # 3. Then reduce settle
        elif settle_hr > MIN_SETTLE:
            old = settle_hr
            settle_hr = max(MIN_SETTLE, settle_hr - 0.25)
            changes.append(f"Settle reduced from {old:.2f} hr to {settle_hr:.2f} hr.")
        else:
            # Cannot reduce further → break
            changes.append("No further cycle reduction possible to reach ≥3 cycles/day.")
            break

        tuned = recompute()
        cycles_per_day = tuned["cycle"]["cycles_per_day"]

    # ---- PART 4: final check: if Ten States still fails, report but keep best solution
    final_checks = tuned["ten_states"]["checks"]
    all_pass = tuned["ten_states"]["overall_pass"]

    notes = (
        "Full cycle auto-tune completed. "
        "All Ten States requirements satisfied."
        if all_pass
        else "Auto-tune reached closest feasible solution; some Ten States checks still unmet."
    )

    return {
        "original_design": original,
        "autotuned_design": tuned,
        "changes_made": changes,
        "notes": notes,
    }
