"""
===========================================================
RealWW – Oxidation Ditch Design Module
===========================================================
"""

from backend.core.constants import (
    DEFAULT_BOD,
    DEFAULT_TKN,
    O2_REQUIRED_PER_LB_BOD,
    O2_REQUIRED_PER_LB_NH3,
    GAL_TO_CUFT,
)

# ============================================================
# Helper Functions
# ============================================================

def bod_load_lb_day(flow_mgd, bod_mgL=DEFAULT_BOD):
    """Return total BOD load as lb/day."""
    return flow_mgd * bod_mgL * 8.34  # mg/L → lb/day


def nh3_load_lb_day(flow_mgd, tkn_mgL=DEFAULT_TKN):
    """Return TKN load as lb/day (approx ammonia load)."""
    return flow_mgd * tkn_mgL * 8.34


def tank_volume_ft3(length_ft, width_ft, depth_ft):
    """Compute oxidation ditch volume (ft³)."""
    return length_ft * width_ft * depth_ft


# ============================================================
# Main Design Function
# ============================================================

def design_oxidation_ditch(
    flow_mgd,
    bod_mgL=DEFAULT_BOD,
    tkn_mgL=DEFAULT_TKN,
    length_ft=131.58,
    width_ft=32.67,
    depth_ft=12.0,
    rotors_total_o2_capacity_lb_day=5088,
):
    """
    Perform a full oxidation ditch design evaluation.
    Uses Hodgenville geometry unless overridden.
    """

    # ---------------------------
    # Loads
    # ---------------------------
    bod_lb_day = bod_load_lb_day(flow_mgd, bod_mgL)
    nh3_lb_day = nh3_load_lb_day(flow_mgd, tkn_mgL)

    # ---------------------------
    # Oxygen demand
    # ---------------------------
    o2_bod = bod_lb_day * O2_REQUIRED_PER_LB_BOD
    o2_nh3 = nh3_lb_day * O2_REQUIRED_PER_LB_NH3
    total_o2_demand = o2_bod + o2_nh3

    # ---------------------------
    # Volume & detention time
    # ---------------------------
    volume_ft3 = tank_volume_ft3(length_ft, width_ft, depth_ft)
    volume_mg = (volume_ft3 * GAL_TO_CUFT) / 1_000_000  # ft³ → MG

    detention_time_hr = (volume_mg / flow_mgd) * 24

    # ---------------------------
    # Volumetric loading
    # ---------------------------
    vol_loading = bod_lb_day / (volume_ft3 / 1000)

    meets_vol_loading = (8 <= vol_loading <= 15)  # Ten States guideline
    o2_deficit = total_o2_demand - rotors_total_o2_capacity_lb_day
    meets_o2 = o2_deficit <= 0

    # ---------------------------
    # Return dictionary (NEVER returns None)
    # ---------------------------
    return {
        "input": {
            "flow_mgd": flow_mgd,
            "bod_mgL": bod_mgL,
            "tkn_mgL": tkn_mgL,
            "geometry_ft": (length_ft, width_ft, depth_ft),
            "rotor_capacity_lb_day": rotors_total_o2_capacity_lb_day,
        },
        "loads": {
            "bod_lb_day": bod_lb_day,
            "nh3_lb_day": nh3_lb_day,
            "o2_bod": o2_bod,
            "o2_nh3": o2_nh3,
            "total_o2_demand": total_o2_demand,
        },
        "hydraulics": {
            "volume_ft3": volume_ft3,
            "volume_mg": volume_mg,
            "detention_time_hr": detention_time_hr,
            "volumetric_loading_lb_1000ft3": vol_loading,
        },
        "standards_check": {
            "meets_volumetric_loading": meets_vol_loading,
            "meets_oxygen_requirements": meets_o2,
            "oxygen_deficit_lb_day": o2_deficit,
        }
    }
