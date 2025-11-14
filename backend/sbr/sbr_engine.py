# ============================================================
# SBR ENGINE — Unified Module 1 + Module 2
# RealWW Backend / Python 3.12 / OOP Architecture
# ------------------------------------------------------------
# This module provides:
#  - Cycle design
#  - Basin geometry
#  - MLSS & SRT tuning
#  - Oxygen demand (BOD + NH3)
#  - Effluent BOD & NH3 predictions
#  - Nitrification kinetic model
#  - Winter temperature checks
#  - Warning & recommendation engine
# Output is nested dict for GUI integration.
# ============================================================

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List
import math


# ============================================================
# INPUT DATA MODEL
# ============================================================

@dataclass
class SBRInputs:
    # ---------------------------------------
    # Hydraulic / Load Inputs
    # ---------------------------------------
    flow_mgd: float
    bod_mgL: float
    tss_mgL: float
    nh3_mgL: float
    phosphorus_mgL: float = 0.0

    # ---------------------------------------
    # Basin / Geometry
    # ---------------------------------------
    n_basins: int = 2
    basin_depth_ft: float = 18.0
    freeboard_ft: float = 2.0

    # Optional future geometry expansion
    basin_shape: str = "rectangular"     # future expansion
    basin_length_ft: float = 100.0
    basin_width_ft: float = 60.0

    # ---------------------------------------
    # Operational Controls
    # ---------------------------------------
    cycles_per_day: float = 3.0
    fill_hr: float = 0.5
    react_hr: float = 4.0
    settle_hr: float = 1.0
    decant_hr: float = 1.5
    idle_hr: float = 0.0

    # MLSS / SRT Controls
    mlss_mgL: float = 3_000.0
    srt_target_days: float = 15.0
    solids_yield_coeff: float = 0.65  # typical for extended aeration

    # Decant Controls
    decant_fraction_of_depth: float = 0.30

    # DO / Aeration Inputs
    target_DO_mgL: float = 2.0
    alpha: float = 0.8
    beta: float = 1.0
    theta_DO: float = 1.024
    temp_C: float = 20.0

    # Mixing / Aeration Efficiency
    ote_standard_lbO2_hp_hr: float = 8.0
    field_aeration_efficiency: float = 0.65

    # Nitrification
    mu_max_20C: float = 0.9    # 1/day
    ks_NH3: float = 1.0        # mg/L
    kd_decay: float = 0.08     # 1/day
    theta_nitrifiers: float = 1.123

    # Winter Scenarios
    winter_temp_C: float = 8.0

    # Ten States
    max_MLSS: float = 10_000.0
    min_MLSS: float = 1_500.0


# ============================================================
# MAIN ENGINE CLASS
# ============================================================

class SBRSystem:

    def __init__(self, inputs: SBRInputs):
        self.inputs = inputs

    # ========================================================
    # PUBLIC API CALL
    # ========================================================
    def run_full_design(self) -> dict[str, Any]:
        """Run entire SBR design sequence and return nested dict."""
        cycle = self._calc_cycle_design()
        basin = self._calc_basin_design(cycle)
        biomass = self._calc_biomass_and_srt(basin)
        oxygen = self._calc_oxygen_demand(basin, biomass)
        eff = self._calc_effluent_quality(basin, biomass, oxygen)

        warnings = self._collect_warnings(basin, biomass, oxygen, eff)
        recommendations = self._collect_recommendations(basin, biomass, oxygen, eff)

        return {
            "inputs": vars(self.inputs),
            "cycle_design": cycle,
            "basin_design": basin,
            "oxygen_demand": oxygen,
            "effluent_quality": eff,
            "warnings": warnings,
            "recommendations": recommendations,
        }
    # ========================================================
    # INTERNAL METHOD 1 — CYCLE DESIGN
    # ========================================================
    def _calc_cycle_design(self) -> dict[str, Any]:
        """Compute full SBR cycle timing and cycles/day validation."""

        i = self.inputs

        cycle_hours = (
            i.fill_hr +
            i.react_hr +
            i.settle_hr +
            i.decant_hr +
            i.idle_hr
        )

        cycles_per_day_actual = 24.0 / cycle_hours

        return {
            "cycle_hours": cycle_hours,
            "cycles_per_day_input": i.cycles_per_day,
            "cycles_per_day_actual": cycles_per_day_actual,
            "fill_hr": i.fill_hr,
            "react_hr": i.react_hr,
            "settle_hr": i.settle_hr,
            "decant_hr": i.decant_hr,
            "idle_hr": i.idle_hr,
        }


    # ========================================================
    # INTERNAL METHOD 2 — BASIN DESIGN (VOLUME, AREA, DEPTH)
    # ========================================================
    def _calc_basin_design(self, cycle: dict[str, Any]) -> dict[str, Any]:
        """Compute basin volume per basin, operating volume, decant volume."""

        i = self.inputs

        # -----------------------------------------
        # Flow per cycle
        # -----------------------------------------
        Q_gpd = i.flow_mgd * 1_000_000.0
        cycles_per_day = cycle["cycles_per_day_input"]

        flow_per_cycle_gal = Q_gpd / max(cycles_per_day, 1e-6)

        # -----------------------------------------
        # Basin dimensions: rectangular for now
        # -----------------------------------------
        surface_area_ft2 = i.basin_length_ft * i.basin_width_ft
        sidewater_depth_ft = i.basin_depth_ft

        total_depth_ft = sidewater_depth_ft + i.freeboard_ft

        # -----------------------------------------
        # Volume (geometric)
        # -----------------------------------------
        geometric_volume_ft3 = surface_area_ft2 * sidewater_depth_ft
        geometric_volume_gal = geometric_volume_ft3 * 7.48

        # -----------------------------------------
        # Operating Volume (per cycle fill)
        # -----------------------------------------
        # "fill depth" is tied to flow per cycle
        fill_depth_ft = (flow_per_cycle_gal / 7.48) / surface_area_ft2

        # -----------------------------------------
        # Decant Volume
        # -----------------------------------------
        decant_depth_ft = i.decant_fraction_of_depth * sidewater_depth_ft
        decant_volume_ft3 = surface_area_ft2 * decant_depth_ft
        decant_volume_gal = decant_volume_ft3 * 7.48

        return {
            "surface_area_ft2": surface_area_ft2,
            "sidewater_depth_ft": sidewater_depth_ft,
            "total_depth_ft": total_depth_ft,
            "geometric_volume_gal": geometric_volume_gal,
            "geometric_volume_ft3": geometric_volume_ft3,
            "flow_per_cycle_gal": flow_per_cycle_gal,
            "fill_depth_ft": fill_depth_ft,
            "decant_depth_ft": decant_depth_ft,
            "decant_volume_gal": decant_volume_gal,
            "n_basins": i.n_basins,
        }
    # ========================================================
    # INTERNAL METHOD 3 — BIOMASS, MLSS, AND SRT MODEL
    # ========================================================
    def _calc_biomass_and_srt(self, basin: dict[str, Any]) -> dict[str, Any]:
        """
        Compute:
          • MLSS inventory (lb)
          • MLVSS (assume 0.8 fraction)
          • SRT (days)
          • WAS rate required to hit target SRT
        """

        i = self.inputs

        # -----------------------------------------
        # Basic volumes
        # -----------------------------------------
        volume_gal_per_basin = basin["geometric_volume_gal"]
        total_volume_MG = (volume_gal_per_basin * i.n_basins) / 1_000_000.0

        # -----------------------------------------
        # MLSS inventory (lb)
        # MLSS (mg/L) × 8.34 × MG of aeration volume
        # -----------------------------------------
        mlss_mgL = i.mlss_mgL
        mlss_lb = mlss_mgL * 8.34 * total_volume_MG

        # VSS fraction
        mlvss_lb = 0.80 * mlss_lb

        # -----------------------------------------
        # Daily solids production (lb/day)
        # Based on influent BOD loading × yield
        # -----------------------------------------
        influent_bod_lbday = i.flow_mgd * i.bod_mgL * 8.34
        daily_biomass_production_lb = influent_bod_lbday * i.solids_yield_coeff

        # -----------------------------------------
        # Target SRT (user sets)
        # -----------------------------------------
        srt_target = i.srt_target_days

        # To maintain SRT, WAS required per day:
        # WAS (lb/day) = (Inventory lb) / SRT(days)
        was_required_lbday = mlss_lb / max(srt_target, 1e-6)

        # -----------------------------------------
        # Actual achieved SRT based on current WAS rules
        # For now, assume WAS_actual = WAS_required (can be upgraded)
        # -----------------------------------------
        srt_actual = mlss_lb / max(was_required_lbday, 1e-6)

        # -----------------------------------------
        # Return block
        # -----------------------------------------
        return {
            "mlss_mgL": mlss_mgL,
            "mlss_lb": mlss_lb,
            "mlvss_lb": mlvss_lb,
            "total_volume_MG": total_volume_MG,
            "daily_biomass_prod_lb": daily_biomass_production_lb,
            "srt_target_days": srt_target,
            "was_required_lbday": was_required_lbday,
            "srt_actual_days": srt_actual,
        }
    # ========================================================
    # INTERNAL METHOD 4 — OXYGEN DEMAND (AOR + SOR)
    # ========================================================
    def _calc_oxygen_demand(
        self,
        basin: dict[str, Any],
        biomass: dict[str, Any]
    ) -> dict[str, Any]:
        """
        Compute full oxygen demand:
          • Carbonaceous BOD oxidation
          • Nitrification oxygen
          • Endogenous respiration
          • AOR (actual oxygen requirement)
          • SOR (standard oxygen requirement)
        """

        i = self.inputs

        # -----------------------------------------
        # Influent loads
        # -----------------------------------------
        Q_mgd = i.flow_mgd
        bod_inf_mgL = i.bod_mgL
        nh3_inf_mgL = i.nh3_mgL

        bod_lbday = Q_mgd * bod_inf_mgL * 8.34
        nh3_lbday = Q_mgd * nh3_inf_mgL * 8.34

        # -----------------------------------------
        # Oxygen for carbonaceous BOD removal
        # O2 ≈ 1.1 lb O2 / lb BOD5 removed
        # -----------------------------------------
        O2_cBOD = 1.1 * bod_lbday

        # -----------------------------------------
        # Nitrification oxygen
        # 4.57 lb O2 / lb NH3-N oxidized (M&E 5th Ed)
        # -----------------------------------------
        O2_N = 4.57 * nh3_lbday

        # -----------------------------------------
        # Endogenous respiration
        # Assume 15% of BOD load (typical extended aeration)
        # -----------------------------------------
        O2_endogenous = 0.15 * bod_lbday

        # -----------------------------------------
        # Total Oxygen Requirement (lb/day)
        # -----------------------------------------
        total_O2_lbday = O2_cBOD + O2_N + O2_endogenous

        # -----------------------------------------
        # Aerated hours per day (from cycle)
        # -----------------------------------------
        aerated_hr_per_day = i.react_hr * i.cycles_per_day

        # AOR = daily oxygen / aerated hours
        AOR_lbO2_per_hr = total_O2_lbday / max(aerated_hr_per_day, 1e-6)

        # -----------------------------------------
        # Transfer Efficiency Correction (Standard → Field)
        # SOR = AOR / (alpha * beta * (ΔC / C*sat) * θ^(T-20))
        # -----------------------------------------
        alpha = i.alpha
        beta = i.beta
        DO_target = i.target_DO_mgL
        C_sat = 9.1  # mg/L (approx at 20°C)
        delta_C = C_sat - DO_target
        theta_factor = i.theta_DO ** (i.temp_C - 20.0)

        transfer_factor = alpha * beta * (delta_C / C_sat) * theta_factor
        transfer_factor = max(transfer_factor, 0.05)  # prevent zero

        SOR_lbO2_per_hr = AOR_lbO2_per_hr / transfer_factor

        # -----------------------------------------
        # Total O2 / basin (use for blower sizing)
        # -----------------------------------------
        AOR_per_basin = AOR_lbO2_per_hr / i.n_basins
        SOR_per_basin = SOR_lbO2_per_hr / i.n_basins

        return {
            "bod_lbday": bod_lbday,
            "nh3_lbday": nh3_lbday,
            "O2_cBOD_lbday": O2_cBOD,
            "O2_nitrification_lbday": O2_N,
            "O2_endogenous_lbday": O2_endogenous,
            "O2_total_lbday": total_O2_lbday,
            "aerated_hr_per_day": aerated_hr_per_day,
            "AOR_lbO2_per_hr": AOR_lbO2_per_hr,
            "SOR_lbO2_per_hr": SOR_lbO2_per_hr,
            "AOR_per_basin_lbO2_per_hr": AOR_per_basin,
            "SOR_per_basin_lbO2_per_hr": SOR_per_basin,
            "transfer_factor": transfer_factor,
        }
    # ========================================================
    # INTERNAL METHOD 5 — EFFLUENT QUALITY & NITRIFICATION
    # ========================================================
    def _calc_effluent_quality(
        self,
        basin: dict[str, Any],
        biomass: dict[str, Any],
        oxygen: dict[str, Any]
    ) -> dict[str, Any]:
        """
        Predict effluent quality using:
          • Nitrifier kinetics (Metcalf & Eddy 5th Ed)
          • Critical SRT for nitrification
          • Temperature corrections
          • DO limitations
          • Simple soluble BOD polishing model
        """

        i = self.inputs

        # ============================================
        # ---- 1. NITRIFICATION KINETICS (M&E) -------
        # ============================================
        temp = i.temp_C
        S = i.nh3_mgL           # ammonia concentration entering SBR
        ks = i.ks_NH3
        mu_max_20 = i.mu_max_20C
        kd = i.kd_decay
        theta_n = i.theta_nitrifiers

        # Temperature-corrected maximum growth rate
        mu_max_T = mu_max_20 * (theta_n ** (temp - 20.0))

        # Actual growth rate based on ammonia availability (Monod)
        mu_actual = mu_max_T * (S / (ks + S))

        # Net specific growth rate
        mu_net = mu_actual - kd
        mu_net = max(mu_net, 1e-6)  # avoid zero or negative values

        # Critical SRT for nitrification (M&E)
        SRT_crit = 1.0 / mu_net

        # Actual SRT from Part 3
        SRT_actual = biomass["srt_actual_days"]

        # Nitrification fraction completed
        nitrification_fraction = min(SRT_actual / SRT_crit, 1.0)

        # ============================================
        # ---- 2. EFFLUENT NH3-N PREDICTION ----------
        # ============================================
        influent_NH3 = i.nh3_mgL

        # Remaining NH3 after nitrification
        eff_NH3_mgL = influent_NH3 * (1.0 - nitrification_fraction)

        # ============================================
        # ---- 3. WINTER TEMPERATURE CHECKS ----------
        # ============================================
        winter_temp = i.winter_temp_C
        mu_max_winter = mu_max_20 * (theta_n ** (winter_temp - 20.0))
        mu_net_winter = mu_max_winter * (S / (ks + S)) - kd
        mu_net_winter = max(mu_net_winter, 1e-6)

        SRT_crit_winter = 1.0 / mu_net_winter
        nitrification_fraction_winter = min(SRT_actual / SRT_crit_winter, 1.0)

        eff_NH3_winter_mgL = influent_NH3 * (1.0 - nitrification_fraction_winter)

        # ============================================
        # ---- 4. EFFLUENT BOD PREDICTION -------------
        # ============================================
        # Simple SBR soluble BOD polish model:
        # Remaining soluble BOD ≈ 3–8 mg/L for SRT 12–20 days
        # Scale by SRT performance:
        #
        # eff_BOD = 30 * exp(-0.18 * SRT_actual)
        #
        eff_BOD_mgL = 30.0 * math.exp(-0.18 * SRT_actual)
        eff_BOD_mgL = max(eff_BOD_mgL, 2.0)   # typical lower bound for SBR

        # ============================================
        # ---- 5. DO LIMITATION (NITRIFIER INHIB.) ---
        # ============================================
        DO_target = i.target_DO_mgL

        if DO_target < 1.5:
            do_limitation_flag = True
        else:
            do_limitation_flag = False

        # ============================================
        # ---- 6. DECANT NH3 SPIKE CHECK --------------
        # ============================================
        # If eff_NH3 very low (<1 mg/L), risk of spike is minimal.
        # Otherwise use a simple heuristic.
        #
        if eff_NH3_mgL > 2.0:
            decant_spike_risk = "Moderate Risk"
        elif eff_NH3_mgL > 5.0:
            decant_spike_risk = "High Risk"
        else:
            decant_spike_risk = "Low Risk"

        return {
            # Nitrification kinetics
            "mu_max_T": mu_max_T,
            "mu_actual": mu_actual,
            "mu_net": mu_net,
            "SRT_crit_days": SRT_crit,
            "SRT_actual_days": SRT_actual,
            "nitrification_fraction": nitrification_fraction,

            # Effluent ammonia
            "effluent_NH3_mgL": eff_NH3_mgL,
            "winter_effluent_NH3_mgL": eff_NH3_winter_mgL,
            "SRT_crit_winter_days": SRT_crit_winter,

            # Effluent BOD
            "effluent_BOD_mgL": eff_BOD_mgL,

            # DO limitation
            "DO_limitation": do_limitation_flag,

            # Decant ammonia spike risk
            "decant_spike_risk": decant_spike_risk,
        }
    # ========================================================
    # INTERNAL METHOD 6 — WARNINGS ENGINE
    # ========================================================
    def _collect_warnings(
        self,
        basin: dict[str, Any],
        biomass: dict[str, Any],
        oxygen: dict[str, Any],
        eff: dict[str, Any]
    ) -> list[str]:

        i = self.inputs
        warnings = []

        # ----------------------------------------------------
        # MLSS checks (Ten States)
        # ----------------------------------------------------
        if i.mlss_mgL > i.max_MLSS:
            warnings.append(
                f"MLSS of {i.mlss_mgL:,.0f} mg/L exceeds Ten States max of {i.max_MLSS:,.0f} mg/L."
            )

        if i.mlss_mgL < i.min_MLSS:
            warnings.append(
                f"MLSS of {i.mlss_mgL:,.0f} mg/L is below typical minimum of {i.min_MLSS:,.0f} mg/L."
            )

        # ----------------------------------------------------
        # Depth checks (settle + decant issues)
        # ----------------------------------------------------
        if basin["fill_depth_ft"] > (i.basin_depth_ft - i.freeboard_ft):
            warnings.append(
                "Fill depth exceeds available basin water depth (risk of overflow)."
            )

        # ----------------------------------------------------
        # Nitrification + SRT warnings
        # ----------------------------------------------------
        SRT_actual = biomass["srt_actual_days"]
        SRT_crit = eff["SRT_crit_days"]

        if SRT_actual < SRT_crit:
            warnings.append(
                f"SRT of {SRT_actual:.1f} days is below the critical SRT for nitrification ({SRT_crit:.1f} days)."
            )

        # Winter check
        SRT_crit_winter = eff["SRT_crit_winter_days"]
        if SRT_actual < SRT_crit_winter:
            warnings.append(
                f"SRT of {SRT_actual:.1f} days is insufficient for winter nitrification (requires {SRT_crit_winter:.1f} days)."
            )

        # ----------------------------------------------------
        # Effluent NH3 exceedance
        # ----------------------------------------------------
        if eff["effluent_NH3_mgL"] > 1.0:
            warnings.append(
                f"Predicted effluent ammonia is elevated at {eff['effluent_NH3_mgL']:.2f} mg/L."
            )

        if eff["winter_effluent_NH3_mgL"] > 2.0:
            warnings.append(
                f"Winter effluent ammonia predicted high at {eff['winter_effluent_NH3_mgL']:.2f} mg/L."
            )

        # ----------------------------------------------------
        # Effluent BOD warning
        # ----------------------------------------------------
        if eff["effluent_BOD_mgL"] > 10.0:
            warnings.append(
                f"Predicted effluent BOD is elevated at {eff['effluent_BOD_mgL']:.1f} mg/L (check aeration and SRT)."
            )

        # ----------------------------------------------------
        # DO limitation warnings
        # ----------------------------------------------------
        if eff["DO_limitation"]:
            warnings.append(
                "DO < 1.5 mg/L detected — nitrification inhibition expected."
            )

        # ----------------------------------------------------
        # Decant ammonia spike
        # ----------------------------------------------------
        if eff["decant_spike_risk"] == "High Risk":
            warnings.append(
                "High risk of ammonia spike during decant — consider longer react or step-feed."
            )

        return warnings


    # ========================================================
    # INTERNAL METHOD 7 — RECOMMENDATIONS ENGINE
    # ========================================================
    def _collect_recommendations(
        self,
        basin: dict[str, Any],
        biomass: dict[str, Any],
        oxygen: dict[str, Any],
        eff: dict[str, Any]
    ) -> list[str]:

        rec = []
        i = self.inputs

        # ----------------------------------------------------
        # Improve nitrification
        # ----------------------------------------------------
        if biomass["srt_actual_days"] < eff["SRT_crit_days"]:
            rec.append("Increase SRT by reducing WAS to improve nitrification reliability.")

        if eff["DO_limitation"]:
            rec.append("Increase aeration or blower output to maintain DO ≥ 2 mg/L.")

        # ----------------------------------------------------
        # Improve effluent BOD
        # ----------------------------------------------------
        if eff["effluent_BOD_mgL"] > 10:
            rec.append("Increase react time or aeration intensity to reduce effluent BOD.")

        # ----------------------------------------------------
        # Winter performance
        # ----------------------------------------------------
        if eff["winter_effluent_NH3_mgL"] > 2.0:
            rec.append("For winter operation, increase SRT, DO, or consider adding a swing zone.")

        # ----------------------------------------------------
        # Decant improvements
        # ----------------------------------------------------
        if eff["decant_spike_risk"] == "Moderate Risk":
            rec.append("Consider step-feed or intermittent react to minimize decant ammonia spike.")

        if eff["decant_spike_risk"] == "High Risk":
            rec.append("Increase overall nitrification capacity to reduce decant ammonia spikes.")

        # ----------------------------------------------------
        # General Optimization
        # ----------------------------------------------------
        rec.append("Verify cycle timing and adjust react period to match diurnal loading.")
        rec.append("Consider installing DO control to reduce energy costs while improving nitrification.")

        return rec
