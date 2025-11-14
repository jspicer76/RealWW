# ============================================================
# TEST — Unified SBR Engine (Module 1 + Module 2)
# ============================================================

from sbr_engine import SBRInputs, SBRSystem
import json

def main():
    # --------------------------------------------------------
    # Sample design input (typical municipal SBR)
    # --------------------------------------------------------
    inputs = SBRInputs(
        flow_mgd=0.43,
        bod_mgL=250,
        tss_mgL=250,
        nh3_mgL=30,
        phosphorus_mgL=7,

        n_basins=2,
        basin_depth_ft=18,
        basin_length_ft=100,
        basin_width_ft=60,

        cycles_per_day=3.0,
        fill_hr=0.5,
        react_hr=4.0,
        settle_hr=1.0,
        decant_hr=1.5,
        idle_hr=0.0,

        mlss_mgL=3000,
        srt_target_days=15,
        decant_fraction_of_depth=0.30,

        target_DO_mgL=2.0,
        alpha=0.8,
        beta=1.0,
        theta_DO=1.024,
        temp_C=20.0,

        winter_temp_C=8.0,
    )

    # --------------------------------------------------------
    # Run full model
    # --------------------------------------------------------
    system = SBRSystem(inputs)
    results = system.run_full_design()

    # --------------------------------------------------------
    # Pretty-print nested results as JSON
    # --------------------------------------------------------
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
