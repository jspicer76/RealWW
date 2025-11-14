# backend/sbr/test_sbr.py

from __future__ import annotations

from backend.sbr.design import design_sbr_autotune_full_cycle


from pprint import pprint

from backend.sbr.design import design_sbr


def print_header(title: str) -> None:
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)


def main() -> None:
    # Example design run
    results = design_sbr(flow_mgd=0.43)

    print_header("REALWW — SBR DESIGN SUMMARY")

    print_header("1. INPUT PARAMETERS")
    pprint(results.get("inputs", {}), sort_dicts=False)

    print_header("2. INFLUENT LOADS (lb/day)")
    pprint(results.get("loads", {}), sort_dicts=False)

    print_header("3. SBR CYCLE BREAKDOWN (hr)")
    pprint(results.get("cycle", {}), sort_dicts=False)

    print_header("4. BASIN SIZING")
    pprint(results.get("basin", {}), sort_dicts=False)

    print_header("5. DECANTER SIZING")
    pprint(results.get("decanter", {}), sort_dicts=False)

    print_header("6. OXYGEN / AERATION REQUIREMENTS")
    pprint(results.get("oxygen", {}), sort_dicts=False)

    print_header("7. TEN STATES COMPLIANCE CHECK + ADVISORIES")
    ten_states = results.get("ten_states", {})
    pprint(ten_states, sort_dicts=False)

    # Print advisories more nicely (if present)
    advisories = ten_states.get("advisories", [])
    if advisories:
        print("\nAdvisories:")
        for idx, msg in enumerate(advisories, start=1):
            print(f"  {idx}. {msg}")

    print("\nDone.\n")


if __name__ == "__main__":
    main()

from backend.sbr.design import design_sbr_autotune

results = design_sbr_autotune(flow_mgd=0.43)

print(results["changes_made"])
pprint(results["autotuned_design"]["ten_states"])

from backend.sbr.design import design_sbr_autotune_full_cycle

results = design_sbr_autotune_full_cycle(flow_mgd=0.43)

print("\nCHANGES MADE:")
for c in results["changes_made"]:
    print(" -", c)

print("\nFINAL TEN STATES CHECK:")
from pprint import pprint
pprint(results["autotuned_design"]["ten_states"])

from backend.sbr.design import design_sbr_autotune_mlss


bio_results = design_sbr_autotune_mlss(
    flow_mgd=0.43,
    target_f_m=0.15,
)

print_header("8. BIOLOGY AUTO-TUNE (F:M TARGET)")
print(bio_results["notes"])
print("Last iteration:")
from pprint import pprint
pprint(bio_results["iterations"][-1], sort_dicts=False)
print("\nFinal biology block:")
pprint(bio_results["final_design"]["biology"], sort_dicts=False)


from backend.sbr.design import design_sbr_autotune_mlss


bio_mlss = design_sbr_autotune_mlss(flow_mgd=0.43, target_f_m=0.15)

print("\nMLSS AUTO-TUNE RESULT:")
print(bio_mlss["notes"])
print("Final MLSS:", bio_mlss["mlss_final"])
print("Final F:M:", bio_mlss["f_m_final"])

from backend.sbr.design import design_sbr_autotune_srt

print_header("9. BIOLOGY AUTO-TUNE (SRT TARGET)")
srt_results = design_sbr_autotune_srt(flow_mgd=0.43, target_srt_days=12)

print(srt_results["notes"])
print("Final SRT (days):", srt_results["srt_final"])
print("Final WAS flow (MGD):", srt_results["was_flow_mgd"])
print("Last iteration:")
pprint(srt_results["iterations"][-1], sort_dicts=False)
