# backend/sbr/test_sbr.py

from __future__ import annotations
from backend.sbr.design import design_sbr_autotune
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
