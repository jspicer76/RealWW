from backend.oxidation_ditch.sizing import design_oxidation_ditch

# Example: Hodgenville design flow = 0.43 MGD
results = design_oxidation_ditch(flow_mgd=0.43)

print("\n=== REALWW OXIDATION DITCH DESIGN ===")
for group, data in results.items():
    print(f"\n--- {group.upper()} ---")
    for k, v in data.items():
        print(f"{k}: {v}")
