import math

class SBRModule4Circular:
    """
    ------------------------------------------------------------
    Module 4 — Circular SBR Geometry (Manual Uplift + Soil Berm)
    ------------------------------------------------------------
    Computes:
      • Basin diameter and geometry
      • Decant depth & volume
      • Concrete quantities
      • Uplift check including:
            - concrete wall weight
            - concrete slab weight
            - concrete ring beam weight
            - soil weight above ring overhang  <<< NEW
      • User manually adjusts slab thickness or overhang as needed
    ------------------------------------------------------------
    """

    def __init__(self, inputs):
        self.inputs = inputs

        # Physical constants
        self.gamma_water = 62.4        # pcf
        self.gamma_conc = 150          # pcf concrete
        self.gamma_soil = inputs.get("soil_unit_weight_pcf", 120)  # pcf soil
        self.FS_required = inputs.get("FS_uplift_required", 1.5)

        # Defaults (user adjustable)
        self.default_ring_ft = inputs.get("ring_extension_ft", 1.5)    # 18" default overhang
        self.default_slab_ft = inputs.get("slab_thickness_ft", 2.0)    # 24" slab
        self.default_soil_cover_ft = inputs.get("soil_cover_depth_ft", 2.0)  # 2 ft soil over ring

    # ------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------
    def compute_diameter(self, V_ft3, depth_ft):
        return math.sqrt(4 * V_ft3 / (math.pi * depth_ft))

    # ------------------------------------------------------------
    # Decant
    # ------------------------------------------------------------
    def compute_decant(self, D, decant_depth_ft):
        area = math.pi * (D / 2) ** 2
        vol_ft3 = area * decant_depth_ft
        return {
            "area_ft2": area,
            "volume_ft3": vol_ft3,
            "volume_gal": vol_ft3 * 7.48052
        }

    # ------------------------------------------------------------
    # Concrete Quantities
    # ------------------------------------------------------------
    def concrete_quantities(self, D, wall_height_ft, t_wall_ft, t_slab_ft):
        r = D / 2

        # Wall
        wall_vol_ft3 = math.pi * D * wall_height_ft * t_wall_ft

        # Slab
        slab_vol_ft3 = math.pi * r ** 2 * t_slab_ft

        return {
            "wall_yd3": wall_vol_ft3 / 27,
            "slab_yd3": slab_vol_ft3 / 27,
            "total_yd3": (wall_vol_ft3 + slab_vol_ft3) / 27
        }

    # ------------------------------------------------------------
    # Uplift Check — NOW INCLUDES SOIL OVER RING
    # ------------------------------------------------------------
    def uplift_check(self, D, wall_height_ft, t_wall_ft, t_slab_ft, ring_ext_ft, soil_cover_ft):
        r = D / 2

        # -------------------------
        # Concrete weights
        # -------------------------
        # Wall
        wall_vol = math.pi * D * wall_height_ft * t_wall_ft
        W_wall = wall_vol * self.gamma_conc

        # Slab
        slab_vol = math.pi * r**2 * t_slab_ft
        W_slab = slab_vol * self.gamma_conc

        # Ring beam (concrete annulus)
        outer_area = math.pi * (r + ring_ext_ft)**2
        inner_area = math.pi * r**2
        ring_area = outer_area - inner_area
        ring_vol = ring_area * t_slab_ft
        W_ring = ring_vol * self.gamma_conc

        # -------------------------
        # Soil Above Ring (NEW)
        # -------------------------
        soil_vol = ring_area * soil_cover_ft
        W_soil = soil_vol * self.gamma_soil

        # -------------------------
        # Uplift Force
        # -------------------------
        uplift_force = math.pi * r**2 * wall_height_ft * self.gamma_water

        # Total resisting load
        resisting = W_wall + W_slab + W_ring + W_soil

        FS = resisting / uplift_force

        return {
            "uplift_lb": uplift_force,
            "dead_weight_lb": resisting,
            "factor_of_safety": FS,
            "ok": FS >= self.FS_required,

            # Components
            "W_wall_lb": W_wall,
            "W_slab_lb": W_slab,
            "W_ring_lb": W_ring,
            "W_soil_lb": W_soil,

            "ring_extension_ft": ring_ext_ft,
            "soil_cover_depth_ft": soil_cover_ft,
            "slab_thickness_ft": t_slab_ft
        }

    # ------------------------------------------------------------
    # MAIN RUN
    # ------------------------------------------------------------
    def run(self):
        inp = self.inputs

        # Geometry
        V_ft3 = inp["required_volume_ft3"]
        depth = inp["working_depth_ft"]
        freeboard = inp.get("freeboard_ft", 2.0)
        H_total = depth + freeboard

        D = self.compute_diameter(V_ft3, depth)

        # Decant
        decant_fraction = inp.get("decant_fraction", 0.30)
        decant_depth_ft = depth * decant_fraction
        decant = self.compute_decant(D, decant_depth_ft)

        # Wall thickness
        t_wall = inp.get("wall_thickness_ft", 0.75)

        # Defaults for this mode
        t_slab = self.default_slab_ft
        ring_ext = self.default_ring_ft
        soil_cover = self.default_soil_cover_ft

        # Uplift Eval
        uplift = self.uplift_check(D, H_total, t_wall, t_slab, ring_ext, soil_cover)

        # Concrete
        concrete = self.concrete_quantities(D, H_total, t_wall, t_slab)

        # Final
        return {
            "geometry": {
                "diameter_ft": D,
                "working_depth_ft": depth,
                "total_wall_height_ft": H_total
            },
            "decant": decant,
            "concrete": concrete,
            "uplift_check": uplift,
            "note": (
                "FS below required; adjust slab thickness, ring extension, or soil cover."
                if not uplift["ok"]
                else "Uplift OK"
            )
        }
