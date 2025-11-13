from core.constants import (
    DEFAULT_BOD,
    O2_REQUIRED_PER_LB_BOD,
)

def daily_bod_load(flow_mgd: float, bod_mgL: float = DEFAULT_BOD) -> float:
    """lb/day of BOD at given flow"""
    return flow_mgd * bod_mgL * 8.34  # mg/L → lb/day

def daily_o2_required_bod(flow_mgd, bod_mgL=DEFAULT_BOD):
    bod_lb_day = daily_bod_load(flow_mgd, bod_mgL)
    return bod_lb_day * O2_REQUIRED_PER_LB_BOD

if __name__ == "__main__":
    flow = 0.43
    print("=== RealWW Test ===")
    print("Flow:", flow, "MGD")
    print("BOD load:", daily_bod_load(flow), "lb/day")
    print("O2 required:", daily_o2_required_bod(flow), "lb/day")

