#include "stl_star/robustness_monitor_general.h"
#include "stl_star/boolean_monitor_general.h"
#include "robustness_common.h"

namespace stl_star {

RobustnessResult robustness_monitor_general(const Formula& f, const Trace& trace,
                                            double epsilon) {
    return robustness_detail::bisect(f, trace, epsilon, boolean_monitor_general_holds_at_start);
}

}  // namespace stl_star
