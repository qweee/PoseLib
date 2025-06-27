#ifndef POSELIB_HC_HC_PROBLEM_BASE_H_
#define POSELIB_HC_HC_PROBLEM_BASE_H_

#include "PoseLib/types.h"

namespace poselib {

template <typename Model = Image>
class HCProblemBase {

    public:
    HCProblemBase() {}
    virtual void simulator(const Model &model) = 0;
    virtual void compute_polys(const Model &model) = 0;
    virtual void compute_jacobian(const Model &model) = 0;
    virtual Model step(const Model &model) = 0;

    size_t num_params;
    size_t num_polys;

    typedef Model param_t;

};

} // namespace poselib

#endif