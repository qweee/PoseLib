#ifndef POSELIB_HOMOTOPY_CONTINUATION_H_
#define POSELIB_HOMOTOPY_CONTINUATION_H_


#include "PoseLib/types.h"
#include "PoseLib/HCsolvers/HCproblems/helper.h"

#include <vector>
#include <iostream>

namespace poselib {

// Notes: Problem (HC problems) is for polynomial functions and Model (Image struct) is for the parameters to be optimized
// problem.x is 2D points and problem.x_simulated is 2D points simulated by the initial solution

template<typename Problem, typename Solution, typename Poly, typename Jacobian>
bool compute_dx(double t, Problem &problem, Solution &sol, Solution &dx) {

    Poly polysF, polysG;
    Jacobian JF, JG;

    problem.compute_PolysandJacobian(sol, problem.x, polysF, JF);
    problem.compute_PolysandJacobian(sol, problem.x_simulated, polysG, JG);

    Jacobian Jx = t * JF + (1 - t) * JG;
    Poly t_grad = polysF - polysG;

    dx = -(Jx.transpose() * Jx).inverse() * Jx.transpose() * t_grad;

    return true;
}


template<typename Problem, typename Solution, typename Poly, typename Jacobian>
HCStats HC_impl(Problem &problem, const HCOptions &opt, Solution &sol)
{

    HCStats stats;

	double t = 0.0;
	
    sol = problem.get_sol_vector();

	for (stats.iterations = 0; stats.iterations < opt.max_iterations+1; ++stats.iterations)
	{   
		// Prediction(Kutta Runge or Euler's method)
		// t = i * step_size;
		if (opt.forth_predictor) {
            
			Solution k1;
            bool k1_flag = compute_dx<Problem, Solution, Poly, Jacobian>(t, problem, sol, k1);
            if (!k1_flag) {
                // std::cout << "Kutta k1 is not full rank" << std::endl;
                stats.success = false;
                break;
            }

            Solution sol_temp1 = sol + opt.step_size/2 * k1;
			Solution k2;
            bool k2_flag = compute_dx<Problem, Solution, Poly, Jacobian>(t + opt.step_size/2, problem, sol_temp1, k2);
            if (!k2_flag) {
                std::cout << "Kutta k2 is not full rank" << std::endl;
                stats.success = false;
                break;
            }

            Solution sol_temp2 = sol + opt.step_size/2 * k2;
			Solution k3;
            bool k3_flag = compute_dx<Problem, Solution, Poly, Jacobian>(t + opt.step_size/2, problem, sol_temp2, k3);
            if (!k3_flag) {
                // std::cout << "Kutta k3 is not full rank" << std::endl;
                stats.success = false;
                break;
            }

            Solution sol_temp3 = sol + opt.step_size * k3;
			Solution k4;
            bool k4_flag = compute_dx<Problem, Solution, Poly, Jacobian>(t + opt.step_size, problem, sol_temp3, k4);
            if (!k4_flag) {
                // std::cout << "Kutta k4 is not full rank" << std::endl;
                stats.success = false;
                break;
            }

			sol += (k1 + 2*k2 + 2*k3 + k4)* opt.step_size/6;

		} else {
            Solution k1;
            bool k1_flag = compute_dx<Problem, Solution, Poly, Jacobian>(t, problem, sol, k1);
            if (!k1_flag) {
                // std::cout << "Euler k1 is not full rank" << std::endl;
                stats.success = false;
                break;
            }
			sol += k1 * opt.step_size;
		}

        t += opt.step_size;


		// Correction(Newton method)
		for (int iter = 0; iter < opt.newton_iter; ++iter)
		{
            Poly Hpolys;
            Jacobian JH;
            problem.compute_HpolysandJacobian(t, sol, Hpolys, JH);

            // Only print debugging info if needed
            if (opt.debug_output) {
                std::cout << "Hpolys: " << Hpolys << std::endl;
            }

            if (!checkValid(JH)) {
                // std::cout << "Newton JH is not full rank" << std::endl;
                stats.success = false;
                break;
            }

			if (opt.adaptive_flag) {
				Solution sol_temp = sol;
				sol = sol_temp - (JH.transpose() * JH).inverse() * JH.transpose() * Hpolys;
                
				if ( (sol - sol_temp).norm() < 1e-6) {
					break;
				}
			} else {
				sol = sol - (JH.transpose() * JH).inverse() * JH.transpose() * Hpolys;
			}
		}

        bool containsNaN = !(sol.array() == sol.array()).all(); // NaN check
        bool containsInf = !sol.array().isFinite().all(); // Inf check
        bool containsLarge = !(sol.array().abs() < 1e3).all(); // large check

        if (containsNaN || containsInf || containsLarge) {
            // break;
            stats.success = false;
            break;
        }

    }

    return stats;

};


} // namespace poselib

#endif