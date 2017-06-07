#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <cppad/cppad.hpp>

using namespace std;

class MPC {
 public:
  
 /* 
  MPC project hyperparameters
  */
 // Number of future points to model, and time step
  const size_t N = 8;
  const double dt = 0.10;
 
 // Weight factors for cost optimization
  struct MPC_Weights {
    const double w_cte = 100;
    const double w_epsi = 100;
    const double w_v = 0.17;
    const double w_angle = 100;
    const double w_accel = 1;
    const double w_angle_jerk = 50000;
    const double w_accel_jerk = 10000;
    const double w_norm = 1 / (w_cte + w_epsi + w_v + w_angle + w_accel + w_angle_jerk + w_accel_jerk);
  } w;
    
  MPC();

  virtual ~MPC();
  
  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
