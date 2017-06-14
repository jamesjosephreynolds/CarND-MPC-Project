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
  const size_t N = 12;
  const double dt = 0.15;
 
 // Weight factors for cost optimization
  struct MPC_Weights {
    const double w_cte = 0.1;          // limit to 0.5 (cte max = 2 meter)
    const double w_epsi = 0.3;      // limit to 0.636 (epsi max = 1.57 rad (pi/2))
    const double w_v = 0.01;        // limit to 0.0179 (v err max = 56 m/s)
    const double w_angle = 1.0;     // limit to 2.29 (angle max = 0.436 rad)
    const double w_accel = 0.25;     // limit to 1.0 (accel max = 1.0 m/s/s)
    const double w_angle_jerk = 1.146;// limit to 1.146 (jerk max = 0.872 rad)
    const double w_accel_jerk = 0.5;// limit to 0.5 (jerk max = 2.0 m/s/s);
    const double w_norm = 1 / (w_cte + w_epsi + w_v + w_angle + w_accel + w_angle_jerk + w_accel_jerk);
  } w;
  
  const int vertex_idx_max = 8;
  vector<double> vertex_tot_cost;
  
  int vertex_idx = 0;
  int vertex_cnt = 0;
  const int vertex_cnt_max = 400;
  const int vertex_cnt_min = 100;
  
  bool restart_sim = false;
    
  MPC();

  virtual ~MPC();
  
  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
