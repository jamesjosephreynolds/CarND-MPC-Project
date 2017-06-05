#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration


MPC m;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
const double v_tgt = 99; // set target speed to 35

// Define the starting indices for different variables
// Variable array - {x, y, psi, v, cte, epsi, angle, accel}
size_t x_st = 0;
size_t y_st = x_st + m.N;
size_t psi_st = y_st + m.N;
size_t v_st = psi_st + m.N;
size_t cte_st = v_st + m.N;
size_t epsi_st = cte_st + m.N;
size_t angle_st = epsi_st + m.N;
size_t accel_st = angle_st + m.N - 1; // use N-1, there is no "initial" angle

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // fg a vector of constraints and cost, vars is a vector of variables.
    
    // Initialize cost
    fg[0] = 0.0;
    
    // Cost associated with the desired state and speed
    // Exclude i = 0, as the current state cannot be changed by control updates
    for (int i = 1; i < m.N; ++i) {
      AD<double> cte = vars[cte_st + i];
      AD<double> epsi = vars[epsi_st + i];
      AD<double> v_err = vars[v_st + i] - v_tgt;
      fg[0] += m.w.w_cte*(cte*cte) + m.w.w_epsi*(epsi*epsi) + m.w.w_v*(v_err*v_err);
    }
    
    // Cost associated with the magnitude of actuators
    for (int i = 0; i < m.N - 1; ++i) {
      AD<double> angle = vars[angle_st + i];
      AD<double> accel = vars[accel_st + i];
      fg[0] += m.w.w_angle*(angle*angle) + m.w.w_accel*(accel*accel);
    }
    
    // Cost associated with harsh actuations (step change)
    for (int i = 0; i < m.N - 2; ++i) {
      AD<double> angle_dt = vars[angle_st + i + 1] - vars[angle_st + i];
      AD<double> accel_dt = vars[accel_st + i + 1] - vars[accel_st + i];
      fg[0] += m.w.w_angle_jerk*(angle_dt*angle_dt) + m.w.w_accel_jerk*(accel_dt*accel_dt);
    }
    
    // Normalize
    fg[0] *= m.w.w_norm;
    
    // Initial constraints are just current state
    fg[x_st + 1] = vars[x_st];
    fg[y_st + 1] = vars[y_st];
    fg[psi_st + 1] = vars[psi_st];
    fg[v_st + 1] = vars[v_st];
    fg[cte_st + 1] = vars[cte_st];
    fg[epsi_st + 1] = vars[epsi_st];
    
    /*std::cout << "state[k] = {" << fg[x_st + 1] << ", ";
    std::cout << fg[y_st + 1] << ", ";
    std::cout << fg[psi_st + 1] << ", ";
    std::cout << fg[v_st + 1] << ", ";
    std::cout << fg[cte_st + 1] << ", ";
    std::cout << fg[epsi_st + 1] << "}" << std::endl;*/
    
    // Define local variables for kinematic model equation readability
    for (int i = 0; i < m.N - 1; ++i) {
      AD<double> x0 = vars[x_st + i];
      AD<double> x1 = vars[x_st + i + 1];
      AD<double> y0 = vars[y_st + i];
      AD<double> y1 = vars[y_st + i + 1];
      AD<double> psi0 = vars[psi_st + i];
      AD<double> psi1 = vars[psi_st + i + 1];
      AD<double> v0 = vars[v_st + i];
      AD<double> v1 = vars[v_st + i + 1];
      AD<double> cte0 = vars[cte_st + i];
      AD<double> cte1 = vars[cte_st + i + 1];
      AD<double> epsi0 = vars[epsi_st + i];
      AD<double> epsi1 = vars[epsi_st + i + 1];
      
      AD<double> angle0 = vars[angle_st + i];
      AD<double> accel0 = vars[accel_st + i];
      
      AD<double> f0 = coeffs[0] + coeffs[1]*x0 + coeffs[2]*x0*x0;// + coeffs[3]*x0*x0*x0;
      AD<double> psi_ref = CppAD::atan(coeffs[1]);//+2*coeffs[2]*x0);
      
      // Kinematic model update equations rewritten as equality constraints
      fg[x_st + i + 1 + 1] = x1 - (x0 + v0*CppAD::cos(psi0)*m.dt); // x(t+dt) = x(t) + v(t)*cos(psi(t))*dt
      fg[y_st + i + 1 + 1] = y1 - (y0 + v0*CppAD::sin(psi0)*m.dt); // y(t+dt) = y(t) + v(t)*sin(psi(t))*dt
      fg[psi_st + i + 1 + 1] = psi1 - (psi0 + v0*angle0*m.dt/Lf); // psi(t+dt) = psi(t) + v(t)*delta(t)*dt/Lf
      fg[v_st + i + 1 + 1] = v1 - (v0 + accel0*m.dt); // v(t+dt) = v(t) + a(t)*dt
      fg[cte_st + i + 1 + 1] = cte1 - ((f0 - y0) + v0*CppAD::sin(epsi0)*m.dt); // cte(t+dt) = cte(t) + v(t)*sin(epsi(t))*dt
      fg[epsi_st + i + 1 + 1] = epsi1 - ((psi0-psi_ref) + v0*angle0*m.dt/Lf); // epsi(t+dt) = epsi(t) + v(t)*delta(t)*dt
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  //size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x, y, psi, v, cte, epsi;
  x = state[0];
  y = state[1];
  psi = state[2];
  v = state[3];
  cte = state[4];
  epsi = state[5];
  
  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  
  int n_st = 6; // number of states: {x, y, psi, v, cte, eps}
  int n_ac = 2; // number of actuators: {angle, accel}
  size_t n_vars = n_st*N + n_ac*(N-1);
  
  // Set the number of constraints
  // state dimension * number of steps
  size_t n_constraints = n_st*N;

  // Initial value of the independent variables.
  // First pass set all to 0.0, then set initial values from input array
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  
  vars[x_st] = x;
  vars[y_st] = y;
  vars[psi_st] = psi;
  vars[v_st] = v;
  vars[cte_st] = cte;
  vars[epsi_st] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  
  // Don't limit state values
  for (int i = 0; i < angle_st; ++i) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  
  // Limit actuator values
  for (int i = angle_st; i < accel_st; ++i) {
    vars_lowerbound[i] = -0.436;
    vars_upperbound[i] = 0.436;
  }
  
  for (int i = accel_st; i < n_vars; ++i) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  constraints_lowerbound[x_st] = x;
  constraints_lowerbound[y_st] = y;
  constraints_lowerbound[psi_st] = psi;
  constraints_lowerbound[v_st] = v;
  constraints_lowerbound[cte_st] = cte;
  constraints_lowerbound[epsi_st] = epsi;
  
  constraints_upperbound[x_st] = x;
  constraints_upperbound[y_st] = y;
  constraints_upperbound[psi_st] = psi;
  constraints_upperbound[v_st] = v;
  constraints_upperbound[cte_st] = cte;
  constraints_upperbound[epsi_st] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> output;
  
  
  /*for (int i = 0; i < n_vars; ++i) {
    std::cout << solution.x[i] << ", ";;
  }
  std::cout << std::endl;*/
  /*
  for (int i = 0; i < n_vars; ++i) {
    output.push_back(vars[i]);
  }*/
  
  output.push_back(solution.x[angle_st]);
  output.push_back(solution.x[accel_st]);
  
  for (int i = 0; i < N; ++i) {
    output.push_back(solution.x[x_st + i]);
    output.push_back(solution.x[y_st + i]);
  }
  
  return output;
}
