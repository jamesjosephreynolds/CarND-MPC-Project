# Model-Predictive Controls for Steering Angle and Vehicle Speed Project

This project is originally forked from https://github.com/udacity/CarND-MPC-Project. This repository includes starter code, that is used herein.

## Controller Structure

For my controller, I implemented the vehicle kinematic model.  I chose not to implement the dynamic model because, for the purposes of running in a simulator, the additional effort would be very high, relative the realizable performance benefit.  Moreover, the dynamic model would require more parameter tuning, and the dimensionality is already high for the kinematic model.

The model update equations are implemented in `MPC.cpp` as follows:

```C++
// Kinematic model update equations rewritten as equality constraints
fg[x_st + i + 1 + 1]    = x1 - (x0 + v0*CppAD::cos(psi0)*m.dt);           // x(t+dt) = x(t) + v(t)*cos(psi(t))*dt
fg[y_st + i + 1 + 1]    = y1 - (y0 + v0*CppAD::sin(psi0)*m.dt);           // y(t+dt) = y(t) + v(t)*sin(psi(t))*dt
fg[psi_st + i + 1 + 1]  = psi1 - (psi0 + v0*angle0*m.dt/Lf);              // psi(t+dt) = psi(t) + v(t)*delta(t)*dt/Lf
fg[v_st + i + 1 + 1]    = v1 - (v0 + accel0*m.dt);                        // v(t+dt) = v(t) + a(t)*dt
fg[cte_st + i + 1 + 1]  = cte1 - ((f0 - y0) + v0*CppAD::sin(epsi0)*m.dt); // cte(t+dt) = cte(t) + v(t)*sin(epsi(t))*dt
fg[epsi_st + i + 1 + 1] = epsi1 - ((psi0-psi_ref) + v0*angle0*m.dt/Lf);   // epsi(t+dt) = epsi(t) + v(t)*delta(t)*dt
```

## Latency

I modeled the actuator latency into my initial value for the state vector.  In order to accomplish this, I took the current state values from the simulator, and projected them forward one latency time step.

```C++
double dt = 0.1; // actuator latency
          
// Put latency into initial state values
px = v*dt;
psi = -v*steer_angle*dt/2.67;
```

### Polynomial Fit with Latent Value

I chose a second-order polynomial for fitting waypoints.  Intuitively, any higher order polynomials seems unreasonable for a very short (1 to 2 second) time horizon.

```C++
// Fit a polynomial to upcoming waypoints          
Eigen::VectorXd coeffs = polyfit(X_w_raw, Y_w_raw, int(2));
```

I further used this polynomial and the latent `x` position, from above, to set the initial value `t=0.0` for the cross-track error (CTE) and heading error.
```C++
// cross track error is distance in y, from the vehicle coordinate systems's perspective
double cte = polyeval(coeffs, px);
          
// epsi is the difference between desired heading and actual
double epsi = atan(coeffs[1]+2*coeffs[2]*px);
```

## Optimality

Optimal control general uses a weighting notion, so that different cost contributors are considered more or less important to the "optimal" solution.  This was not mentioned in the lecture material, but it is the method I've chosen.  Each of the cost terms has its own unique weight.  With 7 weight factors, and the number of model steps and timestep, there are a total of 9 model hyperparameters to tune.

```C++
// Number of future points to model, and time step
const size_t N = 8;
const double dt = 0.10;
 
// Weight factors for cost optimization
struct MPC_Weights {
  const double w_cte = 100;
  const double w_epsi = 100;
  const double w_v = 0.1;
  const double w_angle = 100;
  const double w_accel = 1;
  const double w_angle_jerk = 50000;
  const double w_accel_jerk = 10000;
  const double w_norm = 1 / (w_cte + w_epsi + w_v + w_angle + w_accel + w_angle_jerk + w_accel_jerk);
} w;
```

One nice benefit of using weighting factors, along with a normalization term, is that the model tuning becomes intuitive.  If one aspect of the model performance is undesirable (e.g. jerky steering), then just increase the weight of that cost term `w_angle_jerk`.

In order to tune `N` and `dt`, I took the following approach.  I started with `dt` equal to the system latency, and `N` equal to `20`, for a 2 second time horizon.  I tried to decrease `dt` first, and the control became very oscillatory, so I did not change it.  Then, I decreased `N` until, at high-speed, the resulting prediction had roughly the same horizon as the next waypoints.

Most of my actual model tuning was in the weighting factors.  You can see that `w_v` is very small, because I preferred safe steering over achieving the target speed.  In this way, I could set the target speed to 99 MPH, and just let the controller regulate the safest speed.

## Results

I was able to achieve a very pleasing result, with top speed over 80 MPH, and smooth handling.  The control is stable, and the top speed of 80 MPH is more than 2.5 times higher than I achieved with PID control.

## Reflections
