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

### Polynomial Fit

I chose a second-order polynomial for fitting waypoints.  Intuitively, any higher order polynomials seems unreasonable for a very short (1 to 2 second) time horizon.



## Latency

## Optimality

## Results

## Reflections
