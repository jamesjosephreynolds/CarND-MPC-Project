#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double steer_angle = j[1]["steering_angle"];
          
          double dt = 0.1; // actuator latency
          size_t N = mpc.N;

          // Convert waypoints from map coordinate system to car coordinate system
          int N_way = ptsx.size() <= ptsy.size() ? ptsx.size() : ptsy.size();
          Eigen::VectorXd X_w_raw(N_way); // X waypoints raw
          Eigen::VectorXd Y_w_raw(N_way); // Y waypoints raw
          
          // Convert waypoints into vehicle coordinate system
          for (int i = 0; i < N_way; ++i) {
            double x_pc; // x-coordinate of waypoint in car coordinate system
            double y_pc; // y-coordinate of waypoint in car coordinate system
            
            X_w_raw[i] = (ptsx[i] - px)*cos(psi) + (ptsy[i] - py)*sin(psi);
            Y_w_raw[i] = (ptsy[i] - py)*cos(psi) - (ptsx[i] - px)*sin(psi);
          }
          
          // Fit a polynomial to upcoming waypoints          
          Eigen::VectorXd coeffs = polyfit(X_w_raw, Y_w_raw, int(2));
          
          // Put latency into initial state values
          px = v*dt;
          psi = -v*steer_angle*dt/2.67;
          
          // cross track error is distance in y, from the vehicle coordinate systems's perspective
          double cte = polyeval(coeffs, px);
          
          // epsi is the difference between desired heading and actual
          double epsi = atan(coeffs[1]+2*coeffs[2]*px);
                   
          // Initialize state vector
          Eigen::VectorXd state(6); // {x, y, psi, v, cte, epsi}
          state << px, 0.0, psi, v, cte, epsi;

          // Solve the MPC
          vector<double> output = mpc.Solve(state, coeffs);
          
          // Steering convention is opposite of returned value, negate
          double steer_value = -output[0];
          double throttle_value = output[1];

          //Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          for (int i = 0; i < N; ++i ) {
            mpc_x_vals.push_back(output[2*i + 2]);
            mpc_y_vals.push_back(output[2*i + 3]);
          }

          //Display the waypoints as received (in car coordinate system)
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          for (int i = 0; i < N_way; ++i ) {
            next_x_vals.push_back(X_w_raw[i]);
            next_y_vals.push_back(Y_w_raw[i]);
          }
          
          // Send actuator commands and points to plot back to simulator
          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";

          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
