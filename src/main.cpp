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
using Eigen::MatrixXd;
using Eigen::VectorXd;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

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

  // initialize MPC
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    // cout << sdata << endl;
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
          double cte = 0;
          double epsi = 0;

          // number of waypoints
          size_t n_points = ptsx.size();
          assert(ptsx.size() == ptsy.size());

          // waypoints in world map coordinates
          VectorXd ptsx_world = VectorXd::Map(ptsx.data(), ptsx.size());
          VectorXd ptsy_world = VectorXd::Map(ptsy.data(), ptsy.size());

          // translate to vehicle position (px, py)
          MatrixXd waypoints = MatrixXd(2, n_points);
          waypoints.row(0) = ptsx_world - px * VectorXd::Ones(n_points);
          waypoints.row(1) = ptsy_world - py * VectorXd::Ones(n_points);

          // rotate coordinates about vehicle with angle psi
          MatrixXd rotation_matrix = Eigen::MatrixXd(2, 2);
          rotation_matrix << cos(psi), sin(psi), -sin(psi), cos(psi);
          MatrixXd waypoints_vehicle = rotation_matrix * waypoints;

          // waypoints in vehicle coordinates
          VectorXd ptsx_vehicle = waypoints_vehicle.row(0);
          VectorXd ptsy_vehicle = waypoints_vehicle.row(1);

          // third order polynomial path through waypoints
          VectorXd coeffs = polyfit(ptsx_vehicle, ptsy_vehicle, 3);

          // define state in vehicle coordinates
          Eigen::VectorXd state = Eigen::VectorXd(6);
          px = 0;
          py = 0;
          psi = 0;
          v = v * 1609 / 3600; // meters per second instead of mph
          cte = polyeval(coeffs, px) - py;

          // calculate epsi
          double derivative = 0;
          for (size_t k = 1; k < coeffs.size(); k++) {
            derivative += coeffs[k] * k * pow(px, k-1);
          }
          double psi_desired = atan(derivative);
          epsi = psi - psi_desired;

          size_t latency_mil = 100;
          state << px + v * latency_mil/1000, py, psi, v, cte, epsi;

          vector<double> solution = mpc.Solve(state, coeffs);

          // flip steering sign since left steering should be +ve instead of -ve (simulator convention)
          // denormalize so that 25° (0.436 rad) results in 1.0 actuator value
          double delta = solution[0];
          double steer_value = -delta;
          steer_value /= deg2rad(25);

          // map acceleration to throttle
          double a = solution[1];
          double throttle_value = a;

          // limit steering and throttle absolute value to 1
          if (steer_value < -1) steer_value = -1;
          if (steer_value > 1) steer_value = 1;
          if (throttle_value < -1) throttle_value = -1;
          if (throttle_value > 1) throttle_value = 1;

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          // display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          // 1 delta, 1 a, N x, N y
          size_t N = (solution.size() - 2)/2;
          for (size_t k = 0; k < N; k++) {
            mpc_x_vals.push_back(solution[2 + k]);
            mpc_y_vals.push_back(solution[2 + k + N]);;
          }

          // mpc path points in the simulator are connected by a Green line
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          // display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (size_t k = 0; k < n_points; k++) {
            next_x_vals.push_back(ptsx_vehicle[k]);
            next_y_vals.push_back(ptsy_vehicle[k]);
          }

          // the waypoints in the simulator are connected by a Yellow line
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;

          // latency in milliseconds
          this_thread::sleep_for(chrono::milliseconds(latency_mil));
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
