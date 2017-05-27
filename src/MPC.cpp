#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

const size_t N = 6;
const double dt = 0.250;

// length from front to CoG
const double Lf = 2.67;

// reference errors (0) and speed (converts from mph to m/s)
const double ref_cte = 0;
const double ref_epsi = 0;
const double ref_v = 120 * 1609/3600;

// one vector with all variables
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    // cost
    fg[0] = 0;

    // state cost
    for (size_t i = 0; i < N; i++) {
      fg[0] += 5E1 * CppAD::pow(vars[cte_start + i] - ref_cte, 2);
      fg[0] += 1E4 * CppAD::pow(vars[epsi_start + i] - ref_epsi, 2);
      fg[0] += 1 * CppAD::pow(vars[v_start + i] - ref_v, 2);
    }

    // actuator cost
    for (size_t i = 0; i < N - 1; i++) {
      fg[0] += 1 * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += 1 * CppAD::pow(vars[a_start + i], 2);
    }

    // actuator change cost
    for (size_t i = 0; i < N - 2; i++) {
      fg[0] += 1E5* CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += 1E1 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    // cost is first element and others are pushed one position
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (size_t i = 0; i < N-1; i++) {
      // state at time t
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> cte0 = vars[cte_start + i];
      AD<double> epsi0 = vars[epsi_start + i];

      // state at time t+1
      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i + 1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];

      // actuators at time t
      AD<double> delta0 = vars[delta_start+ i];
      AD<double> a0 = vars[a_start + i];

      // desired location f and heading psi
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      AD<double> derivative0 = coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0;
      AD<double> psi_desired0 = CppAD::atan(derivative0);

      // deviation between desired point and predicted point
      fg[2 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi_start + i] = psi1 - (psi0 + v0/Lf * delta0 * dt);
      fg[2 + v_start + i] = v1 - (v0 + a0 * dt);
      fg[2 + cte_start + i] = cte1 - ((f0 - y0) + v0 * CppAD::sin(epsi0) * dt);
      fg[2 + epsi_start + i] = epsi1 - ((psi0 - psi_desired0) + v0/Lf * delta0 * dt);
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
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // state variables x, y, v, psi
  size_t n_state_vars = 6;
  // actuator variables delta, a
  size_t n_actuator_vars = 2;
  // number of independent variables
  size_t n_vars = n_state_vars * N + n_actuator_vars * (N-1);
  // number of constraints
  size_t n_constraints = n_state_vars * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  // initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // max pos and neg values for variables x, y, v, psi, cte, epsi
  for (size_t i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1e19;
    vars_upperbound[i] =  1e19;
  }

  // minimum/maximum steering angle -25/25° in radians
  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] =  0.436332;
  }

  // minimum/maximum acceleration mph/s
  for (size_t i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -3.0;
    vars_upperbound[i] =  3.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

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

  // return actuator values, and (x, y) positions
  vector <double> result = { solution.x[delta_start], solution.x[a_start] };
  for (size_t i = 0; i < N; i++) {
    result.push_back(solution.x[x_start + i]);
  }

  for (size_t i = 0; i < N; i++) {
    result.push_back(solution.x[y_start + i]);
  }

  return result;
}
