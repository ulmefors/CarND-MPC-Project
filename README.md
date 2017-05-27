# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

---

## Model

### State
The state is defined as [x, y, &psi;, v, cte, e&psi;] where (x, y) is the vehicle position, &psi; is the orientation, v is speed (scalar), cte is the cross track error, and e&psi; is the orientation error.
The actuators are [&delta;, a] and correspond to the steering angle and acceleration.

Steering angle within the MPC is positive when driving to the left (positive y direction) whereas the opposite is true in the simulator. In addition, the optimization calculations use values for steering angle measured in radians whereas the range for the simulator is [-1, 1] (same as brake/accelerator).
The corresponding conversions are performed before sending values from the MPC to the simulator.

The Ipopt solver takes the state vector and upper/lower bounds for all values.
Positions and angles can take very small or large values whereas steering angle and acceleration are more limited (&#177;25° and &#177;3 m/s&#178; respectively).
Max CPU time is set at 50 ms.

### Cost function
The cost function requires careful tuning in order to achieve good control performance. The cte cost forces the vehicle closer to the waypoints but too large values will cause overshoot, particularly when driving with latency.
The e&psi; cost forces the vehicle to steer in parallel with the waypoint path which mitigates overshoot but can allow the vehicle to drive off-center.
The third most heavily weighted cost is rate of change of steering angle which forces the vehicle to take smooth turns but may reduce ability to take sharp corners.

## Time and prediction
Choosing number of timesteps `N` and elapsed time per step `dt` will determine the total prediction time `N * dt`. The prediction time must be large enough to predict a meaningful distance of the road ahead.
If the prediction time is too long there is a risk of too little importance being given to the immediate road ahead which can result in too sharp corners.
In addition, too long prediction time at high speeds can lead to the MPC running out of waypoints to fit which results in instability and undesired paths.

It was found that selecting values `N = 6` and `dt = 0.25 s` gives a good prediction which works well even at higher speeds. Previous values of 25 and 0.06 worked well but required more calculation without much benefit.
 Increasing prediction time further caused problem when fitting the polynomial since the waypoints weren't always visible that far ahead when driving at high speeds.

## Polynomial fit
The waypoints ahead are supplied by the simulator and converted into vehicle coordinates by translation and rotation.

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

The waypoint locations are used to fit a third order polynomial which will serve as the reference path for the vehicle.

      VectorXd coeffs = polyfit(ptsx_vehicle, ptsy_vehicle, 3);


The vehicle state is correspondingly updated by setting x, y, and &psi; to zero. The speed is converted from mph to m/s for consistent unts in the calculations.

## Latency
Latency of 100 ms is introduced which adds a delay between world measurements and engaging actuators. The latency is modelled as

`this_thread::sleep_for(chrono::milliseconds(latency_mil))`

and the vehicle position must be updated by `px += v * time_latency` to take into account the vehicle translation during the latency time until the actuators update.
The model assumes that the vehicle travels in a straight line during the latency time. More accurate results will be achieved by taking the steering angle into account.


## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets) == 0.14, but the master branch will probably work just fine
  * Follow the instructions in the [uWebSockets README](https://github.com/uWebSockets/uWebSockets/blob/master/README.md) to get setup for your platform. You can download the zip of the appropriate version from the [releases page](https://github.com/uWebSockets/uWebSockets/releases). Here's a link to the [v0.14 zip](https://github.com/uWebSockets/uWebSockets/archive/v0.14.0.zip).
  * If you have MacOS and have [Homebrew](https://brew.sh/) installed you can just run the ./install-mac.sh script to install this.
* Fortran Compiler
  * Mac: `brew install gcc` (might not be required)
  * Linux: `sudo apt-get install gfortran`. Additionall you have also have to install gcc and g++, `sudo apt-get install gcc g++`. Look in [this Dockerfile](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/Dockerfile) for more info.
* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt`
  * Linux
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from the Ipopt [releases page](https://www.coin-or.org/download/source/Ipopt/) or the [Github releases](https://github.com/coin-or/Ipopt/releases) page.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `bash install_ipopt.sh Ipopt-3.12.1`. 
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/CarND-MPC-Project/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

