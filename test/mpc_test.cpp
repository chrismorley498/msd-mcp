#include <vector>
#include <iostream>
#include <cstdio>

//For horizon predictor
#include "MPC.hpp"

//For matplotlib
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

//For ode solver
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

//Solver types
using state_type = std::array<double,2>;
using horizon_type = std::vector<state_type>;

int main()
{
    //Define dynamics of the system
    auto mass_spring_damper_dynamics_function = [](const auto &x, auto &dxdt, const double t, const double u){

        dxdt[0] = x[1];    // x' = v
        dxdt[1] = -x[0]-2*x[1] + u;   // v' = -x

    };

    //Define hyper parameters about trajectory
    constexpr double horizon_duration_s{2.0};
    constexpr double horizon_dt{0.1};//Not having a separate dt for now
    constexpr double desired_state_motion_frequency_hz{0.1};
    constexpr std::size_t num_horizon_steps = static_cast<std::size_t>(horizon_duration_s/horizon_dt);

    //Create MPC solver
    MPC<state_type, decltype(mass_spring_damper_dynamics_function)> mpc{mass_spring_damper_dynamics_function};

    //ODE stepper for simulation of plant
    HorizonPrediction<state_type, decltype(mass_spring_damper_dynamics_function)> horizon_predictor{mass_spring_damper_dynamics_function};

    //Define initial sim conditions
    // std::vector<double> u0{1,2,3,-3,-2,0};
    std::vector<double> u0{0,0,0,0,0,0,0,0,0,0};
    state_type x = {0.0, 0.0};

    //Vectors to hold results for plotting
    std::vector<double> pos, vel, desired_pos, desired_vel;//Holds results for plotting
    pos.reserve(1000);
    vel.reserve(1000);
    desired_pos.reserve(1000);
    desired_vel.reserve(1000);

    //Generate desired trajectory
    std::vector<double> eval_time{};
    eval_time.reserve(1000);
    horizon_type desired_state{};
    desired_state.reserve(1000);
    for(int i=0;i<num_horizon_steps;++i){
        const double cur_t = i*horizon_dt;
        eval_time.push_back(cur_t);
        desired_state.push_back({sin(2*M_PI*desired_state_motion_frequency_hz*cur_t), 2*M_PI*desired_state_motion_frequency_hz*cos(2*M_PI*desired_state_motion_frequency_hz*cur_t)});
        desired_pos.push_back(desired_state.back().at(0));
        desired_vel.push_back(desired_state.back().at(1));
    }

    //Solve for optimal control inputs
    mpc.find_optimal_control_inputs(x, desired_state, horizon_duration_s,horizon_dt, u0);

    //Predict expected trajectory using optimal control inputs
    std::vector<double> u_time = mpc.get_control_timepoints();
    horizon_predictor.step_over_horizon(x, horizon_duration_s, horizon_dt, u_time, u0);

    //Parse predicted states for plotting
    const horizon_type predicted_horizon = horizon_predictor.get_state_vec();
    std::vector<double> u{};
    u.reserve(1000);
    for(int i=0;i<num_horizon_steps;++i){
        pos.push_back(predicted_horizon.at(i).at(0));
        vel.push_back(predicted_horizon.at(i).at(1));
        u.push_back(linear_interpolate(u_time, u0,i*horizon_dt));
    }

    //Plot desired and predicted trajectories
    // plt::subplot(2,1,1);
    plt::plot(eval_time, desired_pos);
    plt::plot(eval_time, desired_vel);
    plt::plot(eval_time, pos,{{"linestyle", "--"}});
    plt::plot(eval_time, vel,{{"linestyle", "--"}});
    plt::plot(eval_time, u);
    plt::show();

    return 0;
}
