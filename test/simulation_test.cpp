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
    // //Define dynamics of the system
    auto mass_spring_damper_dynamics_function = [](const auto &x, auto &dxdt, const double t, const double u){

        dxdt[0] = x[1];    // x' = v
        dxdt[1] = -x[0]-2*x[1] + u;   // v' = -x

    };

    auto print_control_input = [](std::vector<double> vec){
        for(const auto val: vec){
            std::cout<<val<<",";
        }
        std::cout<<std::endl;
    };

    //Define hyper parameters about trajectory
    constexpr double horizon_duration_s{10.0};
    constexpr double horizon_dt{0.1};//Not having a separate dt for now
    constexpr double desired_state_motion_frequency_hz{0.1};
    constexpr double total_sim_time{10.0};
    constexpr double sim_dt{0.1};
    constexpr std::size_t num_sim_steps = static_cast<std::size_t>(total_sim_time/sim_dt);
    constexpr std::size_t num_buffer_steps = static_cast<std::size_t>(horizon_duration_s/horizon_dt);

    //Create MPC solver
    MPC<state_type, decltype(mass_spring_damper_dynamics_function)> mpc{mass_spring_damper_dynamics_function};

    //ODE stepper for simulation of plant
    HorizonPrediction<state_type, decltype(mass_spring_damper_dynamics_function)> horizon_predictor{mass_spring_damper_dynamics_function};

    //Define initial sim conditions
    std::vector<double> u0{1,2,3,-3,-2,0};
    state_type x = {0.0, 0.0};
    std::vector<double> pos, vel, desired_pos, desired_vel, sim_time_vec;//Holds results for plotting


    auto step_simulation = [&](){
        std::vector<double> u_time = mpc.get_control_timepoints();
        // print_control_input(u_time);
        horizon_predictor.empty_time_vec();
        horizon_predictor.empty_state_vec();
        std::cout<<"Stepping sim\n";
        std::cout<<"Pos before step: "<<x.at(0);
        // horizon_predictor.step_over_horizon(x, sim_dt, sim_dt, u_time, u0);
        horizon_predictor.single_step(x, sim_dt, sim_dt, u_time, u0);
        // auto calculated_trajectory = horizon_predictor.get_state_vec();
        // x=calculated_trajectory.back();
        std::cout<<"Pos after step: "<<x.at(0)<<'\n';
        // std::cout<<"Size of state is: "<<calculated_trajectory.size()<<std::endl;
        pos.push_back(x.at(0));
        vel.push_back(x.at(1));
    };

    //Loop through simulation steps
    for(int i=0;i<num_sim_steps;++i){

        double cur_t = i*sim_dt;
        sim_time_vec.push_back(cur_t);

        //Create trajectory buffer
        std::vector<double> eval_points{};
        std::vector<double> eval_vels{};
        std::vector<double> eval_time{};
        horizon_type desired_state{};
        for(int j=0;j<num_buffer_steps;++j){
            eval_time.push_back(cur_t);
            desired_state.push_back({sin(2*M_PI*desired_state_motion_frequency_hz*cur_t), 2*M_PI*desired_state_motion_frequency_hz*cos(2*M_PI*desired_state_motion_frequency_hz*cur_t)});
            eval_points.push_back(desired_state.back().at(0));
            eval_vels.push_back(desired_state.back().at(1));
        }

        desired_pos.push_back(eval_points.front());
        desired_vel.push_back(eval_vels.front());

        //Solve for new control inputs
        mpc.find_optimal_control_inputs(x, desired_state, horizon_duration_s,horizon_dt, u0);

        //Step sim
        step_simulation();

    }//End of sim loop










    // std::cout<<"Initial guess\n";
    // print_control_input(u0);
    // std::cout<<std::endl;
    // std::cout<<"Solution\n";
    // print_control_input(u0);
    // std::cout<<std::endl;

    // //Define control inputs
    // std::vector<double> u_time{0,2,4,6,8,10};
    // std::vector<double> u_vals{0,1,2,3,-3,-1};


    // const horizon_type predicted_horizon = horizon_predictor.get_state_vec();
    // const std::vector<double> time_vec = horizon_predictor.get_time_vec();

    // const std::size_t num_prediction_points = predicted_horizon.size();
    // for(auto cur_state: predicted_horizon){
    //     xpos.push_back(cur_state.at(0));
    //     vpos.push_back(cur_state.at(1));
    // }



    plt::plot(sim_time_vec, pos);
    plt::plot(sim_time_vec, vel);
    plt::plot(sim_time_vec, desired_pos,{{"linestyle", "--"}});
    plt::plot(sim_time_vec, desired_vel,{{"linestyle", "--"}});
    plt::show();


    



    // ---------------------------------------------------------
    // GNUPlot (simple pipe-based interface)
    // ---------------------------------------------------------
    // plt::plot(time_vec, xpos);
    // plt::plot(time_vec, vpos);
    // plt::show();

    return 0;
}
