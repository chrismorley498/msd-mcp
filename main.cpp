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
typedef std::array<double, 2> state_type;

int main()
{
    // //Define dynamics of the system
    // std::function<void(const state_type, state_type, const double)> mass_spring_damper_dynamics= [](const state_type &x, state_type &dxdt, const double /* t */){

    // };
    auto mass_spring_damper_dynamics_function = [](const auto &x, auto &dxdt, const double t, const double u){

        dxdt[0] = x[1];    // x' = v
        dxdt[1] = -x[0]-2*x[1] + u;   // v' = -x

    };


    // Initial conditions: x(0) = 1, v(0) = 0
    state_type x = {0.0, 0.0};

    std::vector<double> xpos, vpos;

    runge_kutta_cash_karp54<state_type> stepper;

    double t = 0.0;

    using state_type = std::array<double,2>;
    using horizon_type = std::vector<state_type>;
    HorizonPrediction<state_type, decltype(mass_spring_damper_dynamics_function)> horizon_predictor{mass_spring_damper_dynamics_function};
    // horizon_predictor.set_dynamics_function(mass_spring_damper_dynamics_function);


    //Define hyper parameters about trajectory
    constexpr double horizon_duration_s{10.0};
    constexpr double horizon_dt{0.1};

    //Define target trajectory
    std::vector<double> eval_points{};
    std::vector<double> eval_vels{};
    std::vector<double> eval_time{};
    horizon_type desired_state{};
    const double desired_state_frequency_hz{0.1};
    for(int i=0;i<100;++i){
        double cur_t = i*0.1;
        eval_time.push_back(cur_t);
        desired_state.push_back({sin(2*M_PI*desired_state_frequency_hz*cur_t), 2*M_PI*desired_state_frequency_hz*cos(2*M_PI*desired_state_frequency_hz*cur_t)});
        eval_points.push_back(desired_state.back().at(0));
        eval_vels.push_back(desired_state.back().at(1));
    }

    auto print_control_input = [](std::vector<double> vec){
        for(const auto val: vec){
            std::cout<<val<<",";
        }
        std::cout<<std::endl;
    };


    MPC<state_type, decltype(mass_spring_damper_dynamics_function)> mpc{mass_spring_damper_dynamics_function};
    // mpc.set_dynamics_function(mass_spring_damper_dynamics_function);
    std::vector<double> u0{1,2,3,-3,-2,0};
    std::cout<<"Initial guess\n";
    print_control_input(u0);
    std::cout<<std::endl;
    mpc.find_optimal_control_inputs(x, desired_state, horizon_duration_s,horizon_dt, u0);
    std::cout<<"Solution\n";
    print_control_input(u0);
    std::cout<<std::endl;

    //Define control inputs
    std::vector<double> u_time{0,2,4,6,8,10};
    std::vector<double> u_vals{0,1,2,3,-3,-1};


    horizon_predictor.step_over_horizon(x, horizon_duration_s, horizon_dt, u_time, u0);
    const horizon_type predicted_horizon = horizon_predictor.get_state_vec();
    const std::vector<double> time_vec = horizon_predictor.get_time_vec();

    const std::size_t num_prediction_points = predicted_horizon.size();
    for(auto cur_state: predicted_horizon){
        xpos.push_back(cur_state.at(0));
        vpos.push_back(cur_state.at(1));
    }



    plt::plot(eval_time, eval_points);
    plt::plot(eval_time, eval_vels);
    plt::plot(eval_time, xpos,{{"linestyle", "--"}});
    plt::plot(eval_time, vpos,{{"linestyle", "--"}});
    plt::show();


    



    // ---------------------------------------------------------
    // GNUPlot (simple pipe-based interface)
    // ---------------------------------------------------------
    // plt::plot(time_vec, xpos);
    // plt::plot(time_vec, vpos);
    // plt::show();

    return 0;
}
