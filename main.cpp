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

// Harmonic oscillator system: x' = v, v' = -x
struct harmonic_oscillator
{
    void operator()(const state_type &x, state_type &dxdt, const double /* t */)
    {
        dxdt[0] = x[1];    // x' = v
        dxdt[1] = -x[0]-2*x[1];   // v' = -x
    }

    void set_control_input(const double u){
        u_=u;
    }

    double u_;
};

int main()
{
    // Initial conditions: x(0) = 1, v(0) = 0
    state_type x = {1.0, 0.0};

    std::vector<double> xpos, vpos;

    runge_kutta_cash_karp54<state_type> stepper;

    double t = 0.0, dt = 0.05;

    using state_type = std::array<double,2>;
    using horizon_type = std::vector<state_type>;
    HorizonPrediction<state_type> horizon_predictor{};

    horizon_predictor.step_over_horizon(x, 10.0, 0.1);
    const horizon_type predicted_horizon = horizon_predictor.get_state_vec();
    const std::vector<double> time_vec = horizon_predictor.get_time_vec();

    const std::size_t num_prediction_points = predicted_horizon.size();
    for(auto cur_state: predicted_horizon){
        xpos.push_back(cur_state.at(0));
        vpos.push_back(cur_state.at(1));
    }

    std::cout<<"size of xpos is: "<<xpos.size()<<std::endl;
    std::cout<<"Size of vpoc is: "<<vpos.size()<<std::endl;
    std::cout<<"Size of time vec is: "<<time_vec.size()<<std::endl;


    



    // ---------------------------------------------------------
    // GNUPlot (simple pipe-based interface)
    // ---------------------------------------------------------
    plt::plot(time_vec, xpos);
    plt::plot(time_vec, vpos);
    plt::show();

    return 0;
}
