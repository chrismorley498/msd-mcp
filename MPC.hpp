#pragma once

#include<vector>
#include<array>
#include<functional>

//For ode solver
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
typedef std::array<double, 2> state_type;

template<typename T>
struct SystemDynamics
{
    void operator()(const T &x, T &dxdt, const double /* t */)
    {
        dxdt[0] = x[1];    // x' = v
        dxdt[1] = -x[0]-2*x[1] + u_;   // v' = -x
    }

    // std::function<>

    void set_control_input(const double u){
        u_=u;
    }

    double u_{0.0};
};



template<typename T>
class HorizonPrediction{
    public:

    HorizonPrediction(){


    };

    ~HorizonPrediction(){

    };


    void step_over_horizon(T X, const double tf, const double dt){

        const std::size_t num_steps = tf/dt;
        double t{0.0};
        
        for(int i=0;i<num_steps;++i){
            t_vec_.push_back(t);
            state_vec_.push_back(X);

            stepper_.do_step(SystemDynamics<T>(), X, t, dt);
            t+=dt;
        }


    };

    std::vector<T> get_state_vec(){
        return state_vec_;
    };

    std::vector<double> get_time_vec(){
        return t_vec_;
    }




    private:
    runge_kutta_cash_karp54<state_type> stepper_;
    std::vector<double> t_vec_;
    std::vector<T> state_vec_;
    SystemDynamics<std::array<double,2>> system_;


};