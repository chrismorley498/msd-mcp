#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>
#include<array>
#include<functional>
#include <algorithm>  // for std::lower_bound


//For ode solver
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
typedef std::array<double, 2> state_type;




double linear_interpolate(const std::vector<double>& x,
                          const std::vector<double>& y,
                          double x_eval)
{
    if (x.size() != y.size() || x.empty())
        throw std::invalid_argument("x and y must have the same non-zero size");

    // Clamp x_eval to [x.front(), x.back()]
    if (x_eval <= x.front()) return y.front();
    if (x_eval >= x.back())  return y.back();

    // Find the interval [x[i], x[i+1]] that contains x_eval
    auto it = std::lower_bound(x.begin(), x.end(), x_eval);
    int idx = std::distance(x.begin(), it) - 1;

    // Linear interpolation
    double t = (x_eval - x[idx]) / (x[idx + 1] - x[idx]);
    double y_eval = y[idx] + t * (y[idx + 1] - y[idx]);

    return y_eval;
}

// Fit polynomial of degree 'degree' to points (x[i], y[i])
// Returns polynomial coefficients {c0, c1, ..., c_degree}
Eigen::VectorXd polyfit(const std::vector<double>& x,
                            const std::vector<double>& y,
                            int degree)
{
    const int N = x.size();
    const int M = degree + 1;

    Eigen::MatrixXd A(N, M);
    Eigen::VectorXd b(N);

    for (int i = 0; i < N; ++i) {
        double xi = 1.0;
        for (int j = 0; j < M; ++j) {
            A(i, j) = xi;
            xi *= x[i];
        }
        b(i) = y[i];
    }

    // Solve the least squares problem A * c = b
    Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

    // Copy into std::vector
    std::vector<double> coeffs(M);
    for (int i = 0; i < M; ++i)
        coeffs[i] = c(i);

    return c;
}

double polyeval(const Eigen::VectorXd coeffs, double x){
    const std::size_t num_coeffs = coeffs.size();

    double val{0.0};
    for(int i=0;i<num_coeffs;++i){
        val+=std::pow(x,i)*coeffs(i);
    }

    return val;
}


template<typename T>
struct SystemDynamics
{
    void operator()(const T &x, T &dxdt, const double /* t */)
    {
        dxdt[0] = x[1];    // x' = v
        dxdt[1] = -x[0]-2*x[1] + u_;   // v' = -x
    }


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


    void step_over_horizon(T X, const double tf, const double dt, std::vector<double> u_timepoints, std::vector<double> u_vals){

        assert(u_timepoints.size() == u_vals.size());
        const std::size_t num_steps = tf/dt;
        double t{0.0};
        
        for(int i=0;i<num_steps;++i){
            t_vec_.push_back(t);
            state_vec_.push_back(X);

            system_.set_control_input(linear_interpolate(u_timepoints, u_vals, t));
            stepper_.do_step(system_, X, t, dt);
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


class MPC{
    public:




    private:

    HorizonPrediction horizon_predictor_{};



};