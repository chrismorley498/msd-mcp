#pragma once

// #define BOOST_MPL_LIMIT_LIST_SIZE 30

//For robotics library
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include "pinocchio/algorithm/rnea.hpp"


#include <vector>
#include <eigen3/Eigen/Dense>
#include <array>
#include <algorithm>
#include <nlopt.hpp>
#include <boost/numeric/odeint.hpp>



using namespace boost::numeric::odeint;
typedef std::array<double, 2> state_type;
using control_input_type = double;


/* ============================================================
   Utility functions (unchanged)
   ============================================================ */

double linear_interpolate(const std::vector<double>& x,
                          const std::vector<double>& y,
                          double x_eval)
{
    if (x.size() != y.size() || x.empty())
        throw std::invalid_argument("x and y must have the same non-zero size");

    if (x_eval <= x.front()) return y.front();
    if (x_eval >= x.back())  return y.back();

    auto it = std::lower_bound(x.begin(), x.end(), x_eval);
    int idx = std::distance(x.begin(), it) - 1;

    double t = (x_eval - x[idx]) / (x[idx + 1] - x[idx]);
    return y[idx] + t * (y[idx + 1] - y[idx]);
}

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
    return A.colPivHouseholderQr().solve(b);
}

double polyeval(const Eigen::VectorXd coeffs, double x){
    double val=0;
    for(int i=0;i<coeffs.size();++i)
        val += std::pow(x,i) * coeffs(i);
    return val;
}


/* ============================================================
   1. SystemDynamics — templated, zero overhead
   ============================================================ */

template<typename T, typename Dynamics>
struct SystemDynamics
{
    Dynamics dynamics_function_;   // stored with full type
    double u_{0.0};

    // Constructor injection of the dynamics functor
    SystemDynamics(const Dynamics& dyn) : dynamics_function_(dyn) {}

    void operator()(const T &x, T &dxdt, const double /* t */)
    {
        dynamics_function_(x, dxdt, 0.0, u_);
    }

    void set_control_input(double u) { u_ = u; }

    // Overwrite dynamics (same type)
    void set_dynamics_function(const Dynamics& f) {
        dynamics_function_ = f;
    }
};


/* ============================================================
   2. HorizonPrediction — templated on Dynamics too
   ============================================================ */

template<typename T, typename Dynamics>
class HorizonPrediction {
public:

    HorizonPrediction(const Dynamics& dyn)
        : system_(dyn)
    {
        t_vec_.reserve(1000);
        state_vec_.reserve(1000);
    }

    void set_dynamics_function(const Dynamics& f) {
        system_.set_dynamics_function(f);
    }

    void empty_state_vec(){ state_vec_.clear(); }
    void empty_time_vec(){ t_vec_.clear(); }

    std::vector<T> get_state_vec(){ return state_vec_; }
    std::vector<double> get_time_vec(){ return t_vec_; }

    void step_over_horizon(T X, double tf, double dt,
                           const std::vector<double>& u_timepoints,
                           const std::vector<double>& u_vals)
    {
        const std::size_t num_steps = tf/dt;
        // std::cout<<"Num steps: "<<num_steps<<std::endl;
        double t = 0.0;

        for (std::size_t i = 0; i < num_steps; ++i) {
            // t_vec_.push_back(t);
            // state_vec_.push_back(X);

            double u = linear_interpolate(u_timepoints, u_vals, t);
            system_.set_control_input(u);

            stepper_.do_step(system_, X, t, dt);
            t += dt;


            t_vec_.push_back(t);
            state_vec_.push_back(X);

        }
    }


    void single_step(T& X, double tf, double dt,
                           std::vector<double>& u_timepoints,
                           const std::vector<double>& u_vals)
    {

        const double u0=u_timepoints.front();
        std::for_each(u_timepoints.begin(),u_timepoints.end(), [&](double &val){
            val-=u0;
        });

        double u = linear_interpolate(u_timepoints, u_vals, 0.0);
        // std::cout<<"u: "<<u<<'\n';
        system_.set_control_input(u);
        constexpr double t_place_holder{0.0};
        stepper_.do_step(system_, X, t_place_holder, dt);   
    }


private:

    runge_kutta_cash_karp54<state_type> stepper_;
    std::vector<double> t_vec_;
    std::vector<T> state_vec_;
    SystemDynamics<T, Dynamics> system_;
};


/* ============================================================
   Forward declare IKData for MPC
   ============================================================ */

template<typename T, typename Dynamics>
struct IKData;


/* ============================================================
   3. MPC — also templated on dynamics functor type
   ============================================================ */

template<typename T, typename Dynamics>
class MPC {
public:

    MPC(const Dynamics& dyn)
        : horizon_predictor_(dyn)
    {
        control_point_times_.reserve(1000);
    }

    void set_dynamics_function(const Dynamics& f) {
        horizon_predictor_.set_dynamics_function(f);
    }

    std::vector<double> find_optimal_control_inputs(
            T x0,
            std::vector<T>& desired_trajectory,
            double buffer_duration,
            double dt,
            std::vector<double> &u0)
    {
        const std::size_t num_control_points = u0.size();
        const double seg_dt = buffer_duration/(num_control_points - 1);

        control_point_times_.clear();
        for (std::size_t i = 0; i < num_control_points; ++i)
            control_point_times_.push_back(i * seg_dt);

        buffer_duration_ = buffer_duration;
        dt_ = dt;
        num_control_points_ = num_control_points;

        auto cost_function = [](const std::vector<double> &u,
                                std::vector<double> &grad,
                                void *data)
        {
            IKData<T, Dynamics>* d = static_cast<IKData<T, Dynamics>*>(data);
            MPC<T, Dynamics>* mpc = d->mpc_ptr;

            mpc->horizon_predictor_.empty_state_vec();
            mpc->horizon_predictor_.empty_time_vec();

            mpc->horizon_predictor_.step_over_horizon(
                d->x0, d->buffer_duration, d->dt,
                d->control_point_times, u
            );

            std::vector<T> pred = mpc->horizon_predictor_.get_state_vec();
            const auto& desired = d->desired_state;

            double position_error = 0.0;
            double velocity_error = 0.0;

            for (std::size_t i = 0; i < desired.size(); ++i){
                position_error += std::abs(pred[i][0] - desired[i][0]);
                velocity_error += std::abs(pred[i][1] - desired[i][1]);
            }
            return position_error + 0.2*velocity_error;
        };

        IKData<T, Dynamics> ik_data{
            .mpc_ptr = this,
            .desired_state = desired_trajectory,
            .x0 = x0,
            .buffer_duration = buffer_duration,
            .dt = dt,
            .control_point_times = control_point_times_
        };

        nlopt::opt opt(nlopt::LN_COBYLA, num_control_points);
        opt.set_min_objective(cost_function, &ik_data);
        opt.set_ftol_rel(1e-5);

        double val{};
        opt.optimize(u0, val);

        return u0;
    }

    std::vector<double> get_control_timepoints(){
        return control_point_times_;
    }

private:

    HorizonPrediction<T, Dynamics> horizon_predictor_;

    // Stored for the cost function
    std::vector<double> control_point_times_;
    std::size_t num_control_points_;
    double buffer_duration_;
    double dt_;
};


/* ============================================================
   4. IKData definition
   ============================================================ */

template<typename T, typename Dynamics>
struct IKData {
    MPC<T, Dynamics>* mpc_ptr;

    const std::vector<T>& desired_state;
    const T x0;
    const double buffer_duration;
    const double dt;
    const std::vector<double> control_point_times;
};

