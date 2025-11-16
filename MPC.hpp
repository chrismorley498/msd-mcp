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

//For numerical solver
#include<nlopt.hpp>




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

    void empty_state_vec(){
        state_vec_.clear();
    }

    void empty_time_vec(){
        t_vec_.clear();
    }




    private:
    runge_kutta_cash_karp54<state_type> stepper_;
    std::vector<double> t_vec_;
    std::vector<T> state_vec_;
    SystemDynamics<std::array<double,2>> system_;


};

template<typename T>
struct IKData;//Forward decl

template<typename T>
class MPC{
    public:

    std::vector<double> find_optimal_control_inputs(T x0, std::vector<T>& desired_trajectory, const double buffer_duration, const double dt, std::vector<double> &u0){


        const std::size_t num_control_points = u0.size();
        const double control_point_linear_segment_duration = buffer_duration/(num_control_points-1);
        std::vector<double> control_point_times{};
        for(int i=0;i<num_control_points;++i){
            std::cout<<"Control input time: "<<i*control_point_linear_segment_duration<<std::endl;
            control_point_times.push_back(i*control_point_linear_segment_duration);
        }

        control_point_times_= control_point_times;
        num_control_points_= num_control_points;
        buffer_duration_=buffer_duration;
        dt_=dt;

        auto cost_function = [](const std::vector<double> &u, std::vector<double> &grad, void *data){

            IKData<T> *d = static_cast<IKData<T>*>(data);
            MPC<T>* mpc = d->mpc_ptr;

            mpc->horizon_predictor_.empty_state_vec();
            mpc->horizon_predictor_.empty_time_vec();

            mpc->horizon_predictor_.step_over_horizon(d->x0, d->buffer_duration, d->dt, d->control_point_times, u);

            std::vector<T> predicted_trajectory = mpc->horizon_predictor_.get_state_vec();

            assert(predicted_trajectory.size() == d->desired_state.size());
            const std::size_t num_trajectory_points = d->desired_state.size();
            double position_error{0.0};
            double velocity_error{0.0};
            for(int i=0;i<num_trajectory_points;++i){
                double predicted_position = predicted_trajectory.at(i).at(0);
                double predicted_velocity = predicted_trajectory.at(i).at(1);

                double desired_position = d->desired_state.at(i).at(0);
                double desired_velocity = d->desired_state.at(i).at(1);

                position_error+=std::abs(predicted_position - desired_position);
                velocity_error+=std::abs(predicted_velocity - desired_velocity);
            }

            // std::cout<<"Position error: "<<position_error<<". Velocity error: "<<velocity_error<<std::endl;


            return position_error;

        };



        //Setup cost function
        IKData<T> ik_data{
            .mpc_ptr=this, 
            .desired_state=desired_trajectory, 
            .x0=x0, 
            .buffer_duration=buffer_duration, 
            .dt=dt, 
            .control_point_times=control_point_times
            };
        nlopt::opt optimizer(nlopt::LN_COBYLA, num_control_points);
        optimizer.set_min_objective(cost_function, &ik_data);
        // optimizer.set_xtol_rel(1e-4);
        optimizer.set_ftol_rel(1e-3);
        double val{0.0};
        int result{};
        try{
            result = optimizer.optimize(u0,val);
        } catch (const std::exception &e) {
            std::cerr << "An error occurred: " << e.what() << std::endl;

            // Specific handling for nlopt.RoundoffLimited (if you need special handling)
            if (std::string(e.what()).find("RoundoffLimited") != std::string::npos) {
                std::cerr << "The optimization encountered a roundoff-limited issue." << std::endl;
                const int u_length = u0.size();
                for(int i = 0;i<u_length;++i){
                    u0.at(i) = u0.at(i) + 0.001;
                    result = optimizer.optimize(u0,val);
                }
            }
        }

        return u0;


    };




    private:

    HorizonPrediction<T> horizon_predictor_{};

    //Used within cost function
    std::vector<double> control_point_times_;
    std::size_t num_control_points_;
    double buffer_duration_;
    double dt_;



};

template<typename T>
struct IKData{
    MPC<T>* mpc_ptr;
    const std::vector<T> &desired_state;
    const T x0;
    const double buffer_duration;
    const double dt;
    const std::vector<double> control_point_times;
};