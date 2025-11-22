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


struct RobotModel{

    RobotModel(const std::string urdf_str){

        pinocchio::urdf::buildModel(urdf_str, model);
        data=pinocchio::Data(model);

    };

    void update_state(Eigen::VectorXd& q, Eigen::VectorXd& dq){

        //Update forward kinematics
        pinocchio::forwardKinematics(model, data, q, dq);
        pinocchio::updateFramePlacements(model, data);
        // pinocchio::framesForwardKinematics(model, data, q);//This single line can replace the two lines above.
        // pinocchio::SE3 T_temp = data.oMf.at(ee_frame_id); // T is the transformation matrix (SE3)
    }

    Eigen::VectorXd compute_acceleration(Eigen::VectorXd& q, Eigen::VectorXd& dq, Eigen::VectorXd& tau){

        Eigen::MatrixXd M;
        Eigen::MatrixXd c;

        pinocchio::computeAllTerms(model, data, q, dq);
        M = data.M;       // Mass matrix
        c = data.nle;     // Nonlinear effects (C*v + G)
        // Eigen::VectorXd G = pinocchio::computeGeneralizedGravity(model, data, q);
        constexpr double DAMPING=0.1;
        Eigen::VectorXd a = M.ldlt().solve(tau - c - DAMPING*dq); // joint acceleration: M*vdot = tau - C - G

        return a;

    }





private:

    //Pinocchio Objects
    pinocchio::Model model;
    pinocchio::Data data;


};

int main()
{
    // //Define dynamics of the system
    constexpr double mass{1.0};
    auto mass_spring_damper_dynamics_function = [](const auto &x, auto &dxdt, const double t, const double u){

        dxdt[0] = 1/mass*(x[1]);    // x' = v
        dxdt[1] = 1/mass*(-x[0]-2*x[1] + u);   // v' = -x

    };

    RobotModel robot("/home/chris-morley/code/dev-ws/src/3DofAssemblyWithMotors/3DofAssembly_V2/urdf/3DofAssembly_V2.urdf");

    auto three_link_arm_dynamics_function = [&](auto &x, auto &dxdt, const double t, auto& u){

        Eigen::VectorXd q = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd dq = Eigen::VectorXd::Zero(3);
        for(int i=0;i<3;++i){
            q(i)=x.at(i);
            dq(i)=x.at(i+3);
        }


        // std::cout<<"q: "<<q.transpose()<<std::endl;
        // std::cout<<"dq: "<<dq.transpose()<<std::endl;

        Eigen::VectorXd acc = robot.compute_acceleration(q, dq, u);
        for(int i=0;i<3;++i){
            dxdt.at(i)=x.at(i+3);
            dxdt.at(i+3)=acc(i);
        }

        // std::cout<<"dxdt: [";
        // std::for_each(dxdt.begin(),dxdt.end(),[](const double val){
        //     std::cout<<val<<", ";
        // });
        // std::cout<<"]\n";
        

    };


    Eigen::VectorXd q = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(3);
    std::array<double,6> X = std::array<double,6>{};
    std::array<double,6> dX = std::array<double,6>{};
    Eigen::VectorXd tau = Eigen::VectorXd::Ones(3)*0.01;
    constexpr double temp_t{0.0};

    using ArmStateType = std::array<double,6>;
    using ArmControlInputType = Eigen::VectorXd;

    SystemDynamics<ArmStateType, decltype(three_link_arm_dynamics_function), ArmControlInputType> arm{three_link_arm_dynamics_function};

    runge_kutta_cash_karp54<ArmStateType> stepper;
    for(int i=0;i<1000;++i){
        arm.set_control_input(tau);
        stepper.do_step(arm, X, 0.0, 0.01);
    }


    three_link_arm_dynamics_function(X,dX,temp_t,tau);
    // std::cout<<"Final joint state: "<<dX.transpose()<<std::endl;

    auto print_control_input = [](std::vector<double> vec){
        for(const auto val: vec){
            std::cout<<val<<",";
        }
        std::cout<<std::endl;
    };
    print_control_input(std::vector<double>(dX.begin(), dX.end()));

    //Define hyper parameters about trajectory
    constexpr double horizon_duration_s{0.5};
    constexpr double horizon_dt{0.1};//Not having a separate dt for now
    constexpr double desired_state_motion_frequency_hz{0.1};
    constexpr double total_sim_time{10.0};
    constexpr double sim_dt{0.1};
    constexpr std::size_t num_sim_steps = static_cast<std::size_t>(total_sim_time/sim_dt);
    constexpr std::size_t num_buffer_steps = static_cast<std::size_t>(horizon_duration_s/horizon_dt);

    //Create MPC solver
    MPC<state_type, decltype(mass_spring_damper_dynamics_function), double> mpc{mass_spring_damper_dynamics_function};

    //ODE stepper for simulation of plant
    HorizonPrediction<state_type, decltype(mass_spring_damper_dynamics_function), double> horizon_predictor{mass_spring_damper_dynamics_function};

    //Define initial sim conditions
    // std::vector<double> u0{1,2,3,-3,-2,0};
    std::vector<double> u0{0,0,0,0,0,0,0,0};
    u0.reserve(1000);
    state_type x = {0.0, 0.0};
    std::vector<double> pos, vel, desired_pos, desired_vel, sim_time_vec;//Holds results for plotting
    pos.reserve(1000);
    vel.reserve(1000);
    desired_pos.reserve(1000);
    desired_vel.reserve(1000);
    sim_time_vec.reserve(1000);

    auto step_simulation = [&](){
        std::vector<double> u_time = mpc.get_control_timepoints();
        // print_control_input(u_time);
        horizon_predictor.empty_time_vec();
        horizon_predictor.empty_state_vec();

        //Whole trajectory mode
        horizon_predictor.empty_state_vec();
        horizon_predictor.empty_time_vec();
        horizon_predictor.step_over_horizon(x, sim_dt, sim_dt, u_time, u0);
        auto calculated_trajectory = horizon_predictor.get_state_vec();
        x=calculated_trajectory.back();

        //Update plotting vectors
        pos.push_back(x.at(0));
        vel.push_back(x.at(1));

    };

    //Loop through simulation steps
    for(int i=0;i<num_sim_steps;++i){
        u0 = std::vector<double>{0,0,0,0,0,0,0,0};//Try set initial guess at each iteration

        double cur_t = i*sim_dt;
        sim_time_vec.push_back(cur_t);

        //Create trajectory buffer
        std::vector<double> eval_points{};
        std::vector<double> eval_vels{};
        std::vector<double> eval_time{};
        horizon_type desired_state{};
        eval_points.reserve(1000);
        eval_vels.reserve(1000);
        eval_time.reserve(1000);
        desired_state.reserve(1000);
        for(int j=0;j<num_buffer_steps;++j){
            const double desired_trajectory_buffer_time = cur_t+j*horizon_dt;
            eval_time.push_back(desired_trajectory_buffer_time);
            desired_state.push_back({sin(2*M_PI*desired_state_motion_frequency_hz*desired_trajectory_buffer_time), 2*M_PI*desired_state_motion_frequency_hz*cos(2*M_PI*desired_state_motion_frequency_hz*desired_trajectory_buffer_time)});
            eval_points.push_back(desired_state.back().at(0));
            eval_vels.push_back(desired_state.back().at(1));
        }

        desired_pos.push_back(eval_points.front());
        desired_vel.push_back(eval_vels.front());

        //Solve for new control inputs
        mpc.find_optimal_control_inputs(x, desired_state, horizon_duration_s,horizon_dt, u0);
        std::vector<double> control_timepoints = mpc.get_control_timepoints();

        //Step sim
        step_simulation();



        // std::vector<double> u_time = mpc.get_control_timepoints();
        // horizon_predictor.empty_state_vec();
        // horizon_predictor.empty_time_vec();
        // horizon_predictor.step_over_horizon(x, horizon_duration_s, horizon_dt, u_time, u0);
        // auto calculated_trajectory = horizon_predictor.get_state_vec();
        // x=calculated_trajectory.front();
        // std::vector<double> inter_pos, inter_vel, inter_u;
        // for(int j=0;j<num_buffer_steps;++j){
        //     inter_pos.push_back(calculated_trajectory.at(j).at(0));
        //     inter_vel.push_back(calculated_trajectory.at(j).at(1));
        //     inter_u.push_back(linear_interpolate(u_time, u0,j*horizon_dt));
        // }
        // plt::plot(eval_time, eval_points);
        // plt::plot(eval_time, eval_vels);
        // plt::plot(eval_time, inter_pos,{{"linestyle", "--"}});
        // plt::plot(eval_time, inter_vel,{{"linestyle", "--"}});
        // plt::plot(eval_time, inter_u);
        // plt::show();
        // plt::clear();



    }//End of sim loop



    // // //Plot final results of sim
    // plt::plot(sim_time_vec, desired_pos);
    // plt::plot(sim_time_vec, desired_vel);
    // plt::plot(sim_time_vec, pos,{{"linestyle", "--"}});
    // plt::plot(sim_time_vec, vel,{{"linestyle", "--"}});
    // plt::show();

    return 0;
}
