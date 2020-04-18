/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

//------------------------------------------
//------------------------------------------
//--------  Trim and Linearize  ------------
//------------------------------------------
//------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#include "flight_sim.h"
#include "vehicles.h"

#include <cdyn/integrate.h>
#include <cdyn/simulate.h>

struct TrimSpec
{
    real z_dot;
    real yaw_dot;
    real target_vel;
    
    struct Aircraft * ac;

    real thresh;
};

inline real trim_spec_get_climb_rate(const struct TrimSpec * spec){ return -spec->z_dot; }
inline real trim_spec_get_yaw_rate(const struct TrimSpec * spec){ return spec->yaw_dot; }
inline real trim_spec_get_speed(const struct TrimSpec * spec){ return spec->target_vel; }

double trim_objective(unsigned n, const double * x, double * grad, void * f_data)
{
    (void) n;
    assert (grad == NULL);
    
    struct TrimSpec * data = f_data;

    double ic[12];
    ic[0] = 0.0; // x
    ic[1] = 0.0; // y
    ic[2] = 0.0; // z
    ic[3] = x[0]; // u
    ic[4] = x[1]; // v
    ic[5] = x[2]; // w
    ic[6] = x[3]; // p
    ic[7] = x[4]; // q
    ic[8] = x[5]; // r
    ic[9] = x[6]; // roll
    ic[10] = x[7]; // pitch
    ic[11] = 0.0; //yaw

    double control[4];
    control[0] = x[8];
    control[1] = x[9];
    control[2] = x[10];
    control[3] = x[11];
    
    double sol[12];
    rigid_body_lin_forces(0.0, ic, control, sol, NULL, data->ac);

    double out_trans = pow(sol[3], 2) + pow(sol[4], 2) + pow(sol[5], 2);
    double out_rot = pow(sol[6], 2) + pow(sol[7], 2) + pow(sol[8], 2);
    double out_trim = pow(sol[9], 2) + pow(sol[10], 2);

    double out_zdot = pow(sol[2] - data->z_dot, 2);
    double out_yawdot = pow(sol[11] - data->yaw_dot, 2);

    
    double vel = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
    double out_vel = pow(vel - data->target_vel, 2);

    double out_sideslip = pow(x[1], 2);

    /* printf("%3.2E, %3.2E, %3.2E, %3.2E, %3.2E, %3.2E, %3.2E\n", */
    /*        out_trans, out_rot, out_trim, out_zdot, out_yawdot, out_vel, out_sideslip); */
    
    double out = out_trans + out_rot + out_trim + out_zdot + out_yawdot + out_vel + out_sideslip;

    /* printf("out = %3.5E\n", out); */
    /* double out = out_trans; */
    return out;
}

double trim_objective_g(unsigned n, const double * x, double * grad, void * f_data)
{
    if (grad == NULL){
        return trim_objective(n, x, grad, f_data);
    }
    
    struct TrimSpec * data = f_data;

    double ic[12];
    ic[0] = 0.0; // x
    ic[1] = 0.0; // y
    ic[2] = 0.0; // z
    ic[3] = x[0]; // u
    ic[4] = x[1]; // v
    ic[5] = x[2]; // w
    ic[6] = x[3]; // p
    ic[7] = x[4]; // q
    ic[8] = x[5]; // r
    ic[9] = x[6]; // roll
    ic[10] = x[7]; // pitch
    ic[11] = 0.0; //yaw

    double control[4];
    control[0] = x[8];
    control[1] = x[9];
    control[2] = x[10];
    control[3] = x[11];
    
    double sol[12];
    double jac[144 + 48]; // states and controls
    rigid_body_lin_forces_jac(0.0, ic, control, sol, jac, data->ac);

    double out_trans = pow(sol[3], 2) + pow(sol[4], 2) + pow(sol[5], 2);
    // U, V, W equations
    for (size_t jj = 0; jj < 9; jj++){
        grad[jj] = 2 * sol[3] * jac[(jj+3)*12 + 3] +
                   2 * sol[4] * jac[(jj+3)*12 + 4] +
                   2 * sol[5] * jac[(jj+3)*12 + 5];
    }
    for (size_t jj = 0; jj < 4; jj++){
        grad[jj+8] += 2 * sol[3] * jac[(12+jj)*12 + 3] +
                      2 * sol[4] * jac[(12+jj)*12 + 4] +
                      2 * sol[5] * jac[(12+jj)*12 + 5];
    }
    
    double out_rot = pow(sol[6], 2) + pow(sol[7], 2) + pow(sol[8], 2);
    // P, Q, R equations
    for (size_t jj = 0; jj < 9; jj++){
        grad[jj] += (2 * sol[6] * jac[(jj+3)*12 + 6] +
                     2 * sol[7] * jac[(jj+3)*12 + 7] +
                     2 * sol[8] * jac[(jj+3)*12 + 8]);
    }

    for (size_t jj = 0; jj < 4; jj++){
        grad[jj+8] += (2 * sol[6] * jac[(12+jj)*12 + 6] +
                       2 * sol[7] * jac[(12+jj)*12 + 7] +
                       2 * sol[8] * jac[(12+jj)*12 + 8]);
    }    
    
    double out_trim = pow(sol[9], 2) + pow(sol[10], 2);
    // Roll and pitch equations
    for (size_t jj = 0; jj < 9; jj++){
        grad[jj] += 2 * sol[9] * jac[(jj+3)*12 + 9] +
                    2 * sol[10] * jac[(jj+3)*12 + 10];
    }
    for (size_t jj = 0; jj < 4; jj++){
        grad[jj+8] += 2 * sol[9] * jac[(jj+12)*12 + 9] +
                      2 * sol[10] * jac[(jj+12)*12 + 10];
    }
    
    double out_zdot = pow(sol[2] - data->z_dot, 2);
    for (size_t jj = 0; jj < 9; jj++){
        grad[jj] += 2 * (sol[2] - data->z_dot) * jac[(jj+3)*12 + 2];
    }
    for (size_t jj = 0; jj < 4; jj++){
        grad[jj+8] += 2 * (sol[2] - data->z_dot) * jac[(jj+12)*12 + 2];
    }
    
    double out_yawdot = pow(sol[11] - data->yaw_dot, 2);
    for (size_t jj = 0; jj < 9; jj++){
        grad[jj] += 2 * (sol[11] - data->yaw_dot) * jac[(jj+3)*12 + 11];
    }
    for (size_t jj = 0; jj < 4; jj++){
        grad[jj+8] += 2 * (sol[11] - data->yaw_dot) * jac[(jj+12)*12 + 11];
    }
    
    double vel = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
    double out_vel = pow(vel - data->target_vel, 2);
    for (size_t jj = 0; jj < 3; jj++){
        grad[jj] += 2 * (vel -data->target_vel) * x[jj] / vel;
    }
    
    double out_sideslip = pow(x[1], 2);
    grad[1] += 2 * x[1];
    
    double out = out_trans + out_rot + out_trim + out_zdot + out_yawdot + out_vel + out_sideslip;

    return out;
}


int steady_state_print(FILE * fp, const struct SteadyState * ss)
{
    fprintf(fp, "\n");
    fprintf(fp, "========================================================\n");
    fprintf(fp, "                        TRIM RESULT                     \n");
    fprintf(fp, "========================================================\n");
    /* fprintf(fp, "Optimizer result = %c\nObjective value = %3.5E\n", nlopt_result_to_string(ss->res), ss->obj_val);*/
    fprintf(fp, "Optimizer result = %d\nObjective value = %3.5E\n", ss->res, ss->obj_val);    

    fprintf(fp, "\n\n\n");
    fprintf(fp, "Spec                :        Targets          Achieved   \n");
    fprintf(fp, "---------------------------------------------------------\n");    
    fprintf(fp, "Speed       (ft/s)  :       %3.5E      %3.5E\n", ss->target_speed, ss->achieved_speed);
    fprintf(fp, "Climb Rate  (ft/s)  :       %3.5E      %3.5E\n", ss->target_climb_rate, ss->achieved_climb_rate);
    fprintf(fp, "Yaw Rate    (rad/s) :       %3.5E      %3.5E\n", ss->target_yaw_rate, ss->achieved_yaw_rate);

        
    fprintf(fp, "\n\n\n");
    fprintf(fp, "State             :       x               dx   \t \tx (secondary unit)\n");
    fprintf(fp, "------------------------------------------------------------------------------\n");
    fprintf(fp, "U     (ft/s)      :  %3.5E     %3.5E\n", ss->UVW.v1, ss->dUVW.v1);
    fprintf(fp, "V     (ft/s)      :  %3.5E     %3.5E\n", ss->UVW.v2, ss->dUVW.v2);
    fprintf(fp, "W     (ft/s)      :  %3.5E     %3.5E\n", ss->UVW.v3, ss->dUVW.v3);
    fprintf(fp, "P     (rad,deg/s) :  %3.5E     %3.5E \t %3.5E\n", ss->PQR.v1, ss->dPQR.v1, rad2deg(ss->PQR.v1));
    fprintf(fp, "Q     (rad,deg/s) :  %3.5E     %3.5E \t %3.5E\n", ss->PQR.v2, ss->dPQR.v2, rad2deg(ss->PQR.v2));
    fprintf(fp, "R     (rad,deg/s) :  %3.5E     %3.5E \t %3.5E\n", ss->PQR.v3, ss->dPQR.v3, rad2deg(ss->PQR.v3));
    fprintf(fp, "Roll  (rad,deg)   :  %3.5E     %3.5E \t %3.5E\n", ss->roll, ss->droll, rad2deg(ss->roll));
    fprintf(fp, "Pitch (rad,deg)   :  %3.5E     %3.5E \t %3.5E\n", ss->pitch, ss->dpitch, rad2deg(ss->pitch));

    fprintf(fp, "\n\n\n");
    fprintf(fp, "Input :\n");
    fprintf(fp, "---------------------------------------------------------------\n");
    fprintf(fp, "Elevator (rad,deg) :   %3.5E \t %3.5E\n", ss->aero_con.v1, rad2deg(ss->aero_con.v1));
    fprintf(fp, "Aileron  (rad,deg) :   %3.5E \t %3.5E\n", ss->aero_con.v2, rad2deg(ss->aero_con.v2));
    fprintf(fp, "Rudder   (rad,deg) :   %3.5E \t %3.5E\n", ss->aero_con.v3, rad2deg(ss->aero_con.v3));
    fprintf(fp, "Thrust   (lb-slug) :   %3.5E\n", ss->thrust);

    fprintf(fp, "\n\n\n");
    fprintf(fp, "Derived Quantities :\n");
    fprintf(fp, "---------------------------------------------------------------\n");
    fprintf(fp, "Angle of Attack   (rad,deg) :   %3.5E \t %3.5E\n", ss->aero.aoa, rad2deg(ss->aero.aoa));
    fprintf(fp, "Sideslip Angle    (rad,deg) :   %3.5E \t %3.5E\n", ss->aero.sideslip, rad2deg(ss->aero.sideslip));
    fprintf(fp, "Flight Path Angle (rad,deg) :   %3.5E \t %3.5E\n", ss->flight_path_angle, rad2deg(ss->flight_path_angle));
    fprintf(fp, "Bank Angle        (rad,deg) :   %3.5E \t %3.5E\n", ss->bank_angle, rad2deg(ss->bank_angle));

    fprintf(fp, "\n");
    fprintf(fp, "===============================================================\n");
    fprintf(fp, "\n");
    
    return 0;
}

int trimmer(struct TrimSpec * data, struct SteadyState * ss){


    // Unknowns are in the following order
    // trim for 8 states  (U, V, W, P, Q, R, roll, pitch)    
    // trim for 4 control inputs (elevator aileron rudder thrust)
    double x[12] = {data->target_vel, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0};
    double val;
    nlopt_result res;
    nlopt_opt opt;

    // Bounds
    double lb[12], ub[12];
    for (size_t ii = 0; ii < 12; ii++){
        lb[ii] = -HUGE_VAL;
        ub[ii] = HUGE_VAL;
    }
    lb[11] = 0.1; // lower bound thrust

    // no sideslip
    lb[1] = 0.0;  // lower bound V
    ub[1] = 0.0;  // upper bound V

    x[1] = 0.0; // initial condition V
    x[11] = 10; // initial condition thrust
    
    // run without bounds
    /* opt = nlopt_create(NLOPT_LN_NELDERMEAD, 12); */
    opt = nlopt_create(NLOPT_LN_NEWUOA, 12);
    /* nlopt_set_ftol_rel(opt, -1.0); */
    /* nlopt_set_ftol_abs(opt, 1e-20); */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, trim_objective, data);
    nlopt_set_stopval(opt, 1e-12);
    
    res = nlopt_optimize(opt, x, &val);
    nlopt_destroy(opt);
    for (size_t ii = 0; ii < 12; ii++){
        if (fabs(x[ii]) < data->thresh){
            x[ii] = 0.0;
        }
    }

    // run with bounds (again)
    /* opt = nlopt_create(NLOPT_LN_NELDERMEAD, 12); */
    /* opt = nlopt_create(NLOPT_LN_SBPLX, 12); */
    opt = nlopt_create(NLOPT_LN_NEWUOA, 12);
    /* opt = nlopt_create(NLOPT_LD_LBFGS, 12); */
    nlopt_set_ftol_rel(opt, -1.0);
    nlopt_set_ftol_abs(opt, 0.0);
    nlopt_set_xtol_rel(opt, 1e-19);
    nlopt_set_xtol_abs1(opt, 1e-19);    
    
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, trim_objective, data);
    /* nlopt_set_min_objective(opt, trim_objective_g, data); */

    res = nlopt_optimize(opt, x, &val);
    nlopt_destroy(opt);

    for (size_t ii = 0; ii < 12; ii++){
        if (fabs(x[ii]) < data->thresh){
            x[ii] = 0.0;
        }
    }
    
    double sol[12];
    double ic[12];
    ic[0] = 0.0; ic[1] = 0.0; ic[2] = 0.0; 
    ic[3] = x[0]; ic[4] = x[1]; ic[5] = x[2]; 
    ic[6] = x[3]; ic[7] = x[4]; ic[8] = x[5]; 
    ic[9] = x[6]; ic[10] = x[7]; ic[11] = 0.0;    
    rigid_body_lin_forces(0.0, ic, x+8, sol, NULL, data->ac);

    ss->res = res;
    ss->obj_val = val;

    ss->UVW.v1 = x[0];
    ss->UVW.v2 = x[1];
    ss->UVW.v3 = x[2];
    ss->dUVW.v1 = sol[3];
    ss->dUVW.v2 = sol[4];
    ss->dUVW.v3 = sol[5];
    

    ss->PQR.v1 = x[3];
    ss->PQR.v2 = x[4];
    ss->PQR.v3 = x[5];
    ss->dPQR.v1 = sol[6];
    ss->dPQR.v2 = sol[7];
    ss->dPQR.v3 = sol[8];    
    
    ss->roll = x[6];
    ss->droll = sol[9];
        
    ss->pitch = x[7];
    ss->dpitch = sol[10];

    ss->aero_con.v1 = x[8];
    ss->aero_con.v2 = x[9];
    ss->aero_con.v3 = x[10];    
    
    ss->thrust = x[11];

    ss->target_climb_rate = trim_spec_get_climb_rate(data);
    ss->achieved_climb_rate = -sol[2];
    
    ss->target_yaw_rate = trim_spec_get_yaw_rate(data);
    ss->achieved_yaw_rate = sol[11];
    
    ss->target_speed = trim_spec_get_speed(data);
    ss->achieved_speed = sqrt(pow(x[0],2) + pow(x[1], 2) + pow(x[2], 2));

    aero_angles(&(ss->UVW), ss->achieved_speed, &(ss->aero));

    ss->flight_path_angle = asin(-sol[2] / ss->achieved_speed);


    // 0 1 2 3 4 5 6    7
    // U V W P Q R Roll Pitch
    real ca = cos(ss->aero.aoa);
    real sa = sin(ss->aero.aoa);
    real num = ca * cos(x[7]) * cos(x[6]) + sa * sin(x[7]);
    real den = cos(ss->flight_path_angle);                                                         

    ss->bank_angle = acos(num / den);
    
    return 0;
}


int controller(double time, const double * x, double * u, void * arg)
{
    (void)(time);
    (void)(x);

    struct SteadyState * ss = arg;

    u[0] = ss->aero_con.v1; // elevator
    u[1] = ss->aero_con.v2; // aileron
    u[2] = ss->aero_con.v3; // rudder
    u[3] = ss->thrust;      // thrust

    // step response    
    /* if (time > 0){ */
    /*     u[0] += -0.2 / 180 * M_PI; */
    /*     /\* u[3] = 5.0; *\/ */
    /* } */
    
    return 0;
}

struct Trajectory *
flight_sim_ss(struct Vec3 * xyz, real yaw, struct SteadyState * ss, struct Aircraft * ac,
              double dt_save, size_t nsteps)
{
    /* double dtmin = 1e-16; */
    /* double dtmax = dt_save; */
    /* double tol = 1e-14; */

    struct Integrator * ode = integrator_create_controlled(12, 4, rigid_body_lin_forces, ac, controller, ss);
    integrator_set_type(ode, "rk4");
    integrator_set_dt(ode, 1e-4);
    /* integrator_set_type(ode,"rkf45");     */
    /* integrator_set_adaptive_opts(ode, dtmin, dtmax, tol); */
    integrator_set_verbose(ode, 0);
    

    double start_time = 0.0;
    double control[4];
    double ic[12];
    ic[0] = xyz->v1;     // x
    ic[1] = xyz->v2;     // y
    ic[2] = xyz->v3;     // z
    ic[3] = ss->UVW.v1;  // U
    ic[4] = ss->UVW.v2;  // V
    ic[5] = ss->UVW.v3;  // W
    ic[6] = ss->PQR.v1;  // P
    ic[7] = ss->PQR.v2;  // Q
    ic[8] = ss->PQR.v3;  // R
    ic[9] = ss->roll;    // roll
    ic[10] = ss->pitch;  // pitch
    ic[11] = yaw;        // yaw    
    controller(0.0, ic, control, ss);

    struct Trajectory * traj = NULL;        
    int res = trajectory_add(&traj, 12, 4, start_time, ic, control);

    double dt = dt_save;
    for (size_t ii = 0; ii < nsteps; ii++){
        res = trajectory_step(traj, ode, dt);
    }
    integrator_destroy(ode);
    return traj;
}

int no_controller(double time, const double * x, double * u, void * arg)
{
 
    (void)(x);
    (void)(arg);

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
    // step response    
    /* if (time > 0){ */
    /*     u[0] = -0.2 / 180 * M_PI; */
    /*     /\* u[3] = 5.0;  *\/ */
    /* } */
    
    return 0;
}

struct Trajectory *
flight_sim_lin(struct Vec3 * xyz, real yaw, struct SteadyState * ss, struct Aircraft * ac,
               double dt_save, size_t nsteps, real * AB)
{
    /* double dtmin = 1e-16; */
    /* double dtmax = dt_save; */
    /* double tol = 1e-14; */

    struct Integrator * ode = integrator_create_controlled(12, 4, rigid_body_linearized, AB, no_controller, NULL);
    integrator_set_type(ode, "rk4");
    integrator_set_dt(ode, 1e-4);
    /* integrator_set_type(ode,"rkf45");     */
    /* integrator_set_adaptive_opts(ode, dtmin, dtmax, tol); */
    integrator_set_verbose(ode, 0);
    

    double start_time = 0.0;
    double control[4];
    double ic[12];
    ic[0] = xyz->v1;     // x
    ic[1] = xyz->v2;     // y
    ic[2] = xyz->v3;     // z
    ic[3] = 0.0; //ss->UVW.v1;  // U
    ic[4] = 0.0; //ss->UVW.v2;  // V
    ic[5] = 0.0; //ss->UVW.v3;  // W  
    ic[6] = 0.0; //ss->PQR.v1;  // P
    ic[7] = 0.0; //ss->PQR.v2;  // Q
    ic[8] = 0.0; //ss->PQR.v3;  // R
    ic[9] = 0.0; //ss->roll;    // roll
    ic[10] = 0.0; //ss->pitch;  // pitch
    ic[11] = yaw;        // yaw    
    no_controller(0.0, ic, control, ss);

    struct Trajectory * traj = NULL;        
    int res = trajectory_add(&traj, 12, 4, start_time, ic, control);

    double dt = dt_save;
    for (size_t ii = 0; ii < nsteps; ii++){
        res = trajectory_step(traj, ode, dt);
    }
    integrator_destroy(ode);
    return traj;
}



int print_A_B(FILE * fp, real jac[192])
{

    fprintf(fp, "========================================================\n");
    fprintf(fp, "                  LINEARIZATION RESULT                  \n");
    fprintf(fp, "========================================================\n");
    fprintf(fp, "\n\n");
    fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n",
            "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw");

    for (size_t ii = 0; ii < 12; ii++){
        for (size_t jj = 0; jj < 12; jj++){
            fprintf(fp, "%-11.3E ", jac[jj*12 + ii]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n\n\n");

    fprintf(fp, "%-9s %-9s %-9s %-9s\n", "Elev.", "Ail.", "Rud.", "Thrust");

    for (size_t ii = 0; ii < 12; ii++){
        for (size_t jj = 0; jj < 4; jj++){
            fprintf(fp, "%-9.3f ", jac[(jj+12)*12 + ii]);
        }
        fprintf(fp, "\n");
    }

    
    return 0;
}


static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "Usage: %s <filename> options \n\n", program_name);
    fprintf(stream,
            " Required Arguments\n"
            " ------------------\n"
            " <filename> must be a json file with the vehicle details.\n"
            " \n\n\n "
            " Optional Arguments\n"
            " ------------------\n"            
            " -h --help                Display this usage information.\n"
            " -s --speed      <val>    Desired speed (e.g., 120 --> flight at 120 ft/s, groundspeed), default 120.\n"
            " -c --climb-rate <val>    Desired climb rate (e.g., -5 --> climb at 5 ft/s), default 0.\n"
            " -y --yaw-rate   <val>    Desired turn rate (e.g., 3.14 --> turn 'right' at pi rad/s), default 0.\n"
            " -t --threshold  <val>    Threshold value for setting states to zero default 1e-10.\n"
            " --linearize              Return linear system\n"
            " --simulate      <file>   Simulate the system and print to file\n"
        );
    exit (exit_code);
}


int main(int argc, char* argv[]){

    int next_option;
    const char * const short_options = "hs:c:y:t:";
    const struct option long_options[] = {
        { "help"       ,  0, NULL, 'h' },
        { "speed"      , 1, NULL, 's' },
        { "climb-rate" , 1, NULL, 'c' },
        { "yaw-rate"   , 1, NULL, 'y' },
        { "threshold"  , 1, NULL, 't' },
        { "linearize"  , 0, NULL, 'l' },
        { "simulate"   , 1, NULL, 1 },        
        /* { "verbose"    , 1, NULL, 'v' }, */
        { NULL         , 0, NULL, 0   }
    };
    
    program_name = argv[0];
    
    real speed = 120.0;
    real climb_rate = 0.0;
    real yaw_rate = 0.0;
    real thresh = 1e-10; // threshold for zero
    int linearize = 0;
    char * sim_name = NULL;
    do {
        next_option = getopt_long (argc, argv, short_options, long_options, NULL);
        switch (next_option)
        {
            case 'h': 
                print_code_usage(stdout, 0);
            case 's':
                speed = atof(optarg);
                break;
            case 'c':
                climb_rate = atof(optarg);
                break;
            case 'y':
                yaw_rate = atof(optarg);
                break;
            case 't':
                thresh = atof(optarg);
                break;
            case 'l':
                linearize = 1;
                break;
            case 1:
                sim_name = optarg;
                break;
            case '?': // The user specified an invalid option
                printf("invalid option %s\n\n",optarg);
                print_code_usage (stderr, 1);
            case -1: // Done with options. 
                break;
            default: // Something unexpected
                abort();
        }

    } while (next_option != -1);

    if (argc -  optind != 1){ // one non-optional argument
        fprintf(stderr, "Incorrect input arguments. Vehicle file is required \n");
        fprintf(stderr, "\n\n\n");
        print_code_usage (stderr, 0);
    }

    // Name of the vehicle file
    char * filename = argv[optind];
    struct Aircraft aircraft;
    aircraft_load(&aircraft, filename);
    
    struct TrimSpec trim_spec;
    trim_spec.z_dot = -climb_rate;
    trim_spec.yaw_dot = yaw_rate; ///3.0 * 2.0 * M_PI / 500.0;
    trim_spec.target_vel = speed; // ft/s
    trim_spec.ac = &aircraft;
    trim_spec.thresh = thresh; 

    struct SteadyState ss;
    trimmer(&trim_spec, &ss);
    steady_state_print(stdout, &ss);

    real jac[144 + 48];   
    if (linearize > 0){

        real ic[12];
        ic[0] = 0.0;
        ic[1] = 0.0;
        ic[2] = 0.0;
        ic[3] = steady_state_get_U(&ss);
        ic[4] = steady_state_get_V(&ss);
        ic[5] = steady_state_get_W(&ss);
        ic[6] = steady_state_get_P(&ss);
        ic[7] = steady_state_get_Q(&ss);
        ic[8] = steady_state_get_R(&ss);
        ic[9] = steady_state_get_Roll(&ss);
        ic[10] = steady_state_get_Pitch(&ss);
        ic[11] = 0.0;//steady_state_get_Yaw(&ss);

        real control[4];
        control[0] = steady_state_get_elevator(&ss);
        control[1] = steady_state_get_aileron(&ss);
        control[2] = steady_state_get_rudder(&ss);
        control[3] = steady_state_get_thrust(&ss);

        real rhs[12];
        int res = rigid_body_lin_forces_jac(0.0, ic, control, rhs, jac, &aircraft);
        assert(res == 0);

        print_A_B(stdout, jac);
        
    }

    if (sim_name != NULL){

        fprintf(stdout, "========================================================\n");
        fprintf(stdout, "      Simulating Nonlinear Dynamics at Steady State     \n");
        fprintf(stdout, "========================================================\n");
        struct Vec3 xyz = {0, 0, -5};
        real yaw = M_PI / 4.0;
        double dt_save = 1e-1;
        size_t nsteps = 1000;
        struct Trajectory * traj = flight_sim_ss(&xyz, yaw, &ss, &aircraft, dt_save, nsteps);

        char filename[256];
        sprintf(filename, "nrb_%s",sim_name);
        printf("Saving simulation to %s\n", filename);
        FILE * fp = fopen(filename, "w");

        if (fp == NULL){
            fprintf(stdout, "Cannot open file %s\n", filename);
            return 1;
        }
        fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n",
                "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw");        
        trajectory_print(traj, fp, 10);
        trajectory_free(traj); traj = NULL;
        fclose(fp);
        printf("\n\n");


        if (linearize){
            fprintf(stdout, "========================================================\n");
            fprintf(stdout, "     Simulating Linearized Dynamics at Steady State     \n");
            fprintf(stdout, "========================================================\n");


            sprintf(filename, "lrb_%s", sim_name);
            printf("Saving simulation to %s\n", filename);
            FILE * fp = fopen(filename, "w");
            
            if (fp == NULL){
                fprintf(stdout, "Cannot open file %s\n", filename);
                return 1;                
            }

            struct Vec3 xyz_ss = {0, 0, -5};
            real yaw_ss = M_PI/4.0;
            struct Trajectory * traj_ss = flight_sim_ss(&xyz_ss, yaw_ss, &ss, &aircraft, dt_save, nsteps);
        
            
            struct Vec3 xyz_perturbed = {0, 0, 0};
            real yaw_perturbed = 0.0;
            struct Trajectory * traj_lin = flight_sim_lin(&xyz_perturbed, yaw_perturbed, &ss, &aircraft,
                                  dt_save, nsteps, jac);
            fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n",
                    "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw");
            

            trajectory_print(traj_lin, fp, 10);
            trajectory_print(traj_ss, fp, 10);            
            trajectory_free(traj_lin); traj_lin = NULL;
            trajectory_free(traj_ss); traj_ss = NULL;
            
            fclose(fp);
        }
        
    }

    return 0;
}

