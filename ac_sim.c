/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

//------------------------------------------
//------------------------------------------
//---- Simulate nonlin and linearized  -----
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


struct Trajectory * flight_sim_ss(struct Vec3 * xyz, real yaw, struct SteadyState * ss,
                                  struct Aircraft * ac, double dt_save, size_t nsteps);
/* struct Trajectory * flight_sim_lin(struct Vec3 * xyz, real yaw, struct SteadyState * ss, */
/*                                    struct Aircraft * ac, double dt_save, size_t nsteps, real * AB); */


static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "Usage: %s <vehicle_filename> <initial_condition_filename> \n\n", program_name);
    fprintf(stream,
            "Simulate a 6DOF aircraft\n "

            "\nSample call\n"
            "./ac_sim pioneer_uav.json ic.json\n"
            
            "\n\n\n\n\n\n"
            " Required Arguments\n"
            " ------------------\n"
            " <filename> must be a json file with the vehicle details.\n"
            " \n\n\n "
            " Optional Arguments\n"
            " ------------------\n"            
            " -h --help                Display this usage information.\n"
            /* " --simulate      <file>   Simulate the system and print to file\n" */
        );
    exit (exit_code);
}


int main(int argc, char* argv[]){

    int next_option;
    const char * const short_options = "h:";
    const struct option long_options[] = {
        { "help"       ,  0, NULL, 'h' },
        { NULL         , 0, NULL, 0   }
    };
    
    program_name = argv[0];
    
    do {
        next_option = getopt_long (argc, argv, short_options, long_options, NULL);
        switch (next_option)
        {
            case 'h': 
                print_code_usage(stdout, 0);
            case '?': // The user specified an invalid option
                printf("invalid option %s\n\n",optarg);
                print_code_usage (stderr, 1);
            case -1: // Done with options. 
                break;
            default: // Something unexpected
                abort();
        }

    } while (next_option != -1);

    if (argc - optind != 2){ //two non-optional argument
        fprintf(stderr, "Incorrect input arguments. Vehicle and initial condition files are required \n");
        fprintf(stderr, "\n\n\n");
        print_code_usage (stderr, 0);
    }

    // Name of the vehicle file
    char * filename = argv[optind];
    struct Aircraft aircraft;
    aircraft_load(&aircraft, filename);
    

    fprintf(stdout, "========================================================\n");
    fprintf(stdout, "      Simulating Nonlinear Dynamics at Steady State     \n");
    fprintf(stdout, "========================================================\n");
    /* struct Vec3 xyz = {0, 0, -5}; */
    /* real yaw = M_PI / 4.0; */
    /* double dt_save = 1e-1; */
    /* size_t nsteps = 1000; */
    /* struct Trajectory * traj = flight_sim_ss(&xyz, yaw, &ss, &aircraft, dt_save, nsteps); */

    /* char filename[256]; */
    /* sprintf(filename, "nrb_%s",sim_name); */
    /* printf("Saving simulation to %s\n", filename); */
    /* FILE * fp = fopen(filename, "w"); */

    /* if (fp == NULL){ */
    /*     fprintf(stdout, "Cannot open file %s\n", filename); */
    /*     return 1; */
    /* } */
    /* fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n", */
    /*         "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw"); */
    /* trajectory_print(traj, fp, 10); */
    /* trajectory_free(traj); traj = NULL; */
    /* fclose(fp); */
    /* printf("\n\n"); */


    /*     if (linearize){ */
    /*         fprintf(stdout, "========================================================\n"); */
    /*         fprintf(stdout, "     Simulating Linearized Dynamics at Steady State     \n"); */
    /*         fprintf(stdout, "========================================================\n"); */


    /*         sprintf(filename, "lABmat_%s", sim_name); */
    /*         printf("Saving AB mat to %s\n", filename); */
    /*         fp = fopen(filename, "w"); */

    /*         fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n", */
    /*                 "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw", "Elevator", "Aileron", "Rudder", "Thrust"); */
    /*         for (size_t ii = 0; ii < 12; ii++){ */
    /*             for (size_t jj = 0; jj < 16; jj++){ */
    /*                 fprintf(fp, "%3.15f ", jac[jj*12 + ii]); */
    /*             } */
    /*             fprintf(fp, "\n"); */
    /*         } */
    /*         fclose(fp); */
            
            
    /*         sprintf(filename, "lrb_%s", sim_name); */
    /*         printf("Saving simulation to %s\n", filename); */
    /*         FILE * fp = fopen(filename, "w"); */
            
    /*         if (fp == NULL){ */
    /*             fprintf(stdout, "Cannot open file %s\n", filename); */
    /*             return 1; */
    /*         } */

    /*         struct Vec3 xyz_ss = {0, 0, -5}; */
    /*         real yaw_ss = M_PI/4.0; */
    /*         struct Trajectory * traj_ss = flight_sim_ss(&xyz_ss, yaw_ss, &ss, &aircraft, dt_save, nsteps); */
        
            
    /*         struct Vec3 xyz_perturbed = {0, 0, 0}; */
    /*         real yaw_perturbed = 0.0; */
    /*         struct Trajectory * traj_lin = flight_sim_lin(&xyz_perturbed, yaw_perturbed, &ss, &aircraft, */
    /*                               dt_save, nsteps, jac); */
    /*         fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n", */
    /*                 "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw"); */
            

    /*         trajectory_print(traj_lin, fp, 10); */
    /*         trajectory_print(traj_ss, fp, 10); */
    /*         trajectory_free(traj_lin); traj_lin = NULL; */
    /*         trajectory_free(traj_ss); traj_ss = NULL; */
            
    /*         fclose(fp); */
    /*     } */
        
    /* } */

    return 0;
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
    double ic[12];    
    double control[4];
    steady_state_set_vec(ss, xyz->v1, xyz->v2, xyz->v3, yaw, ic, control);
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




    
/* struct Trajectory * */
/* flight_sim_lin(struct Vec3 * xyz, real yaw, struct SteadyState * ss, struct Aircraft * ac, */
/*                double dt_save, size_t nsteps, real * AB) */
/* { */
/*     /\* double dtmin = 1e-16; *\/ */
/*     /\* double dtmax = dt_save; *\/ */
/*     /\* double tol = 1e-14; *\/ */

/*     struct Integrator * ode = integrator_create_controlled(12, 4, rigid_body_linearized, AB, no_controller, NULL); */
/*     integrator_set_type(ode, "rk4"); */
/*     integrator_set_dt(ode, 1e-4); */
/*     /\* integrator_set_type(ode,"rkf45");     *\/ */
/*     /\* integrator_set_adaptive_opts(ode, dtmin, dtmax, tol); *\/ */
/*     integrator_set_verbose(ode, 0); */
    

/*     double start_time = 0.0; */
/*     double ic[12];     */
/*     double control[4]; */
/*     steady_state_set_vec(ss, xyz->v1, xyz->v2, xyz->v3, yaw, ic, control);     */
/*     no_controller(0.0, ic, control, ss); */

/*     struct Trajectory * traj = NULL;         */
/*     int res = trajectory_add(&traj, 12, 4, start_time, ic, control); */

/*     double dt = dt_save; */
/*     for (size_t ii = 0; ii < nsteps; ii++){ */
/*         res = trajectory_step(traj, ode, dt); */
/*     } */
/*     integrator_destroy(ode); */
/*     return traj; */
/* } */
