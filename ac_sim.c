/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

//------------------------------------------
//------------------------------------------
//---  Simulate Nonlinear Rigid Body  ------
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


int load_ic(char * filename, real * ic, real * ic_control);
struct Trajectory * flight_sim(real * ic, real * control_fixed, struct Aircraft * ac,
                               double dt_save, size_t nsteps);


static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "\nUsage: %s <vehicle_filename> <initial_condition_filename> options \n\n", program_name);
    fprintf(stream,
            "Simulate a 6DOF aircraft\n "

            "\nSample call\n"
            "./ac_sim pioneer_uav.json ic.json\n"
            
            "\n\n\n"
            " Optional Arguments\n"
            " ------------------\n"
            " -h --help                Display this usage information.\n"
            " -t --time                Length of integration period (in seconds) (default 10)\n"
            " -d --dtsave              Size of timestep to save (default 1e-2)\n"
            /* " --simulate      <file>   Simulate the system and print to file\n" */
            "\n\n"
        );
    exit (exit_code);
}


int main(int argc, char* argv[]){

    int next_option;
    const char * const short_options = "ht:d:";
    const struct option long_options[] = {
        { "help"       ,  0, NULL, 'h' },
        { "time"       ,  1, NULL, 't' },
        { "dtsave"     ,  1, NULL, 'd' },        
        { NULL         , 0, NULL, 0   }
    };
    
    program_name = argv[0];

    real dt_save = 1e-2;
    real T = 10.0;
    int num_opts = 0;
    do {
        next_option = getopt_long (argc, argv, short_options, long_options, NULL);
        num_opts++;
        switch (next_option)
        {
            case 'h': 
                print_code_usage(stdout, 0);
            case 'd':
                dt_save = atof(optarg);
                break;
            case 't':
                T = atof(optarg);
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


    if (argc - optind != 2){ 
        fprintf(stderr, "Called as: %s ", argv[0]);
        for (size_t ii = 1; ii < argc; ii++){
            fprintf(stderr, " %s ", argv[ii]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "\n\n\n");

        print_code_usage (stderr, 0);        
    }
    
    // Name of the vehicle file
    char * ac_filename = argv[optind];
    char * ic_filename = argv[optind+1];

    printf("Aircraft filename = %s\n", ac_filename);
    printf("Initial condition filename = %s\n", ic_filename);    

    real ic[12];
    real control_fixed[4];

    int ret = load_ic(ic_filename, ic, control_fixed);
    if (ret == 1){
        fprintf(stderr, "Could not load initial condition\n");
        return 1;
    }

    struct Aircraft aircraft;
    ret = aircraft_load(&aircraft, ac_filename);
    if (ret == 1){
        fprintf(stderr, "Could not load aircraft file\n");
        return 1;
    }
    
    fprintf(stdout, "========================================================\n");
    fprintf(stdout, "             Simulating Nonlinear Dynamics              \n");
    fprintf(stdout, "========================================================\n");


    size_t nsteps = (size_t) floor(T / dt_save);
    printf("Time = %3.4E\n", T);
    printf("save_time = %3.4E\n", dt_save);
    printf("nsteps = %zu\n", nsteps );
    
    
    struct Trajectory * traj = flight_sim(ic, control_fixed, &aircraft, dt_save, nsteps);

    char filename[256];
    sprintf(filename, "%s.run", ic_filename);
    printf("Saving simulation to %s\n", filename);
    FILE * fp = fopen(filename, "w");

    if (fp == NULL){
        fprintf(stdout, "Cannot open file %s\n", filename);
        return 1;
    }
    fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n",
            "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw");
    trajectory_print(traj, fp, 20);

    trajectory_free(traj); traj = NULL;
    fclose(fp);
    printf("\n\n");



    return 0;
}

int controller(double time, const double * x, double * u, void * arg)
{
    (void)(time);
    (void)(x);

    real * u_fixed = arg;

    u[0] = u_fixed[0]; // elevator
    u[1] = u_fixed[1]; // aileron
    u[2] = u_fixed[2]; // rudder
    u[3] = u_fixed[3]; // thrust

    // step response    
    /* if (time > 0){ */
    /*     u[0] += -0.2 / 180 * M_PI; */
    /*     /\* u[3] = 5.0; *\/ */
    /* } */
    
    return 0;
}

struct Trajectory * flight_sim(real * ic, real * control_use, struct Aircraft * ac,
                               double dt_save, size_t nsteps)
{
    /* double dtmin = 1e-16; */
    /* double dtmax = dt_save; */
    /* double tol = 1e-14; */

    struct Integrator * ode = integrator_create_controlled(12, 4, rigid_body_lin_forces, ac,
                                                           controller,
                                                           control_use);
    integrator_set_type(ode, "rk4");
    integrator_set_dt(ode, 1e-4);
    /* integrator_set_type(ode,"rkf45");     */
    /* integrator_set_adaptive_opts(ode, dtmin, dtmax, tol); */
    integrator_set_verbose(ode, 0);
    

    double start_time = 0.0;
    real control[4];
    controller(0.0, ic, control, control_use);

    struct Trajectory * traj = NULL;        
    int res = trajectory_add(&traj, 12, 4, start_time, ic, control);

    double dt = dt_save;
    for (size_t ii = 0; ii < nsteps; ii++){
        res = trajectory_step(traj, ode, dt);
    }
    integrator_destroy(ode);
    return traj;
}

#define JSMN_HEADER
#include "jsmn/jsmn.h"

static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
    if (tok->type == JSMN_STRING && (int)strlen(s) == tok->end - tok->start &&
        strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
        return 0;
    }
    return -1;
}


int load_ic(char * filename, real * ic, real * ic_control)
{

    FILE * fp = fopen(filename, "r");
    if (fp == NULL){
        fprintf(stderr, "Cannot open file %s\n", filename);
    }
    
    // Determine file size
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);

    /* printf("size = %zu\n", size); */
    char* input = malloc(size * sizeof(char));
    
    rewind(fp);
    fread(input, sizeof(char), size, fp);

    /* printf("input = %s\n", input); */
    
    jsmn_parser p;
    jsmntok_t t[256];
    jsmn_init(&p);
    int r = jsmn_parse(&p, input, strlen(input), t, 256);

    if (r < 0) {
        printf("Failed to parse JSON: %d\n", r);
        printf("input = %s\n", input);
        return 1;
    }

    /* Assume the top-level element is an object */
    if (r < 1 || t[0].type != JSMN_OBJECT) {
        printf("Object expected\n");
        return 1;
    }

    /* Loop over all keys of the root object */
    for (size_t i = 1; i < r; i++) {
        /* printf("i = %zu, r = %d \n", i, r); */
           
        if (jsoneq(input, &t[i], "x") == 0) {
            ic[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "y") == 0) {
            ic[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "z") == 0) {
            ic[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "U") == 0) {
            ic[3] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "V") == 0) {
            ic[4] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "W") == 0) {
            ic[5] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "P") == 0) {
            ic[6] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Q") == 0) {
            ic[7] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "R") == 0) {
            ic[8] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "roll") == 0) {
            ic[9] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "pitch") == 0) {
            ic[10] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "yaw") == 0) {
            ic[11] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "elevator") == 0) {
            ic_control[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "aileron") == 0) {
            ic_control[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "rudder") == 0) {
            ic_control[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "thrust") == 0) {
            ic_control[3] = atof(input + t[i + 1].start);
            i++;
        }        
    }
    
    fclose(fp);
    free(input);  input = NULL;

    return 0;
}

    
