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

#include "vehicles.h"
#include "trimming.h"


int steady_state_print(FILE * fp, const struct SteadyState * ss);
int steady_state_print_to_json(FILE * fp, const struct SteadyState * ss, const real * jac);
int print_A_B(FILE * fp, real jac[192]);
    
static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "\nUsage: %s <filename> options \n\n", program_name);
    fprintf(stream,
            "This script trims a six degree of freedom aircraft model to a target speed, climb-rate, and yaw-rate. The script can also return the linear system obtained around the trim condition. Results are printed to the screen by default. However, the trim condition can be written to a json file by specifying an output filename\n "
            "\n\n\n\n\n\n"
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
            " -o --output     <file>   Filename for json to output.\n"
            " --linearize              Return linear system\n"
            "\n\n"
        );
    exit (exit_code);
}

int main(int argc, char* argv[]){

    int next_option;
    const char * const short_options = "hs:c:y:o:t:";
    const struct option long_options[] = {
        { "help"       ,  0, NULL, 'h' },
        { "speed"      , 1, NULL, 's' },
        { "climb-rate" , 1, NULL, 'c' },
        { "yaw-rate"   , 1, NULL, 'y' },
        { "threshold"  , 1, NULL, 't' },
        { "output"  , 1, NULL, 'o' },
        { "linearize"  , 0, NULL, 'l' },
        /* { "verbose"    , 1, NULL, 'v' }, */
        { NULL         , 0, NULL, 0   }
    };
    
    program_name = argv[0];
    
    real speed = 120.0;
    real climb_rate = 0.0;
    real yaw_rate = 0.0;
    real thresh = 1e-10; // threshold for zero
    int linearize = 0;
    /* char * sim_name = NULL; */
    char * output_name = NULL;
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
            case 'o':
                output_name = optarg;
                break;                
            /* case 1: */
            /*     sim_name = optarg; */
            /*     break; */
            case '?': // The user specified an invalid option
                printf("invalid option %s\n\n",optarg);
                print_code_usage (stderr, 1);
            case -1: // Done with options. 
                break;
            default: // Something unexpected
                abort();
        }

    } while (next_option != -1);


    if (argc - optind != 1){ //three non-optional argument
        fprintf(stderr, "Called as: %s ", argv[0]);
        for (size_t ii = 1; ii < argc; ii++){
            fprintf(stderr, " %s ", argv[ii]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "\n\n\n");
        print_code_usage (stderr, 0);        
    }    

    // Name of the vehicle file
    char * filename = argv[optind];
    struct Aircraft aircraft;
    int ret = aircraft_load(&aircraft, filename);
    if (ret == 1){
        fprintf(stderr, "Could not load aircraft file\n");
        return 1;
    }    
    
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
        real control[4];
        real rhs[12];
        
        steady_state_set_vec(&ss, 0.0, 0.0, 0.0, 0.0, ic, control);
        int res = rigid_body_lin_forces_jac(0.0, ic, control, rhs, jac, &aircraft);
        assert(res == 0);

        print_A_B(stdout, jac);
        
    }

    if (output_name != NULL){
        /* char filename[256]; */
        /* sprintf(filename, "nrb_%s",sim_name); */
        printf("Saving to %s\n", output_name);
        FILE * fp = fopen(output_name, "w");
        if (fp == NULL){
            fprintf(stdout, "Cannot open file %s\n", filename);
            return 1;
        }

        if (linearize == 0){
            steady_state_print_to_json(fp, &ss, NULL);
        }
        else{
            steady_state_print_to_json(fp, &ss, jac);
        }

        fclose(fp);
        
        /* struct SteadyState ss_check; */
        /* steady_state_load(output_name, &ss_check);         */
        /* steady_state_print(stdout, &ss_check); */

        /* real jac_check[144 + 48]; */
        /* steady_state_load_jac(output_name, jac_check); */
        /* print_A_B(stdout, jac); */
        
    }

    return 0;
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


int steady_state_print_to_json(FILE * fp, const struct SteadyState * ss, const real * jac)
{
    // Should add accessors for all elements
    fprintf(fp, "{\n");
    
    fprintf(fp, "\"U\":%3.15f,\n", ss->UVW.v1);
    fprintf(fp, "\"V\":%3.15f,\n", ss->UVW.v2);
    fprintf(fp, "\"W\":%3.15f,\n", ss->UVW.v3);

    fprintf(fp, "\"P\":%3.15f,\n", ss->PQR.v1);
    fprintf(fp, "\"Q\":%3.15f,\n", ss->PQR.v2);    
    fprintf(fp, "\"R\":%3.15f,\n", ss->PQR.v3);

    fprintf(fp, "\"roll\":%3.15f,\n", ss->roll);
    fprintf(fp, "\"pitch\":%3.15f,\n", ss->pitch);

    fprintf(fp, "\"elevator\":%3.15f,\n", steady_state_get_elevator(ss));
    fprintf(fp, "\"aileron\":%3.15f,\n",  steady_state_get_aileron(ss));
    fprintf(fp, "\"rudder\":%3.15f,\n",   steady_state_get_rudder(ss));
    fprintf(fp, "\"thrust\":%3.15f,\n",   steady_state_get_thrust(ss));
    
    
    fprintf(fp, "\"Target_Speed\":%3.15f,\n",        ss->target_speed);
    fprintf(fp, "\"Achieved_Speed\":%3.15f,\n",      ss->achieved_speed);
    fprintf(fp, "\"Target_Climb_Rate\":%3.15f,\n",   ss->target_climb_rate);
    fprintf(fp, "\"Achieved_Climb_Rate\":%3.15f,\n", ss->achieved_climb_rate);    
    fprintf(fp, "\"Target_Yaw_Rate\":%3.15f,\n",     ss->target_yaw_rate);
    fprintf(fp, "\"Achieved_Climb_Rate\":%3.15f,\n", ss->achieved_yaw_rate);


    fprintf(fp, "\"Angle_of_Attack\":%3.15f,\n",     ss->aoa);
    fprintf(fp, "\"Sideslip_Angle\":%3.15f,\n",      ss->sideslip);
    fprintf(fp, "\"flight_path_angle\":%3.15f,\n",   ss->flight_path_angle);
    fprintf(fp, "\"bank_angle\":%3.15f,\n",          ss->bank_angle);

    if (jac != NULL){
        fprintf(fp, "\"A\":[");
        for (size_t ii = 0; ii < 143; ii++){ // A matrix column order followed by B matrix column order
            fprintf(fp, "%3.15f", jac[ii]);
            fprintf(fp, ",");
        }
        fprintf(fp, "%3.15f", jac[143]);
        fprintf(fp, "],\n");

        fprintf(fp, "\"B\":[");
        for (size_t ii = 144; ii < 144 + 47; ii++){
            fprintf(fp, "%3.15f", jac[ii]);
            fprintf(fp, ",");
        }
        fprintf(fp, "%3.15f", jac[144+47]);        
        fprintf(fp, "]\n");
    }
    fprintf(fp, "}\n");

    return 0;
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
