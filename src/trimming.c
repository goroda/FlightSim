/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>


#include "trimming.h"


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

    return out;
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

#include "jsmn/jsmn.h"

static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
    if (tok->type == JSMN_STRING && (int)strlen(s) == tok->end - tok->start &&
        strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
        return 0;
    }
    return -1;
}


int steady_state_load(char * filename, struct SteadyState * ss)
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

    jsmn_parser p;
    jsmntok_t t[10000];
    int r = jsmn_parse(&p, input, strlen(input), t, 10000);

    if (r < 0) {
        printf("Failed to parse JSON: %d\n", r);
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
           
        if (jsoneq(input, &t[i], "U") == 0) {
            ss->UVW.v1 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "V") == 0) {
            ss->UVW.v2 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "W") == 0) {
            ss->UVW.v3 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "P") == 0) {
            ss->PQR.v1 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Q") == 0) {
            ss->PQR.v2 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "R") == 0) {
            ss->PQR.v3 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "elevator") == 0) {
            ss->aero_con.v1 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "aileron") == 0) {
            ss->aero_con.v2 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "rudder") == 0) {
            ss->aero_con.v3 = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "thrust") == 0) {
            ss->thrust = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "roll") == 0) {
            ss->roll = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "pitch") == 0) {
            ss->pitch = atof(input + t[i + 1].start);
            i++;
        }  
        /* else { */
        /*     printf("Unused key: %.*s\n", t[i].end - t[i].start, */
        /*            input + t[i].start); */
        /* } */
    }
    
    fclose(fp);
    free(input); 

    return 0;
}

int steady_state_load_jac(char * filename, real * jac)
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

    jsmn_parser p;
    jsmntok_t t[10000];
    int r = jsmn_parse(&p, input, strlen(input), t, 10000);

    if (r < 0) {
        printf("Failed to parse JSON: %d\n", r);
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
        
        if (jsoneq(input, &t[i], "A") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            for (size_t kk = 0; kk < 144; kk++){
                jsmntok_t * aa = &t[i + kk + 2];
                jac[kk] = atof(input + aa->start);
            }
            i += t[i + 1].size + 1;        
        }
        else if (jsoneq(input, &t[i], "A") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            for (size_t kk = 144; kk < 144+48; kk++){
                jsmntok_t * aa = &t[i + kk + 2];
                jac[kk] = atof(input + aa->start);
            }
            i += t[i + 1].size + 1;        
        }

    }

    return 0;
}


//////////////////////////////////
// Old -- not used
//////////////////////////////////
/* double trim_objective_g(unsigned n, const double * x, double * grad, void * f_data) */
/* { */
/*     if (grad == NULL){ */
/*         return trim_objective(n, x, grad, f_data); */
/*     } */
    
/*     struct TrimSpec * data = f_data; */

/*     double ic[12]; */
/*     ic[0] = 0.0; // x */
/*     ic[1] = 0.0; // y */
/*     ic[2] = 0.0; // z */
/*     ic[3] = x[0]; // u */
/*     ic[4] = x[1]; // v */
/*     ic[5] = x[2]; // w */
/*     ic[6] = x[3]; // p */
/*     ic[7] = x[4]; // q */
/*     ic[8] = x[5]; // r */
/*     ic[9] = x[6]; // roll */
/*     ic[10] = x[7]; // pitch */
/*     ic[11] = 0.0; //yaw */

/*     double control[4]; */
/*     control[0] = x[8]; */
/*     control[1] = x[9]; */
/*     control[2] = x[10]; */
/*     control[3] = x[11]; */
    
/*     double sol[12]; */
/*     double jac[144 + 48]; // states and controls */
/*     rigid_body_lin_forces_jac(0.0, ic, control, sol, jac, data->ac); */

/*     double out_trans = pow(sol[3], 2) + pow(sol[4], 2) + pow(sol[5], 2); */
/*     // U, V, W equations */
/*     for (size_t jj = 0; jj < 9; jj++){ */
/*         grad[jj] = 2 * sol[3] * jac[(jj+3)*12 + 3] + */
/*                    2 * sol[4] * jac[(jj+3)*12 + 4] + */
/*                    2 * sol[5] * jac[(jj+3)*12 + 5]; */
/*     } */
/*     for (size_t jj = 0; jj < 4; jj++){ */
/*         grad[jj+8] += 2 * sol[3] * jac[(12+jj)*12 + 3] + */
/*                       2 * sol[4] * jac[(12+jj)*12 + 4] + */
/*                       2 * sol[5] * jac[(12+jj)*12 + 5]; */
/*     } */
    
/*     double out_rot = pow(sol[6], 2) + pow(sol[7], 2) + pow(sol[8], 2); */
/*     // P, Q, R equations */
/*     for (size_t jj = 0; jj < 9; jj++){ */
/*         grad[jj] += (2 * sol[6] * jac[(jj+3)*12 + 6] + */
/*                      2 * sol[7] * jac[(jj+3)*12 + 7] + */
/*                      2 * sol[8] * jac[(jj+3)*12 + 8]); */
/*     } */

/*     for (size_t jj = 0; jj < 4; jj++){ */
/*         grad[jj+8] += (2 * sol[6] * jac[(12+jj)*12 + 6] + */
/*                        2 * sol[7] * jac[(12+jj)*12 + 7] + */
/*                        2 * sol[8] * jac[(12+jj)*12 + 8]); */
/*     }     */
    
/*     double out_trim = pow(sol[9], 2) + pow(sol[10], 2); */
/*     // Roll and pitch equations */
/*     for (size_t jj = 0; jj < 9; jj++){ */
/*         grad[jj] += 2 * sol[9] * jac[(jj+3)*12 + 9] + */
/*                     2 * sol[10] * jac[(jj+3)*12 + 10]; */
/*     } */
/*     for (size_t jj = 0; jj < 4; jj++){ */
/*         grad[jj+8] += 2 * sol[9] * jac[(jj+12)*12 + 9] + */
/*                       2 * sol[10] * jac[(jj+12)*12 + 10]; */
/*     } */
    
/*     double out_zdot = pow(sol[2] - data->z_dot, 2); */
/*     for (size_t jj = 0; jj < 9; jj++){ */
/*         grad[jj] += 2 * (sol[2] - data->z_dot) * jac[(jj+3)*12 + 2]; */
/*     } */
/*     for (size_t jj = 0; jj < 4; jj++){ */
/*         grad[jj+8] += 2 * (sol[2] - data->z_dot) * jac[(jj+12)*12 + 2]; */
/*     } */
    
/*     double out_yawdot = pow(sol[11] - data->yaw_dot, 2); */
/*     for (size_t jj = 0; jj < 9; jj++){ */
/*         grad[jj] += 2 * (sol[11] - data->yaw_dot) * jac[(jj+3)*12 + 11]; */
/*     } */
/*     for (size_t jj = 0; jj < 4; jj++){ */
/*         grad[jj+8] += 2 * (sol[11] - data->yaw_dot) * jac[(jj+12)*12 + 11]; */
/*     } */
    
/*     double vel = sqrt(pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2)); */
/*     double out_vel = pow(vel - data->target_vel, 2); */
/*     for (size_t jj = 0; jj < 3; jj++){ */
/*         grad[jj] += 2 * (vel -data->target_vel) * x[jj] / vel; */
/*     } */
    
/*     double out_sideslip = pow(x[1], 2); */
/*     grad[1] += 2 * x[1]; */
    
/*     double out = out_trans + out_rot + out_trim + out_zdot + out_yawdot + out_vel + out_sideslip; */

/*     return out; */
/* } */
