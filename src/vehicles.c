/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "vehicles.h"

#define JSMN_HEADER
#include "jsmn/jsmn.h"

static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
    if (tok->type == JSMN_STRING && (int)strlen(s) == tok->end - tok->start &&
        strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
        return 0;
    }
    return -1;
}


int aircraft_load(struct Aircraft * ac, char * filename)
{

    FILE * fp = fopen(filename, "r");
    if (fp == NULL){
        fprintf(stderr, "Cannot open file %s\n", filename);
    }
    aircraft_init(ac);
    
    // Determine file size
    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);

    /* printf("size = %zu\n", size); */
    char* input = malloc(size * sizeof(char));

    rewind(fp);
    fread(input, sizeof(char), size, fp);
    /* printf("input = %s\n", input); */
    
    jsmn_parser p;
    jsmntok_t t[128];
    jsmn_init(&p);
    int r = jsmn_parse(&p, input, strlen(input), t, 128);

    if (r < 0) {
        printf("Failed to parse JSON: %d\n", r);
        fprintf(stderr, "Could not parse aircraft file\n");
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
           
        if (jsoneq(input, &t[i], "mass") == 0) {
            ac->m = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "area") == 0) {
            ac->area = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "span") == 0) {
            ac->span = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "chord") == 0) {
            ac->chord = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "mac") == 0) {
            ac->mac = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "AR") == 0) {
            ac->AR = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "e") == 0) {
            ac->e = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "K") == 0) {
            ac->K= atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Ixx") == 0) {
            ac->Ixx = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Ixz") == 0) {
            ac->Ixz = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Iyy") == 0) {
            ac->Iyy = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Izz") == 0) {
            ac->Izz = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "engine_angle") == 0) {
            real ang = atof(input + t[i + 1].start);
            ac->cphit = cos(ang);
            ac->sphit = sin(ang);
            i++;
        }
        else if (jsoneq(input, &t[i], "CD_parasitic") == 0) {
            ac->CD[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CD_alpha") == 0) {
            ac->CD[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CD_mach") == 0) {
            ac->CD[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CL_0") == 0) {
            ac->CL[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CL_alpha") == 0) {
            ac->CL[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CL_alpha_dot") == 0) {
            ac->CL[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CL_mach") == 0) {
            ac->CL[3] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CL_q") == 0) {
            ac->CL[4] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CL_elevator") == 0) {
            ac->CL[5] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CE_beta") == 0) {
            ac->CE[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CE_p") == 0) {
            ac->CE[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CE_rudder") == 0) {
            ac->CE[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CM_0") == 0) {
            ac->Cm[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CM_alpha") == 0) {
            ac->Cm[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CM_alpha_dot") == 0) {
            ac->Cm[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CM_mach") == 0) {
            ac->Cm[3] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CM_q") == 0) {
            ac->Cm[4] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "CM_elevator") == 0) {
            ac->Cm[5] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cl_beta") == 0) {
            ac->Cl[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cl_p") == 0) {
            ac->Cl[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cl_r") == 0) {
            ac->Cl[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cl_aileron") == 0) {
            ac->Cl[3] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cl_rudder") == 0) {
            ac->Cl[4] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cn_beta") == 0) {
            ac->Cn[0] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cn_p") == 0) {
            ac->Cn[1] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cn_r") == 0) {
            ac->Cn[2] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cn_aileron") == 0) {
            ac->Cn[3] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "Cn_rudder") == 0) {
            ac->Cn[4] = atof(input + t[i + 1].start);
            i++;
        }
        else if (jsoneq(input, &t[i], "alfa_limits") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            jsmntok_t * lb = &t[i + 2];
            ac->aoa_bounds[0] = atof(input + lb->start);
            jsmntok_t * ub = &t[i + 3];
            ac->aoa_bounds[1] = atof(input + ub->start); 
            /* printf("ac->bounds = %3.5E, %3.5E\n", ac->aoa_bounds[0], ac->aoa_bounds[1]); */
            i += t[i + 1].size + 1;
        }
        else if (jsoneq(input, &t[i], "beta_limits") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            jsmntok_t * lb = &t[i + 2];
            ac->sideslip_bounds[0] = atof(input + lb->start);
            jsmntok_t * ub = &t[i + 3];
            ac->sideslip_bounds[1] = atof(input + ub->start); 
            i += t[i + 1].size + 1;
        }
        else if (jsoneq(input, &t[i], "elevator_limits") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            jsmntok_t * lb = &t[i + 2];
            ac->elevator_bounds[0] = atof(input + lb->start);
            jsmntok_t * ub = &t[i + 3];
            ac->elevator_bounds[1] = atof(input + ub->start); 
            i += t[i + 1].size + 1;
        }
        else if (jsoneq(input, &t[i], "aileron_limits") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            jsmntok_t * lb = &t[i + 2];
            ac->aileron_bounds[0] = atof(input + lb->start);
            jsmntok_t * ub = &t[i + 3];
            ac->aileron_bounds[1] = atof(input + ub->start); 
            i += t[i + 1].size + 1;
        }
        else if (jsoneq(input, &t[i], "rudder_limits") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            jsmntok_t * lb = &t[i + 2];
            ac->rudder_bounds[0] = atof(input + lb->start);
            jsmntok_t * ub = &t[i + 3];
            ac->rudder_bounds[1] = atof(input + ub->start); 
            i += t[i + 1].size + 1;
        }
        else if (jsoneq(input, &t[i], "thrust_limits") == 0) {
            if (t[i+1].type != JSMN_ARRAY) {
                continue;
            }
            jsmntok_t * lb = &t[i + 2];
            ac->thrust_bounds[0] = atof(input + lb->start);
            jsmntok_t * ub = &t[i + 3];
            ac->thrust_bounds[1] = atof(input + ub->start); 
            i += t[i + 1].size + 1;
        }
        else {
            printf("Unused key: %.*s\n", t[i].end - t[i].start,
                   input + t[i].start);
        }
    }
    
    fclose(fp); 
    free(input);  input = NULL;

    return aircraft_inertia(ac);
}

// Old hardcoded implementation
int pioneer_uav(struct Aircraft * ac)
{

    ac->m = 13.055641902393534;
    
    ac->span = 16.90;
    ac->chord = 1.80;
    ac->area = ac->span * ac->chord;
    ac->mac = 1.80,
    ac->AR = 9.3889;
    ac->e = 0.601427647984054;
    ac->K = 0.056370517893960;
    
    ac->Ixx = 34.832;
    ac->Ixz = -4.902;
    ac->Iyy = 67.08;
    ac->Izz = 82.22;

    ac->cphit = 1;
    ac->sphit = 0;

    //0 aoa aoa_dot mach q elevator
    ac->CL[0] = 0.385;
    ac->CL[1] = 4.78;
    ac->CL[2] = 2.42;
    ac->CL[3] = 0.0;
    ac->CL[4] = 8.05;
    ac->CL[5] = 0.401;

    //para aoa mach
    ac->CD[0] = 0.057579683439860;
    ac->CD[1] = 0.430;
    ac->CD[2] = 0.0;
    
    //beta p rudder
    ac->CE[0] = -0.819;
    ac->CE[1] = 0.0;
    ac->CE[2] = 0.191;
    
    // 0 aoa aoa_dot mach q elev
    ac->Cm[0] = 0.194;
    ac->Cm[1] = -2.12;
    ac->Cm[2] = -11.0;
    ac->Cm[3] = 0.0;
    ac->Cm[4] = -36.6;
    ac->Cm[5] = -1.76;
    
    // beta p r aileron rudder
    ac->Cl[0] = -0.023;
    ac->Cl[1] = -0.450;
    ac->Cl[2] = 0.265;
    ac->Cl[3] = -0.161;
    ac->Cl[4] = -0.00229;
    
    // beta p r aileron rudder
    ac->Cn[0] = 0.109;
    ac->Cn[1] = -0.110;
    ac->Cn[2] = -0.200;
    ac->Cn[3] = 0.020;
    ac->Cn[4] = -0.0917;

    return aircraft_inertia(ac);
}
