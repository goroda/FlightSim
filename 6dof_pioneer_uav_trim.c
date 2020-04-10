/*
  This is a 6DOF Rigid body aircraft simulator

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#include <cdyn/integrate.h>
#include <cdyn/simulate.h>

#include <nlopt.h>

typedef double real;

#define G 32.17

struct Vec3 {
    real v1;
    real v2;
    real v3;
};

inline real vec3_norm(const struct Vec3 * v) {return sqrt(v->v1 * v->v1 + v->v2 * v->v2 + v->v3 * v->v3);};

int vec3_cross(const struct Vec3 * a, const struct Vec3 * b, struct Vec3 * out)
{
    out->v1 = -a->v3 * b->v2 + a->v2 * b->v3;
    out->v2 = a->v3 * b->v1 - a->v1 * b->v3;
    out->v3 = -a->v2 * b->v1 + a->v1 * b->v2;
    return 0;
}

int vec3_add(const struct Vec3 * a, const struct Vec3 * b, struct Vec3 * out)
{
    out->v1 = a->v1 + b->v1;
    out->v2 = a->v2 + b->v2;
    out->v3 = a->v3 + b->v3;    
    return 0;
}

struct AeroAngles {

    real aoa;
    real sideslip;
    
    real aoa_g_w;
    real aoa_g_u;    
    real sideslip_g_v;
    real sideslip_g_vac;
};

int compute_aero_angles(const struct Vec3* UVW, real Vac, struct AeroAngles * ab)
{
    ab->aoa = atan2(UVW->v3, UVW->v1);
    ab->sideslip = asin(UVW->v2/Vac);
    return 0;
}

int compute_aero_angles_g (const struct Vec3* UVW, real Vac, struct AeroAngles * ab)
{
    ab->aoa = atan2(UVW->v3, UVW->v1);
    real den = pow(UVW->v3, 2) + pow(UVW->v1, 2);
    ab->aoa_g_w = UVW->v1 / den;
    ab->aoa_g_u = - UVW->v3 / den;


    
    ab->sideslip = asin(UVW->v2/Vac);
    real rat = Vac * sqrt(1 - pow(UVW->v2, 2) / pow(Vac, 2));
    real ratt = Vac * rat;

    ab->sideslip_g_v =  1.0 / rat;
    ab->sideslip_g_vac = - UVW->v2 / ratt;

    return 0;
}


struct EulerAngles{
    real yaw; // 1 
    real pitch; // 2 
    real roll; // 3

    int precomp;
    
    real cp;
    real sp;
    
    real cr;
    real sr;

    real cy;
    real sy;

    real sr_sp; // sin_roll * sin_pitch;
    real sr_cp; // sin_roll * cos_pitch;
    real sr_sy; // sin_roll * sin yaw
    real sr_cy; // sin_roll * cos yaw    
    real cr_sp; // cos_roll * sin_pitch
    real cr_cp; // cos_roll * cos_pitch
    real cr_cy; // cos_roll * cos_yaw
    real cr_sy; // cos_roll * sin_yaw
    real cp_cy; // cos_pitch * cos_yaw
    real cp_sy; // cos_pitch * sin_yaw

    real sr_sp_cy; // s roll s pitch c yaw
    real sr_sp_sy; // s roll s pitch s yaw

    real cr_sp_cy; // c roll s pitch c yaw
    real cr_sp_sy; // c roll s pitch s yaw


    real sr_sp_g_r, sr_sp_g_p;
    real sr_cp_g_r, sr_cp_g_p; 
    real sr_sy_g_r, sr_sy_g_y; 
    real sr_cy_g_r, sr_cy_g_y; 
    real cr_sp_g_r, cr_sp_g_p; 
    real cr_cp_g_r, cr_cp_g_p; 
    real cr_cy_g_r, cr_cy_g_y; 
    real cr_sy_g_r, cr_sy_g_y; 
    real cp_cy_g_p, cp_cy_g_y; 
    real cp_sy_g_p, cp_sy_g_y; 

    real sr_sp_cy_g_r, sr_sp_cy_g_p, sr_sp_cy_g_y; 
    real sr_sp_sy_g_r, sr_sp_sy_g_p, sr_sp_sy_g_y; 

    real cr_sp_cy_g_r, cr_sp_cy_g_p, cr_sp_cy_g_y; 
    real cr_sp_sy_g_r, cr_sp_sy_g_p, cr_sp_sy_g_y; 
};

int euler_angles_precompute(struct EulerAngles * ea)
{
    ea->precomp = 1;
    ea->cp = cos(ea->pitch);
    ea->sp = sin(ea->pitch);

    ea->cr = cos(ea->roll);
    ea->sr = sin(ea->roll);

    ea->cy = cos(ea->yaw);
    ea->sy = sin(ea->yaw);

    ea->sr_sp = ea->sr * ea->sp; // sin_roll * sin_pitch;
    ea->sr_cp = ea->sr * ea->cp; // sin_roll * cos_pitch;
    ea->sr_sy = ea->sr * ea->sy; // sin_roll * sin yaw
    ea->sr_cy = ea->sr * ea->cy; // sin_roll * cos yaw    
    ea->cr_sp = ea->cr * ea->sp; // cos_roll * sin_pitch
    ea->cr_cp = ea->cr * ea->cp; // cos_roll * cos_pitch
    ea->cr_cy = ea->cr * ea->cy; // cos_roll * cos_yaw
    ea->cr_sy = ea->cr * ea->sy; // cos_roll * sin_yaw
    ea->cp_cy = ea->cp * ea->cy; // cos_pitch * cos_yaw
    ea->cp_sy = ea->cp * ea->sy; // cos_pitch * sin_yaw

    ea->sr_sp_cy = ea->sr_sp * ea->cy; // s roll s pitch c yaw
    ea->sr_sp_sy = ea->sr_sp * ea->sy; // s roll s pitch s yaw

    ea->cr_sp_cy = ea->cr_sp * ea->cy; // c roll s pitch c yaw
    ea->cr_sp_sy = ea->cr_sp * ea->sy; // c roll s pitch s yaw    
    
    return 0;
}


int euler_angles_precompute_g(struct EulerAngles * ea)
{
    euler_angles_precompute(ea);

    ea->sr_sp_g_r = ea->cr_sp; 
    ea->sr_sp_g_p = ea->sr_cp;

    ea->sr_cp_g_r = ea->cr_cp;
    ea->sr_cp_g_p = -ea->sr_sp;

    ea->sr_sy_g_r = ea->cr_sy;
    ea->sr_sy_g_y = ea->sr_cy; 

    ea->sr_cy_g_r = ea->cr_cy;
    ea->sr_cy_g_y = -ea->sr_sy;

    ea->cr_sp_g_r = -ea->sr_sp;
    ea->cr_sp_g_p = ea->cr_cp;

    ea->cr_cp_g_r = -ea->sr_cp;
    ea->cr_cp_g_p = -ea->cr_sp;

    ea->cr_cy_g_r = -ea->sr_cy;
    ea->cr_cy_g_y = -ea->cr_sy; 

    ea->cr_sy_g_r = -ea->sr_sy;
    ea->cr_sy_g_y = ea->cr_cy;

    ea->cp_cy_g_p = -ea->sp * ea->cy; //-ea->sp_cy;
    ea->cp_cy_g_y = -ea->cp_sy;

    ea->cp_sy_g_p = -ea->sp * ea->sy; //-ea->sp_sy,
    ea->cp_sy_g_y = ea->cp_cy;     
    
    ea->sr_sp_cy_g_r = ea->cr_sp_cy;
    ea->sr_sp_cy_g_p = ea->sr_cp * ea->cy;
    ea->sr_sp_cy_g_y = -ea->sr_sp_sy; 

    ea->sr_sp_sy_g_r = ea->cr_sp_sy;
    ea->sr_sp_sy_g_p = ea->sr_cp * ea->sy;
    ea->sr_sp_sy_g_y = ea->sr_sp_cy; 

    ea->cr_sp_cy_g_r = -ea->sr_sp_cy;
    ea->cr_sp_cy_g_p = ea->cr_cp * ea->cy;
    ea->cr_sp_cy_g_y = -ea->cr_sp * ea->sy;

    ea->cr_sp_sy_g_r = -ea->sr_sp_sy;
    ea->cr_sp_sy_g_p = ea->cr_cp * ea->sy;
    ea->cr_sp_sy_g_y = ea->cr_sp_cy;
    
    return 0;
}


int orient_e_to_ac(const struct EulerAngles *ea, const struct Vec3* e, struct Vec3* ac)
{
    // yaw, pitch, roll
    
    ac->v1 = ea->cp_cy * e->v1 +
             ea->cp_sy * e->v2 -
             ea->sp * e->v3;

    ac->v2 = (ea->sr_sp_cy - ea->cr_sy) * e->v1 +
             (ea->sr_sp_sy + ea->cr_cy) * e->v2 +
                                ea->sr_cp * e->v3;

    ac->v3 = (ea->cr_sp_cy + ea->sr_sy) * e->v1 +
             (ea->cr_sp_sy - ea->sr_cy) * e->v2 +
                                ea->cr_cp * e->v3;

    return 0;
}

int orient_ac_to_e(const struct EulerAngles *ea, const struct Vec3* ac, struct Vec3* e)
{
    // -roll, -pitch, -yaw

    e->v1 =              ea->cp_cy * ac->v1 +
        (ea->sr_sp_cy - ea->cr_sy) * ac->v2 + 
        (ea->cr_sp_cy + ea->sr_sy) * ac->v3;

    e->v2 =              ea->cp_sy * ac->v1 +
        (ea->sr_sp_sy + ea->cr_cy) * ac->v2 +  
        (ea->cr_sp_sy - ea->sr_cy) * ac->v3;  

    e->v3 = - ea->sp * ac->v1 +
           ea->sr_cp * ac->v2 +
           ea->cr_cp * ac->v3;

    /* e->v3 = - ea->sp * ac->v1; */
    
    return 0;
}


struct StateGrad
{
    real U_g;
    real V_g;
    real W_g;
    real P_g;
    real Q_g;
    real R_g;
    real Roll_g;
    real Pitch_g;
    real Yaw_g;
};

int orient_ac_to_e_g(const struct EulerAngles *ea, const struct Vec3* ac, struct Vec3* e,
                       struct StateGrad sg[3])
{
    // -roll, -pitch, -yaw

    sg[0].U_g = ea->cp_cy;
    sg[0].V_g = ea->sr_sp_cy - ea->cr_sy;
    sg[0].W_g = ea->cr_sp_cy + ea->sr_sy;
    sg[0].P_g = 0; sg[0].Q_g = 0; sg[0].R_g = 0;

    e->v1 = sg[0].U_g * ac->v1 + sg[0].V_g * ac->v2 + sg[0].W_g * ac->v3;

    sg[0].Roll_g = (ea->sr_sp_cy_g_r - ea->cr_sy_g_r) * ac->v2 +
                    (ea->cr_sp_cy_g_r + ea->sr_sy_g_r) * ac->v3;
    
    sg[0].Pitch_g = ea->cp_cy_g_p * ac->v1 +
        ea->sr_sp_cy_g_p * ac->v2 +
        ea->cr_sp_cy_g_p * ac->v3;
        
    sg[0].Yaw_g = ea->cp_cy_g_y * ac->v1 +
         (ea->sr_sp_cy_g_y - ea->cr_sy_g_y) * ac->v2 +
         (ea->cr_sp_cy_g_y + ea->sr_sy_g_y) * ac->v3;
        
    sg[1].U_g = ea->cp_sy;
    sg[1].V_g = ea->sr_sp_sy + ea->cr_cy; 
    sg[1].W_g = ea->cr_sp_sy - ea->sr_cy; 
    sg[1].P_g = 0; sg[1].Q_g = 0; sg[1].R_g = 0;
    
    e->v2 = sg[1].U_g * ac->v1 +
            sg[1].V_g * ac->v2 +
            sg[1].W_g * ac->v3;

    sg[1].Roll_g = (ea->sr_sp_sy_g_r + ea->cr_cy_g_r) * ac->v2 +
                    (ea->cr_sp_sy_g_r - ea->sr_cy_g_r) * ac->v3;
    
    sg[1].Pitch_g = ea->cp_sy_g_p * ac->v1 +
                     ea->sr_sp_sy_g_p * ac->v2 +
                     ea->cr_sp_sy_g_p * ac->v3;

    sg[1].Yaw_g = ea->cp_sy_g_y * ac->v1 +
                     (ea->sr_sp_sy_g_y + ea->cr_cy_g_y) * ac->v2 +
                     (ea->cr_sp_sy_g_y - ea->sr_cy_g_y) * ac->v3;    
    

    sg[2].U_g = -ea->sp;
    sg[2].V_g = ea->sr_cp;
    sg[2].W_g = ea->cr_cp;
    sg[2].P_g = 0; sg[2].Q_g = 0; sg[2].R_g = 0;
    
    e->v3 = - ea->sp * ac->v1 +
              ea->sr_cp * ac->v2 +
              ea->cr_cp * ac->v3;

    sg[2].Roll_g = ea->sr_cp_g_r * ac->v2 + ea->cr_cp_g_r * ac->v3;
    sg[2].Pitch_g = -ea->cp * ac->v1 +
        ea->sr_cp_g_p * ac->v2 +
        ea->cr_cp_g_p * ac->v3; 
    sg[2].Yaw_g = 0.0;
    
    return 0;
}

int check_grad_ac_to_e(void)
{
    struct EulerAngles ea;
    ea.roll = M_PI/9.0;
    ea.pitch = M_PI/8.0;
    ea.yaw = M_PI/6.0;
    euler_angles_precompute_g(&ea);
    struct Vec3 ac = {120.0, 150.0, 80.0};
    struct Vec3 e;
    
    struct StateGrad sg[3];

    // reference
    orient_ac_to_e(&ea, &ac, &e);
    struct Vec3 e_ref = {e.v1, e.v2, e.v3};
    printf("Checking value computation\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", e.v1, e.v2, e.v3);    

    // compute analytic gradient
    orient_ac_to_e_g(&ea, &ac, &e, sg);
    printf("Analytic %3.5E %3.5E %3.5E\n", e.v1, e.v2, e.v3);
    printf("\n\n\n");
    
    double h = 1e-8;
    
    double grad_U[3];
    ac.v1 += h;
    orient_ac_to_e(&ea, &ac, &e);
    grad_U[0] = (e.v1 - e_ref.v1) / h;
    grad_U[1] = (e.v2 - e_ref.v2) / h;
    grad_U[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to U\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_U[0], grad_U[1], grad_U[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].U_g, sg[1].U_g, sg[2].U_g);
    printf("\n\n");
    
    double grad_V[3];
    ac.v1 -= h;
    ac.v2 += h;
    orient_ac_to_e(&ea, &ac, &e);
    grad_V[0] = (e.v1 - e_ref.v1) / h;
    grad_V[1] = (e.v2 - e_ref.v2) / h;
    grad_V[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to V\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_V[0], grad_V[1], grad_V[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].V_g, sg[1].V_g, sg[2].V_g);
    printf("\n\n");
    
    double grad_W[3];
    ac.v2 -= h;
    ac.v3 += h;
    orient_ac_to_e(&ea, &ac, &e);
    grad_W[0] = (e.v1 - e_ref.v1) / h;
    grad_W[1] = (e.v2 - e_ref.v2) / h;
    grad_W[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to W\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_W[0], grad_W[1], grad_W[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].W_g, sg[1].W_g, sg[2].W_g);
    printf("\n\n");

    double grad_R[3];
    ac.v3 -= h;
    ea.roll += h;
    euler_angles_precompute(&ea);
    orient_ac_to_e(&ea, &ac, &e);
    
    grad_R[0] = (e.v1 - e_ref.v1) / h;
    grad_R[1] = (e.v2 - e_ref.v2) / h;
    grad_R[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Roll\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_R[0], grad_R[1], grad_R[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Roll_g, sg[1].Roll_g, sg[2].Roll_g);
    printf("\n\n");
    
    double grad_P[3];
    ea.roll -= h;
    ea.pitch += h;
    euler_angles_precompute(&ea);
    orient_ac_to_e(&ea, &ac, &e);
    grad_P[0] = (e.v1 - e_ref.v1) / h;
    grad_P[1] = (e.v2 - e_ref.v2) / h;
    grad_P[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Pitch\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_P[0], grad_P[1], grad_P[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Pitch_g, sg[1].Pitch_g, sg[2].Pitch_g);
    printf("\n\n");
    
    double grad_Y[3];
    ea.pitch -= h;
    ea.yaw += h;
    euler_angles_precompute(&ea);
    orient_ac_to_e(&ea, &ac, &e);
    grad_Y[0] = (e.v1 - e_ref.v1) / h;
    grad_Y[1] = (e.v2 - e_ref.v2) / h;
    grad_Y[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Yaw\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Y[0], grad_Y[1], grad_Y[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Yaw_g, sg[1].Yaw_g, sg[2].Yaw_g);
    
    return 0;
}


// Translational kinematics
inline int tkin(const struct EulerAngles * ea, const struct Vec3 *UVW, struct Vec3 * rates){return orient_ac_to_e(ea, UVW, rates);}
inline int tkin_g(const struct EulerAngles * ea, const struct Vec3 *UVW, struct Vec3 * rates, struct StateGrad sg[3]){return orient_ac_to_e_g(ea, UVW, rates, sg);}

// Rotational kinematics
int rkin(const struct EulerAngles * ea, const struct Vec3 * PQR, struct Vec3 * rates)
{
    real secp = 1.0 / ea->cp;
    real tp = ea->sp  * secp;

    rates->v1 = PQR->v1 + ea->sr * tp * PQR->v2 + ea->cr * tp * PQR->v3;
    rates->v2 = ea->cr * PQR->v2 - ea->sr * PQR->v3;
    rates->v3 = ea->sr * secp * PQR->v2 + ea->cr * secp * PQR->v3;

    return 0;
}

// Rotational kinematics
int rkin_g(const struct EulerAngles * ea, const struct Vec3 * PQR, struct Vec3 * rates, struct StateGrad sg[3])
{
    real secp = 1.0 / ea->cp;
    real tp = ea->sp  * secp;

    real secp_g_p = tp * secp;
    real tp_p = secp * secp;

    sg[0].U_g = 0.0; sg[0].V_g = 0.0; sg[0].W_g = 0.0;
    sg[0].P_g = 1.0;
    sg[0].Q_g = ea->sr * tp;
    sg[0].R_g = ea->cr * tp;

    rates->v1 = sg[0].P_g * PQR->v1 + sg[0].Q_g * PQR->v2 + sg[0].R_g * PQR->v3;

    /* sg[0].Roll_g = ea->cr * tp * PQR->v2 - ea->sr * tp * PQR->v3; */
    sg[0].Roll_g = sg[0].R_g * PQR->v2 - sg[0].Q_g * PQR->v3; // turns out to be shortcut for above commented line
    sg[0].Pitch_g = ea->sr * tp_p * PQR->v2 + ea->cr * tp_p * PQR->v3;
    sg[0].Yaw_g = 0.0;

    sg[1].U_g = 0.0; sg[1].V_g = 0.0; sg[1].W_g = 0.0;
    sg[1].P_g = 0.0;
    sg[1].Q_g = ea->cr;
    sg[1].R_g = -ea->sr;

    rates->v2 = ea->cr * PQR->v2 - ea->sr * PQR->v3;

    sg[1].Roll_g = -ea->sr * PQR->v2 - ea->cr * PQR->v3;
    sg[1].Pitch_g = 0.0;
    sg[1].Yaw_g = 0.0;

    
    sg[2].U_g = 0.0; sg[2].V_g = 0.0; sg[2].W_g = 0.0;
    sg[2].P_g = 0.0;
    sg[2].Q_g = ea->sr * secp;
    sg[2].R_g = ea->cr * secp;
    
    rates->v3 = sg[2].Q_g * PQR->v2 + sg[2].R_g * PQR->v3;

    /* sg[2].Roll_g = ea->cr * secp * PQR->v2 - ea->sr * secp * PQR->v3; */
    sg[2].Roll_g = sg[2].R_g * PQR->v2 - sg[2].Q_g * PQR->v3;   //shortand for above 
    sg[2].Pitch_g = ea->sr * secp_g_p * PQR->v2 + ea->cr * secp_g_p * PQR->v3;
    sg[2].Yaw_g = 0.0;
    
    return 0;
}

int check_grad_rkin(void)
{
    struct EulerAngles ea;
    ea.roll = M_PI/9.0;
    ea.pitch = M_PI/8.0;
    ea.yaw = M_PI/6.0;
    euler_angles_precompute_g(&ea);
    struct Vec3 pqr = {120.0, 150.0, 80.0};
    struct Vec3 rates;
    
    struct StateGrad sg[3];

    // reference
    rkin(&ea, &pqr, &rates);
    struct Vec3 e_ref = {rates.v1, rates.v2, rates.v3};
    printf("Checking value computation\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", rates.v1, rates.v2, rates.v3);    

    // compute analytic gradient
    rkin_g(&ea, &pqr, &rates, sg);
    printf("Analytic %3.5E %3.5E %3.5E\n", rates.v1, rates.v2, rates.v3);
    printf("\n\n\n");
    
    double h = 1e-10;
    
    double grad_P[3];
    pqr.v1 += h;
    rkin(&ea, &pqr, &rates);
    grad_P[0] = (rates.v1 - e_ref.v1) / h;
    grad_P[1] = (rates.v2 - e_ref.v2) / h;
    grad_P[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to P\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_P[0], grad_P[1], grad_P[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].P_g, sg[1].P_g, sg[2].P_g);
    printf("\n\n");
    
    double grad_Q[3];
    pqr.v1 -= h;
    pqr.v2 += h;
    rkin(&ea, &pqr, &rates);
    grad_Q[0] = (rates.v1 - e_ref.v1) / h;
    grad_Q[1] = (rates.v2 - e_ref.v2) / h;
    grad_Q[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Q\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Q[0], grad_Q[1], grad_Q[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Q_g, sg[1].Q_g, sg[2].Q_g);
    printf("\n\n");
    
    double grad_R[3];
    pqr.v2 -= h;
    pqr.v3 += h;
    rkin(&ea, &pqr, &rates);
    grad_R[0] = (rates.v1 - e_ref.v1) / h;
    grad_R[1] = (rates.v2 - e_ref.v2) / h;
    grad_R[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to R\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_R[0], grad_R[1], grad_R[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].R_g, sg[1].R_g, sg[2].R_g);
    printf("\n\n");

    double grad_Roll[3];
    pqr.v3 -= h;
    ea.roll += h;
    euler_angles_precompute(&ea);
    rkin(&ea, &pqr, &rates);
    
    grad_Roll[0] = (rates.v1 - e_ref.v1) / h;
    grad_Roll[1] = (rates.v2 - e_ref.v2) / h;
    grad_Roll[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Roll\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Roll[0], grad_Roll[1], grad_Roll[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Roll_g, sg[1].Roll_g, sg[2].Roll_g);
    printf("\n\n");
    
    double grad_Pitch[3];
    ea.roll -= h;
    ea.pitch += h;
    euler_angles_precompute(&ea);
    rkin(&ea, &pqr, &rates);
    grad_Pitch[0] = (rates.v1 - e_ref.v1) / h;
    grad_Pitch[1] = (rates.v2 - e_ref.v2) / h;
    grad_Pitch[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Pitch\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Pitch[0], grad_Pitch[1], grad_Pitch[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Pitch_g, sg[1].Pitch_g, sg[2].Pitch_g);
    printf("\n\n");
    
    double grad_Yaw[3];
    ea.pitch -= h;
    ea.yaw += h;
    euler_angles_precompute(&ea);
    rkin(&ea, &pqr, &rates);
    grad_Yaw[0] = (rates.v1 - e_ref.v1) / h;
    grad_Yaw[1] = (rates.v2 - e_ref.v2) / h;
    grad_Yaw[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Yaw\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Yaw[0], grad_Yaw[1], grad_Yaw[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Yaw_g, sg[1].Yaw_g, sg[2].Yaw_g);
    
    return 0;
}

struct Aircraft
{
    real m; // mass
    real cphit; // cos (phiT) where phiT is the angle of the Thrust (0 is aligne with body)
    real sphit;

    real Ixx;
    real Iyy;
    real Izz;
    real Ixz;

    real rtn1, rtn2, rtn3;
    real rt11, rt12;
    real rt21, rt22;
    real rt31, rt32;

    real span;
    real chord;
    real area;
    real mac;
    real AR;
    real e;
    real K;

    // Force coefficients
    real CL[6]; 
    real CD[3];
    real CE[3];

    real Cm[6];
    real Cl[5];
    real Cn[5];
    
};

int aircraft_precompute_inertia(struct Aircraft * ac)
{

    ac->rtn1 = ac->Ixx * ac->Izz - (ac->Ixz * ac->Ixz);
    ac->rtn2 = ac->Iyy;
    ac->rtn3 = ac->rtn1;

    ac->rt11 = (ac->Iyy * ac->Izz - ac->Izz*ac->Izz - ac->Ixz * ac->Ixz);
    ac->rt12 = ac->Ixz * (ac->Ixx - ac->Iyy + ac->Izz);

    ac->rt21 = ac->Izz - ac->Ixx;
    ac->rt22 = ac->Ixz;

    ac->rt31 = -ac->Ixx * ac->Iyy + ac->Ixx * ac->Ixx + ac->Ixz * ac->Ixz;
    ac->rt32 = ac->Ixz * (-ac->Ixx + ac->Iyy - ac->Izz);

    return 0;
}

int pioneer_uav(struct Aircraft * ac){

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

    return aircraft_precompute_inertia(ac);
}


inline real aircraft_get_mass(const struct Aircraft* ac){ return ac->m; };

inline real aircraft_get_cphit(const struct Aircraft* ac){ return ac->cphit; };
inline real aircraft_get_sphit(const struct Aircraft* ac){ return ac->sphit; };

inline real aircraft_get_Ixx(const struct Aircraft* ac){ return ac->Ixx; };
inline real aircraft_get_Iyy(const struct Aircraft* ac){ return ac->Iyy; };
inline real aircraft_get_Izz(const struct Aircraft* ac){ return ac->Izz; };
inline real aircraft_get_Ixz(const struct Aircraft* ac){ return ac->Ixz; };

inline real aircraft_get_rtn1(const struct Aircraft* ac){ return ac->rtn1; };
inline real aircraft_get_rtn2(const struct Aircraft* ac){ return ac->rtn2; };
inline real aircraft_get_rtn3(const struct Aircraft* ac){ return ac->rtn3; };

inline real aircraft_get_rt11(const struct Aircraft* ac){ return ac->rt11; };
inline real aircraft_get_rt12(const struct Aircraft* ac){ return ac->rt12; };
inline real aircraft_get_rt21(const struct Aircraft* ac){ return ac->rt21; };
inline real aircraft_get_rt22(const struct Aircraft* ac){ return ac->rt22; };
inline real aircraft_get_rt31(const struct Aircraft* ac){ return ac->rt31; };
inline real aircraft_get_rt32(const struct Aircraft* ac){ return ac->rt32; };


// Translational Dynamics
int tdyn(const struct EulerAngles * ea, const struct AeroAngles * aero,
         const struct Vec3 * UVW, const struct Vec3 * PQR, const struct Vec3 * DEL, real Ft,
         const struct Aircraft * ac, struct Vec3 * rates)
{

    real U = UVW->v1;
    real V = UVW->v2;
    real W = UVW->v3;
    
    real P = PQR->v1;
    real Q = PQR->v2;
    real R = PQR->v3;

    real D = DEL->v1;
    /* real E = DEL->v2; */
    real L = DEL->v3;

    double ca = cos(aero->aoa);
    double sa = sin(aero->aoa);
    
    double cb = cos(aero->sideslip);
    double sb = sin(aero->sideslip);

    real m = aircraft_get_mass(ac);
    real cphit = aircraft_get_cphit(ac);
    real sphit = aircraft_get_sphit(ac);
    
    rates->v1 = V * R - W * Q - ea->sp * G - cb * ca * D / m + sa * L / m + cphit * Ft / m;
    rates->v2 = -U * R + W * P + ea->sr * ea->cp * G - sb * D / m;
    rates->v3 = U * Q - V*P + ea->cr * ea->cp * G - cb * sa * D / m - ca * L / m - sphit * Ft / m;

    return 0;
}

// Rotational Dynamics
int rdyn(const struct Vec3 * PQR, const struct Vec3 * LMN, const struct Aircraft * ac,
         struct Vec3 * rates)
{
    real P = PQR->v1;
    real Q = PQR->v2;
    real R = PQR->v3;

    real QR = Q * R;
    real PQ = P * Q;

    real L = LMN->v1;
    real M = LMN->v2;
    real N = LMN->v3;


    real rt11 = aircraft_get_rt11(ac);
    real rt12 = aircraft_get_rt12(ac);
    real rt21 = aircraft_get_rt21(ac);
    real rt22 = aircraft_get_rt22(ac);
    real rt31 = aircraft_get_rt31(ac);
    real rt32 = aircraft_get_rt32(ac);
    real rtn1 = aircraft_get_rtn1(ac);
    real rtn2 = aircraft_get_rtn2(ac);
    real rtn3 = aircraft_get_rtn3(ac);

    real Izz = aircraft_get_Izz(ac);
    real Ixz = aircraft_get_Ixz(ac);
    real Ixx = aircraft_get_Ixx(ac);
    
    // dot{p}
    rates->v1 = (rt11 * QR + rt12 * PQ + Izz * L + Ixz * N) / rtn1;
    rates->v2 = (rt21 * P * R + rt22 * (R * R - P * P) + M) / rtn2;
    rates->v3 = (rt31 * PQ + rt32 * QR + Ixz * L + Ixx * N) / rtn3;
    
    return 0;
}

int compute_aero_forces(const struct AeroAngles * aero,
                        const struct Vec3 * PQR, const struct Vec3 * aero_con,
                        const struct Aircraft * ac,
                        real rho,
                        real vac,
                        struct Vec3 * DEL,
                        struct Vec3 * LMN)
{
    real CL = ac->CL[0] +
        ac->CL[1] * aero->aoa +
        (ac->CL[4] * ac->chord / 2.0 / vac) * PQR->v2 +
        ac->CL[5] * aero_con->v1; // elevator

    real CD = ac->CD[0] + ac->K * CL * CL;
    
    real Cm = ac->Cm[0] +
        ac->Cm[1] * aero->aoa +
        (ac->Cm[4] / 2.0 / vac * ac->chord) * PQR->v2 +
        ac->Cm[5] * aero_con->v1; // elevator

    /* printf("vac = %3.2E\n", vac); */
    real Cl = ac->Cl[0] * aero->sideslip +
        ac->Cl[1] * ac->span / 2.0 / vac * PQR->v1 +
        ac->Cl[2] * ac->span / 2.0 / vac * PQR->v3 +
        ac->Cl[3] * aero_con->v2 + // aileron
        ac->Cl[4] * aero_con->v3; // rudder

    /* printf("ac->cl[1] * p = %3.2E\n", ac->Cl[1] * ac->span / 2.0 / vac * PQR->v1); */
    /* printf("ac->cl[2] * r = %3.2E\n", ac->Cl[2] * ac->span / 2.0 / vac * PQR->v3); */
    /* printf("ac->cl[2] = %3.2E\n", ac->Cl[2] * ac->span / 2.0 / vac);       */
    
    real Cn = ac->Cn[0] * aero->sideslip +
        ac->Cn[1] * ac->span / 2.0 / vac * PQR->v1 +
        ac->Cn[2] * ac->span / 2.0 / vac * PQR->v3 +
        ac->Cn[3] * aero_con->v2 + // aileron
        ac->Cn[4] * aero_con->v3; // rudder    

    
    /* printf("Cl = %3.5E, Cm = %3.5E, Cn = %3.5E\n", Cl, Cm, Cn); */
    
    real vac2 = vac * vac;
    
    real pre = 0.5 * rho * ac->area * vac2;
    real pre_mom = pre * ac->span;

    DEL->v1 = pre * CD; // drag
    DEL->v2 = 0.0;      // side drag
    DEL->v3 = pre * CL; // lift

    LMN->v1 = pre_mom * Cl;         // roll mom
    LMN->v2 = pre * ac->chord * Cm; // pitch mom
    LMN->v3 = pre_mom * Cn;         // yaw mom

    return 0;
}

int rigid_body_lin_forces(double time, const double * state,
                          const double * control,
                          double * out, double * jac,
                          void * arg)
{
    (void) time;
    assert (jac == NULL);

    struct Aircraft * ac = arg;
    
    struct Vec3 UVW = {state[3], state[4], state[5]};
    real vac = vec3_norm(&UVW);

    struct AeroAngles aero;
    compute_aero_angles(&UVW, vac, &aero);

    #ifdef DEBUG
    printf("aoa = %3.2E\n", aero.aoa);
    printf("beta = %3.2E\n", aero.sideslip);
    #endif
    
    struct Vec3 PQR = {state[6], state[7], state[8]};
    struct EulerAngles ea;
    ea.roll = state[9];
    ea.pitch = state[10];
    ea.yaw = state[11];
    euler_angles_precompute(&ea);


    /////////////////////////////
    // Kinematics
    ////////////////////////////
    struct Vec3 xyz_rates;
    tkin(&ea, &UVW, &xyz_rates);

    struct Vec3 euler_rates;
    rkin(&ea, &PQR, &euler_rates);

    /////////////////////////////
    // Dynamics
    ////////////////////////////

    // Force and Moment Computation
    struct Vec3 DEL = {0, 0, 0};
    struct Vec3 LMN = {0, 0, 0};    
    real rho = 0.002376892406675; // slug / ft^3

    // elevator, aileron, rudder
    struct Vec3 aero_con = {control[0], control[1], control[2]};
    // Thrust
    real Ft = control[3]; 
    
    compute_aero_forces(&aero, &PQR, &aero_con, ac, rho, vac, &DEL, &LMN);

    // Momentum equations
    struct Vec3 uvw_rates;
    tdyn(&ea, &aero, &UVW, &PQR, &DEL, Ft, ac, &uvw_rates);

    struct Vec3 pqr_rates;
    rdyn(&PQR, &LMN, ac, &pqr_rates);
    
    out[0] = xyz_rates.v1;
    out[1] = xyz_rates.v2;
    out[2] = xyz_rates.v3;

    out[3] = uvw_rates.v1;
    out[4] = uvw_rates.v2;
    out[5] = uvw_rates.v3;

    out[6] = pqr_rates.v1;
    out[7] = pqr_rates.v2;
    out[8] = pqr_rates.v3;

    out[9]  = euler_rates.v1;
    out[10] = euler_rates.v2;
    out[11] = euler_rates.v3;    

    return 0;
}


struct TrimSpec
{
    real z_dot;
    real yaw_dot;
    real target_vel;
    
    struct Aircraft * ac;

    real thresh;
};

inline real trim_spec_get_climb_rate(const struct TrimSpec * spec){ return spec->z_dot; }
inline real trim_spec_get_yaw_rate(const struct TrimSpec * spec){ return spec->yaw_dot; }
inline real trim_spec_get_speed(const struct TrimSpec * spec){ return spec->target_vel; }

double trim_objective(unsigned n, const double * x, double * grad, void * f_data)
{
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
    return out;
}


struct SteadyState
{
    struct Vec3 UVW, dUVW;
    struct AeroAngles aero;
    struct Vec3 PQR, dPQR;
    struct Vec3 aero_con;    
    real roll, droll, pitch, dpitch;
    real thrust;

    nlopt_result res;
    double obj_val;

    real target_yaw_rate, achieved_yaw_rate;
    real target_speed, achieved_speed;

    real target_climb_rate, achieved_climb_rate;
    real aoa;
    real sideslip;
};

inline real rad2deg(real ang) { return ang / M_PI * 180; }

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
    fprintf(fp, "-Climb Rate (ft/s)  :       %3.5E      %3.5E\n", ss->target_climb_rate, ss->achieved_climb_rate);
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
    fprintf(fp, "Angle of Attack (rad,deg) :   %3.5E \t %3.5E\n", ss->aero.aoa, rad2deg(ss->aero.aoa));
    fprintf(fp, "Sideslip Angle  (rad,deg) :   %3.5E \t %3.5E\n", ss->aero.sideslip, rad2deg(ss->aero.sideslip));

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
    lb[1] = 0.0;  // lower bound V
    ub[1] = 0.0;  // upper bound V

    x[1] = 0.0; // initial condition V
    x[11] = 10; // initial condition thrust
    
    // run without bounds
    /* opt = nlopt_create(NLOPT_LN_NELDERMEAD, 12); */
    opt = nlopt_create(NLOPT_LN_NEWUOA, 12);
    nlopt_set_ftol_rel(opt, -1.0);
    nlopt_set_ftol_abs(opt, 1e-20);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);    
    nlopt_set_min_objective(opt, trim_objective, data);
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
    nlopt_set_ftol_rel(opt, -1.0);
    nlopt_set_ftol_abs(opt, 0.0);
    nlopt_set_xtol_abs1(opt, 1e-20);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, trim_objective, data);

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
    ss->achieved_climb_rate = sol[2];
    
    ss->target_yaw_rate = trim_spec_get_yaw_rate(data);
    ss->achieved_yaw_rate = sol[11];
    
    ss->target_speed = trim_spec_get_speed(data);
    ss->achieved_speed = sqrt(pow(x[0],2) + pow(x[1], 2) + pow(x[2], 2));

    compute_aero_angles(&(ss->UVW), ss->achieved_speed, &(ss->aero));
    
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


static char * program_name;

void print_code_usage (FILE *, int) __attribute__ ((noreturn));
void print_code_usage (FILE * stream, int exit_code)
{

    fprintf(stream, "Usage: %s options \n", program_name);
    fprintf(stream,
            " -h --help                Display this usage information.\n"
            " -s --speed      <val>    Desired speed (e.g., 120 --> flight at 120 ft/s, groundspeed), default 120.\n"
            " -c --climb-rate <val>    Desired climb rate (e.g., -5 --> climb at 5 ft/s), default 0.\n"
            " -y --yaw-rate   <val>    Desired turn rate (e.g., 3.14 --> turn 'right' at pi rad/s), default 0.\n"
            " -t --threshold  <val>    Threshold value for setting states to zero default 1e-10.\n"            
            /* " -v --verbose    <val>      Output words (default 0)\n" */
        );
    exit (exit_code);
}

int main(int argc, char* argv[]){


    /* check_grad_ac_to_e(); */
    check_grad_rkin();    
    return 0;
    
    int next_option;
    const char * const short_options = "hs:c:y:t:";
    const struct option long_options[] = {
        { "help"       ,  0, NULL, 'h' },
        { "speed"      , 1, NULL, 's' },
        { "climb-rate" , 1, NULL, 'c' },
        { "yaw-rate"   , 1, NULL, 'y' },
        { "threshold"  , 1, NULL, 't' },
        /* { "verbose"    , 1, NULL, 'v' }, */
        { NULL         , 0, NULL, 0   }
    };
    

    real speed = 120.0;
    real climb_rate = 0.0;
    real yaw_rate = 0.0;
    real thresh = 1e-10; // threshold for zero
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
            /* case 'v': */
            /*     verbose = strtol(optarg,NULL,10); */
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

    

    
    struct Aircraft aircraft;
    pioneer_uav(&aircraft);
    struct TrimSpec trim_spec;
    trim_spec.z_dot = climb_rate;
    trim_spec.yaw_dot = yaw_rate; ///3.0 * 2.0 * M_PI / 500.0;
    trim_spec.target_vel = speed; // ft/s
    trim_spec.ac = &aircraft;
    trim_spec.thresh = thresh; 

    struct SteadyState ss;
    trimmer(&trim_spec, &ss);
    steady_state_print(stdout, &ss);


    int simulate = 0;    
    if (simulate == 1){
        struct Vec3 xyz = {0, 0, -5};
        real yaw = M_PI / 4.0;
        double dt_save = 1e-1;
        size_t nsteps = 1000;
        struct Trajectory * traj = flight_sim_ss(&xyz, yaw, &ss, &aircraft, dt_save, nsteps);

        trajectory_print(traj,stdout,5);
        trajectory_free(traj); traj = NULL;
    }

    return 0;
    
}
