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
#include <float.h>

#include "flight_sim.h"

typedef double real;
inline real rad2deg(real ang) { return ang / M_PI * 180; }

//-----------------------------------
//-----------------------------------
//--------- Prototypes --------------
//-----------------------------------
//-----------------------------------

int check_grad_ac_to_e(void);
int check_grad_rkin(void);
int check_grad_aero_forces(void);
int check_grad_rdyn(void);
int check_grad_tdyn(void);
int check_grad_rigid_body_dyn(void);
int check_grad_trimmer(void);

#define G 32.17


inline real steady_state_get_U(const struct SteadyState * ss){return ss->UVW.v1;}
inline real steady_state_get_V(const struct SteadyState * ss){return ss->UVW.v2;}
inline real steady_state_get_W(const struct SteadyState * ss){return ss->UVW.v3;}
inline real steady_state_get_P(const struct SteadyState * ss){return ss->PQR.v1;}
inline real steady_state_get_Q(const struct SteadyState * ss){return ss->PQR.v2;}
inline real steady_state_get_R(const struct SteadyState * ss){return ss->PQR.v3;}
inline real steady_state_get_Roll(const struct SteadyState * ss){return ss->roll;}
inline real steady_state_get_Pitch(const struct SteadyState * ss){return ss->pitch;}
inline real steady_state_get_elevator(const struct SteadyState * ss){return ss->aero_con.v1;}
inline real steady_state_get_aileron(const struct SteadyState * ss){return ss->aero_con.v2;}
inline real steady_state_get_rudder(const struct SteadyState * ss){return ss->aero_con.v3;}
inline real steady_state_get_thrust(const struct SteadyState * ss){return ss->thrust;}

int steady_state_set_vec(const struct SteadyState * ss, real x, real y, real z, real yaw,
                         real state[12],
                         real control[4])
{
    state[0] = x;
    state[1] = y;
    state[2] = z;
    state[3] = steady_state_get_U(ss);
    state[4] = steady_state_get_V(ss);
    state[5] = steady_state_get_W(ss);
    state[6] = steady_state_get_P(ss);
    state[7] = steady_state_get_Q(ss);
    state[8] = steady_state_get_R(ss);
    state[9] = steady_state_get_Roll(ss);
    state[10] = steady_state_get_Pitch(ss);
    state[11] = yaw;

    if (control != NULL){
        control[0] = steady_state_get_elevator(ss);
        control[1] = steady_state_get_aileron(ss);
        control[2] = steady_state_get_rudder(ss);
        control[3] = steady_state_get_thrust(ss);
    }

    return 0;
}

/* int  */


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

int aircraft_inertia(struct Aircraft * ac)
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

int aircraft_init(struct Aircraft * ac)
{

    ac->m = 0.0;
    
    ac->span = 0.0;
    ac->chord = 0.0;
    ac->area = ac->span * ac->chord;
    ac->mac = 0.0;
    ac->AR = 0.0;
    ac->e = 0.0;
    ac->K = 0.0;
    
    ac->Ixx = 0.0;
    ac->Ixz = 0.0;
    ac->Iyy = 0.0;
    ac->Izz = 0.0;

    ac->cphit = 1;
    ac->sphit = 0;

    //0 aoa aoa_dot mach q elevator
    ac->CL[0] = 0.0;
    ac->CL[1] = 0.0;
    ac->CL[2] = 0.0;
    ac->CL[3] = 0.0;
    ac->CL[4] = 0.0;
    ac->CL[5] = 0.0;

    //para aoa mach
    ac->CD[0] = 0.0;
    ac->CD[1] = 0.0;
    ac->CD[2] = 0.0;
    
    //beta p rudder
    ac->CE[0] = 0.0;
    ac->CE[1] = 0.0;
    ac->CE[2] = 0.0;
    
    // 0 aoa aoa_dot mach q elev
    ac->Cm[0] = 0.0;
    ac->Cm[1] = 0.0;
    ac->Cm[2] = 0.0;
    ac->Cm[3] = 0.0;
    ac->Cm[4] = 0.0;
    ac->Cm[5] = 0.0;
    
    // beta p r aileron rudder
    ac->Cl[0] = 0.0;
    ac->Cl[1] = 0.0;
    ac->Cl[2] = 0.0;
    ac->Cl[3] = 0.0;
    ac->Cl[4] = 0.0;
    
    // beta p r aileron rudder
    ac->Cn[0] = 0.0;
    ac->Cn[1] = 0.0;
    ac->Cn[2] = 0.0;
    ac->Cn[3] = 0.0;
    ac->Cn[4] = 0.0;

    // degrees
    ac->aoa_bounds[0] = -20.0;
    ac->aoa_bounds[1] = 20.0;

    ac->sideslip_bounds[0] = -20.0;
    ac->sideslip_bounds[1] = 20.0;
    
    ac->elevator_bounds[0] = -20;
    ac->elevator_bounds[1] = 20.0;

    ac->aileron_bounds[0] = -20.0;
    ac->aileron_bounds[1] = 20.0;
    
    ac->rudder_bounds[0] = -20.0;
    ac->rudder_bounds[1] = 20.0;
    
    ac->thrust_bounds[0] = 0.0;
    ac->thrust_bounds[1] = DBL_MAX;

    return aircraft_inertia(ac);
}



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


int aero_angles(const struct Vec3* UVW, real Vac, struct AeroAngles * ab)
{
    ab->aoa = atan2(UVW->v3, UVW->v1);
    ab->sideslip = asin(UVW->v2/Vac);
    return 0;
}

int aero_angles_g (const struct Vec3* UVW, real Vac,
                           const struct StateGrad * vac_g,
                           struct AeroAngles * ab)
{
    ab->aoa = atan2(UVW->v3, UVW->v1);
    real den = pow(UVW->v3, 2) + pow(UVW->v1, 2);
    ab->aoa_g_w = UVW->v1 / den;
    ab->aoa_g_u = - UVW->v3 / den;


    
    ab->sideslip = asin(UVW->v2/Vac);
    real rat = Vac * sqrt(1 - pow(UVW->v2, 2) / pow(Vac, 2));
    real ratt = Vac * rat;

    ab->sideslip_g_v =  (Vac - UVW->v2 * vac_g->V_g) / ratt;
    ab->sideslip_g_vac = - UVW->v2 / ratt;

    return 0;
}


int euler_angles(struct EulerAngles * ea)
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


int euler_angles_g(struct EulerAngles * ea)
{
    euler_angles(ea);

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

// Translational Dynamics
int tdyn_g(const struct EulerAngles * ea, const struct AeroAngles * aero, const struct StateGrad * vac_g,
           const struct Vec3 * UVW, const struct Vec3 * PQR, const struct Vec3 * DEL,
           const struct StateGrad DEL_g[3],
           const struct ControlGrad DEL_cg[3],
           real Ft,
           const struct Aircraft * ac,
           struct Vec3 * rates,
           struct StateGrad sg[3],
           struct ControlGrad cg[3])
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

    real ca_g_u = -sa * aero->aoa_g_u;
    real ca_g_v = 0.0;
    real ca_g_w = -sa * aero->aoa_g_w;

    
    real sa_g_u = ca * aero->aoa_g_u;
    real sa_g_v = 0.0;
    real sa_g_w = ca * aero->aoa_g_w;    

    double cb = cos(aero->sideslip);
    double sb = sin(aero->sideslip);

    real cb_g_u = -sb * aero->sideslip_g_vac * vac_g->U_g;
    real cb_g_v = -sb * aero->sideslip_g_v;
    real cb_g_w = -sb * aero->sideslip_g_vac * vac_g->W_g;    

    real sb_g_u = cb * aero->sideslip_g_vac * vac_g->U_g;
    real sb_g_v = cb * aero->sideslip_g_v;
    real sb_g_w = cb * aero->sideslip_g_vac * vac_g->W_g;    
    
    real cb_ca = cb * ca;
    real cb_ca_g_u = cb * ca_g_u + cb_g_u * ca;
    real cb_ca_g_v = cb * ca_g_v + cb_g_v * ca;
    real cb_ca_g_w = cb * ca_g_w + cb_g_w * ca;

    real cb_sa = cb * sa;
    real cb_sa_g_u = cb_g_u * sa + cb * sa_g_u;
    real cb_sa_g_v = cb_g_v * sa + cb * sa_g_v;
    real cb_sa_g_w = cb_g_w * sa + cb * sa_g_w;
    
    real m = aircraft_get_mass(ac);
    real cphit = aircraft_get_cphit(ac);
    real sphit = aircraft_get_sphit(ac);
    
    rates->v1 = V * R - W * Q - ea->sp * G - cb_ca * D / m + sa * L / m + cphit * Ft / m;

    sg[0].U_g = (1.0/m) * (-cb_ca_g_u * D - cb_ca * DEL_g[0].U_g + sa_g_u * L + sa * DEL_g[2].U_g);
    sg[0].V_g = R + (1.0/m) * (-cb_ca_g_v * D - cb_ca * DEL_g[0].V_g + sa_g_v * L + sa * DEL_g[2].V_g);
    sg[0].W_g = -Q + (1.0/m) * (-cb_ca_g_w * D - cb_ca * DEL_g[0].W_g + sa_g_w * L + sa * DEL_g[2].W_g);
    sg[0].P_g = (1.0/m) * ( - cb_ca * DEL_g[0].P_g + sa * DEL_g[2].P_g);
    sg[0].Q_g = -W + (1.0/m) * ( - cb_ca * DEL_g[0].Q_g + sa * DEL_g[2].Q_g);
    sg[0].R_g = V + (1.0/m) * ( - cb_ca * DEL_g[0].R_g + sa * DEL_g[2].R_g);
    sg[0].Roll_g = (1.0/m) * (- cb_ca * DEL_g[0].Roll_g + sa * DEL_g[2].Roll_g);
    sg[0].Pitch_g = -ea->cp * G + (1.0/m) * (- cb_ca * DEL_g[0].Pitch_g + sa * DEL_g[2].Pitch_g);
    sg[0].Yaw_g = (1.0/m) * (- cb_ca * DEL_g[0].Yaw_g + sa * DEL_g[2].Yaw_g);

    cg[0].elev_g = (1.0/m) * ( - cb_ca * DEL_cg[0].elev_g + sa * DEL_cg[2].elev_g);
    cg[0].thrust_g = cphit / m;
    cg[0].aileron_g = (1.0/m) * ( - cb_ca * DEL_cg[0].aileron_g + sa * DEL_cg[2].aileron_g);
    cg[0].rudder_g = (1.0/m) * ( - cb_ca * DEL_cg[0].rudder_g + sa * DEL_cg[2].rudder_g);

        
    rates->v2 = -U * R + W * P + ea->sr_cp * G - sb * D / m;

    sg[1].U_g = -R + (1.0/m) * (-sb_g_u * D - sb * DEL_g[0].U_g);
    sg[1].V_g = (1.0/m) * (-sb_g_v * D - sb * DEL_g[0].V_g);
    sg[1].W_g = P + (1.0/m) * (-sb_g_w * D - sb * DEL_g[0].W_g);
    sg[1].P_g = W + (1.0/m) * (- sb * DEL_g[0].P_g);
    sg[1].Q_g = (1.0/m) * (- sb * DEL_g[0].Q_g);
    sg[1].R_g = -U + (1.0/m) * (- sb * DEL_g[0].R_g);
    sg[1].Roll_g = ea->sr_cp_g_r * G + (1.0/m) * (- sb * DEL_g[0].Roll_g);
    sg[1].Pitch_g = ea->sr_cp_g_p * G + (1.0/m) * (- sb * DEL_g[0].Pitch_g);
    sg[1].Yaw_g = (1.0/m) * (- sb * DEL_g[0].Yaw_g);    

    cg[1].elev_g = (1.0/m) * (- sb * DEL_cg[0].elev_g);
    cg[1].thrust_g = 0.0;
    cg[1].aileron_g = (1.0/m) * (- sb * DEL_cg[0].aileron_g);
    cg[1].rudder_g = (1.0/m) * (- sb * DEL_cg[0].rudder_g);

    rates->v3 = U * Q - V*P + ea->cr_cp * G - cb_sa * D / m - ca * L / m - sphit * Ft / m;

    sg[2].U_g = Q + (1/m) * (-cb_sa * DEL_g[0].U_g - cb_sa_g_u * D - ca * DEL_g[2].U_g - ca_g_u * L);
    sg[2].V_g = -P + (1/m) * (-cb_sa * DEL_g[0].V_g - cb_sa_g_v * D - ca * DEL_g[2].V_g - ca_g_v * L);
    sg[2].W_g = (1/m) * (-cb_sa * DEL_g[0].W_g - cb_sa_g_w * D - ca * DEL_g[2].W_g - ca_g_w * L);
    sg[2].P_g = -V + (1/m) * (-cb_sa * DEL_g[0].P_g  - ca * DEL_g[2].P_g);
    sg[2].Q_g = U + (1/m) * (-cb_sa * DEL_g[0].Q_g  - ca * DEL_g[2].Q_g);
    sg[2].R_g = (1/m) * (-cb_sa * DEL_g[0].R_g  - ca * DEL_g[2].R_g);
    sg[2].Roll_g = ea->cr_cp_g_r * G + (1/m) * (-cb_sa * DEL_g[0].Roll_g  - ca * DEL_g[2].Roll_g);
    sg[2].Pitch_g = ea->cr_cp_g_p * G + (1/m) * (-cb_sa * DEL_g[0].Pitch_g  - ca * DEL_g[2].Pitch_g);
    sg[2].Yaw_g = (1/m) * (-cb_sa * DEL_g[0].Yaw_g  - ca * DEL_g[2].Yaw_g);    

    cg[2].elev_g = (1/m) * (-cb_sa * DEL_cg[0].elev_g  - ca * DEL_cg[2].elev_g);
    cg[2].thrust_g = -sphit / m;
    cg[2].aileron_g = (1/m) * (-cb_sa * DEL_cg[0].aileron_g  - ca * DEL_cg[2].aileron_g);
    cg[2].rudder_g = (1/m) * (-cb_sa * DEL_cg[0].rudder_g  - ca * DEL_cg[2].rudder_g);
    
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
    
    rates->v1 = (rt11 * QR + rt12 * PQ + Izz * L + Ixz * N) / rtn1;
    rates->v2 = (rt21 * P * R + rt22 * (R * R - P * P) + M) / rtn2;
    rates->v3 = (rt31 * PQ + rt32 * QR + Ixz * L + Ixx * N) / rtn3;
    
    return 0;
}

// Rotational Dynamics
int rdyn_g(const struct Vec3 * PQR,
           const struct Vec3 * LMN,
           const struct StateGrad LMN_g[3],
           const struct ControlGrad LMN_cg[3],
           const struct Aircraft * ac,
           struct Vec3 * rates,
           struct StateGrad sg[3],
           struct ControlGrad cg[3])
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

    sg[0].U_g = (Izz * LMN_g[0].U_g + Ixz * LMN_g[2].U_g) / rtn1;
    sg[0].V_g = (Izz * LMN_g[0].V_g + Ixz * LMN_g[2].V_g) / rtn1;
    sg[0].W_g = (Izz * LMN_g[0].W_g + Ixz * LMN_g[2].W_g) / rtn1;
    sg[0].Roll_g = (Izz * LMN_g[0].Roll_g + Ixz * LMN_g[2].Roll_g) / rtn1;
    sg[0].Pitch_g = (Izz * LMN_g[0].Pitch_g + Ixz * LMN_g[2].Pitch_g) / rtn1;
    sg[0].Yaw_g = (Izz * LMN_g[0].Yaw_g + Ixz * LMN_g[2].Yaw_g) / rtn1;

    cg[0].elev_g = (Izz * LMN_cg[0].elev_g + Ixz * LMN_cg[2].elev_g) / rtn1;
    cg[0].thrust_g = (Izz * LMN_cg[0].thrust_g + Ixz * LMN_cg[2].thrust_g) / rtn1;
    cg[0].aileron_g = (Izz * LMN_cg[0].aileron_g + Ixz * LMN_cg[2].aileron_g) / rtn1;
    cg[0].rudder_g = (Izz * LMN_cg[0].rudder_g + Ixz * LMN_cg[2].rudder_g) / rtn1;    
    
    sg[1].U_g = LMN_g[1].U_g / rtn2;
    sg[1].V_g = LMN_g[1].V_g / rtn2;
    sg[1].W_g = LMN_g[1].W_g / rtn2;
    sg[1].Roll_g = LMN_g[1].Roll_g / rtn2;
    sg[1].Pitch_g = LMN_g[1].Pitch_g / rtn2;
    sg[1].Yaw_g = LMN_g[1].Yaw_g / rtn2;

    cg[1].elev_g = LMN_cg[1].elev_g / rtn2;
    cg[1].thrust_g = LMN_cg[1].thrust_g / rtn2;
    cg[1].aileron_g = LMN_cg[1].aileron_g / rtn2;
    cg[1].rudder_g = LMN_cg[1].rudder_g / rtn2;

    
    sg[2].U_g = (Ixz * LMN_g[0].U_g + Ixx * LMN_g[2].U_g) / rtn3;
    sg[2].V_g = (Ixz * LMN_g[0].V_g + Ixx * LMN_g[2].V_g) / rtn3;
    sg[2].W_g = (Ixz * LMN_g[0].W_g + Ixx * LMN_g[2].W_g) / rtn3;
    sg[2].Roll_g = (Ixz * LMN_g[0].Roll_g + Ixx * LMN_g[2].Roll_g) / rtn3;
    sg[2].Pitch_g = (Ixz * LMN_g[0].Pitch_g + Ixx * LMN_g[2].Pitch_g) / rtn3;
    sg[2].Yaw_g = (Ixz * LMN_g[0].Yaw_g + Ixx * LMN_g[2].Yaw_g) / rtn3;

    cg[2].elev_g = (Ixz * LMN_cg[0].elev_g + Ixx * LMN_cg[2].elev_g) / rtn3;
    cg[2].thrust_g = (Ixz * LMN_cg[0].thrust_g + Ixx * LMN_cg[2].thrust_g) / rtn3;
    cg[2].aileron_g = (Ixz * LMN_cg[0].aileron_g + Ixx * LMN_cg[2].aileron_g) / rtn3;
    cg[2].rudder_g = (Ixz * LMN_cg[0].rudder_g + Ixx * LMN_cg[2].rudder_g) / rtn3;


    sg[0].P_g = (rt12 * Q + Izz * LMN_g[0].P_g + Ixz * LMN_g[2].P_g) / rtn1;
    sg[0].Q_g = (rt11 * R + rt12 * P + Izz * LMN_g[0].Q_g + Ixz * LMN_g[2].Q_g) / rtn1;
    sg[0].R_g = (rt11 * Q + Izz * LMN_g[0].R_g + Ixz * LMN_g[2].R_g) / rtn1;
    
    sg[1].P_g = (rt21 * R - 2 * rt22 * P + LMN_g[1].P_g) / rtn2;
    sg[1].Q_g = LMN_g[1].Q_g / rtn2;
    sg[1].R_g = (rt21 * P + 2 * rt22 * R + LMN_g[1].R_g) / rtn2;
    
    sg[2].P_g = (rt31 * Q + Ixz * LMN_g[0].P_g + Ixx * LMN_g[2].P_g) / rtn3;
    sg[2].Q_g = (rt31 * P + rt32 * R + Ixz * LMN_g[0].Q_g + Ixx * LMN_g[2].Q_g) / rtn3;
    sg[2].R_g = (rt32 * Q + Ixz * LMN_g[0].R_g + Ixx * LMN_g[2].R_g) / rtn3;    
    
    return 0;
}

int aero_forces(const struct AeroAngles * aero,
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

    real Cl = ac->Cl[0] * aero->sideslip +
        ac->Cl[1] * ac->span / 2.0 / vac * PQR->v1 +
        ac->Cl[2] * ac->span / 2.0 / vac * PQR->v3 +
        ac->Cl[3] * aero_con->v2 + // aileron
        ac->Cl[4] * aero_con->v3; // rudder

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

int aero_forces_g(const struct AeroAngles * aero,
                          const struct Vec3 * PQR, const struct Vec3 * aero_con,
                          const struct Aircraft * ac,
                          real rho,
                          real vac,
                          struct StateGrad * vac_g,
                          struct Vec3 * DEL,
                          struct Vec3 * LMN,
                          struct StateGrad sg[6], // gradient of D, E, L, L, M, N
                          struct ControlGrad cg[6]) 
{

    real term_cl_pre_q = ac->CL[4] * ac->chord / 2.0 / vac;
    real term_cl_q = term_cl_pre_q * PQR->v2;
    real term_cl_q_g_vac = -term_cl_q / vac;    
    // Forces
    real CL = ac->CL[0] +
        ac->CL[1] * aero->aoa +
        term_cl_q +
        ac->CL[5] * aero_con->v1;

    sg[2].U_g = ac->CL[1] * aero->aoa_g_u + term_cl_q_g_vac * vac_g->U_g;
    sg[2].V_g = term_cl_q_g_vac * vac_g->V_g;
    sg[2].W_g = ac->CL[1] * aero->aoa_g_w + term_cl_q_g_vac * vac_g->W_g;
    sg[2].P_g = 0.0;
    sg[2].Q_g = term_cl_pre_q; 
    sg[2].R_g = 0.0;
    sg[2].Roll_g = 0.0;
    sg[2].Pitch_g = 0.0;
    sg[2].Yaw_g = 0.0;

    cg[2].elev_g = ac->CL[5];
    cg[2].thrust_g = 0.0;
    cg[2].aileron_g = 0.0;
    cg[2].rudder_g = 0.0;
    
    real CD = ac->CD[0] + ac->K * CL * CL;

    real cd_pre_g = 2.0 * ac->K * CL;
    sg[0].U_g = cd_pre_g * sg[2].U_g;
    sg[0].V_g = cd_pre_g * sg[2].V_g;
    sg[0].W_g = cd_pre_g * sg[2].W_g;
    sg[0].P_g = cd_pre_g * sg[2].P_g;
    sg[0].Q_g = cd_pre_g * sg[2].Q_g;
    sg[0].R_g = cd_pre_g * sg[2].R_g;
    sg[0].Roll_g = cd_pre_g * sg[2].Roll_g;
    sg[0].Pitch_g = cd_pre_g * sg[2].Pitch_g;
    sg[0].Yaw_g = cd_pre_g * sg[2].Yaw_g;    

    cg[0].elev_g = cd_pre_g * cg[2].elev_g;
    cg[0].thrust_g = 0.0;
    cg[0].aileron_g = 0.0;
    cg[0].rudder_g = 0.0;
    

    sg[1].U_g = 0.0;
    sg[1].V_g = 0.0;
    sg[1].W_g = 0.0;
    sg[1].P_g = 0.0;
    sg[1].Q_g = 0.0;
    sg[1].R_g = 0.0;
    sg[1].Roll_g = 0.0;
    sg[1].Pitch_g = 0.0;
    sg[1].Yaw_g = 0.0;

    cg[1].elev_g = 0.0;
    cg[1].thrust_g = 0.0;
    cg[1].aileron_g = 0.0;
    cg[1].rudder_g = 0.0;
    
    real cm_pre_q = ac->Cm[4] / 2.0 / vac * ac->chord;
    real cm_q = cm_pre_q * PQR->v2;
    real cm_q_g_vac = -cm_q / vac;
    
    // Moments
    real Cm = ac->Cm[0] +
        ac->Cm[1] * aero->aoa + cm_q +
        ac->Cm[5] * aero_con->v1; // elevator

    sg[4].U_g = ac->Cm[1] * aero->aoa_g_u + cm_q_g_vac * vac_g->U_g;
    sg[4].V_g = cm_q_g_vac * vac_g->V_g;    
    sg[4].W_g = ac->Cm[1] * aero->aoa_g_w + cm_q_g_vac * vac_g->W_g;
    sg[4].P_g = 0.0;
    sg[4].Q_g = cm_pre_q;
    sg[4].R_g = 0.0;
    sg[4].Roll_g = 0.0;
    sg[4].Pitch_g = 0.0;
    sg[4].Yaw_g = 0.0;

    cg[4].elev_g = ac->Cm[5];
    cg[4].thrust_g = 0.0;
    cg[4].aileron_g = 0.0;
    cg[4].rudder_g = 0.0;    
    
    /* printf("vac = %3.2E\n", vac); */
    real mom_term = ac->span / 2.0 / vac;
    real mom_term_g_vac = - mom_term / vac;

    real pterm = ac->Cl[1] * PQR->v1;
    real rterm = ac->Cl[2] * PQR->v3;
    real prterm = pterm + rterm;
    real Cl = ac->Cl[0] * aero->sideslip +
        pterm * mom_term +
        rterm * mom_term +
        ac->Cl[3] * aero_con->v2 + // aileron
        ac->Cl[4] * aero_con->v3; // rudder

    
    sg[3].U_g = ac->Cl[0] * aero->sideslip_g_vac * vac_g->U_g + prterm * mom_term_g_vac * vac_g->U_g;
    sg[3].V_g = ac->Cl[0] * aero->sideslip_g_v  + prterm * mom_term_g_vac * vac_g->V_g;
    sg[3].W_g = ac->Cl[0] * aero->sideslip_g_vac * vac_g->W_g + prterm * mom_term_g_vac * vac_g->W_g;    
    sg[3].P_g = ac->Cl[1] * mom_term;
    sg[3].Q_g = 0.0;
    sg[3].R_g = ac->Cl[2] * mom_term;
    sg[3].Roll_g = 0.0;
    sg[3].Pitch_g = 0.0;
    sg[3].Yaw_g = 0.0;

    cg[3].elev_g = 0.0;
    cg[3].thrust_g = 0.0;
    cg[3].aileron_g = ac->Cl[3];
    cg[3].rudder_g = ac->Cl[4]; 


    pterm = ac->Cn[1] * PQR->v1;
    rterm = ac->Cn[2] * PQR->v3;
    prterm = pterm + rterm;    
    real Cn = ac->Cn[0] * aero->sideslip +
        pterm * mom_term + 
        rterm * mom_term +
        ac->Cn[3] * aero_con->v2 + // aileron
        ac->Cn[4] * aero_con->v3; // rudder    

    sg[5].U_g = ac->Cn[0] * aero->sideslip_g_vac * vac_g->U_g + prterm * mom_term_g_vac * vac_g->U_g;
    sg[5].V_g = ac->Cn[0] * aero->sideslip_g_v + prterm * mom_term_g_vac * vac_g->V_g;
    sg[5].W_g = ac->Cn[0] * aero->sideslip_g_vac * vac_g->W_g + prterm * mom_term_g_vac * vac_g->W_g;    
    sg[5].P_g = ac->Cn[1] * mom_term;
    sg[5].Q_g = 0.0;
    sg[5].R_g = ac->Cn[2] * mom_term;
    sg[5].Roll_g = 0.0;
    sg[5].Pitch_g = 0.0;
    sg[5].Yaw_g = 0.0;

    cg[5].elev_g = 0.0;
    cg[5].thrust_g = 0.0;
    cg[5].aileron_g = ac->Cn[3];
    cg[5].rudder_g = ac->Cn[4]; 

    
    real vac2 = vac * vac;
    
    real pre = 0.5 * rho * ac->area * vac2;
    real pre_mom = pre * ac->span;

    real pre_g_vac = rho * ac->area * vac;
    real pre_mom_g_vac = pre_g_vac * ac->span;

    DEL->v1 = pre * CD; // drag
    DEL->v2 = 0.0;      // side drag
    DEL->v3 = pre * CL; // lift

    sg[0].U_g = pre * sg[0].U_g + CD * pre_g_vac * vac_g->U_g;
    sg[0].V_g = pre * sg[0].V_g + CD * pre_g_vac * vac_g->V_g;
    sg[0].W_g = pre * sg[0].W_g + CD * pre_g_vac * vac_g->W_g;
    sg[0].P_g = pre * sg[0].P_g;
    sg[0].Q_g = pre * sg[0].Q_g;
    sg[0].R_g = pre * sg[0].R_g;
    sg[0].Roll_g = pre * sg[0].Roll_g;        
    sg[0].Pitch_g = pre * sg[0].Pitch_g;
    sg[0].Yaw_g = pre * sg[0].Yaw_g;

    cg[0].elev_g    = pre * cg[0].elev_g;
    cg[0].thrust_g  = pre * cg[0].thrust_g;
    cg[0].aileron_g = pre * cg[0].aileron_g;
    cg[0].rudder_g  = pre * cg[0].rudder_g;

    sg[2].U_g = pre * sg[2].U_g  + CL * pre_g_vac * vac_g->U_g;
    sg[2].V_g = pre * sg[2].V_g  + CL * pre_g_vac * vac_g->V_g;
    sg[2].W_g = pre * sg[2].W_g  + CL * pre_g_vac * vac_g->W_g;
    sg[2].P_g = pre * sg[2].P_g;
    sg[2].Q_g = pre * sg[2].Q_g;
    sg[2].R_g = pre * sg[2].R_g;
    sg[2].Roll_g  = pre * sg[2].Roll_g;
    sg[2].Pitch_g = pre * sg[2].Pitch_g;
    sg[2].Yaw_g   = pre * sg[2].Yaw_g;

    cg[2].elev_g    = pre * cg[2].elev_g;
    cg[2].thrust_g  = pre * cg[2].thrust_g;
    cg[2].aileron_g = pre * cg[2].aileron_g;
    cg[2].rudder_g  = pre * cg[2].rudder_g;  

    
    LMN->v1 = pre_mom * Cl;         // roll mom
    sg[3].U_g = pre_mom * sg[3].U_g  + Cl * pre_mom_g_vac * vac_g->U_g;
    sg[3].V_g = pre_mom * sg[3].V_g  + Cl * pre_mom_g_vac * vac_g->V_g;
    sg[3].W_g = pre_mom * sg[3].W_g  + Cl * pre_mom_g_vac * vac_g->W_g;
    sg[3].P_g = pre_mom * sg[3].P_g;
    sg[3].Q_g = pre_mom * sg[3].Q_g;
    sg[3].R_g = pre_mom * sg[3].R_g;
    sg[3].Roll_g = pre_mom * sg[3].Roll_g;
    sg[3].Pitch_g = pre_mom * sg[3].Pitch_g;
    sg[3].Yaw_g = pre_mom * sg[3].Yaw_g;

    cg[3].elev_g    = pre_mom * cg[3].elev_g;
    cg[3].thrust_g  = pre_mom * cg[3].thrust_g;
    cg[3].aileron_g = pre_mom * cg[3].aileron_g;
    cg[3].rudder_g  = pre_mom * cg[3].rudder_g;
    

    real pre_chord = pre * ac->chord;
    real pre_chord_g_vac = ac->chord * pre_g_vac;
    LMN->v2 = pre_chord * Cm; // pitch mom

    sg[4].U_g = pre_chord * sg[4].U_g + Cm * pre_chord_g_vac * vac_g->U_g;
    sg[4].V_g = pre_chord * sg[4].V_g + Cm * pre_chord_g_vac * vac_g->V_g;
    sg[4].W_g = pre_chord * sg[4].W_g + Cm * pre_chord_g_vac * vac_g->W_g;
    sg[4].P_g = pre_chord * sg[4].P_g;
    sg[4].Q_g = pre_chord * sg[4].Q_g;
    sg[4].R_g = pre_chord * sg[4].R_g;
    sg[4].Roll_g  = pre_chord * sg[4].Roll_g;
    sg[4].Pitch_g = pre_chord * sg[4].Pitch_g;
    sg[4].Yaw_g   = pre_chord * sg[4].Yaw_g;

    cg[4].elev_g    = pre_chord * cg[4].elev_g;
    cg[4].thrust_g  = pre_chord * cg[4].thrust_g;
    cg[4].aileron_g = pre_chord * cg[4].aileron_g;
    cg[4].rudder_g  = pre_chord * cg[4].rudder_g;

    
    LMN->v3 = pre_mom * Cn;         // yaw mom

    sg[5].U_g = pre_mom * sg[5].U_g  + Cn * pre_mom_g_vac * vac_g->U_g;
    sg[5].V_g = pre_mom * sg[5].V_g  + Cn * pre_mom_g_vac * vac_g->V_g;    
    sg[5].W_g = pre_mom * sg[5].W_g  + Cn * pre_mom_g_vac * vac_g->W_g;
    sg[5].P_g = pre_mom * sg[5].P_g;
    sg[5].Q_g = pre_mom * sg[5].Q_g;
    sg[5].R_g = pre_mom * sg[5].R_g;
    sg[5].Roll_g = pre_mom * sg[5].Roll_g;        
    sg[5].Pitch_g = pre_mom * sg[5].Pitch_g;
    sg[5].Yaw_g = pre_mom * sg[5].Yaw_g;

    cg[5].elev_g    = pre_mom * cg[5].elev_g;
    cg[5].thrust_g  = pre_mom * cg[5].thrust_g;
    cg[5].aileron_g = pre_mom * cg[5].aileron_g;
    cg[5].rudder_g  = pre_mom * cg[5].rudder_g;  
    
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
    aero_angles(&UVW, vac, &aero);

    #ifdef DEBUG
    printf("aoa = %3.2E\n", aero.aoa);
    printf("beta = %3.2E\n", aero.sideslip);
    #endif
    
    struct Vec3 PQR = {state[6], state[7], state[8]};
    struct EulerAngles ea;
    ea.roll = state[9];
    ea.pitch = state[10];
    ea.yaw = state[11];
    euler_angles(&ea);


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
    
    aero_forces(&aero, &PQR, &aero_con, ac, rho, vac, &DEL, &LMN);

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

int rigid_body_lin_forces_jac(double time,
                              const double * state,
                              const double * control,
                              double * out, double * jac,
                              void * arg)
{
    (void) time;
    /* assert (jac == NULL); */

    struct Aircraft * ac = arg;
    
    struct Vec3 UVW = {state[3], state[4], state[5]};
    real vac = vec3_norm(&UVW);
    struct StateGrad vac_g;
    vac_g.U_g = UVW.v1 / vac;
    vac_g.V_g = UVW.v2 / vac;
    vac_g.W_g = UVW.v3 / vac;
    
    struct AeroAngles aero;
    aero_angles_g(&UVW, vac, &vac_g, &aero);

    #ifdef DEBUG
    printf("aoa = %3.2E\n", aero.aoa);
    printf("beta = %3.2E\n", aero.sideslip);
    #endif
    
    struct Vec3 PQR = {state[6], state[7], state[8]};
    struct EulerAngles ea;
    ea.roll = state[9];
    ea.pitch = state[10];
    ea.yaw = state[11];
    euler_angles_g(&ea);


    /////////////////////////////
    // Kinematics
    ////////////////////////////
    struct StateGrad sg_tkin[3];
    struct Vec3 xyz_rates;
    tkin_g(&ea, &UVW, &xyz_rates, sg_tkin);

    struct StateGrad sg_rkin[3];    
    struct Vec3 euler_rates;
    rkin_g(&ea, &PQR, &euler_rates, sg_rkin);

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

    struct StateGrad forces[6];
    struct ControlGrad forces_c[6];
    aero_forces_g(&aero, &PQR, &aero_con, ac, rho, vac,
                          &vac_g, &DEL, &LMN, forces, forces_c);

    // Momentum equations
    struct Vec3 uvw_rates;
    struct StateGrad tmom[3];
    struct ControlGrad tmom_c[3];    
    tdyn_g(&ea, &aero, &vac_g, &UVW, &PQR, &DEL, forces, forces_c, Ft, ac, &uvw_rates, tmom, tmom_c);

    struct Vec3 pqr_rates;
    struct StateGrad rmom[3];
    struct ControlGrad rmom_c[3];
    rdyn_g(&PQR, &LMN, forces + 3, forces_c + 3, ac, &pqr_rates, rmom, rmom_c);

    // recall states are x y z U V W P Q R roll pitch yaw
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


    // jacobian is column major order
    // d / dx, d / dy, d/dz
    for (size_t ii = 0; ii < 36; ii++){
        jac[ii] = 0.0;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // U
        jac[36 + ii] = sg_tkin[ii].U_g;
        jac[39 + ii] = tmom[ii].U_g;
        jac[42 + ii] = rmom[ii].U_g;
        jac[45 + ii] = sg_rkin[ii].U_g;
    }
    
    for (size_t ii = 0; ii < 3; ii++){
        // V
        jac[48 + ii] = sg_tkin[ii].V_g;
        jac[51 + ii] = tmom[ii].V_g;
        jac[54 + ii] = rmom[ii].V_g;
        jac[57 + ii] = sg_rkin[ii].V_g;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // W
        jac[60 + ii] = sg_tkin[ii].W_g;
        jac[63 + ii] = tmom[ii].W_g;
        jac[66 + ii] = rmom[ii].W_g;
        jac[69 + ii] = sg_rkin[ii].W_g;        
    }

    for (size_t ii = 0; ii < 3; ii++){
        // P
        jac[72 + ii] = sg_tkin[ii].P_g;
        jac[75 + ii] = tmom[ii].P_g;
        jac[78 + ii] = rmom[ii].P_g;
        jac[81 + ii] = sg_rkin[ii].P_g;        
    }
    
    for (size_t ii = 0; ii < 3; ii++){
        // Q
        jac[84 + ii] = sg_tkin[ii].Q_g;
        jac[87 + ii] = tmom[ii].Q_g;
        jac[90 + ii] = rmom[ii].Q_g;
        jac[93 + ii] = sg_rkin[ii].Q_g;        
    }

    for (size_t ii = 0; ii < 3; ii++){
        // R
        jac[96 + ii] = sg_tkin[ii].R_g;
        jac[99 + ii] = tmom[ii].R_g;
        jac[102 + ii] = rmom[ii].R_g;
        jac[105 + ii] = sg_rkin[ii].R_g;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // roll
        jac[108 + ii] = sg_tkin[ii].Roll_g;
        jac[111 + ii] = tmom[ii].Roll_g;
        jac[114 + ii] = rmom[ii].Roll_g;
        jac[117 + ii] = sg_rkin[ii].Roll_g;        
    }

    for (size_t ii = 0; ii < 3; ii++){
        // pitch
        jac[120 + ii] = sg_tkin[ii].Pitch_g;
        jac[123 + ii] = tmom[ii].Pitch_g;
        jac[126 + ii] = rmom[ii].Pitch_g;
        jac[129 + ii] = sg_rkin[ii].Pitch_g;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // pitch
        jac[132 + ii] = sg_tkin[ii].Yaw_g;
        jac[135 + ii] = tmom[ii].Yaw_g;
        jac[138 + ii] = rmom[ii].Yaw_g;
        jac[141 + ii] = sg_rkin[ii].Yaw_g;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // elevator
        jac[144 + ii] = 0.0;
        jac[147 + ii] = tmom_c[ii].elev_g;
        jac[150 + ii] = rmom_c[ii].elev_g;
        jac[153 + ii] = 0.0;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // aileron
        jac[156 + ii] = 0.0;
        jac[159 + ii] = tmom_c[ii].aileron_g;
        jac[162 + ii] = rmom_c[ii].aileron_g;
        jac[165 + ii] = 0.0;
    }
    
    for (size_t ii = 0; ii < 3; ii++){
        // rudder
        jac[168 + ii] = 0.0;
        jac[171 + ii] = tmom_c[ii].rudder_g;
        jac[174 + ii] = rmom_c[ii].rudder_g;
        jac[177 + ii] = 0.0;
    }

    for (size_t ii = 0; ii < 3; ii++){
        // thrust
        jac[180 + ii] = 0.0;
        jac[183 + ii] = tmom_c[ii].thrust_g;
        jac[186 + ii] = rmom_c[ii].thrust_g;
        jac[189 + ii] = 0.0;
    }

    return 0;
}

int rigid_body_linearized(double time, const double * state,
               const double * control,
               double * out, double * jac,
               void * arg)
{
    (void) time;
    assert (jac == NULL);

    real * AB = arg;

    for (size_t ii = 0; ii < 12; ii++){
        out[ii] = 0.0;
    }
    
    for (size_t ii = 0; ii < 12; ii++){
        for (size_t jj = 0; jj < 12; jj++){
            out[jj] += (AB[ii*12 + jj] * state[ii]);
        }
    }

    real * B = AB + 144;
    for (size_t ii = 0; ii < 4; ii++){
        for (size_t jj = 0; jj < 12; jj++){
            out[jj] += (B[ii*12 + jj] * control[ii]);
        }
    }

    return 0;
}
