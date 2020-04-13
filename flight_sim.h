/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

#ifndef FLIGHT_SIM_H
#define FLIGHT_SIM_H

#include <nlopt.h>

typedef double real;
inline real rad2deg(real ang);

struct Vec3 {
    real v1;
    real v2;
    real v3;
};

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

struct ControlGrad
{
    real elev_g;
    real aileron_g;
    real rudder_g;
    real thrust_g;    
};

struct AeroAngles {

    real aoa;
    real sideslip;
    
    real aoa_g_w;
    real aoa_g_u;    
    real sideslip_g_v;
    real sideslip_g_vac;
};

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
inline real steady_state_get_U(const struct SteadyState * ss);
inline real steady_state_get_V(const struct SteadyState * ss);
inline real steady_state_get_W(const struct SteadyState * ss);
inline real steady_state_get_P(const struct SteadyState * ss);
inline real steady_state_get_Q(const struct SteadyState * ss);
inline real steady_state_get_R(const struct SteadyState * ss);
inline real steady_state_get_Roll(const struct SteadyState * ss);
inline real steady_state_get_Pitch(const struct SteadyState * ss);
inline real steady_state_get_elevator(const struct SteadyState * ss);
inline real steady_state_get_aileron(const struct SteadyState * ss);
inline real steady_state_get_rudder(const struct SteadyState * ss);
inline real steady_state_get_thrust(const struct SteadyState * ss);

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
int aircraft_inertia(struct Aircraft * ac);

int rigid_body_lin_forces(double time, const double * state,
                          const double * control,
                          double * out, double * jac,
                          void * arg);

int rigid_body_lin_forces_jac(double time,
                              const double * state,
                              const double * control,
                              double * out, double * jac,
                              void * arg);

int aero_angles(const struct Vec3* UVW, real Vac, struct AeroAngles * ab);
int aero_angles_g (const struct Vec3* UVW, real Vac,
                   const struct StateGrad * vac_g,
                   struct AeroAngles * ab);

int euler_angles(struct EulerAngles * ea);
int euler_angles_g(struct EulerAngles * ea);
int orient_ac_to_e(const struct EulerAngles *ea, const struct Vec3* ac, struct Vec3* e);
int orient_ac_to_e_g(const struct EulerAngles *ea, const struct Vec3* ac, struct Vec3* e,
                     struct StateGrad sg[3]);
int rkin(const struct EulerAngles * ea, const struct Vec3 * PQR, struct Vec3 * rates);
int rkin_g(const struct EulerAngles * ea, const struct Vec3 * PQR, struct Vec3 * rates, struct StateGrad sg[3]);
int tdyn(const struct EulerAngles * ea, const struct AeroAngles * aero,
         const struct Vec3 * UVW, const struct Vec3 * PQR, const struct Vec3 * DEL, real Ft,
         const struct Aircraft * ac, struct Vec3 * rates);
int tdyn_g(const struct EulerAngles * ea, const struct AeroAngles * aero, const struct StateGrad * vac_g,
           const struct Vec3 * UVW, const struct Vec3 * PQR, const struct Vec3 * DEL,
           const struct StateGrad DEL_g[3],
           const struct ControlGrad DEL_cg[3],
           real Ft,
           const struct Aircraft * ac,
           struct Vec3 * rates,
           struct StateGrad sg[3],
           struct ControlGrad cg[3]);
int rdyn(const struct Vec3 * PQR, const struct Vec3 * LMN, const struct Aircraft * ac,
         struct Vec3 * rates);
int rdyn_g(const struct Vec3 * PQR,
           const struct Vec3 * LMN,
           const struct StateGrad LMN_g[3],
           const struct ControlGrad LMN_cg[3],
           const struct Aircraft * ac,
           struct Vec3 * rates,
           struct StateGrad sg[3],
           struct ControlGrad cg[3]);
int aero_forces(const struct AeroAngles * aero,
                        const struct Vec3 * PQR, const struct Vec3 * aero_con,
                        const struct Aircraft * ac,
                        real rho,
                        real vac,
                        struct Vec3 * DEL,
                struct Vec3 * LMN);
int aero_forces_g(const struct AeroAngles * aero,
                          const struct Vec3 * PQR, const struct Vec3 * aero_con,
                          const struct Aircraft * ac,
                          real rho,
                          real vac,
                          struct StateGrad * vac_g,
                          struct Vec3 * DEL,
                          struct Vec3 * LMN,
                          struct StateGrad sg[6], // gradient of D, E, L, L, M, N
                  struct ControlGrad cg[6]);
    
#endif 
