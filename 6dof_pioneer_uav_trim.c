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

    return aircraft_inertia(ac);
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
    ss->achieved_climb_rate = sol[2];
    
    ss->target_yaw_rate = trim_spec_get_yaw_rate(data);
    ss->achieved_yaw_rate = sol[11];
    
    ss->target_speed = trim_spec_get_speed(data);
    ss->achieved_speed = sqrt(pow(x[0],2) + pow(x[1], 2) + pow(x[2], 2));

    aero_angles(&(ss->UVW), ss->achieved_speed, &(ss->aero));
    
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

    fprintf(stream, "Usage: %s options \n", program_name);
    fprintf(stream,
            " -h --help                Display this usage information.\n"
            " -s --speed      <val>    Desired speed (e.g., 120 --> flight at 120 ft/s, groundspeed), default 120.\n"
            " -c --climb-rate <val>    Desired climb rate (e.g., -5 --> climb at 5 ft/s), default 0.\n"
            " -y --yaw-rate   <val>    Desired turn rate (e.g., 3.14 --> turn 'right' at pi rad/s), default 0.\n"
            " -t --threshold  <val>    Threshold value for setting states to zero default 1e-10.\n"
            " -l --linearize  <int>    Return linear system. default 0\n"
            /* " -v --verbose    <val>      Output words (default 0)\n" */
        );
    exit (exit_code);
}

int main(int argc, char* argv[]){


    /* check_grad_ac_to_e(); */
    /* check_grad_rkin(); */
    /* check_grad_aero_forces(); */
    /* check_grad_rdyn(); */
    /* check_grad_tdyn(); */
    /* check_grad_rigid_body_dyn(); */
    /* check_grad_trimmer(); */
    /* return 0; */
    
    int next_option;
    const char * const short_options = "hs:c:y:t:l:";
    const struct option long_options[] = {
        { "help"       ,  0, NULL, 'h' },
        { "speed"      , 1, NULL, 's' },
        { "climb-rate" , 1, NULL, 'c' },
        { "yaw-rate"   , 1, NULL, 'y' },
        { "threshold"  , 1, NULL, 't' },
        { "linearize"  , 1, NULL, 'l' },        
        /* { "verbose"    , 1, NULL, 'v' }, */
        { NULL         , 0, NULL, 0   }
    };
    

    real speed = 120.0;
    real climb_rate = 0.0;
    real yaw_rate = 0.0;
    real thresh = 1e-10; // threshold for zero
    int linearize = 0;
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
                linearize = atoi(optarg);
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

    if (linearize > 0){

        real jac[144 + 48];
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

//------------------------------------------
//------------------------------------------
//--------- Gradient Checkers --------------
//------------------------------------------
//------------------------------------------

int check_grad_ac_to_e(void)
{
    struct EulerAngles ea;
    ea.roll = M_PI/9.0;
    ea.pitch = M_PI/8.0;
    ea.yaw = M_PI/6.0;
    euler_angles_g(&ea);
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
    euler_angles(&ea);
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
    euler_angles(&ea);
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
    euler_angles(&ea);
    orient_ac_to_e(&ea, &ac, &e);
    grad_Y[0] = (e.v1 - e_ref.v1) / h;
    grad_Y[1] = (e.v2 - e_ref.v2) / h;
    grad_Y[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Yaw\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Y[0], grad_Y[1], grad_Y[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Yaw_g, sg[1].Yaw_g, sg[2].Yaw_g);
    
    return 0;
}

int check_grad_rkin(void)
{
    struct EulerAngles ea;
    ea.roll = M_PI/9.0;
    ea.pitch = M_PI/8.0;
    ea.yaw = M_PI/6.0;
    euler_angles_g(&ea);
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
    euler_angles(&ea);
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
    euler_angles(&ea);
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
    euler_angles(&ea);
    rkin(&ea, &pqr, &rates);
    grad_Yaw[0] = (rates.v1 - e_ref.v1) / h;
    grad_Yaw[1] = (rates.v2 - e_ref.v2) / h;
    grad_Yaw[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to Yaw\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_Yaw[0], grad_Yaw[1], grad_Yaw[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].Yaw_g, sg[1].Yaw_g, sg[2].Yaw_g);
    
    return 0;
}

int check_grad_trimmer(void)
{
    real input[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                      7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

    // elevator aileron rudder thrust
    real out_ref;
    real out;
    real grad[12];

    struct Aircraft aircraft;
    pioneer_uav(&aircraft);

    struct TrimSpec trim_spec;
    trim_spec.z_dot = 5.0;
    trim_spec.yaw_dot = 3.0; 
    trim_spec.target_vel = 120.0;
    trim_spec.ac = &aircraft;
    trim_spec.thresh = 1e-14;     

    out_ref = trim_objective(12, input, NULL, &trim_spec);
    out = trim_objective_g(12, input, grad, &trim_spec);
    printf("Checking value computation\n");
    printf("%10.5E %10.5E \n", out_ref, out);

    printf("Checking gradient computation\n");
    real h = 1e-8;
    for (size_t ii = 0; ii < 12; ii++){
        input[ii] += h;
        out = trim_objective(12, input, NULL, &trim_spec);
        real val = (out - out_ref) / h;        
        printf("\t %10.5E %10.5E\n", val, grad[ii]);
        input[ii] -= h;
    }


    return 0;
}

int check_grad_aero_forces(void)
{
    struct Vec3 uvw = {120, 90, 60};
    real vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    struct StateGrad vac_g;
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
   
    
    struct Vec3 pqr = {10, -2.0, 8.0};
    struct Vec3 aero_con = {3.0, 1.0, 2.5};
    real rho = 0.002376892406675; // slug / ft^3

    struct Aircraft aircraft;
    pioneer_uav(&aircraft);
    
    struct Vec3 DEL_ref;
    struct Vec3 LMN_ref;

    struct Vec3 DEL;
    struct Vec3 LMN;    

    struct StateGrad sg[6];
    struct ControlGrad cg[6];    

    struct AeroAngles aero;

    // Reference    
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL_ref, &LMN_ref);

    // Analytic
    aero_angles_g(&uvw, vac, &vac_g, &aero);
    aero_forces_g(&aero, &pqr, &aero_con, &aircraft, rho, vac, &vac_g, &DEL, &LMN, sg, cg);

    printf("Checking value computation\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", DEL_ref.v1, DEL_ref.v2, DEL_ref.v3, LMN_ref.v1, LMN_ref.v2, LMN_ref.v3);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", DEL.v1, DEL.v2, DEL.v3, LMN.v1, LMN.v2, LMN.v3);    
    printf("\n\n\n");
    
    double h = 1e-6;
    
    double grad_U[6];
    uvw.v1 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_U[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_U[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_U[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_U[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_U[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_U[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to U\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_U[0], grad_U[1], grad_U[2], grad_U[3], grad_U[4], grad_U[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", sg[0].U_g, sg[1].U_g, sg[2].U_g, sg[3].U_g, sg[4].U_g, sg[5].U_g);
    printf("\n\n");
    uvw.v1 -= h;

        
    double grad_V[6];
    uvw.v2 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);    
    grad_V[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_V[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_V[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_V[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_V[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_V[5] = (LMN.v3 - LMN_ref.v3) / h;    

    printf("Checking Gradient with respect to V\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_V[0], grad_V[1], grad_V[2], grad_V[3], grad_V[4], grad_V[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", sg[0].V_g, sg[1].V_g, sg[2].V_g, sg[3].V_g, sg[4].V_g, sg[5].V_g);
    printf("\n\n");
    uvw.v2 -= h;

    double grad_W[6];
    uvw.v3 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_W[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_W[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_W[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_W[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_W[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_W[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to W\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_W[0], grad_W[1], grad_W[2], grad_W[3], grad_W[4], grad_W[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", sg[0].W_g, sg[1].W_g, sg[2].W_g, sg[3].W_g, sg[4].W_g, sg[5].W_g);
    printf("\n\n");
    uvw.v3 -= h;

    // reset
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);


    double grad_P[6];
    pqr.v1 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_P[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_P[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_P[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_P[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_P[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_P[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to P\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_P[0], grad_P[1], grad_P[2], grad_P[3], grad_P[4], grad_P[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", sg[0].P_g, sg[1].P_g, sg[2].P_g, sg[3].P_g, sg[4].P_g, sg[5].P_g);
    printf("\n\n");
    pqr.v1 -= h;
    
    double grad_Q[6];
    pqr.v2 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_Q[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_Q[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_Q[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_Q[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_Q[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_Q[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to Q\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_Q[0], grad_Q[1], grad_Q[2], grad_Q[3], grad_Q[4], grad_Q[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", sg[0].Q_g, sg[1].Q_g, sg[2].Q_g, sg[3].Q_g, sg[4].Q_g, sg[5].Q_g);
    printf("\n\n");
    pqr.v2 -= h;

    double grad_R[6];
    pqr.v3 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_R[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_R[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_R[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_R[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_R[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_R[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to R\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_R[0], grad_R[1], grad_R[2], grad_R[3], grad_R[4], grad_R[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", sg[0].R_g, sg[1].R_g, sg[2].R_g, sg[3].R_g, sg[4].R_g, sg[5].R_g);
    printf("\n\n");
    pqr.v3 -= h;

    double grad_el[6];
    aero_con.v1 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_el[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_el[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_el[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_el[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_el[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_el[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to Elevator\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_el[0], grad_el[1], grad_el[2], grad_el[3], grad_el[4], grad_el[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", cg[0].elev_g, cg[1].elev_g, cg[2].elev_g, cg[3].elev_g, cg[4].elev_g, cg[5].elev_g);
    printf("\n\n");
    aero_con.v1 -= h;

    double grad_ail[6];
    aero_con.v2 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_ail[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_ail[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_ail[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_ail[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_ail[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_ail[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to Aileron\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_ail[0], grad_ail[1], grad_ail[2], grad_ail[3], grad_ail[4], grad_ail[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", cg[0].aileron_g, cg[1].aileron_g, cg[2].aileron_g, cg[3].aileron_g, cg[4].aileron_g, cg[5].aileron_g);
    printf("\n\n");
    aero_con.v2 -= h;

    double grad_rud[6];
    aero_con.v3 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    grad_rud[0] = (DEL.v1 - DEL_ref.v1) / h;
    grad_rud[1] = (DEL.v2 - DEL_ref.v2) / h;
    grad_rud[2] = (DEL.v3 - DEL_ref.v3) / h;
    grad_rud[3] = (LMN.v1 - LMN_ref.v1) / h;
    grad_rud[4] = (LMN.v2 - LMN_ref.v2) / h;
    grad_rud[5] = (LMN.v3 - LMN_ref.v3) / h;

    printf("Checking Gradient with respect to Rudder\n");
    printf("Numerical %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", grad_rud[0], grad_rud[1], grad_rud[2], grad_rud[3], grad_rud[4], grad_rud[5]);
    printf("Analytic %3.5E %3.5E %3.5E %3.5E %3.5E %3.5E\n", cg[0].rudder_g, cg[1].rudder_g, cg[2].rudder_g, cg[3].rudder_g, cg[4].rudder_g, cg[5].rudder_g);
    printf("\n\n");
    aero_con.v3 -= h;
    
    return 0;
}

int check_grad_tdyn(void)
{
    struct Vec3 uvw = {120, 90, 60};
    real vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));

    struct StateGrad vac_g;
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
   
    
    struct Vec3 pqr = {10, -2.0, 8.0};
    struct Vec3 aero_con = {3.0, 1.0, 2.5};
    real Ft = 48.0; // force
    real rho = 0.002376892406675; // slug / ft^3

    struct Aircraft aircraft;
    pioneer_uav(&aircraft);
    
    struct Vec3 DEL_ref;
    struct Vec3 LMN_ref;
    struct Vec3 rates_ref;

    struct Vec3 DEL;
    struct Vec3 LMN;    
    struct Vec3 rates;
    
    struct StateGrad sg_f[6];
    struct ControlGrad cg_f[6];
    struct StateGrad sg[3];
    struct ControlGrad cg[3];        

    struct AeroAngles aero;

    struct EulerAngles ea;
    ea.roll = M_PI/9.0;
    ea.pitch = M_PI/8.0;
    ea.yaw = M_PI/6.0;
    euler_angles_g(&ea);
    
    // Reference    
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac,
                        &DEL_ref, &LMN_ref);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL_ref, Ft, &aircraft, &rates_ref);

    // Analytic
    aero_angles_g(&uvw, vac, &vac_g, &aero);
    aero_forces_g(&aero, &pqr, &aero_con, &aircraft, rho, vac,
                          &vac_g, &DEL, &LMN, sg_f, cg_f);
    tdyn_g(&ea, &aero, &vac_g, &uvw, &pqr, &DEL, sg_f, cg_f, Ft, &aircraft, &rates, sg, cg);

    printf("Checking value computation\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", rates_ref.v1, rates_ref.v2, rates_ref.v3);
    printf("Analytic %3.5E %3.5E %3.5E\n", rates.v1, rates.v2, rates.v3);
    printf("\n\n\n");
    
    double h = 1e-8;
        
    double grad_U[3];
    uvw.v1 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    
    grad_U[0] = (rates.v1 - rates_ref.v1) / h;
    grad_U[1] = (rates.v2 - rates_ref.v2) / h;
    grad_U[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to U\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_U[0], grad_U[1], grad_U[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].U_g, sg[1].U_g, sg[2].U_g);
    printf("\n\n");
    uvw.v1 -= h;

    double grad_V[3];
    uvw.v2 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);

    grad_V[0] = (rates.v1 - rates_ref.v1) / h;
    grad_V[1] = (rates.v2 - rates_ref.v2) / h;
    grad_V[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to V\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_V[0], grad_V[1], grad_V[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].V_g, sg[1].V_g, sg[2].V_g);
    printf("\n\n");
    uvw.v2 -= h;

    double grad_W[3];
    uvw.v3 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);

    grad_W[0] = (rates.v1 - rates_ref.v1) / h;
    grad_W[1] = (rates.v2 - rates_ref.v2) / h;
    grad_W[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to W\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_W[0], grad_W[1], grad_W[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].W_g, sg[1].W_g, sg[2].W_g);
    printf("\n\n");
    uvw.v3 -= h;

    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    
    double grad_P[3];
    pqr.v1 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);

    grad_P[0] = (rates.v1 - rates_ref.v1) / h;
    grad_P[1] = (rates.v2 - rates_ref.v2) / h;
    grad_P[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to P\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_P[0], grad_P[1], grad_P[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].P_g, sg[1].P_g, sg[2].P_g);
    printf("\n\n");
    pqr.v1 -= h;

    double grad_Q[3];
    pqr.v2 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    grad_Q[0] = (rates.v1 - rates_ref.v1) / h;
    grad_Q[1] = (rates.v2 - rates_ref.v2) / h;
    grad_Q[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to Q\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_Q[0], grad_Q[1], grad_Q[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].Q_g, sg[1].Q_g, sg[2].Q_g);
    printf("\n\n");
    pqr.v2 -= h;

    double grad_R[3];
    pqr.v3 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    grad_R[0] = (rates.v1 - rates_ref.v1) / h;
    grad_R[1] = (rates.v2 - rates_ref.v2) / h;
    grad_R[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to R\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_R[0], grad_R[1], grad_R[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].R_g, sg[1].R_g, sg[2].R_g);
    printf("\n\n");
    pqr.v3 -= h;

    double grad_Roll[3];
    ea.roll += h;
    euler_angles(&ea);
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    
    grad_Roll[0] = (rates.v1 - rates_ref.v1) / h;
    grad_Roll[1] = (rates.v2 - rates_ref.v2) / h;
    grad_Roll[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to Roll\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_Roll[0], grad_Roll[1], grad_Roll[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].Roll_g, sg[1].Roll_g, sg[2].Roll_g);
    printf("\n\n");
    ea.roll -= h;

    double grad_Pitch[3];
    ea.pitch += h;
    euler_angles(&ea);
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    
    grad_Pitch[0] = (rates.v1 - rates_ref.v1) / h;
    grad_Pitch[1] = (rates.v2 - rates_ref.v2) / h;
    grad_Pitch[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to Pitch\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_Pitch[0], grad_Pitch[1], grad_Pitch[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].Pitch_g, sg[1].Pitch_g, sg[2].Pitch_g);
    printf("\n\n");
    ea.pitch -= h;    

    double grad_Yaw[3];
    ea.yaw += h;
    euler_angles(&ea);
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    
    grad_Yaw[0] = (rates.v1 - rates_ref.v1) / h;
    grad_Yaw[1] = (rates.v2 - rates_ref.v2) / h;
    grad_Yaw[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to Yaw\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_Yaw[0], grad_Yaw[1], grad_Yaw[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].Yaw_g, sg[1].Yaw_g, sg[2].Yaw_g);
    printf("\n\n");
    ea.yaw -= h;
    euler_angles(&ea);

    double grad_el[6];
    aero_con.v1 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);
    grad_el[0] = (rates.v1 - rates_ref.v1) / h;
    grad_el[1] = (rates.v2 - rates_ref.v2) / h;
    grad_el[2] = (rates.v3 - rates_ref.v3) / h;

    printf("Checking Gradient with respect to Elevator\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_el[0], grad_el[1], grad_el[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", cg[0].elev_g, cg[1].elev_g, cg[2].elev_g);
    printf("\n\n");
    aero_con.v1 -= h;

    double grad_ail[6];
    aero_con.v2 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);

    grad_ail[0] = (rates.v1 - rates_ref.v1) / h;
    grad_ail[1] = (rates.v2 - rates_ref.v2) / h;
    grad_ail[2] = (rates.v3 - rates_ref.v3) / h;

    printf("Checking Gradient with respect to Aileron\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", grad_ail[0], grad_ail[1], grad_ail[2]);
    printf("Analytic %3.5E %3.5E %3.5E\n", cg[0].aileron_g, cg[1].aileron_g, cg[2].aileron_g);
    printf("\n\n");
    aero_con.v2 -= h;

    double grad_rud[6];
    aero_con.v3 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);

    grad_rud[0] = (rates.v1 - rates_ref.v1) / h;
    grad_rud[1] = (rates.v2 - rates_ref.v2) / h;
    grad_rud[2] = (rates.v3 - rates_ref.v3) / h;

    printf("Checking Gradient with respect to Rudder\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", grad_rud[0], grad_rud[1], grad_rud[2]);
    printf("Analytic %3.5E %3.5E %3.5E\n", cg[0].rudder_g, cg[1].rudder_g, cg[2].rudder_g);
    printf("\n\n");
    aero_con.v3 -= h;

    double grad_thrust[6];
    Ft += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    tdyn(&ea, &aero, &uvw, &pqr, &DEL, Ft, &aircraft, &rates);

    grad_thrust[0] = (rates.v1 - rates_ref.v1) / h;
    grad_thrust[1] = (rates.v2 - rates_ref.v2) / h;
    grad_thrust[2] = (rates.v3 - rates_ref.v3) / h;

    printf("Checking Gradient with respect to Thrust\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", grad_thrust[0], grad_thrust[1], grad_thrust[2]);
    printf("Analytic %3.5E %3.5E %3.5E\n", cg[0].thrust_g, cg[1].thrust_g, cg[2].thrust_g);
    printf("\n\n");
    Ft -= h;
    
    return 0;
}

int check_grad_rdyn(void)
{
    struct Vec3 uvw = {120, 90, 60};
    real vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));

    struct StateGrad vac_g;
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
   
    
    struct Vec3 pqr = {10, -2.0, 8.0};
    struct Vec3 aero_con = {3.0, 1.0, 2.5};
    real rho = 0.002376892406675; // slug / ft^3

    struct Aircraft aircraft;
    pioneer_uav(&aircraft);
    
    struct Vec3 DEL_ref;
    struct Vec3 LMN_ref;
    struct Vec3 rates_ref;

    struct Vec3 DEL;
    struct Vec3 LMN;    
    struct Vec3 rates;
    
    struct StateGrad sg_f[6];
    struct ControlGrad cg_f[6];
    struct StateGrad sg[3];
    struct ControlGrad cg[3];        

    struct AeroAngles aero;

    // Reference
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac,
                        &DEL_ref, &LMN_ref);
    rdyn(&pqr, &LMN_ref, &aircraft, &rates_ref);

    // Analytic
    aero_angles_g(&uvw, vac, &vac_g, &aero);
    aero_forces_g(&aero, &pqr, &aero_con, &aircraft, rho, vac,
                          &vac_g, &DEL, &LMN, sg_f, cg_f);
    rdyn_g(&pqr, &LMN, sg_f + 3, cg_f + 3, &aircraft, &rates, sg, cg);

    printf("Checking value computation\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", rates_ref.v1, rates_ref.v2, rates_ref.v3);
    printf("Analytic %3.5E %3.5E %3.5E\n", rates.v1, rates.v2, rates.v3);
    printf("\n\n\n");
    
    double h = 1e-6;
    
    double grad_U[3];
    uvw.v1 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_U[0] = (rates.v1 - rates_ref.v1) / h;
    grad_U[1] = (rates.v2 - rates_ref.v2) / h;
    grad_U[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to U\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_U[0], grad_U[1], grad_U[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].U_g, sg[1].U_g, sg[2].U_g);
    printf("\n\n");
    uvw.v1 -= h;

    double grad_V[3];
    uvw.v2 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_V[0] = (rates.v1 - rates_ref.v1) / h;
    grad_V[1] = (rates.v2 - rates_ref.v2) / h;
    grad_V[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to V\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_V[0], grad_V[1], grad_V[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].V_g, sg[1].V_g, sg[2].V_g);
    printf("\n\n");
    uvw.v2 -= h;

    double grad_W[3];
    uvw.v3 += h;
    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_W[0] = (rates.v1 - rates_ref.v1) / h;
    grad_W[1] = (rates.v2 - rates_ref.v2) / h;
    grad_W[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to W\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_W[0], grad_W[1], grad_W[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].W_g, sg[1].W_g, sg[2].W_g);
    printf("\n\n");
    uvw.v3 -= h;

    vac = sqrt(pow(uvw.v1, 2) + pow(uvw.v2, 2) + pow(uvw.v3, 2));
    vac_g.U_g = uvw.v1 / vac;
    vac_g.V_g = uvw.v2 / vac;
    vac_g.W_g = uvw.v3 / vac;
    aero_angles(&uvw, vac, &aero);
    
    double grad_P[3];
    pqr.v1 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_P[0] = (rates.v1 - rates_ref.v1) / h;
    grad_P[1] = (rates.v2 - rates_ref.v2) / h;
    grad_P[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to P\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_P[0], grad_P[1], grad_P[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].P_g, sg[1].P_g, sg[2].P_g);
    printf("\n\n");
    pqr.v1 -= h;

    double grad_Q[3];
    pqr.v2 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_Q[0] = (rates.v1 - rates_ref.v1) / h;
    grad_Q[1] = (rates.v2 - rates_ref.v2) / h;
    grad_Q[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to Q\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_Q[0], grad_Q[1], grad_Q[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].Q_g, sg[1].Q_g, sg[2].Q_g);
    printf("\n\n");
    pqr.v2 -= h;

    double grad_R[3];
    pqr.v3 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_R[0] = (rates.v1 - rates_ref.v1) / h;
    grad_R[1] = (rates.v2 - rates_ref.v2) / h;
    grad_R[2] = (rates.v3 - rates_ref.v3) / h;
    
    printf("Checking Gradient with respect to R\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_R[0], grad_R[1], grad_R[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", sg[0].R_g, sg[1].R_g, sg[2].R_g);
    printf("\n\n");
    pqr.v3 -= h;
    
    double grad_el[6];
    aero_con.v1 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_el[0] = (rates.v1 - rates_ref.v1) / h;
    grad_el[1] = (rates.v2 - rates_ref.v2) / h;
    grad_el[2] = (rates.v3 - rates_ref.v3) / h;


    printf("Checking Gradient with respect to Elevator\n");
    printf("Numerical %3.5E %3.5E %3.5E \n", grad_el[0], grad_el[1], grad_el[2]);
    printf("Analytic %3.5E %3.5E %3.5E \n", cg[0].elev_g, cg[1].elev_g, cg[2].elev_g);
    printf("\n\n");
    aero_con.v1 -= h;

    double grad_ail[6];
    aero_con.v2 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);
    grad_ail[0] = (rates.v1 - rates_ref.v1) / h;
    grad_ail[1] = (rates.v2 - rates_ref.v2) / h;
    grad_ail[2] = (rates.v3 - rates_ref.v3) / h;

    printf("Checking Gradient with respect to Aileron\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", grad_ail[0], grad_ail[1], grad_ail[2]);
    printf("Analytic %3.5E %3.5E %3.5E\n", cg[0].aileron_g, cg[1].aileron_g, cg[2].aileron_g);
    printf("\n\n");
    aero_con.v2 -= h;

    double grad_rud[6];
    aero_con.v3 += h;
    aero_forces(&aero, &pqr, &aero_con, &aircraft, rho, vac, &DEL, &LMN);
    rdyn(&pqr, &LMN, &aircraft, &rates);    
    grad_rud[0] = (rates.v1 - rates_ref.v1) / h;
    grad_rud[1] = (rates.v2 - rates_ref.v2) / h;
    grad_rud[2] = (rates.v3 - rates_ref.v3) / h;

    printf("Checking Gradient with respect to Rudder\n");
    printf("Numerical %3.5E %3.5E %3.5E\n", grad_rud[0], grad_rud[1], grad_rud[2]);
    printf("Analytic %3.5E %3.5E %3.5E\n", cg[0].rudder_g, cg[1].rudder_g, cg[2].rudder_g);
    printf("\n\n");
    aero_con.v3 -= h;
    
    return 0;
}

int check_grad_rigid_body_dyn(void)
{
    real time = 0.0;
    real state[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                      7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

    // elevator aileron rudder thrust
    real control[4] = {1.5, 2.5, 3.5, 4.5};

    real out_ref[12];
    real out[12];
    real jac[144+48];

    struct Aircraft aircraft;
    pioneer_uav(&aircraft);

    rigid_body_lin_forces(time, state, control, out_ref, NULL, &aircraft);
    rigid_body_lin_forces_jac(time, state, control, out, jac, &aircraft);

    printf("Checking value computation\n");
    for (size_t ii = 0; ii < 12; ii++){
        printf("%10.5E %10.5E %10.5E \n", out_ref[ii], out[ii],
               out[ii] - out_ref[ii]);
    }

    printf("Checking gradient computation\n");    
    real h = 1e-8;
    for (size_t ii = 0; ii < 12; ii++){
        state[ii] += h;
        rigid_body_lin_forces(time, state, control, out, NULL, &aircraft);
        printf("\n\n\n");
        printf("ii = %zu\n", ii);
        for (size_t jj = 0; jj < 12; jj++){
            real val = (out[jj] - out_ref[jj]) / h;
            printf("\t %10.5E %10.5E %10.5E\n", val, jac[ii*12+jj], (val - jac[ii*12 + jj])/val);
        }
        state[ii] -= h;
    }

    printf("\n");
    printf("--------------------------\n");
    printf("Gradient wrt controls\n");
    printf("--------------------------\n");
    printf("\n");
    for (size_t ii = 0; ii < 4; ii++){
        control[ii] += h;
        rigid_body_lin_forces(time, state, control, out, NULL, &aircraft);
        printf("\n\n\n");
        printf("ii = %zu\n", ii);
        for (size_t jj = 0; jj < 12; jj++){
            real val = (out[jj] - out_ref[jj]) / h;
            printf("\t %10.5E %10.5E %10.5E\n", val, jac[144 + ii*12+jj],
                   (val - jac[144 + ii*12 + jj])/val);
        }
        control[ii] -= h;
    }    
    

    return 0;
}
