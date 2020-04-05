#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

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
};

int compute_aero_angles(const struct Vec3* UVW, real Vac, struct AeroAngles * ab)
{
    ab->aoa = atan2(UVW->v3, UVW->v1);
    ab->sideslip = asin(UVW->v2/Vac);
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

    return 0;
}

int orient_e_to_ac(const struct EulerAngles *ea, const struct Vec3* e, struct Vec3* ac)
{
    // yaw, pitch, roll
    real sr_sp = ea->sr * ea->sp; // sin_roll * sin_pitch;
    real sr_cp = ea->sr * ea->cp; // sin_roll * cos_pitch;
    real sr_sy = ea->sr * ea->sy; // sin_roll * sin yaw
    real sr_cy = ea->sr * ea->cy; // sin_roll * cos yaw    
    real cr_sp = ea->cr * ea->sp; // cos_roll * sin_pitch
    real cr_cp = ea->cr * ea->cp; // cos_roll * cos_pitch
    real cr_cy = ea->cr * ea->cy; // cos_roll * cos_yaw
    real cr_sy = ea->cr * ea->sy; // cos_roll * sin_yaw
    real cp_cy = ea->cp * ea->cy; // cos_pitch * cos_yaw
    real cp_sy = ea->cp * ea->sy; // cos_pitch * sin_yaw
    
    ac->v1 = cp_cy * e->v1 +
             cp_sy * e->v2 -
             ea->sp * e->v3;

    ac->v2 = (sr_sp * ea->cy - cr_sy) * e->v1 +
             (sr_sp * ea->sy + cr_cy) * e->v2 +
                                sr_cp * e->v3;

    ac->v3 = (cr_sp * ea->cy + sr_sy) * e->v1 +
             (cr_sp * ea->sy - sr_cy) * e->v2 +
                                cr_cp * e->v3;

    return 0;
}

int orient_ac_to_e(const struct EulerAngles *ea, const struct Vec3* ac, struct Vec3* e)
{
    // -roll, -pitch, -yaw
    real sr_sp = ea->sr * ea->sp; // sin_roll * sin_pitch;
    real sr_cp = ea->sr * ea->cp; // sin_roll * cos_pitch;
    real sr_sy = ea->sr * ea->sy; // sin_roll * sin yaw
    real sr_cy = ea->sr * ea->cy; // sin_roll * cos yaw    
    real cr_sp = ea->cr * ea->sp; // cos_roll * sin_pitch
    real cr_cp = ea->cr * ea->cp; // cos_roll * cos_pitch
    real cr_cy = ea->cr * ea->cy; // cos_roll * cos_yaw
    real cr_sy = ea->cr * ea->sy; // cos_roll * sin_yaw
    real cp_cy = ea->cp * ea->cy; // cos_pitch * cos_yaw
    real cp_sy = ea->cp * ea->sy; // cos_pitch * sin_yaw

    e->v1 =                cp_cy * ac->v1 +
        (sr_sp * ea->cy - cr_sy) * ac->v2 + //(sr_sp * ea->cy - cr_sy)
        (cr_sp * ea->cy + sr_sy) * ac->v3; //(cr_sp * ea->cy + sr_sy)

    e->v2 =                cp_sy * ac->v1 +
        (sr_sp * ea->sy + cr_cy) * ac->v2 +  //            (sr_sp * ea->sy + cr_cy) * e->v2 +
        (cr_sp * ea->sy - sr_cy) * ac->v3;  //             (cr_sp * ea->sy - sr_cy) * e->v2 +

    e->v3 = - ea->sp * ac->v1 +
               sr_cp * ac->v2 +
               cr_cp * ac->v3;

    return 0;
}

// Translational kinematics
inline int tkin(const struct EulerAngles * ea, const struct Vec3 *UVW, struct Vec3 * rates){return orient_ac_to_e(ea, UVW, rates);}

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
    double z_dot;
    double yaw_dot;
    double target_vel;
    
    struct Aircraft * ac;
};

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
    struct Vec3 UVW;
    struct Vec3 PQR;
    struct Vec3 aero_con;
    double roll, pitch;
    double thrust;

    nlopt_result res;
    double obj_val;
};

int steady_state_print(FILE * fp, const struct SteadyState * ss)
{
    fprintf(fp, "\n");
    fprintf(fp, "========================================================\n");
    fprintf(fp, "                        TRIM RESULT                     \n");
    fprintf(fp, "========================================================\n");
    fprintf(fp, "Optimizer result = %d\nObjective value = %3.5E\n", ss->res, ss->obj_val);
    
    fprintf(fp, "\n\n\n");
    fprintf(fp, "Steady-state state values\n");
    fprintf(fp, "-------------------------\n");    
    fprintf(fp, "(U, V, W)     = (%3.5E \t %3.5E \t %3.5E) \n", ss->UVW.v1, ss->UVW.v2, ss->UVW.v3);
    fprintf(fp, "(P, Q, R)     = (%3.5E \t %3.5E \t %3.5E) \n", ss->PQR.v1, ss->PQR.v2, ss->PQR.v3);    
    fprintf(fp, "(Roll, Pitch) = (%3.5E \t %3.5E) \n", ss->roll, ss->pitch);

    fprintf(fp, "\n\n\n");
    fprintf(fp, "Steady-state control values\n");
    fprintf(fp, "-------------------------\n");    
    fprintf(fp, "(Elev, Ail, Rud) = (%3.5E \t %3.5E \t %3.5E) \n",
            ss->aero_con.v1, ss->aero_con.v2, ss->aero_con.v3);
    fprintf(fp, "Thrust = %3.5E\n", ss->thrust);
    fprintf(fp, "\n");
    fprintf(fp, "========================================================\n");
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

    // run without bounds
    /* opt = nlopt_create(NLOPT_LN_NELDERMEAD, 12); */
    /* opt = nlopt_create(NLOPT_LN_SBPLX, 12);     */
    /* nlopt_set_ftol_rel(opt, -1.0); */
    /* nlopt_set_ftol_abs(opt, 1e-14); */
    /* nlopt_set_min_objective(opt, trim_objective, data); */
    /* res = nlopt_optimize(opt, x, &val); */
    /* nlopt_destroy(opt); */

    // run with bounds (again)
    /* opt = nlopt_create(NLOPT_LN_NELDERMEAD, 12); */
    /* opt = nlopt_create(NLOPT_LN_SBPLX, 12); */
    opt = nlopt_create(NLOPT_LN_NEWUOA, 12);
    nlopt_set_ftol_rel(opt, -1.0);
    nlopt_set_ftol_abs(opt, 0.0);
    nlopt_set_min_objective(opt, trim_objective, data);

    // no sideslip
    double lb[12], ub[12];
    for (size_t ii = 0; ii < 12; ii++){
        lb[ii] = -HUGE_VAL;
        ub[ii] = HUGE_VAL;
    }
    lb[11] = 0.0; // lower bound thrust
    lb[1] = 0.0; // lower bound V
    ub[1] = 0.0; // upper bound V
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    x[1] = 0.0;
    res = nlopt_optimize(opt, x, &val);
    nlopt_destroy(opt);
    
    ss->res = res;
    ss->obj_val = val;


    ss->UVW.v1 = x[0];
    ss->UVW.v2 = x[1];
    ss->UVW.v3 = x[2];

    ss->PQR.v1 = x[3];
    ss->PQR.v2 = x[4];
    ss->PQR.v3 = x[5];    
    
    ss->roll = x[6];
    ss->pitch = x[7];

    ss->aero_con.v1 = x[8];
    ss->aero_con.v2 = x[9];
    ss->aero_con.v3 = x[10];    
    
    ss->thrust = x[11];

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
    double dtmin = 1e-8;
    double dtmax = 1e-2;
    double tol = 5e-12;

    struct Integrator * ode = integrator_create_controlled(12, 4, rigid_body_lin_forces, ac, controller, ss);
    integrator_set_type(ode,"rkf45");
    /* integrator_set_dt(ode, 1e-4); */
    integrator_set_adaptive_opts(ode, dtmin, dtmax, tol);
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

int main(int argc, char* argv[]){

    (void) argc;
    (void) argv;


    int simulate = 1;
    

    struct Aircraft aircraft;
    pioneer_uav(&aircraft);
    struct TrimSpec trim_spec;
    trim_spec.z_dot = 0.0;
    trim_spec.yaw_dot = 3 * 2.0 * M_PI / 500;
    trim_spec.target_vel = 120; // ft/s
    trim_spec.ac = &aircraft;

    struct SteadyState ss;
    trimmer(&trim_spec, &ss);
    steady_state_print(stdout, &ss);

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
