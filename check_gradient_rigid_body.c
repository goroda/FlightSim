/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

//------------------------------------------
//------------------------------------------
//--------- Gradient Checkers --------------
//------------------------------------------
//------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>

#include "flight_sim.h"
#include "vehicles.h"

int check(real tol, real a, real b)
{
    if (pow(a - b, 2) / pow(a, 2) > tol){
        return 1;
    }

    return 0;

}

int check_grad_ac_to_e(real h, real tol)
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
    printf("Numerical %3.5E %3.5E %3.5E\n", e_ref.v1, e_ref.v2, e_ref.v3);

    // compute analytic gradient
    orient_ac_to_e_g(&ea, &ac, &e, sg);
    printf("Analytic %3.5E %3.5E %3.5E\n", e.v1, e.v2, e.v3);
    printf("\n\n\n");
    
    
    real grad_U[3];
    ac.v1 += h;
    orient_ac_to_e(&ea, &ac, &e);
    grad_U[0] = (e.v1 - e_ref.v1) / h;
    grad_U[1] = (e.v2 - e_ref.v2) / h;
    grad_U[2] = (e.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to U\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_U[0], grad_U[1], grad_U[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].U_g, sg[1].U_g, sg[2].U_g);
    printf("\n\n");
    
    real grad_V[3];
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
    
    real grad_W[3];
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

    real grad_R[3];
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
    
    real grad_P[3];
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
    
    real grad_Y[3];
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

int check_grad_rkin(real h, real tol)
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
    
    real grad_P[3];
    pqr.v1 += h;
    rkin(&ea, &pqr, &rates);
    grad_P[0] = (rates.v1 - e_ref.v1) / h;
    grad_P[1] = (rates.v2 - e_ref.v2) / h;
    grad_P[2] = (rates.v3 - e_ref.v3) / h;

    printf("Checking Gradient with respect to P\n");
    printf("Numerical: %3.5E, %3.5E, %3.5E\n", grad_P[0], grad_P[1], grad_P[2]);
    printf("Analytic: %3.5E, %3.5E, %3.5E\n", sg[0].P_g, sg[1].P_g, sg[2].P_g);
    printf("\n\n");
    
    real grad_Q[3];
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
    
    real grad_R[3];
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

    real grad_Roll[3];
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
    
    real grad_Pitch[3];
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
    
    real grad_Yaw[3];
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



int check_grad_aero_forces(real h, real tol)
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
    
    
    real grad_U[6];
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

        
    real grad_V[6];
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

    real grad_W[6];
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


    real grad_P[6];
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
    
    real grad_Q[6];
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

    real grad_R[6];
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

    real grad_el[6];
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

    real grad_ail[6];
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

    real grad_rud[6];
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

int check_grad_tdyn(real h, real tol)
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
    
    real grad_U[3];
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

    real grad_V[3];
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

    real grad_W[3];
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
    
    real grad_P[3];
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

    real grad_Q[3];
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

    real grad_R[3];
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

    real grad_Roll[3];
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

    real grad_Pitch[3];
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

    real grad_Yaw[3];
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

    real grad_el[6];
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

    real grad_ail[6];
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

    real grad_rud[6];
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

    real grad_thrust[6];
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

int check_grad_rdyn(real h, real tol)
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
    
    real grad_U[3];
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

    real grad_V[3];
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

    real grad_W[3];
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
    
    real grad_P[3];
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

    real grad_Q[3];
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

    real grad_R[3];
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
    
    real grad_el[6];
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

    real grad_ail[6];
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

    real grad_rud[6];
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

int check_grad_rigid_body_dyn(real h, real tol)
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

int main(int argc, char* argv[])
{
    real h = 1e-8;
    real tol = 1e-5;
    int res;

    fprintf(stdout, "Checking Gradient Computations \n");
    

    fprintf(stdout, "Aircraft Frame -> Euler Frame...\n");
    res = check_grad_ac_to_e(h, tol);
    /* if (res == 0){ */
    /*     fprintf(stdout, "\t.......... OK\n"); */
    /* } */
    /* else{ */
    /*     fprintf(stdout, "\t.......... failed\n");    */
    /* } */
    res = check_grad_rkin(h, tol);
    res = check_grad_aero_forces(h, tol);
    res = check_grad_rdyn(h, tol);
    res = check_grad_tdyn(h, tol);
    res = check_grad_rigid_body_dyn(h, tol);
    /* res = check_grad_trimmer(); */
    return 0;
}

/* int check_grad_trimmer(real h, real tol) */
/* { */
/*     real input[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, */
/*                       7.0, 8.0, 9.0, 10.0, 11.0, 12.0}; */

/*     // elevator aileron rudder thrust */
/*     real out_ref; */
/*     real out; */
/*     real grad[12]; */

/*     struct Aircraft aircraft; */
/*     pioneer_uav(&aircraft); */

/*     struct TrimSpec trim_spec; */
/*     trim_spec.z_dot = 5.0; */
/*     trim_spec.yaw_dot = 3.0;  */
/*     trim_spec.target_vel = 120.0; */
/*     trim_spec.ac = &aircraft; */
/*     trim_spec.thresh = 1e-14;      */

/*     out_ref = trim_objective(12, input, NULL, &trim_spec); */
/*     out = trim_objective_g(12, input, grad, &trim_spec); */
/*     printf("Checking value computation\n"); */
/*     printf("%10.5E %10.5E \n", out_ref, out); */

/*     printf("Checking gradient computation\n"); */
/*     for (size_t ii = 0; ii < 12; ii++){ */
/*         input[ii] += h; */
/*         out = trim_objective(12, input, NULL, &trim_spec); */
/*         real val = (out - out_ref) / h;         */
/*         printf("\t %10.5E %10.5E\n", val, grad[ii]); */
/*         input[ii] -= h; */
/*     } */


/*     return 0; */
/* } */
