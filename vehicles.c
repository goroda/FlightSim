/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

#include "vehicles.h"


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
