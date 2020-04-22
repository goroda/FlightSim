/*
  This file is part of the FlightSim project

  Author: Alex Gorodetsky
  email: goroda@umich.edu
  Copyright (c) 2020 Alex Gorodetsky 
  License: GPL3

*/

#ifndef TRIMMING_H
#define TRIMMING_H

#include "flight_sim.h"

struct TrimSpec
{
    real z_dot;
    real yaw_dot;
    real target_vel;
    
    struct Aircraft * ac;

    real thresh;
};

int trimmer(struct TrimSpec * data, struct SteadyState * ss);

#endif
