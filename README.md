# Flight Simulator Scripts

Scripts for various aircraft dynamics analysis and simulations

## Preqrequisites
* cdyn https://github.com/goroda/cdyn
* nlopt https://nlopt.readthedocs.io/en/latest/


## Installation
The provided Makefile can be used to compile the scripts if the dependencies are installed

``
make all
``

## Background

Currently, this code implements a 6DOF rigid body model under no wind conditions. Linear aerodynamic coefficient models are used to calculate the forces and moments.

## Example Usage
There is currently only a single script called <kbd> ac_trim </kbd>. This script allows for computing the trim conditions for steady trimmed flight. To see what the script does you can give the help flag <kbd> ./ac_trim -h </kbd>

Here is sample output for steady flight at a speed of 120 ft/s with a 10 ft/s straight climb rate.



``` shell
Alexs-MacBook-Air:flight_sim_c alex$ ./ac_trim -s 120 -c 10 -y 0

========================================================
                        TRIM RESULT                     
========================================================
Optimizer result = 4
Objective value = 2.21131E-18



Spec                :        Targets          Achieved   
---------------------------------------------------------
Speed       (ft/s)  :       1.20000E+02      1.20000E+02
-Climb Rate (ft/s)  :       1.00000E+01      1.00000E+01
Yaw Rate    (rad/s) :       0.00000E+00      -2.15110E-11



State             :       x               dx   	 	x (secondary unit)
------------------------------------------------------------------------------
U     (ft/s)      :  1.19550E+02     -9.20912E-11
V     (ft/s)      :  0.00000E+00     -2.95861E-13
W     (ft/s)      :  1.03883E+01     4.40181E-12
P     (rad,deg/s) :  -7.71925E-13     2.50677E-13 	 -4.42281E-11
Q     (rad,deg/s) :  -1.48382E-09     -2.75975E-12 	 -8.50166E-08
R     (rad,deg/s) :  -2.15108E-11     2.22293E-13 	 -1.23248E-09
Roll  (rad,deg)   :  -7.96985E-11     -8.41775E-13 	 -4.56639E-09
Pitch (rad,deg)   :  3.24720E-03     -1.48382E-09 	 1.86051E-01



Input :
---------------------------------------------------------------
Elevator (rad,deg) :   5.82054E-03 	 3.33493E-01
Aileron  (rad,deg) :   -2.38835E-12 	 -1.36842E-10
Rudder   (rad,deg) :   2.82377E-12 	 1.61790E-10
Thrust   (lb-slug) :   1.38872E+01



Derived Quantities :
---------------------------------------------------------------
Angle of Attack (rad,deg) :   8.66773E-02 	 4.96624E+00
Sideslip Angle  (rad,deg) :   0.00000E+00 	 0.00000E+00

===============================================================
```

## Acknowledgements

The stability and control derivatives are taken from the database provided by the UIUC Applied Aerodynamics Group https://m-selig.ae.illinois.edu/apasim/Aircraft-uiuc.html


## More information

Author: [Alex A. Gorodetsky](https://www.alexgorodetsky.com)  
Contact: [goroda@umich.edu](mailto:goroda@umich.edu)  
Copyright (c) 2020 Alex Gorodetsky  
License: GPL
