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
There is currently only a single script called ~ac_trim~. This script allows for computing the trim conditions for steady trimmed flight. To see what the script does you can give the help flag <kbd> ./ac_trim -h </kbd>

Here is sample output for steady flight at a 10 ft/s straight climb rate.



``` shell
Alexs-MacBook-Air:flight_sim_c alex$ ./ac_trim -s 120 -c -10 -y 0 


========================================================
                        TRIM RESULT                     
========================================================
Optimizer result = 4
Objective value = 8.87054E-20



Spec      :         Targets          Achieved
-----------------------------------------------
Speed       (ft/s)  :       1.20000E+02      1.20000E+02
-Climb-rate (ft/s)  :       -1.00000E+01      -1.00000E+01
Yaw Rate    (rad/s) :       0.00000E+00      2.39908E-12



State      : 	      x 	     dx   	       x  (secondary unit)
-----------------------------------------------------------------------------------
U     (ft/s)      :  1.19576E+02     -2.78932E-11
V     (ft/s)      :  0.00000E+00     4.43692E-14
W     (ft/s)      :  1.00772E+01     8.17124E-14
P     (rad,deg/s) :  -1.47887E-12     7.11748E-14 	 -8.47330E-11
Q     (rad,deg/s) :  -2.96373E-10     8.84997E-14 	 -1.69809E-08
R     (rad,deg/s) :  2.36550E-12     3.40649E-13 	 1.35533E-10
Roll  (rad,deg)   :  9.38861E-12     -1.07889E-12 	 5.37928E-10
Pitch (rad,deg)   :  1.67506E-01     -2.96373E-10 	 9.59739E+00



Input :
---------------------------------------------
Elevator (rad,deg) :   8.95386E-03 	 5.13019E-01
Aileron  (rad,deg) :   5.64446E-13 	 3.23403E-11
Rudder   (rad,deg) :   -1.50416E-13 	 -8.61823E-12
Thrust   (lb-slug) :   8.36081E+01



Derived Quantities :
---------------------------------------------
Angle of Attack (rad,deg) :   8.40760E-02 	 4.81720E+00
Sideslip Angle  (rad,deg) :   0.00000E+00 	 0.00000E+00

========================================================

```

## Acknowledgements

The stability and control derivatives are taken from the database provided by the UIUC Applied Aerodynamics Group https://m-selig.ae.illinois.edu/apasim/Aircraft-uiuc.html


## More information

Author: [Alex A. Gorodetsky](https://www.alexgorodetsky.com)  
Contact: [goroda@umich.edu](mailto:goroda@umich.edu)  
Copyright (c) 2020 Alex Gorodetsky  
License: GPL
