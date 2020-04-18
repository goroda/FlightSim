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
$ ./ac_trim vehicles/pioneer_uav.json -s 120 -c 10 -y 0 --linearize

========================================================
                        TRIM RESULT                     
========================================================
Optimizer result = 4
Objective value = 1.37438E-18



Spec                :        Targets          Achieved   
---------------------------------------------------------
Speed       (ft/s)  :       1.20000E+02      1.20000E+02
-Climb Rate (ft/s)  :       1.00000E+01      1.00000E+01
Yaw Rate    (rad/s) :       0.00000E+00      0.00000E+00



State             :       x               dx   	 	x (secondary unit)
------------------------------------------------------------------------------
U     (ft/s)      :  1.19550E+02     -1.06112E-10
V     (ft/s)      :  0.00000E+00     0.00000E+00
W     (ft/s)      :  1.03883E+01     5.68434E-13
P     (rad,deg/s) :  0.00000E+00     0.00000E+00 	 0.00000E+00
Q     (rad,deg/s) :  -1.16716E-09     -3.53684E-13 	 -6.68732E-08
R     (rad,deg/s) :  0.00000E+00     0.00000E+00 	 0.00000E+00
Roll  (rad,deg)   :  0.00000E+00     0.00000E+00 	 0.00000E+00
Pitch (rad,deg)   :  3.24720E-03     -1.16716E-09 	 1.86051E-01



Input :
---------------------------------------------------------------
Elevator (rad,deg) :   5.82054E-03 	 3.33493E-01
Aileron  (rad,deg) :   0.00000E+00 	 0.00000E+00
Rudder   (rad,deg) :   0.00000E+00 	 0.00000E+00
Thrust   (lb-slug) :   1.38872E+01



Derived Quantities :
---------------------------------------------------------------
Angle of Attack (rad,deg) :   8.66773E-02 	 4.96624E+00
Sideslip Angle  (rad,deg) :   0.00000E+00 	 0.00000E+00

===============================================================

========================================================
                  LINEARIZATION RESULT                  
========================================================


x           y           z           U           V           W           P           Q           R           Roll        Pitch       Yaw        
0.000E+00   0.000E+00   0.000E+00   1.000E+00   0.000E+00   3.247E-03   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+01   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   -1.039E+01  0.000E+00   1.196E+02   
0.000E+00   0.000E+00   0.000E+00   -3.247E-03  0.000E+00   1.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   -1.196E+02  0.000E+00   
0.000E+00   0.000E+00   0.000E+00   -3.866E-02  0.000E+00   2.602E-01   0.000E+00   -1.040E+01  0.000E+00   0.000E+00   -3.217E+01  0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   -3.117E-02  0.000E+00   1.039E+01   -0.000E+00  -1.195E+02  3.217E+01   -0.000E+00  -0.000E+00  
0.000E+00   0.000E+00   0.000E+00   -3.954E-01  -0.000E+00  -1.643E+00  -0.000E+00  1.171E+02   -0.000E+00  -0.000E+00  -1.045E-01  -0.000E+00  
0.000E+00   0.000E+00   0.000E+00   0.000E+00   -6.262E-02  0.000E+00   -7.954E+00  0.000E+00   4.967E+00   0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   2.136E-02   -0.000E+00  -2.459E-01  0.000E+00   -3.835E+00  0.000E+00   0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.009E-01   0.000E+00   -3.546E-01  0.000E+00   -1.803E+00  0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   0.000E+00   3.247E-03   -3.790E-12  0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   -0.000E+00  0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   -1.167E-09  0.000E+00   0.000E+00   



Elev.     Ail.      Rud.      Thrust   
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000     
-0.056    0.000     0.000     0.077     
-0.000    -0.000    -0.000    0.000     
-16.055   -0.000    -0.000    -0.000    
0.000     -41.314   0.809     0.000     
-24.586   0.000     0.000     0.000     
0.000     4.603     -9.861    0.000     
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000   
```

## Acknowledgements

The stability and control derivatives are taken from the database provided by the UIUC Applied Aerodynamics Group https://m-selig.ae.illinois.edu/apasim/Aircraft-uiuc.html


## More information

Author: [Alex A. Gorodetsky](https://www.alexgorodetsky.com)  
Contact: [goroda@umich.edu](mailto:goroda@umich.edu)  
Copyright (c) 2020 Alex Gorodetsky  
License: GPL
