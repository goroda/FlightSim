# Flight Simulator Scripts

Scripts for various aircraft dynamics analysis and simulations

## Preqrequisites
* cdyn https://github.com/goroda/cdyn
* nlopt https://nlopt.readthedocs.io/en/latest/


## Installation
The provided Makefile can be used to compile the scripts if the dependencies are installed

``
make
``

This command will compile three scripts
* <kbd> ac_trim </kbd> for trimming an aircraft in a given flight condition
* <kbd> ac_sim </kbd> for simulating the aircraft
* <kbd> ac_sim_lin </kbd> for simulating a linearized version of the aircraft

## Background

Currently, this code implements a 6DOF rigid body model under no wind conditions. Linear aerodynamic coefficient models are used to calculate the forces and moments.

## Example Usage

There are several analysis scripts available. Here we describe their usage

### Aircraft trimming
Trimming is performed via the script called <kbd> ac_trim </kbd>. This script allows for computing the trim conditions for steady trimmed flight. To see what the script does you can give the help flag <kbd> ./ac_trim -h </kbd> The 

Here is sample output for steady flight at a speed of 120 ft/s with a 10 ft/s straight climb rate. The script also performs linearization around the trim condition because the <kbd> --linearize </kbd> flag is specified. Finally it saves both the trim condition and the linearization (A and B) matrices to the the output file <kbd> trim_conditions.json </kbd> This file can then be read by the other scripts


``` shell
$ ./ac_trim vehicles/pioneer_uav.json -s 120 -c 10 -y 0 --linearize -o trim_conditions.json

========================================================
                        TRIM RESULT                     
========================================================
Optimizer result = 4
Objective value = 3.01356E-19



Spec                :        Targets          Achieved   
---------------------------------------------------------
Speed       (ft/s)  :       1.20000E+02      1.20000E+02
Climb Rate  (ft/s)  :       1.00000E+01      1.00000E+01
Yaw Rate    (rad/s) :       0.00000E+00      0.00000E+00



State             :       x               dx   	 	x (secondary unit)
------------------------------------------------------------------------------
U     (ft/s)      :  1.19576E+02     -3.71481E-11
V     (ft/s)      :  0.00000E+00     0.00000E+00
W     (ft/s)      :  1.00772E+01     1.33937E-12
P     (rad,deg/s) :  0.00000E+00     0.00000E+00 	 0.00000E+00
Q     (rad,deg/s) :  -5.47605E-10     -9.17324E-13 	 -3.13755E-08
R     (rad,deg/s) :  0.00000E+00     0.00000E+00 	 0.00000E+00
Roll  (rad,deg)   :  0.00000E+00     0.00000E+00 	 0.00000E+00
Pitch (rad,deg)   :  1.67506E-01     -5.47605E-10 	 9.59739E+00



Input :
---------------------------------------------------------------
Elevator (rad,deg) :   8.95386E-03 	 5.13019E-01
Aileron  (rad,deg) :   0.00000E+00 	 0.00000E+00
Rudder   (rad,deg) :   0.00000E+00 	 0.00000E+00
Thrust   (lb-slug) :   8.36081E+01



Derived Quantities :
---------------------------------------------------------------
Angle of Attack   (rad,deg) :   8.40760E-02 	 4.81720E+00
Sideslip Angle    (rad,deg) :   0.00000E+00 	 0.00000E+00
Flight Path Angle (rad,deg) :   8.34301E-02 	 4.78019E+00
Bank Angle        (rad,deg) :   1.49012E-08 	 8.53774E-07

===============================================================

========================================================
                  LINEARIZATION RESULT                  
========================================================


x           y           z           U           V           W           P           Q           R           Roll        Pitch       Yaw        
0.000E+00   0.000E+00   0.000E+00   9.860E-01   0.000E+00   1.667E-01   0.000E+00   0.000E+00   0.000E+00   0.000E+00   -1.000E+01  0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   -1.008E+01  0.000E+00   1.196E+02   
0.000E+00   0.000E+00   0.000E+00   -1.667E-01  0.000E+00   9.860E-01   0.000E+00   0.000E+00   0.000E+00   0.000E+00   -1.196E+02  0.000E+00   
0.000E+00   0.000E+00   0.000E+00   -3.883E-02  0.000E+00   2.543E-01   0.000E+00   -1.009E+01  0.000E+00   0.000E+00   -3.172E+01  0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   -3.084E-02  0.000E+00   1.008E+01   -0.000E+00  -1.196E+02  3.172E+01   -0.000E+00  -0.000E+00  
0.000E+00   0.000E+00   0.000E+00   -3.922E-01  -0.000E+00  -1.642E+00  -0.000E+00  1.172E+02   -0.000E+00  -0.000E+00  -5.364E+00  -0.000E+00  
0.000E+00   0.000E+00   0.000E+00   0.000E+00   -6.262E-02  0.000E+00   -7.954E+00  0.000E+00   4.967E+00   0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   2.073E-02   -0.000E+00  -2.459E-01  0.000E+00   -3.835E+00  0.000E+00   0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.009E-01   0.000E+00   -3.546E-01  0.000E+00   -1.803E+00  0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   0.000E+00   1.691E-01   -9.259E-11  0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.000E+00   -0.000E+00  0.000E+00   0.000E+00   0.000E+00   
0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   1.014E+00   -5.554E-10  0.000E+00   0.000E+00   



Elev.     Ail.      Rud.      Thrust   
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000     
-0.077    0.000     0.000     0.077     
-0.000    -0.000    -0.000    0.000     
-16.053   -0.000    -0.000    -0.000    
0.000     -41.314   0.809     0.000     
-24.586   0.000     0.000     0.000     
0.000     4.603     -9.861    0.000     
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000     
0.000     0.000     0.000     0.000  
```


### Aircraft simulation
Simulation is performed via the script called <kbd> ac_sim </kbd>. This script simulates an aircraft for a given initial condition. To see what the script does you can give the help flag <kbd> ./ac_sim -h </kbd> This script requires the specification an initial condition in addition to a vehicle file. A reference initial condition is provided [examples/example_ic.json](examples/example_ic.json). This initial condition is provided below.


``` json
{
    "x":0,
    "y":0,
    "z":0,
    "U":119.556604117257407,
    "V":0.000000000000000,
    "W":10.306231704642098,
    "P":0.000000000000000,
    "Q":-0.000000001679369,
    "R":0.000000000000000,
    "roll":0.000000000000000,
    "pitch":0.085991201830814,
    "yaw":0.0,
    "elevator":0.006646961692774,
    "aileron":0.000000000000000,
    "rudder":0.000000000000000,
    "thrust":48.877175456058318,    
}
```

An example execution that simulates an aircraft for 100 seconds and saves every 0.1 seconds is shown below. The resulting file has the same name as the initial condition, with an addition ".run" extension.


```shell
$ ./ac_sim vehicles/pioneer_uav.json examples/example_ic.json -t 100 -d 0.1 

Aircraft filename = vehicles/pioneer_uav.json
Initial condition filename = examples/example_ic.json
========================================================
             Simulating Nonlinear Dynamics              
========================================================
Time = 1.0000E+02
save_time = 1.0000E-01
nsteps = 1000
Saving simulation to examples/example_ic.json.run
```


## Acknowledgements

The stability and control derivatives are taken from the database provided by the UIUC Applied Aerodynamics Group https://m-selig.ae.illinois.edu/apasim/Aircraft-uiuc.html


## More information

Author: [Alex A. Gorodetsky](https://www.alexgorodetsky.com)  
Contact: [goroda@umich.edu](mailto:goroda@umich.edu)  
Copyright (c) 2020 Alex Gorodetsky  
License: GPL
