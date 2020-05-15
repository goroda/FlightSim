/* struct Trajectory * flight_sim_lin(struct Vec3 * xyz, real yaw, struct SteadyState * ss, */
/*                                    struct Aircraft * ac, double dt_save, size_t nsteps, real * AB); */


/*     if (linearize){ */
    /*         fprintf(stdout, "========================================================\n"); */
    /*         fprintf(stdout, "     Simulating Linearized Dynamics at Steady State     \n"); */
    /*         fprintf(stdout, "========================================================\n"); */


    /*         sprintf(filename, "lABmat_%s", sim_name); */
    /*         printf("Saving AB mat to %s\n", filename); */
    /*         fp = fopen(filename, "w"); */

    /*         fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n", */
    /*                 "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw", "Elevator", "Aileron", "Rudder", "Thrust"); */
    /*         for (size_t ii = 0; ii < 12; ii++){ */
    /*             for (size_t jj = 0; jj < 16; jj++){ */
    /*                 fprintf(fp, "%3.15f ", jac[jj*12 + ii]); */
    /*             } */
    /*             fprintf(fp, "\n"); */
    /*         } */
    /*         fclose(fp); */
            
            
    /*         sprintf(filename, "lrb_%s", sim_name); */
    /*         printf("Saving simulation to %s\n", filename); */
    /*         FILE * fp = fopen(filename, "w"); */
            
    /*         if (fp == NULL){ */
    /*             fprintf(stdout, "Cannot open file %s\n", filename); */
    /*             return 1; */
    /*         } */

    /*         struct Vec3 xyz_ss = {0, 0, -5}; */
    /*         real yaw_ss = M_PI/4.0; */
    /*         struct Trajectory * traj_ss = flight_sim_ss(&xyz_ss, yaw_ss, &ss, &aircraft, dt_save, nsteps); */
        
            
    /*         struct Vec3 xyz_perturbed = {0, 0, 0}; */
    /*         real yaw_perturbed = 0.0; */
    /*         struct Trajectory * traj_lin = flight_sim_lin(&xyz_perturbed, yaw_perturbed, &ss, &aircraft, */
    /*                               dt_save, nsteps, jac); */
    /*         fprintf(fp, "%-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s\n", */
    /*                 "t", "x", "y", "z", "U", "V", "W", "P", "Q", "R", "Roll", "Pitch", "Yaw"); */
            

    /*         trajectory_print(traj_lin, fp, 10); */
    /*         trajectory_print(traj_ss, fp, 10); */
    /*         trajectory_free(traj_lin); traj_lin = NULL; */
    /*         trajectory_free(traj_ss); traj_ss = NULL; */
            
    /*         fclose(fp); */
    /*     } */
        
    /* } */

/* struct Trajectory * */
/* flight_sim_lin(struct Vec3 * xyz, real yaw, struct SteadyState * ss, struct Aircraft * ac, */
/*                double dt_save, size_t nsteps, real * AB) */
/* { */
/*     /\* double dtmin = 1e-16; *\/ */
/*     /\* double dtmax = dt_save; *\/ */
/*     /\* double tol = 1e-14; *\/ */

/*     struct Integrator * ode = integrator_create_controlled(12, 4, rigid_body_linearized, AB, no_controller, NULL); */
/*     integrator_set_type(ode, "rk4"); */
/*     integrator_set_dt(ode, 1e-4); */
/*     /\* integrator_set_type(ode,"rkf45");     *\/ */
/*     /\* integrator_set_adaptive_opts(ode, dtmin, dtmax, tol); *\/ */
/*     integrator_set_verbose(ode, 0); */
    

/*     double start_time = 0.0; */
/*     double ic[12];     */
/*     double control[4]; */
/*     steady_state_set_vec(ss, xyz->v1, xyz->v2, xyz->v3, yaw, ic, control);     */
/*     no_controller(0.0, ic, control, ss); */

/*     struct Trajectory * traj = NULL;         */
/*     int res = trajectory_add(&traj, 12, 4, start_time, ic, control); */

/*     double dt = dt_save; */
/*     for (size_t ii = 0; ii < nsteps; ii++){ */
/*         res = trajectory_step(traj, ode, dt); */
/*     } */
/*     integrator_destroy(ode); */
/*     return traj; */
/* } */
