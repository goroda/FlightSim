all: ac_trim check_grad

ac_trim:
	gcc -O2 6dof_pioneer_uav_trim.c flight_sim.c vehicles.c -lm -lcdyn -lnlopt -o ac_trim

check_grad:
	gcc -O2 check_gradient_rigid_body.c flight_sim.c vehicles.c -lm -lcdyn -lnlopt -o check_grad

clean:
	rm -f *.out
