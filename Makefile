all: main 

main:
	gcc -O2 6dof_pioneer_uav_trim.c -lm -lcdyn -lnlopt -o ac_trim

clean:
	rm -f *.out
