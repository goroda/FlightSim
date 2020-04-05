all: main 

main:
	gcc -O2 main.c -lm -lcdyn -lnlopt -o ac_trim.sh

clean:
	rm -f *.out
