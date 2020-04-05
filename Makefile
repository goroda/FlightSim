all: main 

main:
	gcc -O2 main.c -lm -lcdyn -lnlopt

clean:
	rm -f *.out
