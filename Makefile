IDIR=include
CC=gcc
CFLAGS=-O2 -I$(IDIR) -Itpl

ODIR=obj
LDIR=./

LIBS=-lm -lcdyn -lnlopt

_DEPS = trimming.h flight_sim.h vehicles.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = trimming.o flight_sim.o vehicles.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: ac_trim ac_sim ac_sim_lin check_grad

$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ac_trim: $(OBJ) ac_trim.c
	$(CC) -o ac_trim ac_trim.c $(OBJ) $(CFLAGS) $(LIBS)

ac_sim: $(OBJ) ac_sim.c
	$(CC) -o ac_sim ac_sim.c $(OBJ) $(CFLAGS) $(LIBS)

ac_sim_lin: $(OBJ) ac_sim_lin.c
	$(CC) -o ac_sim_lin ac_sim_lin.c $(OBJ) $(CFLAGS) $(LIBS)

check_grad: $(OBJ) check_gradient_rigid_body.c
	$(CC) -o check_grad check_gradient_rigid_body.c $(OBJ) $(CFLAGS) $(LIBS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *- core $(INCDIR)/*-
	rm -f ac_trim ac_sim ac_sim_lin check_grad
