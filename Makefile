IDIR=include
CC=gcc
CFLAGS=-O2 

ODIR=obj
LDIR=./

LIBS=-lm -lcdyn -lnlopt

_DEPS = flight_sim.h vehicles.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = flight_sim.o vehicles.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ac_trim: $(OBJ) 6dof_pioneer_uav_trim.c
	$(CC) -o ac_trim 6dof_pioneer_uav_trim.c $(OBJ) $(CFLAGS) $(LIBS)

check_grad: $(OBJ) check_gradient_rigid_body.c
	$(CC) -o check_grad check_gradient_rigid_body.c $(OBJ) $(CFLAGS) $(LIBS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *- core $(INCDIR)/*-
