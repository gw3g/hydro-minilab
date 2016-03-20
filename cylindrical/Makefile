CC= gcc 	
DEPS= adv.h
CFLAGS= -O2 -g -pg -I. -lgsl -lgslcblas -lm
ODIR=obj
_OBJ= shasta.o rad_BI.o
OBJ= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

advection: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

