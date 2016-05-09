# basic C/C++ compiler. NO cuda yet
CC	= g++
#CC	= gcc
CFLAGS	= -lgsl -lgslcblas -lm
SDIR	= src
ODIR	= build
TARGET	= bin/runner

# need to make a choice of ONE suffix
SRCEXT	= c

SRC	= $(wildcard $(SDIR)/*.$(SRCEXT))
OBJ	= $(patsubst $(SDIR)/%,$(ODIR)/%,$(SRC:.$(SRCEXT)=.o))
INC	= -I include

$(TARGET): $(OBJ)
	@echo "Compiling : $(CC) $^ -o $(TARGET)"; $(CC) $^ -o $(TARGET)

$(ODIR)/%.o: $(SDIR)/%.$(SRCEXT)
	@mkdir -p $(ODIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

# auxiliary compiles go here:

plotter1:
	gle -o "out/unstable.pdf" -d pdf "out/plotter/unstable.gle"
plotter2:
	python "out/plotter/ene_tr.py"

clean:
	@echo "Cleaning : $(RM) -r $(ODIR) $(TARGET)"; $(RM) -r $(ODIR) $(TARGET)

.PHONY: clean


