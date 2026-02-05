FC = gfortran
FFLAGS = -Wall -Wextra -O2 -pedantic

OBJECTS = grid_module.o main.o

EXECUTABLE = critical_cell

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

grid_module.o: grid_module.f90
	$(FC) $(FFLAGS) -c grid_module.f90

main.o: main.f90 grid_module.o
	$(FC) $(FFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(EXECUTABLE)

run: $(EXECUTABLE)
	./$(EXECUTABLE)

rebuild: clean all

.PHONY: all clean run rebuild
