FC = gfortran
FFLAGS = -Wall -Wextra -O2 -pedantic
OBJECTS = grid_module.o main.o
EXECUTABLE = critical_cell

# Use /usr/bin/time for more detailed output, or just 'time' for the shell default
TIME = /usr/bin/time -p

# Visualization settings
PYTHON = python3
VIS_SCRIPT = visualise.py
VIS_OUTPUT = animation.mp4
VIS_FPS = 15
VIS_GRID_SIZE ?= 20 20

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	@echo "Linking executable..."
	@$(TIME) $(FC) $(FFLAGS) -o $@ $^

grid_module.o: grid_module.f90
	@echo "Compiling grid_module..."
	@$(TIME) $(FC) $(FFLAGS) -c grid_module.f90

main.o: main.f90 grid_module.o
	@echo "Compiling main..."
	@$(TIME) $(FC) $(FFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(EXECUTABLE)

run: $(EXECUTABLE)
	@echo "Running execution timer..."
	@time ./$(EXECUTABLE)

visualise:
	@echo "Generating visualization..."
	@$(TIME) $(PYTHON) $(VIS_SCRIPT) --output $(VIS_OUTPUT) --fps $(VIS_FPS) --grid-size $(VIS_GRID_SIZE)

rebuild: clean all

.PHONY: all clean run rebuild visualise
