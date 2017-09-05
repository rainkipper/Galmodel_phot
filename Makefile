PROG = gm
COMPILER = ifort

PROG_DIR = ~/bin
OBJ_DIR = obj
SRC_DIR = src

_OBJ = 	constants.o \
		yldine_matemaatika.o \
		integration_module.o \
		konvolutsioon.o \
		file_operations.o \
		filters.o \
		profiles_den.o \
		profile_collector.o \
		comp.o \
		all_comp.o \
		images.o \
		read_input.o \
		los_real_integration.o \
		pixel.o \
		comp_image.o \
		adaptive_image_real.o \
		fill_comp_image.o \
		likelihood.o \
		only_output_image.o \
		fitting_multinest.o \
		program.o

FLAGS = -O2 -module $(OBJ_DIR) -lcfitsio -lnest3 -llapack
# FLAGS = -O2 -module $(OBJ_DIR) -lcfitsio -lnest3 -llapack -prof-use -profile-loops=all -profile-loops-report=2 -prof-dir Prof
# FLAGS = -O2 -module $(OBJ_DIR) -lcfitsio -lnest3 -llapack -prof-gen -profile-loops-report=2 -prof-dir Prof
OBJ = $(patsubst %,$(OBJ_DIR)/%,$(_OBJ))
$(PROG).o: $(OBJ)
	$(COMPILER) $(OBJ) $(FLAGS) -o $(PROG_DIR)/$(PROG)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(COMPILER) -c -o $@ $< $(FLAGS)
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod