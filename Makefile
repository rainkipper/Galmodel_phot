PROG = gm
COMPILER = ifort

INSTALL_DIR = ~/bin
OBJ_DIR = obj
SRC_DIR = src
_OBJ = program.o fitting.o fitting_omalooming.o fitting_amoeba.o fitting_multinest.o output.o likelihood.o mass_to_lum_image.o fill_comp_image.o adaptive_image_real.o comp_image.o pixel.o los_real_integration.o approx_Einasto.o read_input.o images.o all_comp.o comp.o profile_collector.o profiles_den.o filters.o psf_rakendamine.o file_operations.o amoeba.o omalooming.o newton_Rhapson.o integration.o populations.o dust.o distributions.o yldine_matemaatika.o constants.o
# _OBJ = constants.o, yldine_matemaatika.o, integration_module.o, newton_Phapson.o, omalooming.o, amoeba.o, file_operations.o, psf_rakendamine.o, filters.o, profiles_den.o, profile_collector.o, comp.o, all_comp.o, images.o, read_input.o, los_real_integration.o, pixel.o, comp_image.o, adaptive_image.o, fill_comp_image.o, likelihood.o, output.o, fitting_multinest.o, fitting_amoeba.o, fitting_omalooming.o, fitting.o , program.o
# OBJ = $(patsubst %,$(OBJ_DIR)/%,$(_OBJ))
OBJ = $(_OBJ)

# -module $(OBJ_DIR)
FLAGS = -O2 -lcfitsio -lnest3 -llapack
# FLAGS = -lcfitsio -lnest3 -llapack -g -O0 -error-limit 10 -warn all
# FLAGS = -O2 -lcfitsio -lnest3 -llapack -prof-use -profile-loops=all -profile-loops-report=2 -prof-dir Prof
# FLAGS = -O2 -lcfitsio -lnest3 -llapack -prof-gen -profile-loops-report=2 -prof-dir Prof


all: $(PROG)
	
$(PROG): $(OBJ)
	$(COMPILER) $(OBJ) -o $(INSTALL_DIR)/$(PROG) $(FLAGS)
	
$(OBJ): %.o: $(SRC_DIR)/%.f90
	$(COMPILER) -c $< $(FLAGS)
clean:
	rm -f *.o *.mod

# all: $(PROG)
# $(PROG).o: $(OBJ)
# 	$(COMPILER) $(OBJ) $(FLAGS) -o $(PROG_DIR)/$(PROG)
# $(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
# 	$(COMPILER) -c -o $@ $< $(FLAGS)
# clean:
# 	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod
#
# #


constants.o:
distributions.o: constants.o
yldine_matemaatika.o: distributions.o constants.o
integration.o: constants.o
newton_Rhapson.o: constants.o
omalooming.o: constants.o
amoeba.o: constants.o
file_operations.o: constants.o
psf_rakendamine.o: constants.o file_operations.o
dust.o: yldine_matemaatika.o
populations.o: yldine_matemaatika.o file_operations.o
filters.o: file_operations.o distributions.o
profiles_den.o: constants.o
profile_collector.o: profiles_den.o
comp.o: profile_collector.o yldine_matemaatika.o populations.o
all_comp.o: comp.o
images.o: filters.o yldine_matemaatika.o psf_rakendamine.o
read_input.o: all_comp.o filters.o images.o file_operations.o populations.o dust.o filters.o
los_real_integration.o: comp.o integration.o
approx_Einasto.o: constants.o
pixel.o: constants.o
comp_image.o: pixel.o yldine_matemaatika.o images.o
adaptive_image_real.o: yldine_matemaatika.o
fill_comp_image.o: comp_image.o all_comp.o los_real_integration.o adaptive_image_real.o approx_Einasto.o
mass_to_lum_image.o: comp_image.o comp.o populations.o dust.o filters.o
likelihood.o: fill_comp_image.o images.o newton_Rhapson.o mass_to_lum_image.o 
output.o: likelihood.o
fitting_multinest.o: likelihood.o output.o
fitting_amoeba.o: amoeba.o likelihood.o
fitting_omalooming.o: omalooming.o likelihood.o
fitting.o: fitting_multinest.o fitting_amoeba.o fitting_omalooming.o
program.o: fitting.o read_input.o
