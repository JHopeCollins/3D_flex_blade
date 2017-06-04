
#########################################


FSRC = constants.f95 \
	options.f95 \
	calls.f95 \
	test.f95

FMOD = constants.f95 \
	options.f95 \
	calls.f95

##########################################


FOBJ = $(FSRC:.f95=.o)

MODS = $(FMOD:.f95=.mod)

FCMP = gfortran
FOPT = -g -Wall -Wextra -fbounds-check

LIBS = -llapack -lblas


##########################################


a.test.out : $(FOBJ)
	$(FCMP) $(FOPT) -o $@ $^ $(LIBS)


##########################################


%.o : %.f95
	$(FCMP) $(FOPT) -c $<


##########################################
