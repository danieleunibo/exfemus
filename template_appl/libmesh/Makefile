######### MAKEFILE FOR FEMUS APPLICATION ##################
# ----------------------------------------------------------
# --------------------- USES -------------------------------
# ----------------------------------------------------------
# make : compile the main application we are running
# make gencase : compile the gencase
# make clean : delete the objects files and the executable just in the app directory
# make src_clean : delete everything compiled in femus
# make resu_clean : delete the solutions in the app RESU dir, excluding mesh and matrixes
# make all_clean : delete everything compiled in femus and clean the app RESU directory


# ----------------------------------------------------------
# -------------------- FLAGS -------------------------------
# ----------------------------------------------------------
# define the compiler to use
CXX = mpicxx
# define the optimized and the debugging flag for compiling
ifeq ($(METHOD),opt)
  CFLAGS = -O2 -w -felide-constructors -fPIC -funroll-loops -fstrict-aliasing 
else
  CFLAGS = -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra \
  -Wno-long-long -Wpointer-arith -Wformat -Wparentheses -fPIC 

endif

#### ------------   LIBMESH -----------------------------
LIBMESH_INCLUDE=$(LIBMESH_PATH)include/
LIBMESH_CONTRIB=$(LIBMESH_PATH)contrib/
LIBMESH_LIB=$(LIBMESH_PATH)lib64/
#### ------------   PETSC -----------------------------
PETSC_PATH=$(PETSC_DIR)
PETSC_LIB=$(PETSC_PATH)/linux-opt/lib/
# ----------------------------------------------------------
# -------------------- INCLUDE -----------------------------
# ----------------------------------------------------------
# define any directories containing header files other than /usr/include
CONTRIB_HDF5_INCPATH =  $(HDF5_PATH)/include
FEMUS_CONTRIB        =  $(FEMUS_DIR)/contrib
MED_ROOT_DIR         =  $(MED_PATH)
med_ROOT_DIR         =  $(med_PATH)


INCLUDES =  -I./
INCLUDES += -I./SRC/
INCLUDES += -I$(FEMUS_DIR)USER_APPL/$(FM_MYAPP)/DATA/
INCLUDES += -I$(FEMUS_DIR)include/
INCLUDES += -I$(FEMUS_DIR)USER_APPL/$(FM_MYAPP)/SRC/
# INCLUDES += -I$(FEMUS_CONTRIB)/matrix/
# INCLUDES += -I$(FEMUS_CONTRIB)/parallel/

INCLUDES += -I$(LIBMESH_INCLUDE)/base/  
INCLUDES += -I$(LIBMESH_INCLUDE)/ 
INCLUDES += -I$(LIBMESH_INCLUDE)/enums 
INCLUDES += -I$(LIBMESH_INCLUDE)/error_estimation
INCLUDES += -I$(LIBMESH_INCLUDE)/fe
INCLUDES += -I$(LIBMESH_INCLUDE)/libmesh
INCLUDES += -I$(LIBMESH_INCLUDE)/mesh
INCLUDES += -I$(LIBMESH_INCLUDE)/numerics
INCLUDES += -I$(LIBMESH_INCLUDE)/parallel
INCLUDES += -I$(LIBMESH_INCLUDE)/partitioning
INCLUDES += -I$(LIBMESH_INCLUDE)/physics
INCLUDES += -I$(LIBMESH_INCLUDE)/quadrature
INCLUDES += -I$(LIBMESH_INCLUDE)/reduced_basis
INCLUDES += -I$(LIBMESH_INCLUDE)/solution_transfer
INCLUDES += -I$(LIBMESH_INCLUDE)/solvers
INCLUDES += -I$(LIBMESH_INCLUDE)/geom
INCLUDES += -I$(LIBMESH_INCLUDE)/systems
INCLUDES += -I$(LIBMESH_INCLUDE)/utils

INCLUDES += -I$(PETSC_DIR)/include
INCLUDES += -I$(PETSC_DIR)/linux-opt/include

INCLUDES += -I$(CONTRIB_HDF5_INCPATH)

INCLUDES += -I${MED_ROOT_DIR}/include/salome/ 
INCLUDES += -I${med_ROOT_DIR}/include/

INCLUDES += -I${FEMUS_CONTRIB}/Libmesh_cpp/
INCLUDES += -I${FEMUS_CONTRIB}/Libmesh_cpp/libmesh_problem/
# the files below are set in order to know when it is necessary to recompile.
# If the INCLUDES_SRC are changed than the SRC files will be recompiled,
# if the INCLUDES_CONTRIB are changed than the SRC_CONTRIB files will be recompiled,
INCLUDES_SRC  =  $(wildcard $(FEMUS_DIR)USER_APPL/$(FM_MYAPP)/DATA/*.h)
# INCLUDES_SRC +=  $(wildcard $(FEMUS_DIR)include/*.h)
# INCLUDES_CONTRIB  =  $(wildcard $(FEMUS_CONTRIB)/matrix/*.h)
# INCLUDES_CONTRIB +=  $(wildcard $(FEMUS_CONTRIB)/parallel/*.h)
INCLUDES_CONTRIB :=  $(wildcard $(FEMUS_CONTRIB)/Libmesh_cpp/*.h)
INCLUDES_CONTRIB +=  $(wildcard $(FEMUS_CONTRIB)/Libmesh_cpp/libmesh_problem/*.h)

# ----------------------------------------------------------
# -------------------- LIBRARIES ---------------------------
# ----------------------------------------------------------
# specify library paths needed by FEMUS (we need PETSC and 
# HDF5 libraries to compile FEMUS)
LFLAGS += -L$(HDF5_PATH)lib64/
LFLAGS += -L$(PETSC_LIB)
LFLAGS += -L$(LIBMESH_LIB)
LFLAGS += -L$(MED_ROOT_DIR)/lib/salome/
LFLAGS += -L$(med_ROOT_DIR)/lib/
# and now specify the libraries' names
LIBS =  -lm   -lmesh_opt -lhdf5 -lpetsc -lmed -lmedloader -lmedcoupling \
	-linterpkernel -lmedcouplingremapper -lmedC 



# ----------------------------------------------------------
# ---------------- SOURCES and OBJECTS ---------------------
# ----------------------------------------------------------
# we divide the source files in two groups, the sources of
# the application and of FEMUS (SRC) and the sources of the
# contrib (SRC_CONTRIB). They are recompiled in different cases.
SRC := $(wildcard $(FEMUS_DIR)USER_APPL/$(FM_MYAPP)/SRC/*.C)
# SRC += $(wildcard $(FEMUS_DIR)src/*.C)
SRC_CONTRIB := $(wildcard $(FEMUS_CONTRIB)/Libmesh_cpp/*.C)
SRC_CONTRIB += $(wildcard $(FEMUS_CONTRIB)/Libmesh_cpp/libmesh_problem/*.C)
# SRC_CONTRIB += $(wildcard $(FEMUS_CONTRIB)/parallel/*.C)

# now define the correspondig objects of the source files simply by
# replacing the suffix .c of all words with the .o suffix.
SRC_OBJS = $(patsubst %.C, %-$(METHOD).o, $(SRC))
CONTRIB_OBJS = $(patsubst %.C, %-$(METHOD).o, $(SRC_CONTRIB))
# finally define the executable file name 
MAIN = $(FM_MYAPP)-$(METHOD)



# --------------------------------------------------------
# -------------- MAKE COMMANDS ---------------------------
# --------------------------------------------------------
# set PHONY variables to distinguish between file names 
# which could be in the directory and make rules
.PHONY: depend clean gencase
# this is the executable file rule, the one executed when you call make
all:    $(MAIN)
	@echo  $(MAIN) succesfully compiled!

# this is the rule for linking the objects to obtain the executable file.
# There are two groups of objects, the SRC and the CONTRIB objects, which
# are recompiled in different cases
$(MAIN): $(SRC_OBJS) $(CONTRIB_OBJS)
	@echo 
	@echo "Linking "$@"..."
	@echo 
	@echo $^          | tr " " "\n"
	@echo $(LFLAGS)   | tr " " "\n"
	@echo $(LIBS)     | tr " " "\n"
	@echo "with compiling flags" $(CFLAGS)
	@echo "=========" $(MAIN) "==========="
	@echo 
	@$(CXX) $(CFLAGS) $(INCLUDES) -o $(MAIN) $^ $(LFLAGS) $(LIBS)

# this is the rule for compiling the sources to obtain the objects.
# The actual rule is the third line below, which is very generic,
# while the first two lines define the dependencies of the two different
# groups of objects. Note that the SRC_OBJS are recompiled if the
# INCLUDES_SRC files change, while the CONTRIB_OBJS are recompiled if the
# INCLUDES_CONTRIB files change.
$(SRC_OBJS): $(INCLUDES_SRC)
$(CONTRIB_OBJS): $(INCLUDES_CONTRIB)
%-$(METHOD).o : %.C
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Remark 1 
# In this rule the automatic variables $< and $@ are used.
# The first is the name of the prerequisite (the .C file) while the
# $@ is the target of the rule (the .o file). By using these variables
# we can compile the two lists of sources without defining anything else 
# because we use an automatic suffix replacement. This rule is very generic
# so be careful by using it.
# Remark 2 
# Writing a command in this way @command the command is not printed on the stdout.
# If you call @echo "Hello" it will be printed only 
# 
# Hello
# 
# if you omit the @ you get
# 
# echo "Hello"
# Hello
# 

# The other useful rules below are explained in the MAKEFILE uses, see above
clean:
	$(RM) ./SRC/*-$(METHOD).o *~ $(MAIN) ../../contrib/Libmesh_cpp/*.o ../../contrib/Libmesh_cpp/libmesh_problem/*.o

gencase:
	cd ../gencase; $(MAKE); cd -;

resu_clean:
	@rm -f RESU/sol.*.h5  RESU/sol.*.xmf RESU/mesh.*.*.h5  RESU/mesh.*.xmf RESU/case.*.h5  RESU/case.*.xmf  RESU/time.*.xmf

src_clean:
	@$(MAKE) clean
	@rm -f ../../src/*.o ../../contrib/matrix/*.o ../../contrib/parallel/*.o

all_clean:
	@$(MAKE) resu_clean
	@$(MAKE) src_clean
	@rm .depend


# -------------   MAKE depend and depend file ------------------------
# include the dependency list 
DEPFILE	= .depend
DEPTOKEN	= '\# MAKEDEPENDS'
DEPFLAGS	= -Y -f $(DEPFILE) -s $(DEPTOKEN) -p ./


 depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)


 $(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CFLAGS) -- $(SRCS) >&/dev/null
	@echo "Updated .depend"




# "sinclude" allows for the Makefile to work when $(DEPFILE) is missing. 
sinclude  $(DEPFILE)

# DO NOT DELETE THIS LINE -- make depend needs it