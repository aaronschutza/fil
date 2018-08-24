F90_DIR := filsrc
OBJ_DIR := filobj
MOD_DIR := filmod
TSY_DIR := tsysrc
FRC_DIR := fricsrc
IRI2007_DIR := iri2007
IRI_FILES := $(wildcard $(IRI2007_DIR)/*.for)
OBJ_IRI_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(IRI_FILES:.for=.o)))
GCPM_DIR := gcpm_v24
GCPM_FILES := $(wildcard $(GCPM_DIR)/*.for)
OBJ_GCPM_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(GCPM_FILES:.for=.o)))
F90_FILES := $(wildcard $(F90_DIR)/*.f90)
F90_MOD_FILES := $(wildcard $(F90_DIR)/*.mod.f90)
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(F90_FILES:.f90=.o)))
OBJ_MOD_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(F90_MOD_FILES:.f90=.o)))
TSY_FILES := $(wildcard $(TSY_DIR)/*.f)
OBJ_TSY_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(TSY_FILES:.f=.o)))
FRC_FILES := $(wildcard $(FRC_DIR)/*.f90)
OBJ_FRC_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(FRC_FILES:.f90=.o)))
FRC_MOD_FILES := $(wildcard $(FRC_DIR)/*.mod.f90)
OBJ_FRC_MOD_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(FRC_MOD_FILES:.f90=.o)))
CC := gfortran
LD_FLAGS := -I$(MOD_DIR) -J$(MOD_DIR)
CC_FLAGS := -I$(MOD_DIR) -J$(MOD_DIR) -fdefault-real-8 
TS_FLAGS :=
FR_FLAGS := -I$(MOD_DIR) -J$(MOD_DIR)
OP_FLAG := -O2
GCPM_FLAGS := -fd-lines-as-comments -fno-automatic -O0
OMP_FLAG := -fopenmp
INI_FILE := default.ini
N := fil

all: $(OBJ_MOD_FILES) $(OBJ_FILES) $(OBJ_TSY_FILES) $(OBJ_FRC_MOD_FILES) $(OBJ_FRC_FILES) $(OBJ_IRI_FILES) $(OBJ_GCPM_FILES)
	@$(RM) $(F90_DIR)/wd.mod.f90
	@touch $(F90_DIR)/wd.mod.f90
	@echo "module wdir" >> $(F90_DIR)/wd.mod.f90
	@echo "	CHARACTER(*), PARAMETER :: curdir = \"${CURDIR}\"" >> $(F90_DIR)/wd.mod.f90
	@echo "	CHARACTER(*), PARAMETER :: exenam = \"${N}\"" >> $(F90_DIR)/wd.mod.f90
	@echo "	CHARACTER(*), PARAMETER :: ininam = \"${N}.ini\"" >> $(F90_DIR)/wd.mod.f90
	@echo "end module wdir" >> $(F90_DIR)/wd.mod.f90
	@$(CC) $(CC_FLAGS) $(OP_FLAG) $(OMP_FLAG) -c -o $(OBJ_DIR)/wd.mod.o $(F90_DIR)/wd.mod.f90
	@$(CC) $(LD_FLAGS) $(OMP_FLAG) -o $(N) $^
	@if [ ! -e ./$(N).ini ]; then cp ./$(INI_FILE) ./$(N).ini; fi

gcpm: $(OBJ_IRI_FILES) $(OBJ_GCPM_FILES)
	$(CC) $(GCPM_FLAGS) -c -o $(OBJ_DIR)/gcpm_test.o gcpm_test.for
	$(CC) -o gcpm_test $(OBJ_DIR)/gcpm_test.o $^
	
$(OBJ_DIR)/%.mod.o: $(F90_DIR)/%.mod.f90
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)
	$(CC) $(CC_FLAGS) $(OP_FLAG) $(OMP_FLAG) -c -o $@ $<

$(OBJ_DIR)/%.o: $(F90_DIR)/%.f90
	$(CC) $(CC_FLAGS) $(OP_FLAG) $(OMP_FLAG) -c -o $@ $<

$(OBJ_DIR)/%.o: $(TSY_DIR)/%.f
	$(CC) $(TS_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(FRC_DIR)/%.f90
	$(CC) $(FR_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(FRC_DIR)/%.mod.f90
	$(CC) $(FR_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(IRI2007_DIR)/%.for
	$(CC) $(GCPM_FLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(GCPM_DIR)/%.for
	$(CC) $(GCPM_FLAGS) -c -o $@ $<

revini:
	@cp ./$(INI_FILE) ./$(N).ini

runiniscript:
	@python inireader.py $(F90_DIR)

clean: 
	@$(RM) $(N) $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod
	
cleaner: 
	@$(RM) $(N) $(N).ini $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod Options.ini
	
cleanmods:
	@$(RM) $(OBJ_DIR)/*.mod.o $(MOD_DIR)/*.mod
