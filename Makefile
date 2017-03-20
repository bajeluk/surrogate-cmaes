# Makefile for Matlab Compiler
# 
# It makes binary executable file which afterwards does not need 
# any Matlab license, but it requires Matlab Compiler Runtime (MCR)
# to be installed on the destination system
#
# Final binary file $(FNAME) is copied into $(DESTDIR) directory

FNAME = metacentrum_task_matlab
MODEL = metacentrum_testmodels
DESTDIR = exp

MATLAB_COMPILER = mcc
MC_FLAGS= -R -singleCompThread -R -nojvm -R -nodisplay
MC_INCLUDE= -a exp/opt_cmaes.m -a exp/opt_s_cmaes.m -a exp/util -a exp/vendor/bbob -a src -a exp/log
SRC = exp/$(FNAME).m
OTHERS = exp/*.m exp/pproc/*.m exp/log/*.m src/ src/cmaes/* src/data/* src/model/* src/sample/* src/surrogate/* src/util/* src/surrogateManager.m
OUT = $(DESTDIR)/$(FNAME)

SRC_MODEL = exp/$(MODEL).m
OUT_MODEL = $(DESTDIR)/$(MODEL)

$(OUT):	$(SRC) $(OTHERS)
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $(FNAME) $<
	mv $(FNAME) $(DESTDIR)

$(OUT_MODEL): 	$(SRC_MODEL) $(OTHERS)
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $(MODEL) $<
	mv $(MODEL) $(DESTDIR)

model:	$(OUT_MODEL)

all:	$(OUT)

clean:
	rm mccExcludedFiles.log readme.txt requiredMCRProducts.txt run_$(FNAME).sh
