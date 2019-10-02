# Makefile for Matlab Compiler
# 
# It makes binary executable file which afterwards does not need 
# any Matlab license, but it requires Matlab Compiler Runtime (MCR)
# to be installed on the destination system
#
# Final binary file $(FNAME) is copied into $(DESTDIR) directory

FNAME = metacentrum_task_matlab
DESTDIR = exp

MATLAB_COMPILER = mcc
MC_FLAGS= -R -singleCompThread -R -nojvm -R -nodisplay
MC_INCLUDE= -a exp/opt_cmaes.m -a exp/opt_s_cmaes.m -a exp/util -a exp/vendor/bbob -a src -a exp/log -a exp/vendor/bbob_extra
SRC = exp/$(FNAME).m
OTHERS = exp/*.m exp/pproc/*.m exp/log/*.m src/ src/cmaes/* src/data/* src/model/* src/sample/* src/surrogate/* src/util/* src/surrogateManager.m
OUT = $(DESTDIR)/$(FNAME)

$(OUT):	$(SRC) $(OTHERS)
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $(FNAME) $<
	mv $(FNAME) $(DESTDIR)

all:	$(OUT)

clean:
	rm mccExcludedFiles.log readme.txt requiredMCRProducts.txt run_$(FNAME).sh
