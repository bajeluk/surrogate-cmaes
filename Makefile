# Makefile for Matlab Compiler

FNAME = metacentrum_task_matlab

MATLAB_COMPILER = mcc
MC_FLAGS= -R -singleCompThread -R -nojvm -R -nodisplay
MC_INCLUDE= -a exp/opt_cmaes.m -a exp/opt_s_cmaes.m -a exp/util -a exp/vendor/bbob -a src
OUT = $(FNAME)
SRC = exp/$(FNAME).m
OTHERS = exp/bbob_test_01.m exp/metacentrum_task_matlab.m src/surrogateManager.m

$(OUT):	$(SRC) $(OTHERS)
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $@ $<

all:	$(OUT)

clean:
	rm mccExcludedFiles.log readme.txt requiredMCRProducts.txt run_$(FNAME).sh
