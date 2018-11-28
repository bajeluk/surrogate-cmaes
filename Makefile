# Makefile for Matlab Compiler
#
# It makes binary executable file which afterwards does not need
# any Matlab license, but it requires Matlab Compiler Runtime (MCR)
# to be installed on the destination system
#
# Final binary file $(FNAME) is copied into $(DESTDIR) directory

FNAME = metacentrum_task_matlab
MODEL = metacentrum_testmodels
METALEARN = metacentrum_metalearn
DESTDIR = exp

MATLAB_COMPILER = mcc
MC_FLAGS= -R -nodisplay -R -singleCompThread
MC_INCLUDE= -a exp/opt_cmaes.m -a exp/opt_s_cmaes.m -a exp/util -a exp/vendor/bbob -a src -a exp/log
SRC = exp/$(FNAME).m
OTHERS = exp/*.m exp/pproc/*.m exp/log/*.m src/ src/cmaes/* src/data/* src/model/* src/sample/* src/surrogate/* src/util/* src/surrogateManager.m
OUT = $(DESTDIR)/$(FNAME)

SRC_MODEL = exp/$(MODEL).m
OUT_MODEL = $(DESTDIR)/$(MODEL)

SRC_METALEARN = exp/$(METALEARN).m
OUT_METALEARN = $(DESTDIR)/$(METALEARN)

$(OUT):	$(SRC) $(OTHERS)
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $(FNAME) $<
	mv $(FNAME) $(DESTDIR)

$(OUT_MODEL): 	$(SRC_MODEL) $(OTHERS) exp/experiments/exp_*.m
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $(MODEL) $<
	mv $(MODEL) $(DESTDIR)

$(OUT_METALEARN):  $(SRC_METALEARN) $(OTHERS) exp/experiments/exp_metaLearn*.m
	$(MATLAB_COMPILER) -m $(MC_FLAGS) $(MC_INCLUDE) -o $(METALEARN) $<
	mv $(METALEARN) $(DESTDIR)

model:	$(OUT_MODEL)

metalearn: $(OUT_METALEARN)

jq:
	cd exp/vendor/jq-1.6 &&	git submodule update --init && autoreconf -fi && ./configure --with-oniguruma=builtin && $(MAKE) -j8

all:	$(OUT)

clean:
	rm mccExcludedFiles.log readme.txt requiredMCRProducts.txt run_$(FNAME).sh
