Surrogate CMA-ES
================

Surrogate CMA-ES (S-CMA-ES) is a surrogate-based optimizing evolution strategy. It is based on the N. Hansen's CMA-ES algorithm which is interconnected with Gaussian processes or random forests.

## Log and News ##

2017/03/22: Metacentrum interface migrated to the newer [PBS Pro](https://wiki.metacentrum.cz/wiki/PBS_Professional), default memory limit raised to 1.5GB

## Introduction ##

The optimizer can be called via the similar function interface as the original Matlab CMA-ES code:

```matlab
function [xmin, fmin, counteval, stopflag, out, bestever, y_eval] = ...
    s_cmaes( fitfun, xstart, insigma, inopts, varargin )
```

**Arguments**:

- `fitfun` -- name of the objective function with the interface
```matlab
y = fitness(x, varargin)
```
- `xstart` -- objective variables initial point, determines dimension
- `insigma` -- initial coordinate wise standard deviation(s)
- `inopts` -- CMA-ES options struct
- `varargin` -- variable arguments
    - if the 5th argument is the string `'SurrogateOptions'`, then the **surrogate modelling** is switched **on**. In this case, the next (6th) argument should be a structure array with the options and settings for the second argument of `surrogateManager()`. It will be an instance of parameters defined in **surrogateParams**, see section [Experiment definition](#experiment-definition)
    - the rest of the variable arguments is then always passed as `varargin` into the objective function (see `fitfun` argument above)

**Returns**:

- `xmin` -- search point with the minimum fitness from the last iteration
- `fmin` -- the function value of **xmin**: _fitfun(xmin) = fmin_
- `counteval` -- the number of function evaluations used
- `out` -- various output statistics, histories and solutions
- `bestever` -- struct containing overall best solution (for convenience)
- `y_eval` -- numeric matrix with the best function values (1st column) for respective function evaluations (2nd column) and RMSE and Kendall's tau coefficients (3rd and 4th column)


## COCO/BBOB experiments ##

This repository provides an interface for running experiments within COCO/BBOB framework. The experiment is identified with its *experiment ID* (`exp_id`).

### Experiment definition ###

Each experiment is specified in the m-file with the filename `[exp_id].m` placed in the directory `exp/experiments/`. Example of such file can be found in the file `exp/experiments/exp_restrEC_quicktest.m` (new version), or in the file `exp/experiments/exp_geneEC_08_10D.m` (old version).

The following settings *have to be* defined -- recommended values for BBOB/COCO experiments are in curly brackets, types of the values in normal brackets:

- **exp_id** -- experiment ID -- unique name of each experiment, without spaces (string)
- **exp_description** -- short description of the experiment (string)
- **bbobParams**
    - **dimensions** -- list of dimensions of input space which should be tested { `2,3,5,10,20` } (integer, integer...)
    - **functions** -- list of BBOB function numbers to be tested { `cell2num(1:24)` } (integer, integer...)
    - **opt_function** -- handle to a BBOB wrapper for the optimizer { `@opt_s_cmaes` } (function handle)
    - **instances** -- vector of instances to be tested; for BBOB final results, use the default value { `[1:5 41:50]` } (integer vector)
    - **maxfunevals** -- maximal number of (original) function evaluations -- string to be `eval`-ed in `surrogateManager()`; particularly, the `dim` parameter can be used { `'250 * dim'` } (string)
- **surrogateParams** -- structure array defining behaviour and settings of the surrogate modelling, i.e. for the function `surrogateManager()` and functions called from within there. The following settings is recommended to be set:
    - **evoControl** -- type of evolution control to be used { `'none'` } ( `'doubletrained'` | `'individual'` | `'generation'` | `'none'` )
    - **modelType** -- the type of a surrogate model to be used, Gaussian processes and random forests are implemented so far. The special model `'bbob'` stands for a virtual model which in fact returns exact values from the respective BBOB functions. The default value `''` causes an error, so it really should be set to any of the three valid models. { `''` } ( `'gp'` | `'rf'` | `'bbob'` )
    - **modelOpts** -- structure array defining parameters for the respective surrogate model, see comments in the respective models { `[]` }

### The experiment initialization ###

The experiment specified in the file `[exp_id].m` has to be **initialized first**. After running `startup` from the surrogate-cmaes root folder (it is called automatically if Matlab started from this directory), call in Matlab (in the new version)
```matlab
expInit('[exp_id]')
```
or (in the old version)
```matlab
[exp_id]
```
This generates the definition of the experiment into the folder `exp/experiments/[exp_id]/`, especially into the file `exp/experiments/[exp_id]/scmaes_params.mat` from where this experiment settings are loaded during experiment submission and during running of each experiment's task.

### Jobs, tasks and parameter ordering ###

The initialization of an experiment defines a whole _job_ of _tasks_, where each task can be run as a separate Matlab process on a separate machine. The whole experiment performs a _grid search_ through all the possible combinations of parameter-values, each task tests one such combination.

Each of the *N* parameters *p_i* takes a value from its corresponding value-set **V**(*i*). Thus, the total number of tasks *N_t* is given as a product of the numbers of values of respective parameters. The parameters and its corresponding value-sets are specified in the experiment definition file `[exp_id].m`.

**Example.** Consider the following settings of an experiment. The first section _bbobParams_ defines parameters of the BBOB/COCO framework, the _surrogateParams_ section defines parameters for the surrogate modelling and _modelParams_ specify parameters of the models. Each parameter *p_i* has its value-set **V**(*i*) defined in the curly brackets followed by its name, where only the parameters _dimensions_, _functions_ and _covFcn_ have more than one value in its sets.

```matlab
exp_id = 'exp_example';
exp_description = 'Example of the experiment definition file';

% BBOB/COCO framework settings
bbobParams = { ...
  'dimensions',         { 2, 5 }, ...          % 2 different values
  'functions',          { 1, 2, 3 }, ...       % 3 different values
  'opt_function',       { @opt_s_cmaes }, ...
  'instances',          { [1:10] }, ...
  'maxfunevals',        { '75 * dim' }
};

% Surrogate manager parameters
surrogateParams = { ...
  'evoControl',         { 'doubletrained' }, ...
  'modelType',          { 'gp' }, ...
  'evoControlTrainRange', { 10 }, ...
  'evoControlTrainNArchivePoints', { '15*dim' },...
  'evoControlRestrictedParam', { 0.2 }
};

% Model parameters
modelParams = { ...
  'predictionType',     { 'sd2' }, ...
  'covFcn',             { '@covSEiso', '{@covMaterniso, 5}' }, ...  % 2 different values
  'hyp',                { struct('lik', log(0.01), 'cov', log([0.5; 2])) }, ...
  'normalizeY',         { true }
};
```

The number of tasks is therefore equal to *N_t* = 12:

```matlab
N_t   =   2 * 3 * 1 * 1 * 1   *   1 * 1 * 1 * 1 * 1   *   1 * 2 * 1 * 1   =   12
```
_(end of example)_

Each of the *N_t* different parameter-sets *(p_1, p_2,..., p_N)* is numbered with its unique *task_id*, a number between 1 and *N_t*. The *task_id* = 1 always corresponds to the parameter-set having all the first values from its respective value-sets. Then, other values are selected for the next *task_id*'s such that the later the parameter is defined in `[exp_id].m`, the more frequently is its value changed. Thus, for example, the same settings, but with the last parameter having its second value belongs to *task_id* = 2. The last *task_id* = *N_t* is defined as all the parameters having its last values.

The conversion between the *task_id* number and the indices *(k_1, k_2,..., k_N)* of values of the respective parameters *(p_1, p_2,..., p_N)* in fact corresponds to conversion of a multiple-base number *(k_1, k_2,..., k_N)* to a decimal number where the base *b(i)* of the *i*-th order is equal to products of cardinalities of the value-sets **V**(*j*) for each *j = (i+1),...,N*:

```matlab
% for i = 1,...,(N-1):
b(i) = length(V(i+1)) * length(V(i+2)) * ... * length(V(N))
% for i = N
b(N) = 1
```

where `length(V(j))` stands for the number of different values in *j*-th value-set **V**(*j*). The *task_id* is therefore equal to

1 + 𝛴 (*k_i* - 1) *b*(*i*)

(summation is for *i* = 1,...,*N*), or in Matlab pseudo-code

```matlab
task_id = (k_1 - 1)*b(1) + (k_2 - 1)*b(2) + ... + (k_N - 1)*1 + 1
```

Vice versa, the value-indices *k_1, k_2,...,k_N* can be obtained from *task_id* as

```matlab
% for i = 1:
k_1 = floor( (task_id - 1) / b(1) ) + 1
% for i = 2,...,N:
k_i = floor( mod(task_id - 1, b(i-1)) / b(i) ) + 1
```

**Example.** In the example above, the following settings belongs to *task_id* = 11.

```matlab
bbobParams = struct(...
  'dimensions',         5, ...
  'functions',          3, ...
  'opt_function',       @opt_s_cmaes, ...
  'instances',          [1:10], ...
  'maxfunevals',        '75 * dim');
surrogateParams = struct(...
  'evoControl',         'doubletrained', ...
  'modelType',          'gp', ...
  'evoControlTrainRange', 10, ...
  'evoControlTrainNArchivePoints', '15*dim',...
  'evoControlRestrictedParam', 0.2);
modelParams = struct(...
  'predictionType',     'sd2', ...
  'covFcn',             '@covSEiso', ...
  'hyp',                struct('lik', log(0.01), 'cov', log([0.5; 2])), ...
  'normalizeY',         true);
```

As there are only three multiple-valued parameters _dimensions_, _functions_ and _covFcn_, the mapping to integers 1,...,12 can be simplified to 3-element-basis *b* = *(b(1), b(2), b(3))* = (6, 2, 1), and the coefficients are chosen *k_1* from {1, 2}, *k_2* from {1, 2, 3} and *k_3* from {1, 2}. Each *task_id* can be therefore expressed as

```matlab
task_id = (k_1 - 1)*6 + (k_2 - 1)*2 + (k_3 - 1)*1 + 1
```

which equals to *task_id* = 11 for *(k_1, k_2, k_3)* = (2, 3, 1).

_(end of example)_

### Running an individual task ###

Having the experiment initialized, the individual task (with the chosen `task_id`) can be called with the function `metacentrum_task_matlab()`. It has the following syntax

```matlab
metacentrum_task_matlab(exp_id, exppath_short, id, varargin)
```
where 

- **exp_id** is the _experiment ID_ (`exp_id`) (string)
- **exppath_short** is the absolute path to the directory `exp/experiments` (string)
- **id** is the *task_id* of the task to be run (integer)
- *varargin* is an optional argument specifying additional settings (structure array)

The function should be called from Matlab from the Surrogate CMA-ES root folder, after running the `startup.m` script. Note that only the integer **id** specifies the concrete parameter-values being used for this particular run.

### Submitting tasks on Czech Metacentrum ###

After the initialization of the experiment (and hence creation of the directory `exp/experiments/[exp_id]` with the file `scmaes_params.mat`), all or subset of the *task_id*'s can be submitted for computation on the Czech national computational cluster *Metacentrum* via calling (from shell)

```shell
exp/metacentrum_master_template.sh EXPID META_QUEUE [ID1] [ID2]...
```
where

- **EXPID** is the _experiment ID_ (`exp_id`) (string)
- **META_QUEUE** is the name of the Metacentrum queue to which the tasks should be submitted, i.e. one of the following { `2h` | `4h` | `1d` | `2d` | ... } (string)
- *ID1*, *ID2*,... are the optional integer arguments defining the concrete *task_id*'s to be submitted; if no such *ID*'s are given, all the possible *task_id*'s of the experiment will be submitted (integer, integer...)

Part of this tasks submission is using Matlab Compile Runtime (MCR) to compile the Matlab sources of all the project into a single binary `exp/metacentrum_task_matlab` which will be eventually called within each individual job. This considerably saves the available Matlab licenses since no Matlab license is needed for running this MCR compiled binary.

See the `Makefile` in the project's root folder for the details regarding the MCR compilation, and adjust the `MC_INCLUDE` variable therein if additional files or folders needs to be used for computation. The MCR compilation is not performed if the binary `exp/metacentrum_binary_task` is up to date, which means newer than most of the m-files in directories `src/` and `exp/` -- see the dependencies specified in the variable `OTHERS` in the `Makefile`.

In addition to the `metacentrum_binary_task`, several other core source files are packed into a `tar` archive `deploy/[exp_id]_src.tar` which is at the beginning of each task copied to each computational node -- adjust the variable `FILES_TO_DEPLOY` in the file `exp/bash_settings.sh` if more files are needed. If the archive `deploy/[exp_id]_src.tar` already exists before calling `exp/metacentrum_master_template.sh`, it is used for task submission as it is (i.e. without packing new sources of without MCR re-compilation).

During the second part of the `metacentrum_master_template.sh` run, the script `exp/experiment/[exp_id]/binary_task.sh` (which is in fact a copy of `exp/metacentrum_binary_task_template.m`) will be submitted as individual tasks for the Metacentrum's PBS/Torque scheduler. See the file for the details.

**Warning**. Especially **setting the environment variable `$MCR_CACHE_ROOT`** is essential to be set into the particular computing node's local drive (`$SCRATCHDIR` in Metacentrum) -- it is done in the `exp/bash_settings.sh` script. Otherwise, **extraordinary overuse** of the **shared NFS filesystem** will practically disallow all the particular Metacentrum storage from being used!!!

### Experiment post-processing ###

The results of the computed experiment can be reported using Matlab function `generateReport`

```matlab
generateReport(expFolder, settings)
```
where 

- **expFolder** is the name of the folder (string) containing the experiment results,
- **settings** are pairs of property (string) and value, or struct with properties as fields: 
    - **Description** is description of the report (string),
    - **LegendOption** is an option how legends should be organized (settings similar to function `relativeFValues`, recommended settings for `generateReport` are `out` (legend is in one separate figure) and `manyout` (legend is in multiple separated figures).
    - **OmitEmptyPlots** is an option if plots with no data available will be omitted in report (boolean).
    - **Publish** is the format of the resulting file similar to function `publish`, i.e. `html`, `pdf`; set to `off` to disable publishing (string).

The function generates Matlab script `[expFolder]/pproc/[exp_id]_report.m` which will be published using Matlab function `publish` according to **Publish** option in **settings**. Published file can be found in `[expFolder]/pproc/html` folder.

To report multiple experiment results put all experiment folder names into cell-array **expFolder**. The script `exp_[n]report_[hash].m` will be generated in all `pproc` folders of experiment folders, where `n` denotes number of experiment folders and `hash` is a sequence of characters generated according to experiment folder names.

Generated report contains plots of dependences of minimal function values on function evaluations divided by dimension for individual functions. First, only comparison of algorithm runs are plotted in each tested function and dimension. Second, former algorithms are compared to important algorithms in continuous black-box optimization field.

**Recommendation**. Use `clear all` command before running `generateReport` to ensure that no static variables from the previous run will be used.

**Warning**. Report generating requires a few minutes of rendering Matlab figures. This can cause some troubles with computer usability because of pop-up figure windows in older versions of Matlab. To avoid these problems run `generateReport` in remote desktop or take a break.
