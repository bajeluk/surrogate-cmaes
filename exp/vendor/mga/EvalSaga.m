function [cost_value] = EvalSaga(t)
global MGADSMproblem
[J,DVvec] = mga_dsm(t,MGADSMproblem);
% cost_value = 1000000.0 - J;
cost_value = -J;
