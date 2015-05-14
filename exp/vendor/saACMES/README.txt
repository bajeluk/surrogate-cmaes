1. make .mex files from .cpp code in MatLab/Octave:
mex RankSVMLean.cpp
mex RankSVMFunc.cpp
mex EstimError.cpp
2. launch TestExampleRosenbrock.m, an example of saACM on Rosenbrock 10-dimensional function from BBOB framework
3. to optimize your own function, please change MyFunc.m, opts.StopFitness and all feval('fgeneric',..), fgeneric(...) you can find in the code
4. if you have any questions or suggestions, please write us an email:
ilya.loshchilov@gmail.com, marc.schoenauer@inria.fr, michele.sebag@lri.fr