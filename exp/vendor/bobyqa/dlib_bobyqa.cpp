// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt

#include "call_matlab.h"
#include "dlib/matrix.h"
#include "dlib/optimization.h"

using namespace dlib;
using namespace std;

/*
    This mex function takes a MATLAB function handle, calls it, and
    returns the results.

    For example, you can call this function in MATLAB like so:
        A = magic(3)
        y = example_mex_callback(A, @(x)x+x)

    This will result in y containing the value 2*A.
*/

typedef dlib::matrix<double,0,1> column_vector;

double optfunc( const function_handle& f, const column_vector& x )
{
    double result;

    call_matlab(f, x, returns(result));
    return result;
}

/* Objective function */
class cObj
{
public:
    cObj ( const std::string sFunName )
    {
        mfunc = sFunName;
    }

    double operator() ( const column_vector& x ) const
    {
        double result;

        call_matlab(mfunc, x, returns(result));
        return result;
    }

private:
    std::string mfunc;
};


void mex_function (
    const std::string fname,
    const column_vector& xstart,
    const int nPoints,
    const column_vector& lb,
    const column_vector& ub,
    const double rho_beg,
    const double rho_end,
    const long maxfunevals,
    // const function_handle& f,
    double& fopt,
    column_vector& xopt) 
{

    // The f argument to this function is a function handle passed from MATLAB.  To
    // call it we use the following syntax:
    // call_matlab(f, A, returns(result));
    // This is equivalent to result = f(A). Therefore, the returns(variable) syntax 
    // is used to indicate which variables are outputs of the function.
  
    const cObj fObj = cObj(fname);

    column_vector my_start = xstart;

    fopt = dlib::find_min_bobyqa(fObj,  // objective function
                    my_start, 
                    nPoints,    // number of interpolation points
                    lb,   // lower bound constraint
                    ub,   // upper bound constraint
                    rho_beg,    // initial trust region radius
                    rho_end, // stopping trust region radius
                    maxfunevals   // max number of objective function evaluations
                    );

    // fopt = fObj(xstart);
    xopt = my_start;

    cout << "bestvalue: \n" << fopt << endl;
}

// #including this brings in all the mex boiler plate needed by MATLAB.
#include "mex_wrapper.cpp"

