#include "nonlinear.h"

void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[]
){
    
    if (nrhs != 3 || nlhs != 1) {
        mexWarnMsgTxt("Usage: y = nonlinear_help(x,p,G)");
        return;
    }
    
    int N = mxGetM(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    double * y = mxGetPr(plhs[0]);
    double * x = mxGetPr(prhs[0]);
    
    MEM mem;
    if (!prepareMem(&mem,prhs[1],prhs[2])) return;
    
    nonlinear(N,y,x,&mem);
}
