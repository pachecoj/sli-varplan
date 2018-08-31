#include "nonlinear.h"

void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[]
){
    
    if (nrhs != 8 || nlhs != 2) {
        mexWarnMsgTxt("Usage: [x,traj] = nonlinearsimulate(x0,p,G,nstep,rand,dt,u,sigma_noise)");
        return;
    }
    
    int N = mxGetM(prhs[0]);
    double * x0 = mxGetPr(prhs[0]);
    MEM mem;
    if (!prepareMem(&mem,prhs[1],prhs[2])) return;
    int nstep = (int)mxGetScalar(prhs[3]);
    double * rand = mxGetPr(prhs[4]);
    double dt = mxGetScalar(prhs[5]);
    double sqrtdt = sqrt(dt);
    double * u = mxGetPr(prhs[6]);
    double sigma_noise = mxGetScalar(prhs[7]);
    
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    double * x = mxGetPr(plhs[0]);
    for (int i= 0; i < N; i++) x[i] = x0[i];
    plhs[1] = mxCreateDoubleMatrix(N,nstep,mxREAL);
    double * traj = mxGetPr(plhs[1]);
    
    double * buf = new double[N];
    
    for (int step = 0; step < nstep; step++){
        nonlinear(N,buf,x,&mem);
        for (int i = 0; i < N; i++){
            x[i] += (buf[i] + u[i]) * dt + sigma_noise*rand[i + step*N]*sqrtdt;            
            if (x[i] < 1e-6) x[i] = 1e-6;
            traj[i + step*N] = x[i];
        }
    }
    delete [] buf;
}
