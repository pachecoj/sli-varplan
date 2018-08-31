#ifndef NONLINEAR_H__
#define NONLINEAR_H__

#include "mex.h"
#include "math.h"

typedef struct{
    mxArray * Ais, * Iis, * AiIis;
    double * kappa, * n;
    double *Vd, *Vs, *AA;
} MEM;

#define ACT(y) (pow(y,nu)/(pow(kappa,nu)+pow(y,nu)))
#define REP(y) (1-ACT(y))


void nonlinear(int N, double * y, double * x, MEM * mem){
    for (int i = 0; i < N; i++){
        int nAi = mxGetN(mxGetCell(mem->Ais,i));
        int * Ai = (int *)mxGetPr(mxGetCell(mem->Ais,i));
        int nIi = mxGetN(mxGetCell(mem->Iis,i));
        int * Ii = (int *)mxGetPr(mxGetCell(mem->Iis,i));
        int nAiIi = mxGetN(mxGetCell(mem->AiIis,i));
        int * AiIi = (int *)mxGetPr(mxGetCell(mem->AiIis,i));
        
        // decay part
        y[i] = -mem->Vd[i] * x[i] / (mem->kappa[i + i*N] + x[i]);
        
        double h = 1.;
        for (int j = 0; j < nIi; j++){
            double h2 = pow(x[Ii[j]]/mem->kappa[i + Ii[j]*N],mem->n[i + Ii[j]*N]);
            h *= 1. / (1. + h2);
        }
        for (int j = 0; j < nAi; j++){
            double h2 = pow(x[Ai[j]]/mem->kappa[i + Ai[j]*N],mem->n[i + Ai[j]*N]);
            h *= (1. + mem->AA[i + Ai[j]*N] * h2) / (1. + h2);
        }
        y[i] += mem->Vs[i] * h;
    }
}


bool prepareMem(MEM * mem, const mxArray * p, const mxArray * G){
    mem->kappa = mxGetPr(mxGetField(G,0,"kappa"));
    mem->n = mxGetPr(mxGetField(G,0,"n"));
    mem->AA = mxGetPr(mxGetField(G,0,"AA"));
    mem->Vs = mxGetPr(mxGetField(G,0,"Vs"));
    mem->Vd = mxGetPr(mxGetField(G,0,"Vd"));
    
    mem->Ais = mxGetField(p,0,"Ai");
    mem->Iis = mxGetField(p,0,"Ii");
    mem->AiIis = mxGetField(p,0,"AiIi");
    return true;
}

#endif
