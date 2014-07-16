// ----------------------- //
// -- SolvingParabola.c -- //
// ----------------------- //

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "ia32intrin.h" // intel chrono

#include "PrHit.h"

//=========================================================================
// Solve parabola using Cramer's rule 
//========================================================================
#define N1 1000
#define N2 1000
#define N3 1000

PrHit PrHit1[N1];
PrHit PrHit2[N2];
PrHit PrHit3[N3];

float X1[N1], Z1[N1];
float X2[N2], Z2[N2];
float X3[N3], Z3[N3];

float S[N1];

// ----------------------------------------------------------------------------------------------------------------------------
void solveParabola_AOS(const PrHit* hit1, const PrHit* hit2, const PrHit* hit3, float zReference, float* a, float* b, float* c)
// ----------------------------------------------------------------------------------------------------------------------------
{
    
    const float z1 = hit1->z - zReference;
    const float z2 = hit2->z - zReference;
    const float z3 = hit3->z - zReference;
    
    const float x1 = hit1->x;
    const float x2 = hit2->x;
    const float x3 = hit3->x;
    
    const float det = (z1*z1)*z2 + z1*(z3*z3) + (z2*z2)*z3 - z2*(z3*z3) - z1*(z2*z2) - z3*(z1*z1);
    
    if( fabs(det) < 1e-8 ){
        *a = 0.0f;
        *b = 0.0f;
        *c = 0.0f;
        return;
    }
    
    const float det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
    const float det2 = (z1*z1)*x2 + x1*(z3*z3) + (z2*z2)*x3 - x2*(z3*z3) - x1*(z2*z2) - x3*(z1*z1);
    const float det3 = (z1*z1)*z2*x3 + z1*(z3*z3)*x2 + (z2*z2)*z3*x1 - z2*(z3*z3)*x1 - z1*(z2*z2)*x3 - z3*(z1*z1)*x2;
    
    *a = det1/det;
    *b = det2/det;
    *c = det3/det;
}
// ----------------------------------------------------------------------------------------------------------------------------------------
void solveParabola_SOA(float x1, float z1, float x2, float z2, float x3, float z3, float* restrict a, float* restrict b, float* restrict c)
// ----------------------------------------------------------------------------------------------------------------------------------------
{
    
    /*const float z1 = hit1->z - Reference;
    const float z2 = hit2->z - Reference;
    const float z3 = hit3->z - Reference;
    
    const float x1 = hit1->x;
    const float x2 = hit2->x;
    const float x3 = hit3->x;*/
    
    const float det = (z1*z1)*z2 + z1*(z3*z3) + (z2*z2)*z3 - z2*(z3*z3) - z1*(z2*z2) - z3*(z1*z1);
    
    if( fabs(det) < 1e-8 ){
        *a = 0.0f;
        *b = 0.0f;
        *c = 0.0f;
        return;
    }
    
    const float det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
    const float det2 = (z1*z1)*x2 + x1*(z3*z3) + (z2*z2)*x3 - x2*(z3*z3) - x1*(z2*z2) - x3*(z1*z1);
    const float det3 = (z1*z1)*z2*x3 + z1*(z3*z3)*x2 + (z2*z2)*z3*x1 - z2*(z3*z3)*x1 - z1*(z2*z2)*x3 - z3*(z1*z1)*x2;
    
    *a = det1/det;
    *b = det2/det;
    *c = det3/det;
}
// ----------------------------------------------------------------------------------
void init_AOS(PrHit* restrict PrHit1, PrHit* restrict PrHit2, PrHit* restrict PrHit3)
// ----------------------------------------------------------------------------------
{
    for(int i1=0; i1<N1; i1++) { PrHit1[i1].x = rand(); PrHit1[i1].z = rand(); }
    for(int i2=0; i2<N2; i2++) { PrHit2[i2].x = rand(); PrHit2[i2].z = rand(); }
    for(int i3=0; i3<N3; i3++) { PrHit3[i3].x = rand(); PrHit3[i3].z = rand(); }
}
// ----------------
void test_AOS(void)
// ----------------
{
    float s=0.0f; // sommation pour forcer le calcul (sinon le compilateur supprime le code)
    float a, b, c;
    float zReference = 1.1f;
    double t0, t1, dt;
    
    init_AOS(PrHit1, PrHit2, PrHit3);
    
    t0 = _rdtsc();
    
    #ifdef OPENMP
    #pragma omp parallel for reduction(+:s) 
    //#pragma omp parallel for
    #endif
    for(int i1=0; i1<N1; i1++) {
        for(int i2=0; i2<N2; i2++) {
            #pragma simd
            for(int i3=0; i3<N3; i3++) { // loop vectorized
                                
                solveParabola_AOS(&PrHit1[i1], &PrHit2[i2], &PrHit3[i3], zReference, &a, &b, &c);
                
                s += a+b+c; // reduction pour forcer le calcul
                //S[i1] += a+b+c;
            }
        }
    }
    t1 = _rdtsc();
    dt = t1-t0; printf("cycles = %.2f\n", dt/(N1*N2*N3));
    
    printf("s = %f\n", s);
}
// ----------------------------------------------------------------------------------------------------------------------------------
void init_SOA(float* restrict X1, float* restrict Z1, float* restrict X2, float* restrict Z2, float* restrict X3, float* restrict Z3)
// ----------------------------------------------------------------------------------------------------------------------------------
{
    for(int i1=0; i1<N1; i1++) { X1[i1] = rand(); Z1[i1] = rand(); }
    for(int i2=0; i2<N2; i2++) { X2[i2] = rand(); Z2[i2] = rand(); }
    for(int i3=0; i3<N3; i3++) { X3[i3] = rand(); Z3[i3] = rand(); }
}
// ----------------
void test_SOA(void)
// ----------------
{
    
    
    float s=0.0f, sa=0.0f, sb=0.0f, sc=0.0f; // sommation pour forcer le calcul (sinon le compilateur supprime le code)
    float a, b, c;
    float zReference = 1.1f;
    double t0, t1, dt;
    
    init_SOA(X1, Z1, X2, Z2, X3, Z3);
    
    t0 = _rdtsc();
    
    #ifdef OPENMP
    #pragma omp parallel for reduction(+:s) 
    //#pragma omp parallel for
    #endif
    for(int i1=0; i1<N1; i1++) {
        for(int i2=0; i2<N2; i2++) {
            //#pragma ivdep
            #pragma simd
            for(int i3=0; i3<N3; i3++) { // loop vectorized
                
                float x1 = X1[i1]; float z1 = Z1[i1] - zReference;
                float x2 = X2[i2]; float z2 = Z2[i2] - zReference;
                float x3 = X3[i3]; float z3 = Z3[i3] - zReference;
                
                solveParabola_SOA(x1, z1, x2, z2, x3, z3, &a, &b, &c);
                
                s += a+b+c;
                //S[i1] += a+b+c;
            }
        }
    }
    t1 = _rdtsc();
    dt = t1-t0; printf("cycles = %.2f\n", dt/(N1*N2*N3));
    
    printf("s = %f\n", s);
}
// ---------------------------
int main(int argc, char *argv)
// ---------------------------
{
    test_AOS();
    test_SOA();
    return 0;
}