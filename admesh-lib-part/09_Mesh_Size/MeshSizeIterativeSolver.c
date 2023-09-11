/*-----------------------------------
 * MeshSizeIterativeSolver
 * MEX file
 *
 * h = MeshSizeIter(hn,h0,g,delta);
 *
 * input: 
 * output:
 *
 * -----------------------------------*/

// Include Libraries
#include "mex.h"
#include <math.h>
#include <matrix.h>
#include <limits.h>

// Define variables and functions
#define eps 2.2204e-16
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


// Mex Function
void mexFunction(int nOUT, mxArray *pOUT[], int nIN, const mxArray *pIN[])
{
        
    // Define variable types
    mwSize LY,LX,N;
    int i, j, k, ki;
    double *h, *D, *h0, delta, g, hmax, hmin;
    double tol, R, xfordiff, xbackdiff, yfordiff, ybackdiff, deltat, Delta;
    double inf,hn;
            
    // Read in inputs
    LY = mxGetM(pIN[0]);
    LX = mxGetN(pIN[0]);
    N  = LY*LX;
    
    h0      = mxGetPr(pIN[0]);
    D       = mxGetPr(pIN[1]);
    hmax    = mxGetScalar(pIN[2]);
    hmin    = mxGetScalar(pIN[3]);
    g       = mxGetScalar(pIN[4]);
    delta   = mxGetScalar(pIN[5]);
    
    // Get infinity variable
    inf = mxGetInf();

    // Prepare output
    pOUT[0] = mxCreateDoubleMatrix(LY, LX, mxREAL);
    h       = mxGetPr(pOUT[0]);
    
    // Define variables
    deltat  = delta/2;
    tol     = (10e-6);
    //tol     = (10e-6);
    R       = 0;
    j = 1;
    
    
    // While residual is less than specified tolerance
     while(1)
     {
         
         R = 0;
         
         for (i = 1; i<LX-1; i++)
         {
             for (j = 1; j<LY-1; j++)
             {
                  
                 k = LY * i + j; // 0-based
                 
                 //mexPrintf("\n k=%i",k);
                 
                 if( D[k] > hmax)
                 {
                     continue;
                 }
                                 
                 // Compute upwind differences
                 ki = LY * (i+1) + j; // 0-based
                 xfordiff    = min((h0[ki] - h0[k])/delta,0);
                 xfordiff    = xfordiff*xfordiff;
                 
                 ki = LY * (i-1) + j; // 0-based
                 xbackdiff   = max((h0[k] - h0[ki])/delta,0);
                 xbackdiff   = xbackdiff*xbackdiff;
                 
                 ki = LY * (i) + j+1; // 0-based
                 yfordiff    = min((h0[ki] - h0[k])/delta,0);
                 yfordiff    = yfordiff*yfordiff;
                 
                 ki = LY * (i) + j-1; // 0-based
                 ybackdiff   = max((h0[k] - h0[ki])/delta,0);
                 ybackdiff   = ybackdiff*ybackdiff;
                 
                 // Compute Delta
                 Delta = sqrt(xfordiff + xbackdiff + yfordiff + ybackdiff);
                 
                 // Compute next time step
                 hn = h0[k] + deltat*(min(Delta,g) - Delta);
                 
                 // Compute max Residual
                 R = fabs( (hn-h0[k]) ) + R;
                 
                 // Next time-step
                 h0[k] = hn;
                 
             } // j loop
         } // i loop
         
         // Check tolerance
         if ( R <= tol)
         {
             break;
         }
         
     } // while loop
    
    for (i = 0; i<LX; i++)
    {
        for (j = 0; j<LY; j++)
        {
            k = LY * i + j; // 0-based

            h[k] = h0[k];
        }
    }
        
}
