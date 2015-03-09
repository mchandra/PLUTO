#include "pluto.h"

/* ============= TODO ============================

- Extension to all geometries (pressure term)
- MULTID Shock Flattening --> HLL hybridization
- computation of cmax
- extension to divergence cleaning
- background field
- Powell ?
================================================== */

/* ********************************************************************* */
void MUSTA_Solver (const State_1D *state, int beg, int end, 
                   double *cmax, Grid *grid)
/*!
 *  Implementation of the MUSTA Riemann solver using either a
 *  Lax-Friedrichs solver of FORCE flux.
 *  Reference: Eq. (6) and (7) of
 *
 *  "MUSTA Schemes for multi-dimensional hyperbolic-systems:
 *   analysis and improvements"
 *  Titarev & Toro, Int J. Numer. Meth. Fluids (2005), 49:117
 *
 *
 *********************************************************************** */
{
  int      nv, i, k;
  double scrh, cRL;
  double dtdx, dxdt;
  static unsigned char *flag;
  static double **UL, **VL, **fL, *pL, *a2L;
  static double **UM, **VM, **fM, *pM, *a2M;
  static double **UR, **VR, **fR, *pR, *a2R;
  double smin, smax, **flux, *press, *SL, *SR;

  if (UL == NULL){
    print1 (" >> USING MUSTA: OLD SOLVER\n");
    UL = ARRAY_2D(NMAX_POINT, NVAR, double);
    UM = ARRAY_2D(NMAX_POINT, NVAR, double);
    UR = ARRAY_2D(NMAX_POINT, NVAR, double);

    VL = ARRAY_2D(NMAX_POINT, NVAR, double);
    VM = ARRAY_2D(NMAX_POINT, NVAR, double);
    VR = ARRAY_2D(NMAX_POINT, NVAR, double);

    fL = ARRAY_2D(NMAX_POINT, NVAR, double);
    fM = ARRAY_2D(NMAX_POINT, NVAR, double);
    fR = ARRAY_2D(NMAX_POINT, NVAR, double);

    pL = ARRAY_1D(NMAX_POINT, double);
    pM = ARRAY_1D(NMAX_POINT, double);
    pR = ARRAY_1D(NMAX_POINT, double);

    a2L = ARRAY_1D(NMAX_POINT, double);
    a2M = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);

    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }
  
  SL    = state->SL;
  SR    = state->SR;
  flux  = state->flux;
  press = state->press;  

/* ------------------------------------------------------
            Copy main arrays
   ------------------------------------------------------ */

  for (i = beg; i <= end; i++) {
  for (nv = 0; nv < NFLX; nv++) {
    VL[i][nv] = state->vL[i][nv];
    VR[i][nv] = state->vR[i][nv];
    UL[i][nv] = state->uL[i][nv];
    UR[i][nv] = state->uR[i][nv];
  }}


/* --------------------------------------------------------
          Start MUSTA iteration cycle.
   -------------------------------------------------------- */

  dtdx = g_dt/grid[g_dir].dx[beg];
  dxdt = 1.0/dtdx;
   
  for (k = 0; k <= 1; k++){

  /* ---- compute the maximum propagation speed ----- */

    SoundSpeed2 (VL, a2L, NULL, beg, end, FACE_CENTER, grid);
    SoundSpeed2 (VR, a2R, NULL, beg, end, FACE_CENTER, grid);

    HLL_Speed (state->vL, state->vR, a2L, a2R, NULL, SL, SR, beg, end);
    for (i = beg; i <= end; i++)  cmax[i] = MAX(fabs(SL[i]), fabs(SR[i]));
        
    Flux (UL, VL, a2L, NULL, fL, pL, beg, end);
    Flux (UR, VR, a2R, NULL, fR, pR, beg, end);

  /* -- predictor flux using Lax-Friedrichs -- */
/*
    for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NFLX; nv++) {
        flux[i][nv] =   0.5*(fR[i][nv] + fL[i][nv]) 
                      - 0.5*(UR[i][nv] - UL[i][nv])*cmax[i];
      }
      press[i] = 0.5*(pL[i] + pR[i]);
    }
*/
  /* -- predictor flux using HLL -- */

    for (i = beg; i <= end; i++) {
      smin = MIN(0.0, SL[i]);
      smax = MAX(0.0, SR[i]);
      scrh = 1.0/(smax-smin);
      for (nv = 0; nv < NFLX; nv++) {
        flux[i][nv] = smin*smax*(UR[i][nv] - UL[i][nv]) +
                             smax*fL[i][nv] - smin*fR[i][nv];
        flux[i][nv] *= scrh;
      }
      press[i] = (smax*pL[i] - smin*pR[i])*scrh;
    }
 
  /* -- predictor flux using FORCE -- */
/*
    for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NFLX; nv++) {
        UM[i][nv] = 0.5*(UL[i][nv] + UR[i][nv]) -
                    0.5*(fR[i][nv] - fL[i][nv])*dtdx; 
      }
      UM[i][MXn] -= 0.5*(pR[i] - pL[i])*dtdx;
    }
    ConsToPrim  (UM, VM, beg, end, flag);
    SoundSpeed2 (VM, a2M, NULL, beg, end, FACE_CENTER, grid);
    Flux (UM, VM, a2M, NULL, fM, pM, beg, end);
      		     
    for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = 0.25*(fL[i][nv] + 2.0*fM[i][nv] + fR[i][nv] 
                             - dxdt*(UR[i][nv] - UL[i][nv]));
      }
      state->press[i] = 0.25*(pL[i] + 2.0*pM[i] + pR[i]);
    }
*/
    if (k == 1) break;

  /* -- open Riemann fan -- */

    for (i = beg; i <= end; i++) {
{
double cfl, M, cL, cR; 

M   = 2.0*fabs(pL[i] - pR[i])/(pL[i] + pR[i]);
cfl = (0.8 + M)/(1.0 + 2.0*M);

/*
M   = fabs(pL[i] - pR[i])/(pL[i] + pR[i]);
cfl = (0.8 + M*M)/(1.0 + 2.5*M*M);
*/
cfl = 0.8;
dtdx = cfl/cmax[i];
}

      for (nv = 0; nv < NFLX; nv++) {
        UL[i][nv] -= dtdx*(flux[i][nv] - fL[i][nv]);
        UR[i][nv] -= dtdx*(fR[i][nv] - flux[i][nv]);
      }
      UL[i][MXn] -= dtdx*(press[i] - pL[i]);
      UR[i][MXn] -= dtdx*(pR[i] - press[i]);
    }
    ConsToPrim(UL, VL, beg, end, flag);
    ConsToPrim(UR, VR, beg, end, flag);
  }

}

