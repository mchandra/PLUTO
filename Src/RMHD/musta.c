#include "pluto.h"

/* ********************************************************************* */
void MUSTA_Solver (const State_1D *state, int beg, int end, double *cmax,
                 Grid *grid)
/*
 *
 *
 **************************************************************************** */
{
  int      nv, i, k;
  double scrh1, scrh2, cRL;
  double dtdx, dxdt;
  static unsigned char *flag;
  static double **UL, **VL, **fL, *pL, *a2L, *hL;
  static double **UM, **VM, **fM, *pM, *a2M, *hM;
  static double **UR, **VR, **fR, *pR, *a2R, *hR;
  static double *cmin_L, *cmax_L;
  static double *cmin_R, *cmax_R;

  if (UL == NULL){
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

    hL = ARRAY_1D(NMAX_POINT, double);
    hM = ARRAY_1D(NMAX_POINT, double);
    hR = ARRAY_1D(NMAX_POINT, double);

    a2L = ARRAY_1D(NMAX_POINT, double);
    a2M = ARRAY_1D(NMAX_POINT, double);
    a2R = ARRAY_1D(NMAX_POINT, double);

    cmin_L = ARRAY_1D(NMAX_POINT, double);
    cmax_L = ARRAY_1D(NMAX_POINT, double);
    cmin_R = ARRAY_1D(NMAX_POINT, double);
    cmax_R = ARRAY_1D(NMAX_POINT, double);

    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }

/* -- copy main arrays -- */

  for (i = beg; i <= end; i++) {
  for (nv = 0; nv < NFLX; nv++) {
    VL[i][nv] = state->vL[i][nv];
    VR[i][nv] = state->vR[i][nv];
    UL[i][nv] = state->uL[i][nv];
    UR[i][nv] = state->uR[i][nv];
  }}
        
  dtdx = g_dt/grid[g_dir].dx[beg];
  dxdt = 1.0/dtdx;

  for (k = 1; k <= 2; k++){
    SoundSpeed2 (VL, a2L, hL, beg, end, FACE_CENTER, grid);
    SoundSpeed2 (VR, a2R, hR, beg, end, FACE_CENTER, grid);

    Flux (UL, VL, hL, fL, pL, beg, end);
    Flux (UR, VR, hR, fR, pR, beg, end);

    MaxSignalSpeed (VL, a2L, hL, cmin_L, cmax_L, beg, end);
    MaxSignalSpeed (VR, a2R, hR, cmin_R, cmax_R, beg, end);

    for (i = beg; i <= end; i++) {
      scrh1 = MAX(fabs(cmax_L[i]), fabs(cmax_R[i]));
      scrh2 = MAX(fabs(cmin_L[i]), fabs(cmin_R[i]));
      cRL   = MAX(scrh1, scrh2);
      state->SL[i] = -cRL, state->SR[i] = cRL;
      cmax[i] = cRL;
    }

    for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NFLX; nv++) {
        UM[i][nv] = 0.5*(UL[i][nv] + UR[i][nv]) -
                    0.5*(fR[i][nv] - fL[i][nv])/cmax[i];
      }
      UM[i][MXn] -= 0.5*(pR[i] - pL[i])/cmax[i];
    }
    ConsToPrim  (UM, VM, beg, end, flag);
    SoundSpeed2 (VM, a2M, hM, beg, end, FACE_CENTER, grid);
    Flux (UM, VM, hM, fM, pM, beg, end);
      		     
    for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NFLX; nv++) {
        state->flux[i][nv] = 0.25*(fL[i][nv] + 2.0*fM[i][nv] + fR[i][nv] 
                             - cmax[i]*(UR[i][nv] - UL[i][nv]));
      }
      state->press[i] = 0.25*(pL[i] + 2.0*pM[i] + pR[i]);
    }

    for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NFLX; nv++) {
        UL[i][nv] -= dtdx*(state->flux[i][nv] - fL[i][nv]);
        UR[i][nv] -= dtdx*(fR[i][nv] - state->flux[i][nv]);
      }
      UL[i][MXn] -= dtdx*(state->press[i] - pL[i]);
      UR[i][MXn] -= dtdx*(pR[i] - state->press[i]);
    }
    ConsToPrim(UL, VL, beg, end, flag);
    ConsToPrim(UR, VR, beg, end, flag);

  }
}
