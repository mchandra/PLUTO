#include "pluto.h"

void InitFASTran(int argc, char *argv[], const Data *d, const Grid *grid,
                 struct fastranDataStruct *fastranData)
{
  PetscInitialize(&argc, &argv, NULL, NULL);
  PetscPrintf(MPI_COMM_WORLD, "\n> Initializing FASTran...");

  SNESCreate(MPI_COMM_WORLD, &fastranData->snes);

  #if (DIMENSIONS==2)
    if (grid[0].nghost != grid[1].nghost)
    {
      PetscPrintf(MPI_COMM_WORLD,
                  "\n ERRROR: Unequal ghost zones, ng1 = %d, ng2 = %d.\n",
                  grid[0].nghost, grid[1].nghost);
      exit(1);

    }

    /* Only coded periodic boundary conditions so far */
    DMDACreate2d(MPI_COMM_WORLD,
                 DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                 DMDA_STENCIL_BOX,
                 NX1, NX2,
                 PETSC_DECIDE, PETSC_DECIDE,
                 1, grid[0].nghost, PETSC_NULL, PETSC_NULL,
                 &fastranData->dmda);
  #else
    PetscPrintf(MPI_COMM_WORLD, "\n Dimension %d not coded yet.\n", DIMENSIONS);
    exit(1);
  #endif

  SNESSetDM(fastranData->snes, fastranData->dmda);

  DMCreateGlobalVector(fastranData->dmda, &fastranData->temperatureVec);
  DMCreateGlobalVector(fastranData->dmda, &fastranData->residualVec);
  VecSet(fastranData->temperatureVec, 0.);
  VecSet(fastranData->residualVec, 0.);

  fastranData->d    = d;
  fastranData->grid = grid;
  SNESSetFunction(fastranData->snes, fastranData->residualVec,
                  ComputeResidual, fastranData);

  SNESSetFromOptions(fastranData->snes);
  PetscViewer monviewer;

  SNESMonitorSet(fastranData->snes, SNESMonitorDefault, monviewer,
                 (PetscErrorCode (*)(void**))PetscViewerDestroy);

  /* Set initial dt */
  fastranData->dt = g_dt;
  PetscPrintf(MPI_COMM_WORLD, "done\n\n");
  
  int n;
  int dumpCounter=0;
    /* First copy the temperature from PLUTO to FASTran */
    double **T;
    DMDAVecGetArray(fastranData->dmda, fastranData->temperatureVec, &T);

    int jPluto, iPluto, jPetsc, iPetsc;
    for (jPluto=JBEG; jPluto <= JEND; jPluto++)
    {
      for (iPluto=IBEG; iPluto <= IEND; iPluto++)
      {
        jPetsc = jPluto - grid[1].nghost;
        iPetsc = iPluto - grid[0].nghost;

        T[jPetsc][iPetsc] = 
          d->Vc[PRS][0][jPluto][iPluto]/d->Vc[RHO][0][jPluto][iPluto];

        //printf("i = %d, j = %d, T = %f\n", iPetsc, jPetsc, T[jPetsc][iPetsc]);
      }
    }

    DMDAVecRestoreArray(fastranData->dmda, fastranData->temperatureVec, &T);
    char primVarsFileName[50];
    sprintf(primVarsFileName, "%s%06d.h5", "data", dumpCounter);

    PetscViewer viewer;
    PetscViewerHDF5Open(MPI_COMM_WORLD, primVarsFileName,
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) fastranData->temperatureVec, "temperature");
    VecView(fastranData->temperatureVec, viewer);
    PetscViewerDestroy(&viewer);
    dumpCounter++;
  
  for (n=0; n<100; n++)
  {
    TimeStepSourceTermsUsingFASTran(d, NULL, grid, fastranData);

    sprintf(primVarsFileName, "%s%06d.h5", "data", dumpCounter);

    PetscViewerHDF5Open(MPI_COMM_WORLD, primVarsFileName,
                        FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject) fastranData->temperatureVec, "temperature");
    VecView(fastranData->temperatureVec, viewer);
    PetscViewerDestroy(&viewer);

    PetscPrintf(MPI_COMM_WORLD, "\nDumped temperature\n", primVarsFileName);
    dumpCounter++;
  }
  
  exit(1);
}

void TimeStepSourceTermsUsingFASTran(const Data *d, 
                                     Time_Step *Dts,
                                     Grid *grid,
                                     struct fastranDataStruct *fastranData)
{
  fastranData->d    = d;
  fastranData->grid = grid;

  #if (DIMENSIONS==2)

    /* First copy the temperature from PLUTO to FASTran */
    double **T;
    DMDAVecGetArray(fastranData->dmda, fastranData->temperatureVec, &T);

    int jPluto, iPluto, jPetsc, iPetsc;
    for (jPluto=JBEG; jPluto <= JEND; jPluto++)
    {
      for (iPluto=IBEG; iPluto <= IEND; iPluto++)
      {
        jPetsc = jPluto - grid[1].nghost;
        iPetsc = iPluto - grid[0].nghost;

        T[jPetsc][iPetsc] = 
          d->Vc[PRS][0][jPluto][iPluto]/d->Vc[RHO][0][jPluto][iPluto];
      }
    }

    DMDAVecRestoreArray(fastranData->dmda, fastranData->temperatureVec, &T);

    /* Solve */
    SNESSolve(fastranData->snes, NULL, fastranData->temperatureVec);

    /* Finally copy the new temperature from FASTran to PLUTO */
    DMDAVecGetArray(fastranData->dmda, fastranData->temperatureVec, &T);

    for (jPluto=JBEG; jPluto <= JEND; jPluto++)
    {
      for (iPluto=IBEG; iPluto <= IEND; iPluto++)
      {
        jPetsc = jPluto - grid[1].nghost;
        iPetsc = iPluto - grid[0].nghost;

//        d->Vc[PRS][0][jPluto][iPluto] = 
//          T[jPetsc][iPetsc] * d->Vc[RHO][0][jPluto][iPluto];

        d->Vc[RHO][0][jPluto][iPluto] = 
          d->Vc[PRS][0][jPluto][iPluto] / T[jPetsc][iPetsc] ;
      }
    }

    DMDAVecRestoreArray(fastranData->dmda, fastranData->temperatureVec, &T);

  #endif
}

PetscErrorCode ComputeResidual(SNES snes, 
                               Vec temperatureVec, 
                               Vec residualVec,
                               void *ptr)
{
  struct fastranDataStruct *fastranData = (struct fastranDataStruct*)ptr;

  Vec temperatureVecLocal;
  DMGetLocalVector(fastranData->dmda, &temperatureVecLocal);

  /* Exchange ghost zone data and handles periodic boundary conditions
   * automatically. */
  DMGlobalToLocalBegin(fastranData->dmda,
                       temperatureVec,
                       INSERT_VALUES,
                       temperatureVecLocal);
  DMGlobalToLocalEnd(fastranData->dmda,
                     temperatureVec,
                     INSERT_VALUES,
                     temperatureVecLocal);

  double **T;
  double **residual;

  DMDAVecGetArray(fastranData->dmda, temperatureVecLocal, &T);
  DMDAVecGetArray(fastranData->dmda, residualVec, &residual);

  int jPluto, iPluto, jPetsc, iPetsc;
  for (jPluto=JBEG; jPluto <= JEND; jPluto++)
  {
    for (iPluto=IBEG; iPluto <= IEND; iPluto++)
    {
      jPetsc = jPluto - fastranData->grid[1].nghost;
      iPetsc = iPluto - fastranData->grid[0].nghost;

      double T_i_j              = T[jPetsc][iPetsc];
      double T_iPlus1_j         = T[jPetsc][iPetsc+1];
      double T_iMinus1_j        = T[jPetsc][iPetsc-1];
      double T_i_jPlus1         = T[jPetsc+1][iPetsc];
      double T_i_jMinus1        = T[jPetsc-1][iPetsc];

      double T_old_i_j             =    fastranData->d->Vc[PRS][0][jPluto][iPluto]
                                      / fastranData->d->Vc[RHO][0][jPluto][iPluto];

      double T_old_iPlus1_j        =    fastranData->d->Vc[PRS][0][jPluto][iPluto+1]
                                      / fastranData->d->Vc[RHO][0][jPluto][iPluto+1];

      double T_old_iMinus1_j       =    fastranData->d->Vc[PRS][0][jPluto][iPluto-1]
                                      / fastranData->d->Vc[RHO][0][jPluto][iPluto-1];

      double T_old_i_jPlus1        =    fastranData->d->Vc[PRS][0][jPluto+1][iPluto]
                                      / fastranData->d->Vc[RHO][0][jPluto+1][iPluto];

      double T_old_i_jMinus1       =    fastranData->d->Vc[PRS][0][jPluto-1][iPluto]
                                      / fastranData->d->Vc[RHO][0][jPluto-1][iPluto];

      double T_old_iPlus1_jPlus1   =    fastranData->d->Vc[PRS][0][jPluto+1][iPluto+1]
                                      / fastranData->d->Vc[RHO][0][jPluto+1][iPluto+1];

      double T_old_iMinus1_jPlus1  =    fastranData->d->Vc[PRS][0][jPluto+1][iPluto-1]
                                      / fastranData->d->Vc[RHO][0][jPluto+1][iPluto-1];

      double T_old_iPlus1_jMinus1  =    fastranData->d->Vc[PRS][0][jPluto-1][iPluto+1]
                                      / fastranData->d->Vc[RHO][0][jPluto-1][iPluto+1];

      double T_old_iMinus1_jMinus1 =   fastranData->d->Vc[PRS][0][jPluto-1][iPluto-1]
                                     / fastranData->d->Vc[RHO][0][jPluto-1][iPluto-1];


      /* Bx at (i+1/2, j) */
      double BxRightEdge = fastranData->d->Vs[BX1s][0][jPluto][iPluto];

      /* By at (i+1/2, j): 
       * By(i,  j)   = 0.5*(By(i,  j-1/2) + By(i,   j+1/2) )
       * By(i+1,j)   = 0.5*(By(i+1,j-1/2) + By(i+1, j+1/2) ) 
       * By(i+1/2,j) = 0.5*(By(i,  j)     + By(i+1, j    ) ) */
      double By_i_j      = 0.5*(  fastranData->d->Vs[BX2s][0][jPluto-1][iPluto]
                                + fastranData->d->Vs[BX2s][0][jPluto][iPluto]
                               );
      double By_iPlus1_j = 0.5*(  fastranData->d->Vs[BX2s][0][jPluto-1][iPluto+1]
                                + fastranData->d->Vs[BX2s][0][jPluto][iPluto+1]
                               );
      double ByRightEdge = 0.5*(By_i_j + By_iPlus1_j);

      /* Bx at (i-1/2, j) */
      double BxLeftEdge = fastranData->d->Vs[BX1s][0][jPluto][iPluto-1];

      /* By at (i-1/2, j):
       * By(i-1,j)   = 0.5*(By(i-1,j-1/2) + By(i-1, j+1/2) )
       * By(i,  j)   = 0.5*(By(i  ,j-1/2) + By(i,   j+1/2) ) 
       * By(i-1/2,j) = 0.5*(By(i-1,  j)   + By(i,   j    ) ) */
      double By_iMinus1_j = 0.5*(  fastranData->d->Vs[BX2s][0][jPluto-1][iPluto-1]
                                 + fastranData->d->Vs[BX2s][0][jPluto][iPluto-1]
                               );
      double ByLeftEdge   = 0.5*(By_i_j + By_iMinus1_j);

      /* Bx at (i, j+1/2):
       * Bx(i,j)     = 0.5*(Bx(i-1/2,j)   + Bx(i+1/2, j)   )
       * Bx(i,j+1)   = 0.5*(Bx(i-1/2,j+1) + Bx(i+1/2, j+1) ) 
       * Bx(i,j+1/2) = 0.5*(Bx(i,  j)     + Bx(i,     j+1) ) */
      double Bx_i_j      = 0.5*(  fastranData->d->Vs[BX2s][0][jPluto][iPluto-1]
                                + fastranData->d->Vs[BX2s][0][jPluto][iPluto]
                               );
      double Bx_i_jPlus1 = 0.5*(  fastranData->d->Vs[BX2s][0][jPluto+1][iPluto-1]
                                + fastranData->d->Vs[BX2s][0][jPluto+1][iPluto]
                               );
      double BxTopEdge   = 0.5*(Bx_i_j + Bx_i_jPlus1);

      /* By at (i, j+1/2) */
      double ByTopEdge = fastranData->d->Vs[BX2s][0][jPluto][iPluto];

      /* Bx at (i, j-1/2):
       * Bx(i,j-1)   = 0.5*(Bx(i-1/2,j-1) + Bx(i+1/2, j-1) ) 
       * Bx(i,j)     = 0.5*(Bx(i-1/2,j)   + Bx(i+1/2, j)   )
       * Bx(i,j-1/2) = 0.5*(Bx(i,  j)     + Bx(i,     j-1) ) */
      double Bx_i_jMinus1 = 0.5*(  fastranData->d->Vs[BX2s][0][jPluto-1][iPluto-1]
                                 + fastranData->d->Vs[BX2s][0][jPluto-1][iPluto]
                                );
      double BxBottomEdge = 0.5*(Bx_i_j + Bx_i_jMinus1);

      /* By at (i, j-1/2) */
      double ByBottomEdge = fastranData->d->Vs[BX2s][0][jPluto-1][iPluto];

      double BMagRightEdge  = sqrt(  BxRightEdge*BxRightEdge 
                                   + ByRightEdge*ByRightEdge 
                                  );

      double BMagLeftEdge   = sqrt(  BxLeftEdge*BxLeftEdge 
                                   + ByLeftEdge*ByLeftEdge 
                                  );

      double BMagTopEdge    = sqrt(  BxTopEdge*BxTopEdge 
                                   + ByTopEdge*ByTopEdge 
                                  );

      double BMagBottomEdge = sqrt(  BxBottomEdge*BxBottomEdge 
                                   + ByBottomEdge*ByBottomEdge 
                                  );

      /* Finally get the unit vectors */
      double bxRightEdge  =  BxRightEdge  / BMagRightEdge;
      double bxLeftEdge   =  BxLeftEdge   / BMagLeftEdge;
      double bxTopEdge    =  BxTopEdge    / BMagTopEdge;
      double bxBottomEdge =  BxBottomEdge / BMagBottomEdge;

      double byRightEdge  =  ByRightEdge  / BMagRightEdge;
      double byLeftEdge   =  ByLeftEdge   / BMagLeftEdge;
      double byTopEdge    =  ByTopEdge    / BMagTopEdge;
      double byBottomEdge =  ByBottomEdge / BMagBottomEdge;

      double dx1 = fastranData->grid[0].dx_glob[iPluto];
      double dx2 = fastranData->grid[1].dx_glob[jPluto];

      double x1RightEdge  = fastranData->grid[0].xr_glob[iPluto];
      double x2RightEdge  = fastranData->grid[1].x_glob[iPluto];

      double x1LeftEdge   = fastranData->grid[0].xl_glob[iPluto];
      double x2LeftEdge   = fastranData->grid[1].x_glob[iPluto];

      double x1TopEdge    = fastranData->grid[0].x_glob[iPluto];
      double x2TopEdge    = fastranData->grid[1].xr_glob[iPluto];

      double x1BottomEdge = fastranData->grid[0].x_glob[iPluto];
      double x2BottomEdge = fastranData->grid[1].xl_glob[iPluto];

      double primVarsRightEdge[NVAR];
      double primVarsLeftEdge[NVAR];
      double primVarsTopEdge[NVAR];
      double primVarsBottomEdge[NVAR];

      int var;
      for (var=0; var<NVAR; var++)
      {
        primVarsRightEdge[var] = 
          0.5*(  fastranData->d->Vc[var][0][jPluto][iPluto]
               + fastranData->d->Vc[var][0][jPluto][iPluto+1]
              );
        
        primVarsLeftEdge[var] = 
          0.5*(  fastranData->d->Vc[var][0][jPluto][iPluto-1]
               + fastranData->d->Vc[var][0][jPluto][iPluto]
              );

        primVarsTopEdge[var] = 
          0.5*(  fastranData->d->Vc[var][0][jPluto+1][iPluto]
               + fastranData->d->Vc[var][0][jPluto][iPluto]
              );

        primVarsBottomEdge[var] = 
          0.5*(  fastranData->d->Vc[var][0][jPluto][iPluto]
               + fastranData->d->Vc[var][0][jPluto-1][iPluto]
              );
      }

      double kappaParallelRightEdge, kappaParallelLeftEdge;
      double kappaParallelTopEdge  , kappaParallelBottomEdge;

      double kappaPerpRightEdge, kappaPerpLeftEdge;
      double kappaPerpTopEdge  , kappaPerpBottomEdge;

      double phiRightEdge, phiLeftEdge;
      double phiTopEdge  , phiBottomEdge;

      TC_kappa(primVarsRightEdge,
               x1RightEdge, x2RightEdge, 0, 
               &kappaParallelRightEdge, &kappaPerpRightEdge, &phiRightEdge);

      TC_kappa(primVarsLeftEdge,
               x1LeftEdge, x2LeftEdge, 0, 
               &kappaParallelLeftEdge, &kappaPerpLeftEdge, &phiLeftEdge);

      TC_kappa(primVarsTopEdge,
               x1TopEdge, x2TopEdge, 0, 
               &kappaParallelTopEdge, &kappaPerpTopEdge, &phiTopEdge);

      TC_kappa(primVarsBottomEdge,
               x1BottomEdge, x2BottomEdge, 0, 
               &kappaParallelBottomEdge, &kappaPerpBottomEdge, &phiBottomEdge);

      residual[jPetsc][iPetsc] =
        (T_i_j - T_old_i_j)/fastranData->dt
       -(
          (1./dx1)
        * (   kappaParallelRightEdge * bxRightEdge
            * (
                 (  bxRightEdge * (T_iPlus1_j - T_i_j)/dx1)
               + (  byRightEdge / dx2 
                  * limiter4( T_old_i_jPlus1 - T_old_i_j, T_old_iPlus1_jPlus1 - T_old_iPlus1_j,
                              T_old_i_j - T_old_i_jMinus1, T_old_iPlus1_j - T_old_iPlus1_jMinus1
                            )
                 )
              )
           -  kappaParallelLeftEdge * bxLeftEdge
            * (
                 (  bxLeftEdge * (T_i_j - T_iMinus1_j)/dx1)
               + (  byLeftEdge / dx2
                  * limiter4(T_old_iMinus1_jPlus1 - T_old_iMinus1_j, T_old_i_jPlus1 - T_old_i_j,
                             T_old_iMinus1_j - T_old_iMinus1_jMinus1, T_old_i_j - T_old_i_jMinus1
                            )
                 )   
              )
          )
        + 
          (1./dx2)
        * (   kappaParallelTopEdge * byTopEdge
            * (
                 (  byTopEdge * (T_i_jPlus1 - T_i_j)/dx2)
               + (  bxTopEdge / dx1
                  * limiter4(T_old_iPlus1_j - T_old_i_j, T_old_iPlus1_jPlus1 - T_old_i_jPlus1,
                             T_old_i_j - T_old_iMinus1_j, T_old_i_jPlus1 - T_old_iMinus1_jPlus1
                            )
                 )
              )
           -  kappaParallelBottomEdge * byBottomEdge
            * (
                 (  byBottomEdge * (T_i_j - T_i_jMinus1)/dx2)
               + (  bxBottomEdge / dx1
                  * limiter4(T_old_iPlus1_jMinus1 - T_old_i_jMinus1, T_old_iPlus1_j - T_old_i_j,
                             T_old_i_jMinus1 - T_old_iMinus1_jMinus1, T_old_i_j - T_old_iMinus1_j
                            )
                 )
              )
          )
        );

//      residual[jPetsc][iPetsc] = 
//        (T_i_j - T_old_i_j)/fastranData->dt -
//        100.*(T_iPlus1_j - T_iMinus1_j)/(2.*(1./NX1)) -
//        100.*(T_i_jPlus1 - T_i_jMinus1)/(2.*(1./NX2));

    }
  }


  DMDAVecRestoreArray(fastranData->dmda, temperatureVecLocal, &T);
  DMDAVecRestoreArray(fastranData->dmda, residualVec, &residual);

  DMRestoreLocalVector(fastranData->dmda, &temperatureVecLocal);

}


double minMod(double a, double b)
{
  return ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0);
}

double limiter2(double a, double b)
{
  return minMod(2*minMod(a, b), (a+b)/2.);
}

double limiter4(double a, double b, double c, double d)
{
  return limiter2(limiter2(a, b), limiter2(c, d));
}
