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

  DMCreateGlobalVector(fastranData->dmda, &fastranData->temperatureVec);
  DMCreateGlobalVector(fastranData->dmda, &fastranData->residualVec);
  VecSet(fastranData->temperatureVec, 0.);
  VecSet(fastranData->residualVec, 0.);

  fastranData->d    = d;
  fastranData->grid = grid;
  SNESSetFunction(fastranData->snes, fastranData->residualVec,
                  ComputeResidual, fastranData);

  SNESSetFromOptions(fastranData->snes);

  PetscPrintf(MPI_COMM_WORLD, "done\n\n");
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
    for (jPluto=JBEG; jPluto < JEND; jPluto++)
    {
      for (iPluto=IBEG; iPluto < IEND; iPluto++)
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

    for (jPluto=JBEG; jPluto < JEND; jPluto++)
    {
      for (iPluto=IBEG; iPluto < IEND; iPluto++)
      {
        jPetsc = jPluto - grid[1].nghost;
        iPetsc = iPluto - grid[0].nghost;

        d->Vc[PRS][0][jPluto][iPluto] = 
          T[jPetsc][iPetsc] * d->Vc[RHO][0][jPluto][iPluto];

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
  for (jPluto=JBEG; jPluto < JEND; jPluto++)
  {
    for (iPluto=IBEG; iPluto < IEND; iPluto++)
    {
      jPetsc = jPluto - fastranData->grid[1].nghost;
      iPetsc = iPluto - fastranData->grid[0].nghost;

      double T_i_j              = T[jPetsc][iPetsc];
      double T_iPlus1_j         = T[jPetsc][iPetsc+1];
      double T_iMinus1_j        = T[jPetsc][iPetsc-1];
      double T_i_jPlus1         = T[jPetsc+1][iPetsc];
      double T_i_jMinus1        = T[jPetsc-1][iPetsc];
      double T_iPlus1_jPlus1    = T[jPetsc+1][iPetsc+1];
      double T_iMinus1_jPlus1   = T[jPetsc+1][iPetsc-1];
      double T_iPlus1_jMinus1   = T[jPetsc-1][iPetsc+1];
      double T_iMinus1_jMinus1  = T[jPetsc-1][iPetsc-1];

      double T_old_i_j             =    fastranData->d->Vc[PRS][0][jPluto][iPluto]
                                      / fastranData->d->Vc[PRS][0][jPluto][iPluto];

      double T_old_iPlus1_j        =    fastranData->d->Vc[PRS][0][jPluto][iPluto+1]
                                      / fastranData->d->Vc[PRS][0][jPluto][iPluto+1];

      double T_old_iMinus1_j       =    fastranData->d->Vc[PRS][0][jPluto][iPluto-1]
                                      / fastranData->d->Vc[PRS][0][jPluto][iPluto-1];

      double T_old_i_jPlus1        =    fastranData->d->Vc[PRS][0][jPluto+1][iPluto]
                                      / fastranData->d->Vc[PRS][0][jPluto+1][iPluto];

      double T_old_i_jMinus1       =    fastranData->d->Vc[PRS][0][jPluto-1][iPluto]
                                      / fastranData->d->Vc[PRS][0][jPluto-1][iPluto];

      double T_old_iPlus1_jPlus1   =    fastranData->d->Vc[PRS][0][jPluto+1][iPluto+1]
                                      / fastranData->d->Vc[PRS][0][jPluto+1][iPluto+1];

      double T_old_iMinus1_jPlus1  =    fastranData->d->Vc[PRS][0][jPluto+1][iPluto-1]
                                      / fastranData->d->Vc[PRS][0][jPluto+1][iPluto-1];

      double T_old_iPlus1_jMinus1  =    fastranData->d->Vc[PRS][0][jPluto-1][iPluto+1]
                                      / fastranData->d->Vc[PRS][0][jPluto-1][iPluto+1];

      double T_old_iMinus1_jMinus1 =   fastranData->d->Vc[PRS][0][jPluto-1][iPluto-1]
                                     / fastranData->d->Vc[PRS][0][jPluto-1][iPluto-1];


      /* Bx at (i+1/2, j) */
      double BxRightEdge = fastranData->d->Vs[BX1s][0][jPluto][iPluto];

      /* By at (i+1/2, j) */
      double ByRightEdge = 0.5*(  fastranData->d->Vc[BX2][0][jPluto][iPluto]
                                + fastranData->d->Vc[BX2][0][jPluto][iPluto+1]
                               );

      /* Bx at (i-1/2, j) */
      double BxLeftEdge = fastranData->d->Vs[BX1s][0][jPluto][iPluto-1];

      /* By at (i-1/2, j) */
      double ByLeftEdge = 0.5*(  fastranData->d->Vc[BX2][0][jPluto][iPluto-1]
                               + fastranData->d->Vc[BX2][0][jPluto][iPluto]
                              );

      /* Bx at (i, j+1/2) */
      double BxTopEdge = 0.5*(  fastranData->d->Vc[BX1][0][jPluto][iPluto]
                              + fastranData->d->Vc[BX1][0][jPluto+1][iPluto]
                             );

      /* By at (i, j+1/2) */
      double ByTopEdge = fastranData->d->Vs[BX2s][0][jPluto][iPluto];

      /* Bx at (i, j-1/2) */
      double BxBottomEdge = 0.5*(  fastranData->d->Vc[BX1][0][jPluto-1][iPluto]
                                 + fastranData->d->Vc[BX1][0][jPluto][iPluto]
                                );

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

//      residual[jPetsc][iPetsc] =
//        (T_i_j - T_old_i_j)/fastranData->dt
//       -(
//          (1./DX1)
//        * (   D(X1RightEdge, X2Center)*bxRightEdge
//            * (
//                 (  bxRightEdge * (T_iPlus1_j - T_i_j)/DX1)
//               + (  byRightEdge / DX2 
//                  * limiter4( T_old_i_jPlus1 - T_old_i_j, T_old_iPlus1_jPlus1 - T_old_iPlus1_j,
//                              T_old_i_j - T_old_i_jMinus1, T_old_iPlus1_j - T_old_iPlus1_jMinus1
//                            )
//                 )
//              )
//           -  D(X1LeftEdge, X2Center)*bxLeftEdge
//            * (
//                 (  bxLeftEdge * (T_i_j - T_iMinus1_j)/DX1)
//               + (  byLeftEdge / DX2
//                  * limiter4(T_old_iMinus1_jPlus1 - T_old_iMinus1_j, T_old_i_jPlus1 - T_old_i_j,
//                             T_old_iMinus1_j - T_old_iMinus1_jMinus1, T_old_i_j - T_old_i_jMinus1
//                            )
//                 )   
//              )
//          )
//        + 
//          (1./DX2)
//        * (   D(X1Center, X2TopEdge)*byTopEdge
//            * (
//                 (  byTopEdge * (T_i_jPlus1 - T_i_j)/DX2)
//               + (  bxTopEdge / DX1
//                  * limiter4(T_old_iPlus1_j - T_old_i_j, T_old_iPlus1_jPlus1 - T_old_i_jPlus1,
//                             T_old_i_j - T_old_iMinus1_j, T_old_i_jPlus1 - T_old_iMinus1_jPlus1
//                            )
//                 )
//              )
//          -   D(X1Center, X2BottomEdge)*byBottomEdge
//            * (
//                 (  byBottomEdge * (T_i_j - T_i_jMinus1)/DX2)
//               + (  bxBottomEdge / DX1
//                  * limiter4(T_old_iPlus1_jMinus1 - T_old_i_jMinus1, T_old_iPlus1_j - T_old_i_j,
//                             T_old_i_jMinus1 - T_old_iMinus1_jMinus1, T_old_i_j - T_old_iMinus1_j
//                            )
//                 )
//              )
//          )
//        );
      
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
