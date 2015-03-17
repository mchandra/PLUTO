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

  printf("IBEG = %ld, IEND = %ld\n", IBEG-grid[0].nghost, IEND-grid[0].nghost);

  PetscPrintf(MPI_COMM_WORLD, "done\n\n");
}

void TimeStepUsingFASTran(const Data *d, 
                          Time_Step *Dts,
                          Grid *grid)
                          //struct fastranDataStruct *fastranData)
{
//  fastranData->d    = d;
//  fastranData->grid = grid;

  #if (DIMENSIONS==2)
  


  #endif

}

PetscErrorCode ComputeResidual(SNES snes, 
                               Vec temperatureVec, 
                               Vec residualVec,
                               void *ptr)
{

}
