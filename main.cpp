/* ----------------------------------------------------------------------
 *| This is a Parallel Nonlineat Finite Element Method Pacakage          |
 *| The Codes are written by WangJianHua                                 |
 *|                                                                      |
 *| Reference Book:                                                      |
 *| (1) Non-Linear Finite Element Analysis of Solids and Structures,     |
 *|     R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel   |
 *|     John Wiley and Sons, 2012.                                       |
 *| (2) Computational  Methods for Plasticity: Theory and Applications,  |
 *|     EA de Souza Neto, D Peric and DRJ Owen. Wiley, Chichester, 2008. |
 *| All rights reserved!                                                 |
  ----------------------------------------------------------------------- */

#include <mpi.h>
#include <petsc.h>
#include <iostream>
#include <io/InputReader.h>
#include <solvers/Solver.h>
#include <io/OutputManager.h>

int main(int argc, char **args)
{
  PetscErrorCode ierr;
  int rank;
  int size;

  ierr = PetscInitialize(&argc, &args, NULL, NULL); CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  GlobalData *globalDat = NONLINEARFEMIO::InputReader(rank, args);

  std::shared_ptr<Solver> solver = std::make_shared<Solver>();

  std::shared_ptr<OutputManager> output = std::make_shared<OutputManager>();

  while(GlobalData::GetInstance()->m_active)
  {
    if(nullptr != solver) solver->Run();
    if(nullptr != output) output->Run();
  }

  std::cout << "PyFem analysis terminated successfully." << std::endl;

  GlobalData::DestoryInstance();
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return ierr;
}
