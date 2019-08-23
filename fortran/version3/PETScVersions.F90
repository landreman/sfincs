! Here are some definitions that are required because the syntax for several PETSc objects
! has changed from version to version.
  
#include <petscversion.h>

! For PETSc versions prior to 3.3, the MatCreateAIJ subroutine was called MatCreateMPIAIJ.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 3))
#define MatCreateAIJ MatCreateMPIAIJ
#endif
! Hereafter in this code, use MatCreateAIJ.

! For PETSc versions prior to 3.4, the PetscTime subroutine was called PetscGetTime.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 4))
#define PetscTime PetscGetTime
#endif
!Hereafter in this code, use PetscTime.

! For PETSc versions prior to 3.4, the constant SNESNEWTONLS was called SNESLS.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 4))
#define SNESNEWTONLS SNESLS
#endif
!Hereafter in this code, use SNESNEWTONLS

! For PETSc versions prior to 3.5, PETSC_DEFAULT_DOUBLE_PRECISION was used in place of PETSC_DEFAULT_REAL.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
#define PETSC_DEFAULT_REAL PETSC_DEFAULT_DOUBLE_PRECISION
#endif
!Hereafter in this code, use PETSC_DEFAULT_REAL.

! For PETSc versions prior to 3.5, DMDA_BOUNDARY_NONE was used in place of DM_BOUNDARY_NONE.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 5))
#define DM_BOUNDARY_NONE DMDA_BOUNDARY_NONE
#endif
!Hereafter in this code, use DM_BOUNDARY_NONE.

! For PETSc versions prior to 3.9, MatSolverPackage was called MatSolverType.
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 9))
#define MatSolverType MatSolverPackage
#define PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#define PCFactorGetMatSolverType PCFactorGetMatSolverPackage
#define PCFactorSetUpMatSolverType PCFactorSetUpMatSolverPackage
#endif
!Hereafter in this code, use MatSolverType.

! The remaining code ensures that the proper PETSc include files and modules are included in every subroutine by including just one line:
! #include "PETScVersions.F90"
#if (PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 6))
  ! Version <= 3.5
#include <finclude/petscsnesdef.h>
#elif (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR < 8)
  ! Version 3.6-3.7
#include <petsc/finclude/petscsnesdef.h>
#else
  ! Version 3.8
#include <petsc/finclude/petscsnes.h>
#endif
  use petscsnes
