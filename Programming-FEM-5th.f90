!  ProgrammingFEM5th.f90 
!
!  Author  :  lcy
!  Start   :  2019/12/01
!  Email   :  cunyicom@outlook.com
!  GitHub  :  https://github.com/cunyizju/Programming-FEM-5th-Fortran
!  
!  Acknowledgement :  Prof. I M Smith, Prof. D V Griffiths ( The authors of Progrmming the finite element method )
!                     Prof. X Y Shang ( My supervisor )
!****************************************************************************************************************
!
!  PROGRAM :  ProgrammingFEM5th
!
!  PURPOSE :  Progrmming the finite element method using Microsoft Visual Studio and Intel Fortran.
!             Contains the following subroutines
!               Static Equilibrium Structures()
!               Static Equilibrium Linear Elastic Solids()  
!               Material Nonlinearity()
!               Steady State Flow()
!               Transient Problems()
!               Coupled Problems()
!               Eigenvalue Problems()
!               Forced Vibrations()
!****************************************************************************************************************

program ProgrammingFEM5th
!---------------------------------- declaration statement --------------------------------------    
    IMPLICIT NONE
    
    CALL StaticEquilibriumStructures()
    CALL StaticEquilibriumLinearElasticSolids()  
    !CALL MaterialNonlinearity() !Require more than half an hour. It is 
    CALL SteadyStateFlow()
    CALL TransientProblems()
    CALL CoupledProblems()
    CALL EigenvalueProblems()
    CALL ForcedVibrations()
end program ProgrammingFEM5th
    
   