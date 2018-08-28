module global_variables

  implicit none
  save

!########################constants#########################!
  real*8, parameter :: pi=3.141592653589793D0     !Circumference ratio pi
  real*8, parameter :: gamma=.5772156649015329D0  !Euler constants gamma
!########################constants#########################!

!####################systems coefficient###################!
  integer :: Ngl      !Number of linear chains grafted on plate
  integer :: Nx       !Chains in x direction
  integer :: Ny       !Chains in y direction
  integer :: Nml      !Number of monomers in each chain
  integer :: Nq       !Total charge in the system
  integer :: Npe      !Total monomers in Polyelectrolytes(PE)
  integer :: NN       !Total particles in the system
  integer :: man      !Manning effect, each man particle have one charge
  real*8  :: qq       !Charge of charged monomers
  real*8  :: sigmag   !Grafted density in xy plane
  real*8  :: ratio_xy !Rotio of length x and width y of the box
  real*8  :: Lx       !Length of cell in x direction
  real*8  :: Ly       !Length of cell in y direction
  real*8  :: Lz       !Distance of two plate
  real*8  :: Z_empty  !Empty space ratio of height and length in slab geometry
  real*8  :: Beta     !Beta=1/(kB*T), T is temperature, 
                      !kB is Boltzmann constant
!##################end systems coefficient#################!


!##################running and Histogram###################!
  integer :: restart_or_continue  !restart or continue after breaking off 
  integer :: random_or_uniform
  integer :: StepNum0             !steps of preheating
  integer :: StepNum              !steps of running
  integer :: DeltaStep            !steps of every calculation of physical
                                  !quantities
  integer :: step                 !steps of calculate the physical quantities
  integer :: dstep                !interval of each calculation of the 
                                  !physical quantities
  integer :: multistep            !each longstep recalculate the coulomb force
  real*8  :: dr                   !length of each moving
!
!timing
  real*8  :: started    = 0       !time at starting
  real*8  :: finished   = 0       !time at finishing
  real*8  :: total_time = 0       !total time of the simulation
!
!histogram
  integer :: SizeHist             !number of histogram which is equally divided
!################end running and Histogram#################!

!##########################arrays##########################!
 real*8, allocatable, dimension(:,:) :: pos       !old position array
 real*8, dimension(4) :: pos_ip0                  !old position of ip
 real*8, dimension(4) :: pos_ip1                  !new position of ip
 integer :: ip                                    !The particle that is choosed
!########################end arrays########################!

end module global_variables

