module initialize_update
  !------------------------------------!
  !Use system parameters to initialize
  !the positions and update the positions.
  !------------------------------------!
  implicit none

  contains

subroutine Initialize_position
  !------------------------------------!
  !Initialize position
  !   This program is used to initialize the position of
  !   Polyelectrolytes and ions, and parameters, energy of
  !   the potential.
  !Input
  !   pos, random_or_uniform
  !Output
  !   pos
  !External Variables
  !   pos, random_or_uniform
  !Routine Referenced:
  !1.subroutine random_grafted
  !   initialize chains by randomly grafting on the plate
  !2.subroutine uniform_grafted
  !   initialize chains by uniformly grafting on the plate
  !3.subroutine initialize_ions
  !   initialize ions in the system
  !------------------------------------!
  use global_variables
  implicit none

  pos=0
  
  if ( random_or_uniform == 0 ) then
    call random_grafted     ! Don't forget periodic condition
  else 
    call uniform_grafted
  end if

  if ( qq /= 0 ) then
    call initialize_ions
  end if

end subroutine Initialize_position

subroutine random_grafted
  !--------------------------------------!
  !Initialize position with random grafted in the box
  !of (-Lx/2, Lx/2), (-Ly/2, Ly/2), (0, Lz)
  !
  !Input
  !  pos
  !Output
  !  pos 
  !External Variables
  !  pos
  !  Lx, Ly, Nx, Ny, Nml
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  real*8  :: R_bond, rij(3), rsqr, rnd1, rnd2
  integer :: i, j, k, l, m, n
  R_bond = 0.97D0

!   !position of PE, random grafted but the chains are straight
!   l = 1                   
!   do i = 1, Nx
!     do j = 1, Ny
! !       !
! !       !Random grafted.
! !       m = 1
! !       do while ( m == 1 )
! !         m = 0
! !         !
! !         !Random Grafted
! !         call random_number(rnd1)
! !         call random_number(rnd2)
! !         pos(l,1) = rnd1*Lx - Lx/2
! !         pos(l,2) = rnd2*Ly - Ly/2
! !         pos(l,3) = 0
! !         !
! !         !Keep the paritcles are not very closed.
! !         do n=1,l-1
! !           call rij_and_rr(rij,rsqr,n,l)
! !           if ( rsqr < 0.7 ) then
! !             m = 1
! !             cycle
! !           end if
! !         end do
! !       end do
!      !
!      !Grafted by cubic crystal lattice
!       pos(l,1)=(i-0.5D0)*Lx/Nx-Lx/2
!       pos(l,2)=(j-0.5D0)*Ly/Ny-Ly/2
!       pos(l,3)=0
!       l = l + 1
!       do k = 2, Nml
!         pos(l,1) = pos(l-k+1,1)
!         pos(l,2) = pos(l-k+1,2)
!         pos(l,3) = (k-1)*R_bond
!         l = l + 1
!       end do
!     end do
!   end do


  l = 1                   
  do i = 1, Ngl
    !
    !Random grafted.
    m = 1
    do while ( m == 1 )
      m = 0
      !
      !Random Grafted
      call random_number(rnd1)
      call random_number(rnd2)
      pos(l,1) = rnd1*Lx - Lx/2
      pos(l,2) = rnd2*Ly - Ly/2
      pos(l,3) = 0
      !
      !Keep the paritcles are not very closed.
      do n=1,l-1
        call rij_and_rr(rij,rsqr,n,l)
        if ( rsqr < 0.7 ) then
          m = 1
          cycle
        end if
      end do
    end do
    l = l + 1
    do k = 2, Nml
      pos(l,1) = pos(l-k+1,1)
      pos(l,2) = pos(l-k+1,2)
      pos(l,3) = (k-1)*R_bond
      l = l + 1
    end do
  end do

end subroutine random_grafted


subroutine uniform_grafted
  !--------------------------------------!
  !The end of chains are uniform grafted 
  !and the chains are also random in the box
  !of (-Lx/2, Lx/2), (-Ly/2, Ly/2), (0, Lz)
  !   
  !Input
  !   pos
  !Output
  !   pos
  !External Variables
  !   pos
  !   Lx, Ly, Lz, Nx, Ny, Nml, pi, 
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  real*8  :: R_bond, rij(3), rsqr, rnd1, rnd2
  integer :: i, j, k, l, m, n
  R_bond = 0.97D0
  !
  !position of PE, uniform grafted and the chains are also random
  l=0
  do i=1, Nx
    do j=1, Ny
      l=l+1
  !      !
  !      !Grafted by cubic crystal lattice
  !       pos(l,1)=(i-0.5)*Lx/Nx-Lx/2
  !       pos(l,2)=(j-0.5)*Ly/Ny-Ly/2
  !       pos(l,3)=0
      !
      !Grafted by face cubic crystal lattice
      pos(l,1) = ( i - 1 + 0.1D0 + mod(j,2) / 2 ) * Lx / Nx - Lx/2
      pos(l,2) = ( j-0.5D0 ) * Ly / Ny - Ly/2
      pos(l,3) = 0
      do k=2, Nml
        l=l+1
        m=1
        do while ( m == 1 )
          m = 0
          !
          !New particle is uniform distributed on the surface of the sphere of 
          !former particle with the radius R_bond
          call random_number(rnd1)
          call random_number(rnd2)
          pos(l,1) = pos(l-1,1) + R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
          pos(l,2) = pos(l-1,2) + R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
          pos(l,3) = pos(l-1,3) + R_bond*cos(pi*rnd1**2) !rnd1**2 not rnd1 
                                   !means the chains will tend to striaght
          !
          !Periodic condition
          call periodic_condition( pos(l,1:2) )
          !
          !Keep the particle will not too close to the former paritcles
          do n = 1, l-1
            call rij_and_rr(rij,rsqr,n,l)
            if ( rsqr < 0.8 .or. pos(l,3) < 0.9 ) then
              m = 1
              cycle
            end if
            if ( pos(l,3) < 0.8 .or. ( Lz-pos(l,3) ) < 0.8 ) then
              m = 1
              cycle
            end if              
          end do
        end do
      end do
    end do
  end do

end subroutine uniform_grafted


subroutine initialize_ions
  !--------------------------------------!
  !Initialize position of anions of PE and ions on PE 
  !in the box of (-Lx/2, Lx/2), (-Ly/2, Ly/2), (0, Lz)
  !
  !Input
  !   pos
  !Output
  !   pos
  !External Variables
  !   pos
  !   Lx, Ly, Lz, NN, Npe, qq
  !Routine Referenced:
  !1. 
  !--------------------------------------!
  use global_variables
  implicit none
  real*8  :: rnd1, rnd2, rnd3, rij(3), rsqr
  integer :: i, j, m
  !
  !position of anions of PE
  do i=Npe+1, NN
    m=1
    do while (m==1)
      m=0
      call random_number(rnd1)
      call random_number(rnd2)
      call random_number(rnd3)
      pos(i,1)=rnd1*Lx-Lx/2
      pos(i,2)=rnd2*Ly-Ly/2
      pos(i,3)=rnd3**3*(Lz-1.8D0)+0.9D0
      do j=1,i-1
        call rij_and_rr(rij,rsqr,i,j)
        if (rsqr<0.8) then
          m=1
          cycle
        end if
      end do
    end do
    pos(i,4)=-qq/abs(qq)
  end do
  !
  !ions on PE
  do i = 1, Npe
    if ( man == 1 ) then
      if ( mod(i,Nml) /= 1 ) then
        pos(i,4) = qq
      end if
    else
      if ( mod(i,man) == 0 ) then
        pos(i,4) = qq
      end if
    end if
  end do

end subroutine initialize_ions


subroutine Monte_Carlo_Move( EE, DeltaE )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out)   :: DeltaE
  integer :: j
  real*8 :: EE1, EE2

<<<<<<< HEAD
  do j = 1, NN-Ngl
    call total_energy(EE1)
=======
  do j = 1, NN
!     call total_energy(EE1)
>>>>>>> parent of 6e031f6... 2018083003

    call Choose_Particle
    call New_Position
    call Delta_Energy(DeltaE)
    call Move_or_not(EE, DeltaE)

    !
    !test EE2-EE1 = DeltaE
!     call total_energy(EE2)
!     write(*,*) EE2 - EE1, DeltaE, EE2, EE1  
  end do

end subroutine Monte_Carlo_Move


subroutine Monte_Carlo_Move_and_Time( EE, DeltaE, time )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(out) :: DeltaE
  real*8, intent(inout)  :: EE
  real*8, dimension(3), intent(out) :: time
  integer :: j
  
  time = 0
  do j = 1, NN-Ngl
    call Choose_Particle
    call New_Position
    call Delta_Energy_time(DeltaE,time)
    call Move_or_not(EE, DeltaE)
  end do

end subroutine Monte_Carlo_Move_and_Time


subroutine choose_particle
  !------------------------------------!
  !This subroutine is used to choose a particle ip to move.
  !   
  !Input
  !   
  !Output
  !   ip
  !External Variables
  !   NN, Nm, Npe, ip
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd                 

  call random_number(rnd)
  ip = int(rnd*NN) + 1
  !
  !The monomer anchored on the plate can't move, so we need to choose again.
  do while( mod(ip,Nml) == 1 .and. ip <= Npe )
    call random_number(rnd)
    ip = int(rnd*NN) + 1
  end do

end subroutine choose_particle


subroutine New_Position
  !--------------------------------------!
  !This program is used to generate new position.
  !   
  !Input
  !   ip
  !Output
  !   pos1
  !External Variables
  !   pos, pos1, dr
  !Routine Referenced:
  !1. Periodic_condition( rr(2) )
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd(3)

  pos_ip0 = pos(ip,:)
  call random_number(rnd)
  pos_ip1(1:3) = pos_ip0(1:3) + (rnd - 0.5D0) * dr
  pos_ip1(4)   = pos_ip0(4)
  call periodic_condition( pos_ip1(1:2) )

end subroutine New_Position


subroutine Move_or_not(EE, DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   pos, pos_ip0, pos_ip1, ip, Beta
  !Routine Referenced:
  !1.
  !Reference:
  !
  !--------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8,  intent(in)   :: DeltaE
  real*8,  intent(inout) :: EE
  real*8  :: rnd
  !
  !Judge whether move or not
  if ( DeltaE < 0 ) then
    pos(ip,1:3) = pos(ip,1:3) + pos_ip1(1:3) - pos_ip0(1:3)
    EE = EE + DeltaE
    if ( pos_ip0(4) /= 0 ) then
      call update_rhok
    end if
  else 
    call random_number(rnd)
    if ( rnd < Exp(-Beta*DeltaE) ) then
      pos(ip,1:3) = pos(ip,1:3) + pos_ip1(1:3) - pos_ip0(1:3)
      EE = EE + DeltaE
      if ( pos_ip0(4) /= 0 ) then
        call update_rhok
      end if
    endif
  endif
end subroutine Move_or_not


subroutine periodic_condition(rr)
  !--------------------------------------!
  !Peridodic condition of position vector
  !rr(2) in slab geometry.
  !   
  !Input
  !   rr
  !Output
  !   rr
  !External Variables
  !   Lx, Ly
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: rr(2)

  if ( rr(1) > Lx/2 ) then
    rr(1) = rr(1) - Lx
  elseif( rr(1) <= -Lx/2 ) then
    rr(1) = rr(1) + Lx
  end if
  if ( rr(2) > Ly/2 ) then
    rr(2) = rr(2) - Ly
  elseif( rr(2) <= -Ly/2 ) then
    rr(2) = rr(2) + Ly
  end if

end subroutine periodic_condition


subroutine rij_and_rr(rij, rsqr, i, j)
  !-----------------------------------------!
  !compute displacement vector and displacement of two particles
  !input:
  !  post(pos or pos1), i, j(particle number) 
  !output:
  !  rij(displacement vecter), rr(square of displacement)
  !External Variant:
  !  Lz(used in period condition)
  !note:
  !  including period condition
  !-----------------------------------------!
  use global_variables
  implicit none
  real*8, dimension(3), intent(out) :: rij
  real*8, intent(out) :: rsqr
  integer, intent(in) :: i
  integer, intent(in) :: j

  rij = pos(i,1:3) - pos(j,1:3)

  if ( rij(1) > Lx/2 ) then
    rij(1) = rij(1) - Lx
  elseif( rij(1) <= -Lx/2 ) then
    rij(1) = rij(1) + Lx
  end if
  if ( rij(2) > Ly/2 ) then
    rij(2) = rij(2) - Ly
  elseif( rij(2) <= -Ly/2 ) then
    rij(2) = rij(2) + Ly
  end if

 !   rij(1) = rij(1) - floor(rij(1)/Lx+0.5)*Lx
 !   rij(2) = rij(2) - floor(rij(2)/Ly+0.5)*Ly

  rsqr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

end subroutine rij_and_rr

end module initialize_update



