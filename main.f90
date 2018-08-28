program main
use global_variables
use input_output
use initialize_update
use compute_energy
implicit none

!##########Data Dictionary############!
  integer :: i, j, k 
  real*8  :: EE, EE1=0, st, fn
  real*8  :: DeltaE, time(3)
!#####################################!

!#############Initialize##############!  
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
  if (restart_or_continue <= 1 ) then
    !
    !initialize position
    call Initialize_position
    call write_pos
    call write_pos1(1)
    call write_hist
    !
    !initialize energy and parameters of potential
    call initialize_energy_parameters
    !
    !Error analysis of Ewald sum
    if ( qq /=0 ) then
      call error_analysis
    end if
    !
    !Compute total energy
    call total_energy(EE)
    i=1
  else
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !initialize energy and parameters of potential
    call initialize_energy_parameters
    !
    !Error analysis of Ewald sum
    if ( qq /= 0 ) then
      call error_analysis
    end if
    !
    !Compute total energy
    call total_energy(EE)
  end if
!#####################################!

!##############Preheation#############!
 if ( i <= StepNum0 ) then
  do step = i, StepNum0
    do dstep = 1, DeltaStep
      if (mod(dstep,50)==0) then 
        call Monte_Carlo_Move_and_Time(EE, DeltaE, time)
        call compute_physical_quantities
        call write_physical_quantities( (step-1)*DeltaStep+dstep, &
                                        EE, EE1, DeltaE, time     )
      else
        call Monte_Carlo_Move(EE, DeltaE)
      end if    
      call update_verlet_list
    end do
    if ( mod(step,10) == 0 ) then
      call total_energy(EE1)
      call write_pos1(step)
      call error_analysis
    end if
  end do
  i = step
end if
!#####################################!

!###############Running###############!
  do step=i, StepNum+StepNum0
    do dstep = 1, DeltaStep
      if (mod(dstep,50)==0) then 
        call Monte_Carlo_Move_and_Time(EE, DeltaE, time)
        call histogram
        call compute_physical_quantities
        call write_physical_quantities( (step-1)*DeltaStep+dstep, &
                                        EE, EE1, DeltaE, time     )
      else
        call Monte_Carlo_Move(EE, DeltaE)
      end if
      call update_verlet_list
    end do
    if ( mod(step,10) == 0 ) then
      call total_energy(EE1)
      call write_pos1(step)
      call write_hist
    end if
  end do
!#####################################!

  call cpu_time(finished)
  total_time=finished-started+total_time
  call write_time(total_time)
  write(*,*) 'finished!'

end program main








