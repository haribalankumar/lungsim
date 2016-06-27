module tree_wave_propagation

!*Brief Description:* Simulating wave propagation in a 1D tree structure
!
!*LICENSE:*
!
!
!
!*Full Description:*
!Simulating wave propagation in a 1D tree structure
!
  use arrays, only: dp
  use other_consts, only: PI
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public evaluate_wave_propagation


contains
!
!##############################################################################
!
subroutine evaluate_wave_propagation(a0,no_freq,a,b,n_bcparams,&
    bc_params)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SOLVE_WAVE_PROPAGATION: SOLVE_WAVE_PROPAGATION
  use arrays, only: num_elems,num_nodes,all_admit_bcs
  use diagnostics, only: enter_exit

  integer, intent(in):: no_freq
  real(dp), intent(in) :: a0
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  integer, intent(in) :: n_bcparams
  real(dp), intent(in) :: bc_params(n_bcparams)

  type(all_admit_bcs) :: bc

  !integer, intent(out) :: you
  !real(dp), intent(inout)  :: pass
  !local variables
  integer :: nf
  real(dp) :: omega
  complex(dp), allocatable :: eff_admit(:,:)
  complex(dp), allocatable :: char_admit(:,:)
  real(dp), allocatable :: reflect(:)
  integer :: AllocateStatus

  character(len=60) :: sub_name

  sub_name = 'evalulate_wave_propagation'
  call enter_exit(sub_name,1)
!! SET UP PARAMETERS DEFINING OUTLET BOUNDARY CONDITIONS
   if(bc_params(1)==1.0_dp)then !note we need to check that the right number of parameters have been input
      bc%bc_type='two_unit_wk'
      bc%two_parameter%admit_P1=bc_params(2)
      bc%two_parameter%admit_P2=bc_params(3)
    elseif(bc_params(1)==2.0_dp)then
      bc%bc_type='three_unit_wk'
      bc%three_parameter%admit_P1=bc_params(2)
      bc%three_parameter%admit_P1=bc_params(3)
      bc%three_parameter%admit_P1=bc_params(4)
    else
    !NOT YET IMPLEMENTED - CALL AN ERROR
    endif
!! SET UP PARAMETERS DEFINING COMPLIANCE MODEL

!!ALLOCATE MEMORY
    allocate (eff_admit(no_freq,num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for eff_admit array ***"
        allocate (char_admit(no_freq,num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for char_admit array ***"
    allocate (reflect(num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for reflect array ***"
    !Apply boundary conditions to terminal units
    call terminal_admittance(no_freq,eff_admit,char_admit,bc)


    ! calculate effective admittance through the tree
    do nf=1,no_freq
      omega=nf*2*PI
    enddo
      write(*,*) bc

    deallocate (eff_admit, STAT = AllocateStatus)
    deallocate (char_admit, STAT = AllocateStatus)
    deallocate (reflect, STAT = AllocateStatus)
  call enter_exit(sub_name,2)
end subroutine evaluate_wave_propagation
!
!##############################################################################
!
!*terminal_admittance* applies chosen admittance boundary conditions at the terminal units
subroutine terminal_admittance(no_freq,eff_admit,char_admit,bc)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_terminal_admittance: terminal_admittance
  use arrays,only: num_elems,all_admit_bcs,units,num_units
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(no_freq,num_elems)
  type(all_admit_bcs) :: bc
  !local variables
  integer :: nf,ne,nunit
  real(dp) :: omega,R1,R2,C
  complex(dp) :: wolm

  character(len=60) :: sub_name

  sub_name = 'terminal_admittance'
  call enter_exit(sub_name,1)
    if(bc%bc_type.eq.'two_unit_wk')then
      R1=bc%two_parameter%admit_P1
      C=bc%two_parameter%admit_P2
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI
        do nunit=1,num_units
          ne=units(nunit)
          eff_admit(nf,ne)=(1+cmplx(0,1)*omega*R1*C)/R1
          char_admit(nf,ne)=eff_admit(nf,ne)
        enddo
      enddo
    elseif(bc%bc_type.eq.'three_unit_wk')then
      R1=bc%three_parameter%admit_P1
      R2=bc%three_parameter%admit_P2
      C=bc%three_parameter%admit_P3
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI
        do nunit=1,num_units
          ne=units(nunit)
          eff_admit(nf,ne)=(1+cmplx(0,1)*omega*R2*C)/(R1+R2+cmplx(0,1)*omega*R1*R2*C)
          char_admit(nf,ne)=eff_admit(nf,ne)
        enddo
      enddo

    endif

      call enter_exit(sub_name,2)
end subroutine terminal_admittance

end module tree_wave_propagation
