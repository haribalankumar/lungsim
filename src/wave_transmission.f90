module wave_transmission

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
  public evaluate_wave_transmission


contains
!
!##############################################################################
!
subroutine evaluate_wave_transmission(n_time,heartrate,a0,no_freq,a,b,n_bcparams,bc_params)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_WAVE_TRANSMISSION: EVALUATE_WAVE_PROPAGATION
  use indices
  use arrays, only: dp,all_admit_bcs,num_elems,elem_field
  use diagnostics, only: enter_exit

  integer, intent(in) :: n_time
  real(dp), intent(in) :: heartrate
  real(dp), intent(in) :: a0
  integer, intent(in):: no_freq
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  integer, intent(in) :: n_bcparams
  real(dp), intent(in) :: bc_params(n_bcparams)

  type(all_admit_bcs) :: bc

  real(dp) :: harmonic_scale
  real(dp) :: steady_flow
  complex(dp), allocatable :: eff_admit(:,:)
  complex(dp), allocatable :: char_admit(:,:)
  complex(dp), allocatable :: reflect(:,:)
  complex(dp), allocatable :: prop_const(:,:)
  complex(dp), allocatable :: p_factor(:,:)
  integer :: AllocateStatus
  character(len=60) :: sub_name

  sub_name = 'evalulate_wave_transmission'
  call enter_exit(sub_name,1)

  !! SET UP PARAMETERS DEFINING OUTLET BOUNDARY CONDITIONS
  if(bc_params(1).eq.1.0_dp)then !note we need to check that the right number of parameters have been input
    bc%bc_type='two_unit_wk'
    bc%two_parameter%admit_P1=bc_params(2)
    bc%two_parameter%admit_P2=bc_params(3)
  elseif(bc_params(1).eq.2.0_dp)then
    bc%bc_type='three_unit_wk'
    bc%three_parameter%admit_P1=bc_params(2)
    bc%three_parameter%admit_P2=bc_params(3)
    bc%three_parameter%admit_P3=bc_params(4)
  elseif(bc_params(1).eq.4.0_dp)then
    bc%bc_type='two_wk_plus'
    bc%four_parameter%admit_P1=bc_params(2)
    bc%four_parameter%admit_P2=bc_params(3)
    bc%four_parameter%admit_P3=bc_params(4)
    bc%four_parameter%admit_P4=bc_params(5)
  elseif(bc_params(1).eq.5.0_dp)then
    bc%bc_type='zero_reflection'
  else
    print *, 'ERROR: Your boundary condition choice has not yet been implemented'
    call exit(0)
  endif

  !!Determine steady component of flow
  if(a0.eq.0.0_dp)then !Using steady flow solution at inlet as a0
    steady_flow=elem_field(ne_flow,1)!ASSUMING FIRST ELEMENT
  else!otherwise input a0 is used
    steady_flow=a0
  endif

  !! SET UP PARAMETERS DEFINING COMPLIANCE MODEL
  harmonic_scale=heartrate/60.0_dp !frequency of first harmonic (Hz)
  !!ALLOCATE MEMORY
  allocate (eff_admit(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for eff_admit array ***"
  allocate (char_admit(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for char_admit array ***"
  allocate (reflect(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflect array ***"
  allocate (prop_const(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const array ***"
  allocate (p_factor(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for p_factor array ***"
  !initialise admittance
  char_admit=0.0_dp
  eff_admit=0.0_dp

  !!DEALLOCATE MEMORY
  deallocate (eff_admit, STAT = AllocateStatus)
  deallocate (char_admit, STAT = AllocateStatus)
  deallocate (reflect, STAT = AllocateStatus)
  deallocate (prop_const, STAT=AllocateStatus)
  deallocate (p_factor, STAT=AllocateStatus)
  call enter_exit(sub_name,2)
end subroutine evaluate_wave_transmission
end module wave_transmission
