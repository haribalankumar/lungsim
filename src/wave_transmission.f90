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
subroutine evaluate_wave_transmission()
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_WAVE_TRANSMISSION: EVALUATE_WAVE_PROPAGATION
  use indices
  use arrays, only: dp
  use diagnostics, only: enter_exit
  !n_time,a0,no_freq,a,b,n_bcparams,&
  !  bc_params)
  !integer, intent(in):: no_freq
  !character(len=MAX_STRING_LEN) :: admittance_model
  !integer, intent(in) :: n_time
  !real(dp), intent(in) :: a0
  !real(dp), intent(in) :: a(no_freq)
  !real(dp), intent(in) :: b(no_freq)
  !integer, intent(in) :: n_bcparams
  !real(dp), intent(in) :: bc_params(n_bcparams)

  !type(all_admit_bcs) :: bc

  character(len=60) :: sub_name

  sub_name = 'evalulate_wave_transmission'
  call enter_exit(sub_name,1)

  call enter_exit(sub_name,2)
end subroutine evaluate_wave_transmission
end module wave_transmission
