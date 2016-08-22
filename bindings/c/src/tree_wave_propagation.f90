module tree_wave_propagation_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_wave_propagation_c(n_time,a0,no_freq,a,b,n_bcparams,bc_params) bind(C, name="evaluate_wave_propagation_c")

use tree_wave_propagation, only: evaluate_wave_propagation
use arrays,only:dp
implicit none
integer, intent(in) :: n_time
integer,intent(in) :: no_freq
real(dp),intent(in) :: a0
real(dp),intent(in) :: a(no_freq),b(no_freq)
integer, intent(in) :: n_bcparams
real(dp), intent(in) :: bc_params(n_bcparams)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_wave_propagation(n_time,a0,no_freq,a,b,n_bcparams,bc_params)
#else
call evaluate_wave_propagation(n_time,a0,no_freq,a,b,n_bcparams,bc_params)
#endif

end subroutine evaluate_wave_propagation_c

!###################################################################################
end module tree_wave_propagation_c
