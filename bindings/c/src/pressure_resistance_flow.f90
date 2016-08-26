module pressure_resistance_flow_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_prq_c(bcinlet,bcoutlet,targetflow) bind(C, name="evaluate_prq_c")

use pressure_resistance_flow, only: evaluate_prq
use arrays,only: dp
implicit none
real(dp),intent(in) :: bcinlet
real(dp),intent(in) :: bcoutlet
real(dp),intent(in) :: targetflow

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_prq(bcinlet,bcoutlet,targetflow)
#else
call evaluate_prq(bcinlet,bcoutlet,targetflow)
#endif

end subroutine evaluate_prq_c

!###################################################################################
end module pressure_resistance_flow_c
