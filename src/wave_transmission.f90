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

  !Apply boundary conditions to terminal units
  call boundary_admittance(no_freq,eff_admit,char_admit,bc,harmonic_scale)

  !!DEALLOCATE MEMORY
  deallocate (eff_admit, STAT = AllocateStatus)
  deallocate (char_admit, STAT = AllocateStatus)
  deallocate (reflect, STAT = AllocateStatus)
  deallocate (prop_const, STAT=AllocateStatus)
  deallocate (p_factor, STAT=AllocateStatus)
  call enter_exit(sub_name,2)
end subroutine evaluate_wave_transmission
!
!##############################################################################
!
!*boundary_admittance* applies chosen admittance boundary conditions at the terminal units
subroutine boundary_admittance(no_freq,eff_admit,char_admit,bc,harmonic_scale)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_boundary_admittance: boundary_admittance
  use arrays,only: num_elems,all_admit_bcs,units,num_units
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  type(all_admit_bcs) :: bc
  real(dp), intent(in) :: harmonic_scale
  !local variables
  integer :: nf,ne,nunit
  real(dp) :: omega,R1,R2,C,length,radius,C_term,E,h_bar,density
  real(dp) ::  viscosity,h,L_term,R_term,vein_res
  complex(dp) :: term_admit

  character(len=60) :: sub_name
  sub_name = 'boundary_admittance'
  call enter_exit(sub_name,1)
    if(bc%bc_type.eq.'two_unit_wk')then
      R1=bc%two_parameter%admit_P1
      C=bc%two_parameter%admit_P2
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        do nunit=1,num_units
          ne=units(nunit)
         !temporarily store in eff_admit, to be added to the char admit
          eff_admit(nf,ne)=(1.0_dp+cmplx(0,1.0_dp)*omega*R1*C)/R1
        enddo
      enddo
    elseif(bc%bc_type.eq.'three_unit_wk')then
      R1=bc%three_parameter%admit_P1
      R2=bc%three_parameter%admit_P2
      C=bc%three_parameter%admit_P3
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        do nunit=1,num_units
          ne=units(nunit)
          !temporarily store in eff_admit, to be added to the char admit
          eff_admit(nf,ne)=(1+cmplx(0,1)*omega*R2*C)/(R1+R2+cmplx(0,1)*omega*R1*R2*C)
        enddo
      enddo
    elseif(bc%bc_type.eq.'two_wk_plus')then
      !special case for uterine arteries which are in parallel with shunts
      E=1.5e6_dp !Pa
      h_bar=0.1_dp!this is a fraction of the radius so is unitless
      density=0.10500e-02_dp !g/mm^3
      viscosity=0.3500e-02_dp !pa.s=kg/m.s =.3e-2 = equivalent in g/(mm.s)
      vein_res=0.45_dp
      R1=bc%four_parameter%admit_P1
      C=bc%four_parameter%admit_P2
      length=bc%four_parameter%admit_P3
      radius=bc%four_parameter%admit_P4
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        do nunit=1,num_units
          ne=units(nunit)
         !temporarily store in eff_admit, to be added to the char admit
         !ADMITTANCE DUE TO THE TERMINAL LOAD
          eff_admit(nf,ne)=(1.0_dp+cmplx(0,1.0_dp)*omega*R1*C)/R1
          ! A SECOND ADMITTANCE IN PARALLEL REPRESENTING SHUNTS
          h=h_bar*radius
          C_term=3.0_dp*PI*radius**3/(2.0_dp*h*E)!
          L_term=9.0_dp*density&
             /(4.0_dp*PI*radius**2)!per unit length
          R_term=81.0_dp*viscosity/ &
             (8.0_dp*PI*radius**4) !laminar resistance per unit length
         !G=0.0_dp
          term_admit=sqrt(cmplx(0.0_dp,1.0_dp,8)*omega*C_term)&
            /sqrt(R_term+cmplx(0.0_dp,1.0_dp,8)*omega*L_term)*50.0_dp*1.0_dp
                        term_admit=term_admit/(1+term_admit*vein_res)
          eff_admit(nf,ne)=term_admit+eff_admit(nf,ne)
        enddo
      enddo
    elseif(bc%bc_type.eq.'zero_reflection')then
      !At this stage no need to do anything, need zero reflection coefficient in terminal
    endif
      call enter_exit(sub_name,2)
end subroutine boundary_admittance


end module wave_transmission
