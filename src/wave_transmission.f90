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
subroutine evaluate_wave_transmission(n_time,heartrate,a0,no_freq,a,b,n_adparams,admittance_param,&
  n_model,model_definition)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_WAVE_TRANSMISSION: EVALUATE_WAVE_PROPAGATION
  use indices
  use arrays, only: dp,all_admit_param,num_elems,elem_field,fluid_properties,elasticity_param
  use diagnostics, only: enter_exit

  integer, intent(in) :: n_time
  real(dp), intent(in) :: heartrate
  real(dp), intent(in) :: a0
  integer, intent(in):: no_freq
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  integer, intent(in) :: n_adparams
  real(dp), intent(in) :: admittance_param(n_adparams)
  integer, intent(in) :: n_model
  real(dp), intent(in) :: model_definition(n_model)

  type(all_admit_param) :: admit_param
  type(fluid_properties) :: fluid
  type(elasticity_param) :: elast_param

  character(len=60) :: mesh_type
  real(dp) :: viscosity
  real(dp) :: density
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
  !!MODEL TYPE AND FLUID PROPERTIES
  !mesh_type: can be simple_tree, full_plus_ladder, full_sheet, full_tube The first can be airways, arteries, veins but no special features at the terminal level, the last one has arteries and veins connected by capillary units of some type (lung ladder acinus, lung sheet capillary bed, capillaries are just tubes represented by an element)
  if(model_definition(1).eq.1.0_dp)then
    mesh_type='simple_tree'
  elseif(model_definition(1).eq.2.0_dp)then
    mesh_type='full_plus_ladder'
  !elseif(model_definition(1).eq.3.0_dp)then
  !  full_sheet
  !elseif(model_definition(1).eq.4.0_dp)then
  ! full_tube
  else
    print *, 'ERROR: Your geometry choice has not yet been implemented'
    call exit(0)
  endif
  !viscosity and density of fluid
  if(model_definition(2).eq.1.0_dp)then !BLOOD
    viscosity=fluid%blood_viscosity
    density=fluid%blood_density
  elseif(model_definition(2).eq.2.0_dp)then !AIR
    viscosity=fluid%air_viscosity
    density=fluid%air_density
  else
    viscosity=model_definition(3)
    density=model_definition(4)
  endif

  !!SET UP ADMITTANCE MODEL
  if(admittance_param(1).eq.1.0_dp)then
    admit_param%admittance_type='lachase_standard'
    elast_param%vessel_type='elastic_hooke'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  elseif(admittance_param(1).eq.2.0_dp)then
    admit_param%admittance_type='lachase_modified'
    elast_param%vessel_type='elastic_hooke'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  elseif(admittance_param(1).eq.3.0_dp)then
    admit_param%admittance_type='zhu_chesler'
    elast_param%vessel_type='elastic_hooke'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  elseif(admittance_param(1).eq.4.0_dp)then
    admit_param%admittance_type='duan_zamir'
    elast_param%vessel_type='elastic_alpha'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  else
    print *, 'ERROR: Your admittance model choice has not yet been implemented'
    call exit(0)
  endif

  !! SET UP PARAMETERS DEFINING OUTLET BOUNDARY CONDITIONS
  if(admittance_param(5).eq.1.0_dp)then !note we need to check that the right number of parameters have been input
    admit_param%bc_type='two_unit_wk'
    admit_param%two_parameter%admit_P1=admittance_param(6)
    admit_param%two_parameter%admit_P2=admittance_param(7)
  elseif(admittance_param(5).eq.2.0_dp)then
    admit_param%bc_type='three_unit_wk'
    admit_param%three_parameter%admit_P1=admittance_param(6)
    admit_param%three_parameter%admit_P2=admittance_param(7)
    admit_param%three_parameter%admit_P3=admittance_param(8)
  elseif(admittance_param(5).eq.4.0_dp)then
    admit_param%bc_type='two_wk_plus'
    admit_param%four_parameter%admit_P1=admittance_param(6)
    admit_param%four_parameter%admit_P2=admittance_param(7)
    admit_param%four_parameter%admit_P3=admittance_param(8)
    admit_param%four_parameter%admit_P4=admittance_param(9)
  elseif(admittance_param(5).eq.5.0_dp)then
    admit_param%bc_type='zero_reflection'
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
  call boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
    density,viscosity,elast_param)
   write(*,*) 'admittance_model',admit_param%admittance_type
  !calculate characteristic admittance of each branch
  call characteristic_admittance(no_freq,char_admit,prop_const,harmonic_scale, &
    density,viscosity,admit_param,elast_param)

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
subroutine boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
  density,viscosity,elast_param)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_boundary_admittance: boundary_admittance
  use arrays,only: num_elems,all_admit_param,units,num_units,elasticity_param
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  type(all_admit_param) :: admit_param
  real(dp), intent(in) :: harmonic_scale
  real(dp), intent(in) :: density
  real(dp), intent(in) :: viscosity

  type(elasticity_param) :: elast_param

  !local variables
  integer :: nf,ne,nunit
  real(dp) :: omega,R1,R2,C,length,radius,C_term,E,h_bar
  real(dp) ::  h,L_term,R_term,vein_res
  complex(dp) :: term_admit

  character(len=60) :: sub_name
  sub_name = 'boundary_admittance'
  call enter_exit(sub_name,1)
    if(admit_param%bc_type.eq.'two_unit_wk')then
      R1=admit_param%two_parameter%admit_P1
      C=admit_param%two_parameter%admit_P2
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        do nunit=1,num_units
          ne=units(nunit)
         !temporarily store in eff_admit, to be added to the char admit
          eff_admit(nf,ne)=(1.0_dp+cmplx(0,1.0_dp)*omega*R1*C)/R1
        enddo
      enddo
    elseif(admit_param%bc_type.eq.'three_unit_wk')then
      R1=admit_param%three_parameter%admit_P1
      R2=admit_param%three_parameter%admit_P2
      C=admit_param%three_parameter%admit_P3
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        do nunit=1,num_units
          ne=units(nunit)
          !temporarily store in eff_admit, to be added to the char admit
          eff_admit(nf,ne)=(1+cmplx(0,1)*omega*R2*C)/(R1+R2+cmplx(0,1)*omega*R1*R2*C)
        enddo
      enddo
    elseif(admit_param%bc_type.eq.'two_wk_plus')then
      !special case for uterine arteries which are in parallel with shunts
      E=elast_param%elasticity_parameters(1) !Pa
      h_bar=elast_param%elasticity_parameters(1)!this is a fraction of the radius so is unitless
      vein_res=0.45_dp
      R1=admit_param%four_parameter%admit_P1
      C=admit_param%four_parameter%admit_P2
      length=admit_param%four_parameter%admit_P3
      radius=admit_param%four_parameter%admit_P4
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
    elseif(admit_param%bc_type.eq.'zero_reflection')then
      !At this stage no need to do anything, need zero reflection coefficient in terminal
    endif
      call enter_exit(sub_name,2)
end subroutine boundary_admittance

!
!##############################################################################
!
!*characteristic_admittance* calculates the characteristic admittance of each
subroutine characteristic_admittance(no_freq,char_admit,prop_const,harmonic_scale,&
  density,viscosity,admit_param,elast_param)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAD:"SO_characteristic_admittance: characteristic_admittance
  use other_consts, only: MAX_STRING_LEN
  use indices
  use arrays, only: num_elems,elem_field,elasticity_param,all_admit_param
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: prop_const(1:no_freq,num_elems)
  real(dp), intent(in) :: harmonic_scale
  real(dp), intent(in) :: density
  real(dp), intent(in) :: viscosity

  type(elasticity_param) :: elast_param
  type(all_admit_param) :: admit_param

  !local variables
  real(dp) :: L,C,R, G,omega,gen_factor
  real(dp) :: E,h_bar,h !should be global - maybe express as alpha (i.e. pre multiply)
  integer :: ne,nf
  integer :: exit_status=0
  character(len=60) :: sub_name

  sub_name = 'characteristic_admittance'
  call enter_exit(sub_name,1)
  E=elast_param%elasticity_parameters(1) !Pa
  h_bar=elast_param%elasticity_parameters(2)!this is a fraction of the radius so is unitless

  write(*,*) 'admittance_model',admit_param%admittance_type,E,h_bar
  do ne=1,num_elems
    if(admit_param%admittance_type.eq.'lachase_standard')then
      h=h_bar*elem_field(ne_radius_out0,ne)
      C=3.0_dp*PI*elem_field(ne_radius_out0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
      L=density*elem_field(ne_length,ne)/(4*PI*elem_field(ne_radius_out0,ne)**2)
      R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
          (PI*elem_field(ne_radius_out0,ne)**4) !laminar resistance
      G=0.0_dp
    elseif(admit_param%admittance_type.eq.'lachase_modified')then
      h=h_bar*elem_field(ne_radius_out0,ne)
      C=3.0_dp*PI*elem_field(ne_radius_out0,ne)**3/(2.0_dp*h*E)!
      L=9.0_dp*density&
         /(4.0_dp*PI*elem_field(ne_radius_out0,ne)**2)!per unit length
      R=81.0_dp*viscosity/ &
           (8.0_dp*PI*elem_field(ne_radius_out0,ne)**4) !laminar resistance per unit length
      G=0.0_dp
    elseif(admit_param%admittance_type.eq.'zhu_chesler')then
      h=h_bar*elem_field(ne_radius_out0,ne)
      C=3.0_dp*PI*elem_field(ne_radius_out0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
      L=9.0_dp*density*elem_field(ne_length,ne)/(4.0_dp*PI*elem_field(ne_radius_out0,ne)**2)
      R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius_out0,ne)**4) !laminar resistance
      G=0.0_dp
    elseif(admit_param%admittance_type.eq.'duan_zamir')then
   !like arthurs thesis to be implemented
      print *, "This admitance model isnt yet implemented"
      call exit(exit_status)
    else !Unrecognised admittance model
      print *, "EXITING"
      print *, "Unrecognised admittance model, please check inputs"
      call exit(exit_status)
    endif
    do nf=1,no_freq
      omega=nf*2*PI*harmonic_scale
      char_admit(nf,ne)=sqrt(G+cmplx(0.0_dp,1.0_dp,8)*omega*C)/sqrt(R+cmplx(0.0_dp,1.0_dp,8)*omega*L)
      prop_const(nf,ne)=sqrt((G+cmplx(0.0_dp,1.0_dp,8)*omega*C)*(R+cmplx(0.0_dp,1.0_dp,8)*omega*L))
      !write(*,*) 'TESTING: char_admit',ne,nf,  char_admit(nf,ne)
      !write(*,*) 'TESTING: prop_const', ne,nf, prop_const(nf,ne)
    enddo!nf
  enddo!ne

  call enter_exit(sub_name,2)
end subroutine characteristic_admittance


end module wave_transmission
