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
subroutine evaluate_wave_propagation(n_time,a0,no_freq,a,b,n_bcparams,&
    bc_params)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SOLVE_WAVE_PROPAGATION: SOLVE_WAVE_PROPAGATION
  use other_consts, only: MAX_STRING_LEN
  use arrays, only: num_elems,num_nodes,all_admit_bcs
  use diagnostics, only: enter_exit

  integer, intent(in):: no_freq
  character(len=MAX_STRING_LEN) :: admittance_model
  integer, intent(in) :: n_time
  real(dp), intent(in) :: a0
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  integer, intent(in) :: n_bcparams
  real(dp), intent(in) :: bc_params(n_bcparams)

  type(all_admit_bcs) :: bc
  !local variables
  integer :: nf,tt,time_step
  real(dp) :: omega
  complex(dp), allocatable :: eff_admit(:,:)
  complex(dp), allocatable :: char_admit(:,:)
  complex(dp), allocatable :: reflect(:,:)
  complex(dp), allocatable :: prop_const(:,:)
  complex(dp), allocatable :: p_factor(:,:)
  complex(dp) :: pressure(n_time,no_freq+1),flow(n_time,no_freq+1)
  real(dp) :: flow_phase,press_phase
  real(dp) :: flow_amp,press_amp
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
    call terminal_admittance(no_freq,eff_admit,char_admit,bc)
    admittance_model='lachase_modified'
    !calculate characteristic admittance of each branch
    call characteristic_admittance(admittance_model,no_freq,char_admit,prop_const)
    ! calculate effective admittance through the tree
    call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const)

    call pressure_factor(no_freq,p_factor,reflect,prop_const)
   ! !Define the pressure flow waveform change from upstream vessels

    deallocate (eff_admit, STAT = AllocateStatus)
    deallocate (char_admit, STAT = AllocateStatus)
    deallocate (reflect, STAT = AllocateStatus)
    deallocate (prop_const, STAT=AllocateStatus)
    deallocate (p_factor, STAT=AllocateStatus)
  call enter_exit(sub_name,2)
end subroutine evaluate_wave_propagation

!
!##############################################################################
!
!*characteristic_admittance* calculates the characteristic admittance of each
subroutine characteristic_admittance(admittance_model,no_freq,char_admit,prop_const)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAD:"SO_characteristic_admittance: characteristic_admittance
  use other_consts, only: MAX_STRING_LEN
  use indices
  use arrays, only: num_elems,elem_field
  use diagnostics, only: enter_exit

  character(len=MAX_STRING_LEN), intent(in) :: admittance_model
  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: prop_const(1:no_freq,num_elems)
  !local variables
  real(dp) :: L,C,R, G,omega
  real(dp) :: E,h_bar,h,density,viscosity !should be global - maybe express as alpha (i.e. pre multiply)
  integer :: ne,nf
  integer :: exit_status=0
  character(len=60) :: sub_name

  sub_name = 'characteristic_admittance'
  call enter_exit(sub_name,1)
  E=1.5e6_dp !Pa
  h_bar=0.1_dp!this is a fraction of the radius so is unitless
  density=0.10500e-02_dp !g/mm^3
  viscosity=0.3500e-02_dp !pa.s=kg/m.s =.3e-2 = equivalent in g/(mm.s)

  !!write(*,*) 'admittance_model',admittance_model
  do ne=1,num_elems
      if(admittance_model.eq.'lachase_standard')then
        h=h_bar*elem_field(ne_radius_in0,ne)
        C=3.0_dp*PI*elem_field(ne_radius_in0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=density*elem_field(ne_length,ne)/(4*PI*elem_field(ne_radius_in0,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius_in0,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admittance_model.eq.'lachase_modified')then
        h=h_bar*elem_field(ne_radius_in0,ne)
        C=3.0_dp*PI*elem_field(ne_radius_in0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)!
        L=9.0_dp*density*elem_field(ne_length,ne)&
           /(4.0_dp*PI*elem_field(ne_radius_in0,ne)**2)!per unit length
        R=81.0_dp*viscosity*elem_field(ne_length,ne)/ &
             (8.0_dp*PI*elem_field(ne_radius_in0,ne)**4) !laminar resistance per unit length
        G=0.0_dp
        !!write(*,*) 'ne', ne, h, C, R,L,viscosity
      elseif(admittance_model.eq.'zhu_chesler')then
        h=h_bar*elem_field(ne_radius_in0,ne)
        C=3.0_dp*PI*elem_field(ne_radius_in0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=9.0_dp*density*elem_field(ne_length,ne)/(4.0_dp*PI*elem_field(ne_radius_in0,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
              (PI*elem_field(ne_radius_in0,ne)**4) !laminar resistance
        G=0.0_dp
       ! !!write(*,*) 'ne', ne, h, C, R,L
      elseif(admittance_model.eq.'duan_zamir')then
     !like arthurs thesis to be implemented
        print *, "This admitance model isnt yet implemented"
        call exit(exit_status)
      else !Unrecognised admittance model
        print *, "EXITING"
        print *, "Unrecognised admittance model, please check inputs"
        call exit(exit_status)
      endif
    do nf=1,no_freq
      omega=nf*2*PI
      char_admit(nf,ne)=sqrt(G+cmplx(0.0_dp,1.0_dp,8)*omega*C)/sqrt(R+cmplx(0.0_dp,1.0_dp,8)*omega*L)
      prop_const(nf,ne)=sqrt((G+cmplx(0.0_dp,1.0_dp,8)*omega*C)*(R+cmplx(0.0_dp,1.0_dp,8)*omega*L))/elem_field(ne_length,ne)
    ! if(ne.eq.10)then
      !!write(*,*) 'char_admit',ne,nf,  char_admit(nf,ne)
      !!write(*,*) 'prop_const', ne,nf, prop_const(nf,ne)
    ! endif
    enddo!nf
  enddo!ne
  call enter_exit(sub_name,2)
end subroutine characteristic_admittance
!
!##############################################################################
!
!*terminal_admittance* applies chosen admittance boundary conditions at the terminal units
subroutine terminal_admittance(no_freq,eff_admit,char_admit,bc)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_terminal_admittance: terminal_admittance
  use arrays,only: num_elems,all_admit_bcs,units,num_units
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
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
         !temporarily store in eff_admit, to be added to the char admit
          eff_admit(nf,ne)=(1.0_dp+cmplx(0,1.0_dp)*omega*R1*C)/R1
          !!write(*,*) 'ne', ne, 'nf', nf, eff_admit(nf,ne)
          !char_admit(nf,ne)=eff_admit(nf,ne)
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
          !temporarily store in eff_admit, to be added to the char admit
          eff_admit(nf,ne)=(1+cmplx(0,1)*omega*R2*C)/(R1+R2+cmplx(0,1)*omega*R1*R2*C)
        enddo
      enddo

    endif

      call enter_exit(sub_name,2)
end subroutine terminal_admittance
!
!##################################################################
!
!*tree_resistance:* Calculates the total admittance of a tree
  subroutine tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const)
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: no_freq
    complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
    complex(dp), intent(in) :: char_admit(1:no_freq,num_elems)
    complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
    complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)

    character(len=60) :: sub_name
!local variables

    real(dp) :: invres,elem_res(num_elems),omega
    integer :: num2,ne,ne2,nf
    complex(dp) :: daughter_admit

    sub_name = 'tree_admittance'
    call enter_exit(sub_name,1)
    reflect(:,:)=cmplx(0.0_dp,0.0_dp,8)
    do nf=1,no_freq
      omega=nf*2*PI
      do ne=num_elems,1,-1
        daughter_admit=cmplx(0.0_dp,0.0_dp,8)!
        do num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals will add one daughter if no branching
           ne2=elem_cnct(1,num2,ne)
           daughter_admit=daughter_admit+eff_admit(nf,ne2)
           !!write(*,*) 'daughter', ne,ne2,nf,daughter_admit
        enddo
        if(elem_cnct(1,0,ne).gt.0)then !not a terminal
           reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
            (char_admit(nf,ne)+daughter_admit)!double checked
                     !  if(ne.ge.9) !write(*,*) ne, nf, reflect(nf,ne), daughter_admit
           eff_admit(nf,ne)=char_admit(nf,ne)*(1&
            -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
            (1&
            +reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))!double checked
         else!a terminal
           daughter_admit=eff_admit(nf,ne)
           !write(*,*) eff_admit(nf,ne),nf,ne
           reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
            (char_admit(nf,ne)+daughter_admit)
            !if(ne.ge.9) !write(*,*) ne, nf, reflect(nf,ne), daughter_admit
           ! !write(*,*) 'term reflect',nf,daughter_admit, char_admit(nf,ne),reflect(nf,ne)
            !now we overwrite the effective admittance of the terminal to include reflection from the daughter.
           eff_admit(nf,ne)=char_admit(nf,ne)*(1&
            -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
            (1&
            +reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
         endif
         if(ne.eq.9)then
         write(*,*) 'Element 9 rf',nf, abs(reflect(nf,ne)),&
           atan2(imagpart(reflect(nf,ne)),realpart(reflect(nf,ne)))/pi*180
         endif
      enddo
    enddo!nf

    call enter_exit(sub_name,2)
  end subroutine tree_admittance
!
!##################################################################
!
!*pressure_factor:* Calculates change in pressure through tree
  subroutine pressure_factor(no_freq,p_factor,reflect,prop_const)
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: no_freq
    complex(dp), intent(inout) :: p_factor(1:no_freq,num_elems)
    complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
    complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)


    character(len=60) :: sub_name
!local variables
    integer :: ne, nf,ne_up
    real(dp) :: omega

    sub_name = 'pressure_factor'
    call enter_exit(sub_name,1)
    do nf=1,no_freq
      omega=nf*2*PI
      do ne=1,num_elems
        !look for upstram element
        if(elem_cnct(-1,0,ne).eq.0)then !no upstream elements, inlet, ignore
        ne_up=1
          p_factor(nf,ne)=(1)* &!assumes input admittance is the same as characteristic admittance for this vessel
            exp(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))/&
            (1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
        else
          ne_up=elem_cnct(-1,1,ne)
          !!write(*,*) 'ne_up', ne, ne_up
          p_factor(nf,ne)=p_factor(nf,ne_up)*(1+reflect(nf,ne_up))* &
            exp(-1.0_dp*prop_const(nf,ne_up)*elem_field(ne_length,ne_up))/&
            (1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))

        endif!neup
        !if(nf.eq.1) !write(*,*) ne, prop_const(nf,ne_up),prop_const(nf,ne), p_factor(nf,ne)
         !reflect(nf,ne)!*exp(cmplx(0,-1.0_dp)*prop_const(nf,ne)*elem_field(ne_length,ne))
      enddo
    enddo!nf       !call pressure_factor(no_freq,p_factor,reflect,prop_const)

    call enter_exit(sub_name,2)
  end subroutine pressure_factor
end module tree_wave_propagation
