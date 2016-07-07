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
  use other_consts, only: MAX_STRING_LEN
  use arrays, only: num_elems,num_nodes,all_admit_bcs
  use diagnostics, only: enter_exit

  integer, intent(in):: no_freq
  character(len=MAX_STRING_LEN) :: admittance_model
  real(dp), intent(in) :: a0
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  integer, intent(in) :: n_bcparams
  real(dp), intent(in) :: bc_params(n_bcparams)


  type(all_admit_bcs) :: bc

  !integer, intent(out) :: you
  !real(dp), intent(inout)  :: pass
  !local variables
  integer :: nf,tt,time_step
  real(dp) :: omega
  complex(dp), allocatable :: eff_admit(:,:)
  complex(dp), allocatable :: char_admit(:,:)
  complex(dp), allocatable :: reflect(:,:)
  complex(dp), allocatable :: prop_const(:,:)
  complex(dp) :: pressure(10,no_freq+1),flow(10,no_freq+1)
  complex(dp) :: flow_phase
  complex(dp) :: flow_amp
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
    allocate (reflect(no_freq,num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for reflect array ***"
        allocate (prop_const(no_freq,num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const array ***"
    !initialise admittance
    char_admit=0.0_dp
    eff_admit=0.0_dp
    write(*,*) no_freq
    !Apply boundary conditions to terminal units
    call terminal_admittance(no_freq,eff_admit,char_admit,bc)
    admittance_model='lachase_standard'
    !calculate characteristic admittance of each branch
    call characteristic_admittance(admittance_model,no_freq,char_admit,prop_const)
    ! calculate effective admittance through the tree
    call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const)
    !Define the pressure flow waveform change from upstream vessels
    time_step=10
    do tt=1,time_step!timesteps
      pressure(tt,:)=0.0_dp
      flow(tt,:)=0.0_dp
      do nf=1,no_freq
        omega=nf*2*PI
        flow_amp=cmplx(a(nf),b(nf))*eff_admit(nf,1)!as this stands with zeros this should give a flow that follows pressure
        pressure(tt,nf)=pressure(tt,nf)+a(nf)*cos(omega*(tt-1)/time_step)&
          +b(nf)*sin(omega*(tt-1)/time_step)
        flow(tt,nf)=flow(tt,nf)+realpart(flow_amp)*cos(omega*(tt-1)/time_step)&
          +imagpart(flow_amp)*sin(omega*(tt-1)/time_step)
      enddo
    enddo
    do tt=1,time_step
      pressure(tt,no_freq+1)=sum(pressure(tt,0:no_freq))+a0
      flow(tt,no_freq+1)=sum(flow(tt,0:no_freq))
    enddo



    write(*,*) abs(pressure(:,no_freq+1))
    write(*,*) abs(flow(:,no_freq+1))

    deallocate (eff_admit, STAT = AllocateStatus)
    deallocate (char_admit, STAT = AllocateStatus)
    deallocate (reflect, STAT = AllocateStatus)
    deallocate (prop_const, STAT=AllocateStatus)
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
  complex(dp), intent(inout) :: char_admit(no_freq,num_elems)
  complex(dp), intent(inout) :: prop_const(no_freq,num_elems)
  !local variables
  real(dp) :: L,C,R, G,omega
  real(dp) :: E,h_bar,h,density,viscosity !should be global - maybe express as alpha (i.e. pre multiply)
  integer :: ne,nf
  integer :: exit_status=0
  character(len=60) :: sub_name

  sub_name = 'characteristic_admittance'
  call enter_exit(sub_name,1)
  E=1.5e6_dp !Pa
  h_bar=0.1_dp!mm - would potentially be a property of each branch
  density=0.10500e-02_dp
  viscosity=0.33600e-02_dp

  write(*,*) 'admittance_model',admittance_model
  do ne=1,num_elems
    do nf=1,no_freq
      omega=nf*2*PI
      if(admittance_model.eq.'lachase_standard')then
        h=h_bar/elem_field(ne_radius_in0,ne)
        C=3.0_dp*PI*elem_field(ne_radius,ne)**3/(2.0_dp*h*E)
        L=density/(4*PI*elem_field(ne_radius,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admittance_model.eq.'lachase_modified')then
        h=h_bar/elem_field(ne_radius_in0,ne)
        C=3.0_dp*PI*elem_field(ne_radius,ne)**3/(2.0_dp*h*E)
        L=9.0_dp*density/(4.0_dp*PI*elem_field(ne_radius,ne)**2)
        R=81.0_dp*viscosity*elem_field(ne_length,ne)/ &
              (8.0_dp*PI*elem_field(ne_radius,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admittance_model.eq.'zhu_chesler')then
        h=h_bar/elem_field(ne_radius_in0,ne)
        C=3.0_dp*PI*elem_field(ne_radius,ne)**3/(2.0_dp*h*E)
        L=9.0_dp*density/(4.0_dp*PI*elem_field(ne_radius,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
              (PI*elem_field(ne_radius,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admittance_model.eq.'duan_zamir')then
     !like arthurs thesis to be implemented
        print *, "This admitance model isnt yet implemented"
        call exit(exit_status)
      else !Unrecognised admittance model
        print *, "EXITING"
        print *, "Unrecognised admittance model, please check inputs"
        call exit(exit_status)
      endif
      char_admit(nf,ne)=sqrt(G+cmplx(0,1)*omega*C)/sqrt(R+cmplx(0,1)*omega*L)
      prop_const(nf,ne)=sqrt((G+cmplx(0,1)*omega*C)*(R+cmplx(0,1)*omega*L))
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
    complex(dp), intent(inout) :: eff_admit(no_freq,num_elems)
    complex(dp), intent(in) :: char_admit(no_freq,num_elems)
    complex(dp), intent(inout) :: reflect(no_freq,num_elems)
    complex(dp), intent(in) :: prop_const(no_freq,num_elems)

    character(len=60) :: sub_name
!local variables

    real(dp) :: invres,elem_res(num_elems),omega
    integer :: num2,ne,ne2,nf
    complex(dp) :: daughter_admit

    sub_name = 'tree_admittance'
    call enter_exit(sub_name,1)
    do nf=1,no_freq
      omega=nf*2*PI
      do ne=num_elems,1,-1
        daughter_admit=0.0_dp
        do num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals will add one daughter if no branching
           ne2=elem_cnct(1,num2,ne)
           daughter_admit=daughter_admit+eff_admit(nf,ne2)
        enddo
        if(elem_cnct(1,0,ne).gt.0)then !not a terminal
          reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
            (char_admit(nf,ne)+daughter_admit)
          eff_admit(nf,ne)=char_admit(nf,ne)*(exp(cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne))&
           -reflect(nf,ne)*exp(-1.0_dp*cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
           (exp(cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne))&
           +reflect(nf,ne)*exp(-1.0_dp*cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne)))
         else!a terminal
           daughter_admit=eff_admit(nf,ne)
           reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
            (char_admit(nf,ne)+daughter_admit)
           eff_admit(nf,ne)=char_admit(nf,ne)*(exp(cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne))&
            -reflect(nf,ne)*exp(-1.0_dp*cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
            (exp(cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne))&
            +reflect(nf,ne)*exp(-1.0_dp*cmplx(0,1)*prop_const(nf,ne)*elem_field(ne_length,ne)))
         endif
      enddo
    enddo!nf

    call enter_exit(sub_name,2)
  end subroutine tree_admittance


end module tree_wave_propagation
