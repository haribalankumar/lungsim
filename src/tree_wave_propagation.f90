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
  use indices
  use filenames, only:AIRWAY_EXNODEFILE,AIRWAY_EXELEMFILE,AIRWAY_ELEMFILE
  use other_consts, only: MAX_STRING_LEN
  use arrays, only: num_elems,num_nodes,all_admit_bcs,elem_field
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
  integer,parameter :: fid = 10
  integer,parameter :: fid2 = 20
  integer,parameter :: fid3 = 30
  integer :: nf,tt,time_step
  real(dp) :: omega
  complex(dp), allocatable :: eff_admit(:,:)
  complex(dp), allocatable :: char_admit(:,:)
  complex(dp), allocatable :: reflect(:,:)
  complex(dp), allocatable :: prop_const(:,:)
  complex(dp), allocatable :: p_factor(:,:)
  real(dp) :: pressure(n_time,no_freq+1),flow(n_time),inlet_pressure(n_time,no_freq+1),&
  inlet_flow(n_time,no_freq+1),incident_flow(n_time,no_freq+1),reflected_flow(n_time,no_freq+1)&
  ,velocity(n_time)
  real(dp) :: flow_phase,press_phase,press_phase2,E,h,h_bar,C
  real(dp) :: flow_amp,press_amp,press_amp2,harmonic_scale,time
  integer :: AllocateStatus
  integer :: point1,point2

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
      bc%three_parameter%admit_P2=bc_params(3)
      bc%three_parameter%admit_P3=bc_params(4)
    elseif(bc_params(1)==4.0_dp)then
      bc%bc_type='two_wk_plus'
      bc%four_parameter%admit_P1=bc_params(2)
      bc%four_parameter%admit_P2=bc_params(3)
      bc%four_parameter%admit_P3=bc_params(4)
      bc%four_parameter%admit_P4=bc_params(5)
    else
    !NOT YET IMPLEMENTED - CALL AN ERROR
    endif
!! SET UP PARAMETERS DEFINING COMPLIANCE MODEL
    harmonic_scale=72.0_dp/60.0_dp !frequency of first harmonic (Hz)
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
    call terminal_admittance(no_freq,eff_admit,char_admit,bc,harmonic_scale)
    admittance_model='lachase_modified'
    !calculate characteristic admittance of each branch
    call characteristic_admittance(admittance_model,no_freq,char_admit,prop_const,harmonic_scale)
    ! calculate effective admittance through the tree
    call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale)

    call pressure_factor(no_freq,p_factor,reflect,prop_const,harmonic_scale)
   ! !Define the pressure flow waveform change from upstream vessels
  time_step=n_time

  pressure=0.0_dp
  flow=0.0_dp
  inlet_pressure=0.0_dp
  inlet_flow=0.0_dp
  incident_flow=0.0_dp
  reflected_flow=0.0_dp


      !write(*,*) 'qfactor'
 ! do nf=1,no_freq
  !   omega=nf*2*PI*harmonic_scale
    !PRessure at insonation site
    !as we are using flow boundary conditions we need to divide flow in by the effective admittance of the first branch
    !to get pressure then we use p_factor to move down the tree. This provides p0 for each branch
     !write(*,*) abs(p_factor(nf,9)/eff_admit(nf,1)),',',&
     !atan2(imagpart(p_factor(nf,9)/eff_admit(nf,1)),&
     !realpart(p_factor(nf,9)/eff_admit(nf,1))),','
    !enddo
    write(*,*) 'insonationsite'
   ! do nf=1,no_freq
    ! omega=nf*2*PI*harmonic_scale
    !    write(*,*) abs(eff_admit(nf,9)*p_factor(nf,9)/eff_admit(nf,1)),',',&
    !   atan2(imagpart(eff_admit(nf,9)*p_factor(nf,9)/eff_admit(nf,1)),&
    !   realpart(eff_admit(nf,9)*p_factor(nf,9)/eff_admit(nf,1))),','
 !        write(*,*) abs(char_admit(nf,9)*p_factor(nf,9)/char_admit(nf,1)),',',&
 !      atan2(imagpart(char_admit(nf,9)*p_factor(nf,9)/char_admit(nf,1)),&
 !      realpart(char_admit(nf,9)*p_factor(nf,9)/char_admit(nf,1))),','




   ! enddo
   open(fid, file = AIRWAY_EXNODEFILE,status='old',action='write',position="append")
   write(fid,fmt=*) 'radius',bc%four_parameter%admit_P4
        write(*,*) 'reflect_coeff'
    do nf=1,no_freq
     omega=nf*2*PI*harmonic_scale
        write(*,*) abs(reflect(nf,1)),atan2(imagpart(reflect(nf,1)),&
       realpart(reflect(nf,1)))
       write(fid,fmt=*) abs(reflect(nf,1)),atan2(imagpart(reflect(nf,1)),&
       realpart(reflect(nf,1)))
    enddo

    write(*,*) 'prop_const'
    do nf=1,no_freq
     omega=nf*2*PI*harmonic_scale
        write(*,*) realpart(prop_const(nf,1)),imagpart(prop_const(nf,1))
       ! write(*,*) abs(exp(-2.0_dp*prop_const(nf,9)*elem_field(ne_length,10))),',',&
       ! atan2(imagpart(exp(-2.0_dp*prop_const(nf,9)*elem_field(ne_length,10))),&
       !realpart(exp(-2.0_dp*prop_const(nf,9)*elem_field(ne_length,10)))),','
    enddo

    close(fid)
  !  write(*,*) 'inlet_ref_coeff'
  !  do nf=1,no_freq
  !   omega=nf*2*PI*harmonic_scale
  !      write(*,*) abs(reflect(nf,1)),',',atan2(imagpart(reflect(nf,1)),&
  !     realpart(reflect(nf,1))),','
  !  enddo

   ! write(*,*) 'inlet reflect_offsett'
   ! do nf=1,no_freq
   ! omega=nf*2*PI*harmonic_scale
   !    write(*,*) abs(exp(-2.0_dp*prop_const(nf,1)*elem_field(ne_length,1))),',',&
   !    atan2(imagpart(exp(-2.0_dp*prop_const(nf,1)*elem_field(ne_length,1))),&
   !    realpart(exp(-2.0_dp*prop_const(nf,1)*elem_field(ne_length,1)))),','
   ! enddo


    open(fid2, file = AIRWAY_EXELEMFILE,status='old',action='write',position="append")
       write(fid2,fmt=*) 'radius',bc%four_parameter%admit_P4
    time=0.0_dp
    do tt=1,time_step
      do nf=1,no_freq
      omega=nf*2*PI*harmonic_scale
        inlet_flow(tt,nf)=a(nf)*cos(omega*time+b(nf))
        incident_flow(tt,nf)=a(nf)*exp(-10.0_dp*realpart(prop_const(nf,1)))&
          *COS(omega*time+b(nf)-10.0_dp*imagpart(prop_const(nf,1)))
        reflected_flow(tt,nf)=a(nf)*abs(reflect(nf,1))*exp((10.0_dp-2*100.0_dp)*realpart(prop_const(nf,1)))&
          * COS(omega*time+b(nf)+(10.0_dp-2*100.0_dp)*imagpart(prop_const(nf,1))+atan2(imagpart(reflect(nf,1)),&
         realpart(reflect(nf,1))))
        pressure(tt,nf)=a(nf)*1000.0_dp/60.0_dp*abs(reflect(nf,1))&
          *COS(omega*time+b(nf)+atan2(imagpart(reflect(nf,1)),realpart(reflect(nf,1))))
        !=$A3*PArameters!$B27*EXP(($G$31-2*$I$31)*$A14)*COS(2*PI()*$E32*F$2+PArameters!$C27+$G$31*$B14-2*$I$31*$B14+$B3)
      enddo
      time=time+0.01_dp
    enddo
    time=0.0_dp
    do tt=1,time_step
      inlet_flow(tt,no_freq+1)=sum(inlet_flow(tt,1:no_freq))+a0
      incident_flow(tt,no_freq+1)=sum(incident_flow(tt,1:no_freq))
      reflected_flow(tt,no_freq+1)=-1.0_dp*sum(reflected_flow(tt,1:no_freq))
      pressure(tt,no_freq+1)=80.0_dp*133.0_dp-sum(pressure(tt,1:no_freq))
      flow(tt)=incident_flow(tt,no_freq+1)+reflected_flow(tt,no_freq+1)+a0
      E=1.5e6_dp
      h_bar=0.1_dp
      h=h_bar*elem_field(ne_radius_out0,1)
      C=3.0_dp*PI*elem_field(ne_radius_out0,1)**3/(2.0_dp*h*E)
      velocity(tt)=(flow(tt)*1000.0_dp/60.0_dp)/(pi*elem_field(ne_radius_out0,1)**2&
        +C*(pressure(tt,no_freq+1)-80.0_dp*133.0_dp))/10.0_dp !cm/s
      write(fid2,fmt=*) time,pressure(tt,no_freq+1),flow(tt),velocity(tt)
      time=time+0.01_dp
    enddo
    close(fid2)

   open(fid3, file = AIRWAY_ELEMFILE,status='old',action='write',position="append")
   write(fid3,fmt=*) 'radius',bc%four_parameter%admit_P4
   write(fid3,fmt=*) 's/d',maxval(velocity(1:time_step))/minval(velocity(1:time_step))
   write(fid3,fmt=*) 'RI',(maxval(velocity(1:time_step))-minval(velocity(1:time_step)))/maxval(velocity(1:time_step))
   write(fid3,fmt=*) 'PI',(maxval(velocity(1:time_step))-minval(velocity(1:time_step)))*time_step/sum(velocity(1:time_step))
   point1=1
   do while(velocity(point1+1).gt.velocity(point1).and.point1.lt.time_step)
     point1=point1+1
   enddo
   do while(velocity(point1+1).lt.velocity(point1).and.point1.lt.time_step)
        point1=point1+1
   enddo
   !point1=minloc(velocity(floor(dble(time_step/4)):ceiling(dble(time_step/2))),1)+floor(dble(time_step/4))
   point2=maxloc(velocity(floor(dble(time_step/3)):ceiling(dble(3*time_step/4))),1)+floor(dble(time_step/3))
   write(*,*) point1,point2
   if(point1.lt.point2)then
     write(fid3,fmt=*) 'notch_height',velocity(point1)&
       -minval(velocity(floor(dble(time_step/4)):ceiling(dble(time_step/2))))
     write(fid3,fmt=*) 'notch_ratio',(maxval(velocity(floor(dble(time_step/3)):ceiling(dble(3*time_step/4))))&
       -minval(velocity(floor(dble(time_step/4)):ceiling(dble(time_step/2)))))&
        /(maxval(velocity(1:time_step))-minval(velocity(1:time_step)))
    else
    write(fid3,fmt=*) point1,point2
   endif

    close(fid3)

    !do tt=1,time_step
   ! write(*,*) 'inlet_pressure'
   ! write(*,*) real(inlet_pressure(:,no_freq+1)),real(inlet_pressure(:,no_freq+1))
   ! write(*,*) 'pressure at insonation site'
   ! write(*,*) real(pressure(:,no_freq+1)),real(pressure(:,no_freq+1))
   !write(*,*) 'inlet flow'
   ! write(*,*) real(inlet_flow(:,no_freq+1)),real(inlet_flow(:,no_freq+1))
   !  write(*,*) 'flow at insonation site'
   ! write(*,*) real(flow(:,no_freq+1)),real(flow(:,no_freq+1))
   ! write(*,*) 'steady comp',a0

    !write(*,*) p_factor(1,1)
    !enddo
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
subroutine characteristic_admittance(admittance_model,no_freq,char_admit,prop_const,harmonic_scale)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAD:"SO_characteristic_admittance: characteristic_admittance
  use other_consts, only: MAX_STRING_LEN
  use indices
  use arrays, only: num_elems,elem_field
  use diagnostics, only: enter_exit

  character(len=MAX_STRING_LEN), intent(in) :: admittance_model
  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: prop_const(1:no_freq,num_elems)
    real(dp), intent(in) :: harmonic_scale
  !local variables
  real(dp) :: L,C,R, G,omega,gen_factor
  real(dp) :: E,h_bar,h,density,viscosity !should be global - maybe express as alpha (i.e. pre multiply)
  integer :: ne,nf
  integer :: exit_status=0
  character(len=60) :: sub_name

  sub_name = 'characteristic_admittance'
  call enter_exit(sub_name,1)
  E=1.5e6_dp !Pa !default
  !E=2.0e6_dp
  h_bar=0.1_dp!this is a fraction of the radius so is unitless
  density=0.10500e-02_dp !g/mm^3
  viscosity=0.3500e-02_dp !pa.s=kg/m.s =.3e-2 = equivalent in g/(mm.s)

  !!write(*,*) 'admittance_model',admittance_model
  do ne=1,num_elems
      if(admittance_model.eq.'lachase_standard')then
        h=h_bar*elem_field(ne_radius_out0,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=density*elem_field(ne_length,ne)/(4*PI*elem_field(ne_radius_out0,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius_out0,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admittance_model.eq.'lachase_modified')then
        !write(*,*) 'rad', ne, elem_field(ne_radius_out0,ne)
        h=h_bar*elem_field(ne_radius_out0,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out0,ne)**3/(2.0_dp*h*E)!
        L=9.0_dp*density&
           /(4.0_dp*PI*elem_field(ne_radius_out0,ne)**2)!per unit length
        R=81.0_dp*viscosity/ &
             (8.0_dp*PI*elem_field(ne_radius_out0,ne)**4) !laminar resistance per unit length
        G=0.0_dp
        !write(*,*) 'ne', ne, h, C, R,L,viscosity
      elseif(admittance_model.eq.'zhu_chesler')then
        h=h_bar*elem_field(ne_radius_out0,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out0,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=9.0_dp*density*elem_field(ne_length,ne)/(4.0_dp*PI*elem_field(ne_radius_out0,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
              (PI*elem_field(ne_radius_out0,ne)**4) !laminar resistance
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
      omega=nf*2*PI*harmonic_scale
      char_admit(nf,ne)=sqrt(G+cmplx(0.0_dp,1.0_dp,8)*omega*C)/sqrt(R+cmplx(0.0_dp,1.0_dp,8)*omega*L)
      prop_const(nf,ne)=sqrt((G+cmplx(0.0_dp,1.0_dp,8)*omega*C)*(R+cmplx(0.0_dp,1.0_dp,8)*omega*L))
     if(ne.eq.3)then
        char_admit(nf,ne)=char_admit(nf,ne)*50.0_dp !40 radials in the placenta
     elseif(ne.eq.2)then
        !char_admit(nf,ne)=char_admit(nf,ne)*20.0_dp
     endif
      !!write(*,*) 'char_admit',ne,nf,  char_admit(nf,ne)
      !!write(*,*) 'prop_const', ne,nf, prop_const(nf,ne)
    ! endif
    enddo!nf
      write(*,*) ne,elem_field(ne_radius_out0,ne)
  enddo!ne

  call enter_exit(sub_name,2)
end subroutine characteristic_admittance
!
!##############################################################################
!
!*terminal_admittance* applies chosen admittance boundary conditions at the terminal units
subroutine terminal_admittance(no_freq,eff_admit,char_admit,bc,harmonic_scale)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_terminal_admittance: terminal_admittance
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
  sub_name = 'terminal_admittance'
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
          !!write(*,*) 'ne', ne, 'nf', nf, eff_admit(nf,ne)
          !char_admit(nf,ne)=eff_admit(nf,ne)
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
          eff_admit(nf,ne)=term_admit!+eff_admit(nf,ne)
        enddo
      enddo

    endif

      call enter_exit(sub_name,2)
end subroutine terminal_admittance
!
!##################################################################
!
!*tree_resistance:* Calculates the total admittance of a tree
  subroutine tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale)
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: no_freq
    complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
    complex(dp), intent(in) :: char_admit(1:no_freq,num_elems)
    complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
    complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
    real(dp), intent(in) :: harmonic_scale

    character(len=60) :: sub_name
!local variables

    real(dp) :: invres,elem_res(num_elems),omega
    integer :: num2,ne,ne2,nf
    complex(dp) :: daughter_admit

    sub_name = 'tree_admittance'
    call enter_exit(sub_name,1)
    reflect(:,:)=cmplx(0.0_dp,0.0_dp,8)
    do nf=1,no_freq
      omega=nf*2*PI*harmonic_scale
      do ne=num_elems,1,-1
        daughter_admit=cmplx(0.0_dp,0.0_dp,8)!
        do num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals will add one daughter if no branching
           ne2=elem_cnct(1,num2,ne)
           daughter_admit=daughter_admit+eff_admit(nf,ne2)
        enddo
        if(elem_cnct(1,0,ne).gt.0)then !not a terminal
           reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
            (char_admit(nf,ne)+daughter_admit)!double checked
                      ! if(ne.eq.1) write(*,*) ne, nf, char_admit(nf,ne), daughter_admit
           eff_admit(nf,ne)=char_admit(nf,ne)*(1&
            -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
            (1&
            +reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))!double checked
         else!a terminal
           daughter_admit=eff_admit(nf,ne)
           !eff_admit(nf,ne) !temp just make it the same
           !write(*,*) char_admit(nf,ne),nf,ne
           reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
            (char_admit(nf,ne)+daughter_admit)
          ! if(ne.eq.2)write(*,*) ne, nf, reflect(nf,ne), daughter_admit,char_admit(nf,ne)
           ! !write(*,*) 'term reflect',nf,daughter_admit, char_admit(nf,ne),reflect(nf,ne)
            !now we overwrite the effective admittance of the terminal to include reflection from the daughter.
           eff_admit(nf,ne)=char_admit(nf,ne)*(1&
            -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
            (1&
            +reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
         endif
         !if(ne.eq.9)then
         !write(*,*) 'Element 9 rf',nf, abs(reflect(nf,ne)),&
         !  atan2(imagpart(reflect(nf,ne)),realpart(reflect(nf,ne)))/pi*180
         !endif
      enddo
    enddo!nf

    call enter_exit(sub_name,2)
  end subroutine tree_admittance
!
!##################################################################
!
!*pressure_factor:* Calculates change in pressure through tree
  subroutine pressure_factor(no_freq,p_factor,reflect,prop_const,harmonic_scale)
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: no_freq
    complex(dp), intent(inout) :: p_factor(1:no_freq,num_elems)
    complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
    complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
    real(dp), intent(in) :: harmonic_scale



    character(len=60) :: sub_name
!local variables
    integer :: ne, nf,ne_up
    real(dp) :: omega

    sub_name = 'pressure_factor'
    call enter_exit(sub_name,1)
    do nf=1,no_freq
      omega=nf*2*PI*harmonic_scale
      do ne=1,num_elems
        !look for upstram element
        if(elem_cnct(-1,0,ne).eq.0)then !no upstream elements, inlet, ignore
        ne_up=1
          p_factor(nf,ne)=(1.0_dp)* &!assumes input admittance is the same as characteristic admittance for this vessel
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
