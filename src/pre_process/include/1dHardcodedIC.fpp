#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here

#:enddef

#:def Hardcoded1D()

    select case (patch_icpp(patch_id)%hcid)
    case (100)
        ! Put your variable assignments here
        !
        q_prim_vf(momxb)%sf(i,0,0)=8*exp(-200*((x_cc(i)-0.0054/2)/0.0054)**2)
       q_prim_vf(E_idx)%sf(i,0,0)=1.01325*10**5+8.0d0*0.3d0*800.74d0*exp(-200*((x_cc(i)-0.0054/2)/0.0054)**2)
       ! q_prim_vf(1)%sf(i,0,0)=(1.01325*10**5+2.2d0*0.227d0*800.74d0*exp(-400*((x_cc(i)-0.002/2)/0.002)**2))/(8.314*1000*((0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32) *(300+30*exp(-400*((x_cc(i)-0.002/2)/0.002)**2)))
      ! q_prim_vf(15)%sf(i,0,0)=300+50*exp(-300*((x_cc(i)-0.002/2)/0.002)**2)
    !  q_prim_vf(1)%sf(i,0,0)=(1.01325*10**5)/(8.314*1000*((0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32) *(300+50*exp(-400*((x_cc(i)-0.002/2)/0.002)**2)))
       q_prim_vf(1)%sf(i,0,0)=0.227+0.227/800.74*(8*exp(-200*((x_cc(i)-0.0054/2)/0.0054)**2))
     !  q_prim_vf(5)%sf(i,0,0)=0.4+x_cc(i)/0.004*0.2
           q_prim_vf(7)%sf(i,0,0)=1.0
   !  q_prim_vf(8)%sf(i,0,0)=0.6-x_cc(i)/0.004*0.2
    ! q_prim_vf(10)%sf(i,0,0)=0.1
    !q_prim_vf(18)%sf(i,0,0)=0.1
    !q_prim_vf(19)%sf(i,0,0)=0.1
    !q_prim_vf(20)%sf(i,0,0)=0.1
    !q_prim_vf(52)%sf(i,0,0)=x_cc(i)/0.0054*0.4
    !q_prim_vf(5)%sf(i,0,0)=0.4
    ! q_prim_vf(9)%sf(i,0,0)=0.6-x_cc(i)/0.04*0.2
    !   q_prim_vf(9)%sf(i,0,0)=1.0d0
    !   q_prim_vf(15)%sf(i,0,0)=300.0d0
    !q_prim_vf(1)%sf(1,0,0)=(1.01325*10**5+8.0d0*0.227d0*800.74d0*exp(-400*((x_cc(i)-0.002/2)/0.002)**2))/(8314/17.008*300)

    ! ! (0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32)!
    ! Put your variable assignments here
    ! (0.4+x_cc(i)/0.002*0.2)/2.016+(0.6-x_cc(i)/0.002*0.2)/32)
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
