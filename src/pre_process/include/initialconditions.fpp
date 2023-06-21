#:def analytical_1D()

    if (patch_id == 2) then
         q_prim_vf(contxb)%sf(i, 0, 0) = 1 + 0.2*sin(5*x_cc(i)) 
    endif
       
    if (patch_id == 3) then
         q_prim_vf(contxb)%sf(i, 0, 0) = 1 + 0.2*sin(5*x_cc(i))
    endif

#:enddef



#:def analytical_2D()

#:enddef



#:def analytical_3D()

#:enddef



#:def bubble_pulse_1D()

#:enddef
