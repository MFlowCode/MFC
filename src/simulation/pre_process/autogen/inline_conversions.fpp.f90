# 1 "/Users/anand/MFC/src/common/inline_conversions.fpp"
!>  This procedure calculates the speed of sound
        !! @param pres Pressure
        !! @param rho Cell averaged density
        !! @param pi_inf Cell averaged liquid stiffness
        !! @param gamma Cell averaged specific heat ratio
        !! @param H Cell averaged enthalpy
        !! @param adv Advection Variables
        !! @param vel_sum Sum of all velocities
        !! @param q_prim_vf Primitive vars in 1 direction
        !! @param flg Helps determine which conditionals to be called.
            ! flg >= 2: Check all conditionals
            ! flg =  1: Check for alt_soundspeed, otherwise run the 3rd conditional block
            ! flg = 0: Check for alt_soundspeed, otherwise use enthalpy
        !! @param c Speed of sound
        
# 57 "/Users/anand/MFC/src/common/inline_conversions.fpp"
