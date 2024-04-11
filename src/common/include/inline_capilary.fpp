#:def compute_capilary_stress_tensor()

    if (normW > 1d-6) then
        Omega(1, 1) = -sigma*(w2*w2 + w3*w3)/normW

        Omega(2, 1) = sigma*w1*w2/normW
        Omega(1, 2) = Omega(2, 1)

        Omega(2, 2) = -sigma*(w1*w1 + w3*w3)/normW

        if (p > 0) then

            Omega(3, 1) = sigma*w1*w3/normW
            Omega(1, 3) = Omega(3, 1)

            Omega(3, 2) = sigma*w2*w3/normW
            Omega(2, 3) = Omega(3, 2)

            Omega(3, 3) = -sigma*(w1*w1 + w2*w2)/normW

        end if
    else
        Omega(1, 1) = 0d0
        Omega(1, 2) = 0d0
        Omega(2, 1) = 0d0
        Omega(2, 2) = 0d0
        if (p > 0) then
            Omega(1, 3) = 0d0
            Omega(3, 1) = 0d0
            Omega(2, 3) = 0d0
            Omega(3, 2) = 0d0
            Omega(3, 3) = 0d0
        end if
    end if

#:enddef compute_capilary_stress_tensor
