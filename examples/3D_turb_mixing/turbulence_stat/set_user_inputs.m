function set_user_inputs()

    %% Fluid property (Air)
    gamma_a = 1.4;
    pi_inf_a = 0;
    Gamma = 1/(gamma_a - 1);
    Pi_inf = gamma_a*pi_inf_a/(gamma_a - 1);

    %% Reynolds number
    Re = 320;   % Reynolds number

    %% Configuration
    num_fluids = 1; % Number of fluids
    num_dims = 3;   % Number of dimensions

    %% Grids
    Lx = 160.0; % Domain size in x-direction
    Ly = 160.0; % Domain size in y-direction
    Lz = 80.0;  % Domain size in z-direction

    m = 1023;   % Number of grids in x-direction (0, ..., m)
    n = 1023;   % Number of grids in y-direction (0, ..., n)
    p = 511;    % Number of grids in z-direction (0, ..., p)

    % Stretched grid in y-direction
    stretch_y = "T"; 
    a_y = 2; 
    y_a = -0.3*Ly; 
    y_b = 0.3*Ly; 
    loops_y = 2;

    %% Time
    Nt_beg = 0;
    Nt_end = 384;
    Nt_save = 384;
    dt = 0.02609262883235486;
    Nfiles  = (Nt_end - Nt_beg)/Nt_save + 1;
    timesteps = Nt_beg:Nt_save:Nt_end;
    time = timesteps*dt;

    %% Index
    contxb = 1;
    contxe = num_fluids;
    momxb  = contxe + 1;
    momxe  = contxe + num_dims;
    E_idx  = momxe + 1;
    advxb  = E_idx + 1;
    advxe  = E_idx + num_fluids;
    sys_size = advxe;

    save variables/user_inputs.mat;

end