! AUTO-GENERATED - do not edit directly. Regenerate: ./mfc.sh generate
!
namelist /user_inputs/ Bx0, Ca, R0ref, Re_inv, Web, acoustic, acoustic_source, adap_dt, adap_dt_max_iters, adap_dt_tol, adv_n, &
    & alf_factor, alpha_bar, alt_soundspeed, avg_state, bc_x, bc_y, bc_z, bf_x, bf_y, bf_z, bub_pp, bubble_model, bubbles_euler, &
    & bubbles_lagrange, case_dir, cfl_adap_dt, cfl_const_dt, cfl_target, chem_params, coefficient_of_restitution, &
    & collision_model, collision_time, cont_damage, cont_damage_s, cyl_coord, down_sample, dt, fd_order, fft_wrt, &
    & file_per_process, fluid_pp, g_x, g_y, g_z, hyper_cleaning, hyper_cleaning_speed, hyper_cleaning_tau, hyperelasticity, &
    & hypoelasticity, ib, ib_coefficient_of_friction, ib_state_wrt, ic_beta, ic_eps, int_comp, integral, integral_wrt, k_x, k_y, &
    & k_z, lag_params, low_Mach, m, mixture_err, model_eqns, mp_weno, mpp_lim, muscl_eps, n, n_start, null_weights, &
    & num_bc_patches, num_ibs, num_igr_iters, num_igr_warm_start_iters, num_integrals, num_probes, num_source, &
    & nv_uvm_igr_temps_on_gpu, nv_uvm_out_of_core, nv_uvm_pref_gpu, p, p_x, p_y, p_z, palpha_eps, parallel_io, patch_ib, pi_fac, &
    & poly_sigma, polydisperse, polytropic, precision, pref, prim_vars_wrt, probe, probe_wrt, ptgalpha_eps, qbmm, rdma_mpi, &
    & relax, relax_model, rhoref, riemann_solver, run_time_info, sigma, surface_tension, t_save, t_step_old, t_step_print, &
    & t_step_save, t_step_start, t_step_stop, t_stop, tau_star, teno_CT, thermal, time_stepper, w_x, w_y, w_z, wave_speeds, &
    & weno_Re_flux, weno_avg, weno_eps, x_a, x_b, x_domain, y_a, y_b, y_domain, z_a, z_b, z_domain, &
#:if not MFC_CASE_OPTIMIZATION
    & igr, igr_iter_solver, igr_order, igr_pres_lim, mapped_weno, mhd, muscl_lim, muscl_order, nb, num_fluids, recon_type, &
        & relativity, teno, viscous, weno_order, wenoz, wenoz_q
#:endif
