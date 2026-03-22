close all; clear all;

%% Initialization
disp("Initialize run_turbulence ..."); tic;
% Create directories
create_directory();
% Read user inputs
set_user_inputs(); load variables/user_inputs.mat;
% Set grids
create_grid(); load variables/user_inputs.mat;
% Allocate array for momentum thickness
mth = zeros(1,Nfiles);
toc;

%% Loop over timesteps
for i = 1:Nfiles

    % Read data in conservative form
    disp(" "); disp("[1/6] Read data ..."); tic;
    qp = f_read_data(strcat("../restart_data/lustre_",int2str(timesteps(i)),".dat")); toc;

    % Compute mean quantities
    disp("[2/6] Compute mean and fluctuating quantities ..."); tic;
    [qp_mean qp_fluc] = f_compute_qp_mean_fluc(qp); toc;

    % Compute derivatives
    disp("[3/6] Compute velocity derivatives ..."); tic;
    [dvel_ds dvelmean_dy dvelfluc_ds] = ...
        f_compute_vel_derivatives(qp(momxb:momxe,:,:,:), qp_mean(momxb:momxe,:), qp_fluc(momxb:momxe,:,:,:)); toc;

    % Compute mixing layer thickness
    disp("[4/6] Compute mixing layer thickness ..."); tic;
    [vth mth(i) y_norm_mth y_norm_vth] = f_compute_mixing_layer_thickness(qp_mean,dvelmean_dy); toc;

    % Compute Reynolds stress
    disp("[5/6] Compute Reynolds stress ..."); tic;
    if (Reynolds_stress)
        f_compute_Reynolds_stress(qp(1,:,:,:), qp_fluc(momxb:momxe,:,:,:), y_norm_vth, timesteps(i)); 
    else
        disp("skipped");
    end
    toc;

    % Compute TKE budget
    disp("[6/6] Compute TKE budget ..."); tic;
    if (tke_budget)
        f_compute_tke_budget(dvel_ds, dvelmean_dy, dvelfluc_ds, ...
                        squeeze(qp(1,:,:,:)), qp_fluc(momxb:momxe,:,:,:), squeeze(qp(E_idx,:,:,:)), ...
                        y_norm_mth, mth(i), timesteps(i)); 
    else
        disp("skipped");
    end
    toc;
end
disp("End of loop");

% Compute mixing layer growth rate
disp("Compute growth rate ..."); tic;
plot_growth_rate(time,mth); toc;

disp("End of program");


%% FUNCTIONS
% Create directory for saving results
function create_directory()
    if exist("variables", 'dir')
        delete variables/*.mat;
    else
        mkdir(strcat("variables"));
    end
    if ~exist(strcat("results"), "dir")
        mkdir(strcat("results"));
    end
    if ~exist(strcat("results/tke_budget"), "dir")
        mkdir(strcat("results/tke_budget"));
    end
    if ~exist(strcat("results/tke_budget_data"), "dir")
        mkdir(strcat("results/tke_budget_data"));
    end
    if ~exist(strcat("results/Reynolds_stress"), "dir")
        mkdir(strcat("results/Reynolds_stress"));
    end
    if ~exist(strcat("results/Reynolds_stress_data"), "dir")
        mkdir(strcat("results/Reynolds_stress_data"));
    end
end

% Create grids
function create_grid()

    load variables/user_inputs.mat;

    % Grid-related parameters
    mp = m + 1; mpp = mp + 1;
    np = n + 1; npp = np + 1;
    pp = p + 1; ppp = pp + 1;
    r = mp*np*pp;

    % Allocate memory
    x_cc = zeros(1,mp); x_cb = zeros(1,mpp);
    y_cc = zeros(1,np); y_cb = zeros(1,npp); 
    z_cc = zeros(1,pp); z_cb = zeros(1,ppp);

    % x-dir
    for i = 1:mpp
        x_cb(i) = (i - 1)*Lx/mp;
    end
    dx = x_cb(2:mpp) - x_cb(1:mp);
    x_cc = (x_cb(2:mpp) + x_cb(1:mp))/2;

    % y-dir
    for i = 1:npp
        y_cb(i) = (i - 1)*Ly/np - 0.5*Ly;
    end
    if (stretch_y == "T")
        y_cb = y_cb / Ly;
        y_a = y_a / Ly;
        y_b = y_b / Ly;
        for j = 1:loops_y
            for i = 1:npp
                y_cb(i) = y_cb(i)/a_y* ...
                                 (a_y + log(cosh(a_y*(y_cb(i) - y_a))) ...
                                      + log(cosh(a_y*(y_cb(i) - y_b))) ...
                                    - 2*log(cosh(a_y*(y_b - y_a)/2)));
            end
        end        
        y_cb = y_cb*Ly;
    end
    dy = y_cb(2:npp) - y_cb(1:np);
    y_cc = (y_cb(2:npp) + y_cb(1:np))/2;

    % z-dir
    for i = 1:ppp
        z_cb(i) = (i - 1)*Lz/pp;
    end
    dz = z_cb(2:ppp) - z_cb(1:pp);
    z_cc = (z_cb(2:ppp) + z_cb(1:pp))/2;

    % Save variables
    save variables/user_inputs.mat
end

% Read lustre files
function qp = f_read_data(filename)

    load variables/user_inputs.mat;

    % Read data
    disp(filename);
    fileID = fopen(filename,'r');
    A = fread(fileID,'double');
    fclose(fileID);

    % Reassign density & velocity components
    qc = permute(reshape(A, mp, np, pp, sys_size),[4 1 2 3]);

    % Convert conservative variables to primitive variables
    qp = f_convert_cons_to_prim(qc);
end

% Convert conservative variables to primitive variables
function qp = f_convert_cons_to_prim(qc)

    load variables/user_inputs.mat;

    qp = zeros(size(qc));

    % density
    qp(contxb:contxe,:,:,:) = qc(contxb:contxe,:,:,:);

    dyn_p = 0;
    for i = momxb:momxe
        qp(i,:,:,:) = qc(i,:,:,:) ./ qc(1,:,:,:); % Velocity
        dyn_p = dyn_p + 0.5*qc(i,:,:,:).*qp(i,:,:,:); % Dynamic pressure
    end

    % pressure
    qp(E_idx,:,:,:) = (1/Gamma)*(qc(E_idx,:,:,:) - dyn_p - Pi_inf);

end

% Compute mean and fluctuating flow field data
function [qp_mean qp_fluc] = f_compute_qp_mean_fluc(qp)

    load variables/user_inputs.mat;

    % Compute mean quantities
    qp_mean = zeros(E_idx,np);

    % mean rho
    qp_mean(contxb:contxe,:) = squeeze(mean(qp(contxb:contxe,:,:,:),[2 4]));

    % favre-averaged velocity
    for i = momxb:momxe
        qp_mean(i,:) = squeeze(mean(qp(1,:,:,:).*qp(i,:,:,:),[2 4])./mean(qp(1,:,:,:),[2 4]));
    end

    % mean pressure
    qp_mean(E_idx,:) = squeeze(mean(qp(E_idx,:,:,:),[2 4]));

    % fluctuation
    qp_fluc = qp(1:E_idx,:,:,:) - permute(repmat(qp_mean, [1 1 mp pp]), [1 3 2 4]);
end

% Compute velocity derivatives
function [dvel_ds dvelmean_dy dvelfluc_ds] = f_compute_vel_derivatives(vel, vel_mean, vel_fluc)

    load variables/user_inputs.mat;

    dvelmean_dy = zeros(3,np);
    dvel_ds = zeros(3,3,mp,np,pp);
    dvelfluc_ds = zeros(3,3,mp,np,pp);

    % favre-averaged velocity derivatives
    dvelmean_dy(1,:) = f_compute_derivative_1d(vel_mean(1,:),y_cc);
    dvelmean_dy(2,:) = f_compute_derivative_1d(vel_mean(2,:),y_cc);
    dvelmean_dy(3,:) = f_compute_derivative_1d(vel_mean(3,:),y_cc);

    if (tke_budget)
        % Compute velocity derivatives
        vel1 = squeeze(vel(1,:,:,:));
        vel2 = squeeze(vel(2,:,:,:));
        vel3 = squeeze(vel(3,:,:,:));

        dvel_ds(1,1,:,:,:) = f_compute_derivative_3d(vel1,x_cc,1);
        dvel_ds(2,1,:,:,:) = f_compute_derivative_3d(vel2,x_cc,1);
        dvel_ds(3,1,:,:,:) = f_compute_derivative_3d(vel3,x_cc,1);

        dvel_ds(1,2,:,:,:) = f_compute_derivative_3d(vel1,y_cc,2);
        dvel_ds(2,2,:,:,:) = f_compute_derivative_3d(vel2,y_cc,2);
        dvel_ds(3,2,:,:,:) = f_compute_derivative_3d(vel3,y_cc,2);

        dvel_ds(1,3,:,:,:) = f_compute_derivative_3d(vel1,z_cc,3);
        dvel_ds(2,3,:,:,:) = f_compute_derivative_3d(vel2,z_cc,3);
        dvel_ds(3,3,:,:,:) = f_compute_derivative_3d(vel3,z_cc,3);

        % fluctuating velocity derivatives
        vel_fluc1 = squeeze(vel_fluc(1,:,:,:));
        vel_fluc2 = squeeze(vel_fluc(2,:,:,:));
        vel_fluc3 = squeeze(vel_fluc(3,:,:,:));
        
        dvelfluc_ds(1,1,:,:,:) = f_compute_derivative_3d(vel_fluc1,x_cc,1);
        dvelfluc_ds(2,1,:,:,:) = f_compute_derivative_3d(vel_fluc2,x_cc,1);
        dvelfluc_ds(3,1,:,:,:) = f_compute_derivative_3d(vel_fluc3,x_cc,1);

        dvelfluc_ds(1,2,:,:,:) = f_compute_derivative_3d(vel_fluc1,y_cc,2);
        dvelfluc_ds(2,2,:,:,:) = f_compute_derivative_3d(vel_fluc2,y_cc,2);
        dvelfluc_ds(3,2,:,:,:) = f_compute_derivative_3d(vel_fluc3,y_cc,2);

        dvelfluc_ds(1,3,:,:,:) = f_compute_derivative_3d(vel_fluc1,z_cc,3);
        dvelfluc_ds(2,3,:,:,:) = f_compute_derivative_3d(vel_fluc2,z_cc,3);
        dvelfluc_ds(3,3,:,:,:) = f_compute_derivative_3d(vel_fluc3,z_cc,3);
    end
end

% Compute mixing layer thickness
function [vth mth y_norm_mth y_norm_vth] = f_compute_mixing_layer_thickness(qp_mean, dvelmean_dy)

    load variables/user_inputs.mat;

    % Compute vorticity thickness
    vth = 2 / max(abs(dvelmean_dy(1,:)),[],"all");
    disp("vth: "+num2str(vth));
    
    % Compute momentum thickness
    f = qp_mean(1,:) .* (1 - qp_mean(momxb, :)) .* (1 + qp_mean(momxb, :)) / 4;
    mth = trapz(y_cc,f);
    disp("mth: "+num2str(mth));

    % Compute mth-normalized y value
    y_norm_mth = y_cc/mth;
    y_norm_vth = y_cc/vth;

end

% Compute Reynolds stress
function f_compute_Reynolds_stress(rho, vel_fluc, y_norm_vth, timestep)

    rho = squeeze(rho);
    vel_fluc1 = squeeze(vel_fluc(1, :, :, :));
    vel_fluc2 = squeeze(vel_fluc(2, :, :, :));
    vel_fluc3 = squeeze(vel_fluc(3, :, :, :));

    % Compute Reynolds stress
    ruu = mean(rho .* vel_fluc1.^2, [1 3]);
    rvv = mean(rho .* vel_fluc2.^2, [1 3]);
    rww = mean(rho .* vel_fluc3.^2, [1 3]);
    ruv = mean(rho .* vel_fluc1.*vel_fluc2, [1 3]);

    % Plot Reynolds stress
    plot_Reynolds_stress(ruu, rvv, rww, ruv, y_norm_vth, timestep);

    % Save data
    save("results/Reynolds_stress_data/tstep_"+string(timestep)+".mat","y_norm_vth","ruu","rvv","rww","ruv");

end

% Compute strain rate tensor
function [SS S11 S12 S13 S22 S23 S33] = f_compute_strain_rate_tensor(dvel_ds);
    S11 = squeeze(dvel_ds(1,1,:,:,:));
    S12 = squeeze(0.5*(dvel_ds(1,2,:,:,:) + dvel_ds(2,1,:,:,:)));
    S13 = squeeze(0.5*(dvel_ds(1,3,:,:,:) + dvel_ds(3,1,:,:,:)));
    S22 = squeeze(dvel_ds(2,2,:,:,:));
    S23 = squeeze(0.5*(dvel_ds(2,3,:,:,:) + dvel_ds(3,2,:,:,:)));
    S33 = squeeze(dvel_ds(3,3,:,:,:));

    SS = S11.^2 + 2*S12.^2 + 2*S13.^2 + S22.^2 + 2*S23.^2 + S33.^2;
end

% Compute TKE budget
function f_compute_tke_budget(dvel_ds, dvelmean_dy, dvelfluc_ds, rho, vel_fluc, pres_fluc, y_norm_mth, mth, timestep)

    load variables/user_inputs.mat;

    % Compute strain-rate tensor
    [SS S11 S12 S13 S22 S23 S33] = f_compute_strain_rate_tensor(dvel_ds);
    S21 = S12; S31 = S13; S32 = S23;

    % Viscous stress tensor
    divU = squeeze(dvel_ds(1,1,:,:,:) + dvel_ds(2,2,:,:,:) + dvel_ds(3,3,:,:,:));
    sig = zeros(3,3,mp,np,pp);
    sig(1,1,:,:,:) = 2/Re*(S11 - divU/3);
    sig(1,2,:,:,:) = 2/Re*(S12);
    sig(1,3,:,:,:) = 2/Re*(S13);
    sig(2,1,:,:,:) = 2/Re*(S21);
    sig(2,2,:,:,:) = 2/Re*(S22 - divU/3);
    sig(2,3,:,:,:) = 2/Re*(S23);
    sig(3,1,:,:,:) = 2/Re*(S31);
    sig(3,2,:,:,:) = 2/Re*(S32);
    sig(3,3,:,:,:) = 2/Re*(S33 - divU/3);

    sigmean = zeros(3,3,np);
    for i = 1:3
        for j = 1:3
            sigmean(i,j,:) = squeeze(mean(sig(i,j,:,:,:),[3 5]));
        end
    end
    sigfluc = sig - permute(repmat(sigmean, [1 1 1 mp pp]), [1 2 4 3 5]);
    
    % Transport (T)
    T1 = -0.5*mean(...
         rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)) ...
       + rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:)) ...
       + rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:)),[1 3]);
    T2 = -mean(pres_fluc.*squeeze(vel_fluc(2,:,:,:)),[1 3]);
    T3 = mean(...
         squeeze(sigfluc(1,2,:,:,:)).*squeeze(vel_fluc(1,:,:,:)) ...
       + squeeze(sigfluc(2,2,:,:,:)).*squeeze(vel_fluc(2,:,:,:)) ...
       + squeeze(sigfluc(3,2,:,:,:)).*squeeze(vel_fluc(3,:,:,:)),[1 3]);
    T0 = T1 + T2 + T3;
    T = f_compute_derivative_1d(T0,y_cc); 

    % Turbulent production (P)
    P1 = -mean(rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(1,:));
    P2 = -mean(rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(2,:));
    P3 = -mean(rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(3,:));
    P = P1 + P2 + P3;

    % Dissipation (D)
    D = -mean(squeeze(sigfluc(1,1,:,:,:).*dvelfluc_ds(1,1,:,:,:) ...
                    + sigfluc(2,1,:,:,:).*dvelfluc_ds(2,1,:,:,:) ... 
                    + sigfluc(3,1,:,:,:).*dvelfluc_ds(3,1,:,:,:) ...
                    + sigfluc(1,2,:,:,:).*dvelfluc_ds(1,2,:,:,:) ...
                    + sigfluc(2,2,:,:,:).*dvelfluc_ds(2,2,:,:,:) ...
                    + sigfluc(3,2,:,:,:).*dvelfluc_ds(3,2,:,:,:) ...
                    + sigfluc(1,3,:,:,:).*dvelfluc_ds(1,3,:,:,:) ...
                    + sigfluc(2,3,:,:,:).*dvelfluc_ds(2,3,:,:,:)  ...
                    + sigfluc(3,3,:,:,:).*dvelfluc_ds(3,3,:,:,:)), [1 3]);

    % Plot TKE budget
    plot_tke_budget(T, P, D, y_norm_mth, mth, timestep);

    % Save data
    save("results/tke_budget_data/tstep_"+string(timestep)+".mat","y_norm_mth","mth","T0","P","D");

end

% Compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = f_compute_derivative_1d(fun,s)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector

    for j = 1:length(s) % for each index in s
        % if at bottom boundary, use 1-sided derivative
        dfunds(1) = (fun(2) - fun(1)) / (s(2) - s(1));
        % if at top boundary, use 1-sided derivative
        dfunds(end) = (fun(end) - fun(end-1)) / (s(end) - s(end-1));
        % otherwise, if in the interior of the domain, 2-sided derivative
        for i = 2:length(s)-1
            dfunds(i) = (fun(i+1) - fun(i-1)) / (s(i+1) - s(i-1));
        end
    end
end

% Compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = f_compute_derivative_3d(fun,s,dir)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector
    
    if (dir == 1)
        % Forward difference at start
        dfunds(1,:,:) = (fun(2,:,:) - fun(1,:,:)) / (s(2) - s(1));
        % Backward difference at end
        dfunds(end,:,:) = (fun(end,:,:) - fun(end-1,:,:)) / (s(end) - s(end-1));
        % Central difference for interior
        for i = 2:length(s)-1
            dfunds(i,:,:) = (fun(i+1,:,:) - fun(i-1,:,:)) / (s(i+1) - s(i-1));
        end
    elseif (dir == 2)
        dfunds(:,1,:) = (fun(:,2,:) - fun(:,1,:)) / (s(2) - s(1));
        dfunds(:,end,:) = (fun(:,end,:) - fun(:,end-1,:)) / (s(end) - s(end-1));
        for i = 2:length(s)-1
            dfunds(:,i,:) = (fun(:,i+1,:) - fun(:,i-1,:)) / (s(i+1) - s(i-1));
        end
    elseif (dir == 3)
        dfunds(:,:,1) = (fun(:,:,2) - fun(:,:,1)) / (s(2) - s(1));
        dfunds(:,:,end) = (fun(:,:,end) - fun(:,:,end-1)) / (s(end) - s(end-1));
        for i = 2:length(s)-1
            dfunds(:,:,i) = (fun(:,:,i+1) - fun(:,:,i-1)) / (s(i+1) - s(i-1));
        end
    else
        disp("Wrong dir input");
    end
end

% Plot Reynolds stress
function plot_Reynolds_stress(ruu, rvv, rww, ruv, y_norm_vth, timestep)

    load variables/user_inputs.mat;
    load reference_data/reference.mat;

    % Sqrt & Normalization
    ruu = sqrt(ruu) / 2;    % sqrt(ruu) / Delta U
    rvv = sqrt(rvv) / 2;    % sqrt(rvv) / Delta U
    rww = sqrt(rww) / 2;    % sqrt(rww) / Delta U
    ruv(ruv > 0) = 0;       
    ruv = sqrt(-ruv) / 2;   % sqrt(-ruv) / Delta U

    % Plot
    f1 = figure("DefaultAxesFontSize", 18); 
    set(f1,"Position",[200 200 1000 700]);

    % sqrt(ruu)/Delta U
    subplot(2,2,1);
    plot(y_norm_vth,ruu,'k-','LineWidth',1); hold on; grid on;
    plot(b1990_rs_ruu(:,1),b1990_rs_ruu(:,2),'g+','LineWidth',1,'MarkerSize',8);
    plot(v2014_rs_ruu(:,1),v2014_rs_ruu(:,2),'bo','LineWidth',1,'MarkerSize',8);
    plot(w2022_rs_ruu(:,1),w2022_rs_ruu(:,2),'r^','LineWidth',1,'MarkerSize',8);
    axis([-1.5 1.5 0 0.25]);
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} u^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(rvv)/Delta U
    subplot(2,2,2); 
    plot(y_norm_vth,rvv,'k-','LineWidth',1); hold on; grid on;
    plot(b1990_rs_rvv(:,1),b1990_rs_rvv(:,2),'g+','LineWidth',1,'MarkerSize',8);
    plot(v2014_rs_rvv(:,1),v2014_rs_rvv(:,2),'bo','LineWidth',1,'MarkerSize',8);
    plot(w2022_rs_rvv(:,1),w2022_rs_rvv(:,2),'r^','LineWidth',1,'MarkerSize',8);
    axis([-1.5 1.5 0 0.25]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} v^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(rww)/Delta U
    subplot(2,2,3); 
    plot(y_norm_vth,rww,'k-','LineWidth',1); hold on; grid on;
    plot(b1990_rs_rww(:,1),b1990_rs_rww(:,2),'g+','LineWidth',1,'MarkerSize',8);
    plot(v2014_rs_rww(:,1),v2014_rs_rww(:,2),'bo','LineWidth',1,'MarkerSize',8);
    plot(w2022_rs_rww(:,1),w2022_rs_rww(:,2),'r^','LineWidth',1,'MarkerSize',8);
    axis([-1.5 1.5 0 0.25]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} w^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(-rvu)/Delta U
    subplot(2,2,4); 
    plot(y_norm_vth,ruv,'k-','LineWidth',1); hold on; grid on;
    plot(b1990_rs_ruv(:,1),b1990_rs_ruv(:,2),'g+','LineWidth',1,'MarkerSize',8);
    plot(v2014_rs_ruv(:,1),v2014_rs_ruv(:,2),'bo','LineWidth',1,'MarkerSize',8);
    plot(w2022_rs_ruv(:,1),w2022_rs_ruv(:,2),'r^','LineWidth',1,'MarkerSize',8);
    axis([-1.5 1.5 0 0.25]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{-\left< \overline{\rho} u^{\prime} v^{\prime}\right>}  / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,2);
    legend("$\mbox{Present}$",...
            "$\mbox{Bell \& Mehta (1990)}$",...
            "$\mbox{Vaghefi (2014)}$",...
            "$\mbox{Wang et al. (2022)}$",...
            'Interpreter','latex','location','northeast');

    saveas(f1, "results/Reynolds_stress/tstep_"+string(timestep),"png");
    close(f1);

end

% Plot TKE budget
function plot_tke_budget(T, P, D, y_norm_mth, mth, timestep)

    % Normalization
    T = T / (8/mth);    % T / (Delta U^3 / mth)
    P = P / (8/mth);    % P / (Delta U^3 / mth)
    D = D / (8/mth);    % D / (Delta U^3 / mth)

    load variables/user_inputs.mat;
    load reference_data/reference.mat;

    % Plot
    f1 = figure("DefaultAxesFontSize",18);
    set(f1,"Position",[200 200 1000 700]);

    % Present
    h1 = plot([-100 -100],[-100 -100],'-k','LineWidth',2); hold on; grid on;
    plot(y_norm_mth,T,'-b','LineWidth',2); 
    plot(y_norm_mth,P,'-g','LineWidth',2);
    plot(y_norm_mth,D,'-r','LineWidth',2);
    xlim([-5 5]); xticks([-5:1:5]);
    ylim([-0.002 0.003]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    % Pantano & Sarkar (2002)
    h2 = plot([-100 -100],[-100 -100],'ko','LineWidth',2,'MarkerSize',8);
    plot(p2002_tke_transport(:,1),p2002_tke_transport(:,2),'bo','LineWidth',2,'MarkerSize',8);
    plot(p2002_tke_production(:,1),p2002_tke_production(:,2),'go','LineWidth',2,'MarkerSize',8);
    plot(p2002_tke_dissipation(:,1),p2002_tke_dissipation(:,2),'ro','LineWidth',2,'MarkerSize',8);
    
    % Rogers & Moser (1994)
    h3 = plot([-100 -100],[-100 -100],'k^','LineWidth',2,'MarkerSize',8);
    plot(r1994_tke_transport(:,1),r1994_tke_transport(:,2),'b^','LineWidth',2,'MarkerSize',8);
    plot(r1994_tke_production(:,1),r1994_tke_production(:,2),'g^','LineWidth',2,'MarkerSize',8);
    plot(r1994_tke_dissipation(:,1),r1994_tke_dissipation(:,2),'r^','LineWidth',2,'MarkerSize',8);
    
    % Vaghefi (2014)
    h4 = plot([-100 -100],[-100 -100],'k+','LineWidth',2,'MarkerSize',8);
    plot(v2014_tke_transport(:,1),v2014_tke_transport(:,2),'b+','LineWidth',2,'MarkerSize',8);
    plot(v2014_tke_production(:,1),v2014_tke_production(:,2),'g+','LineWidth',2,'MarkerSize',8);
    plot(v2014_tke_dissipation(:,1),v2014_tke_dissipation(:,2),'r+','LineWidth',2,'MarkerSize',8);
    
    % Wang et al. (2022)
    h5 = plot([-100 -100],[-100 -100],'k*','LineWidth',2,'MarkerSize',8);
    plot(w2022_tke_transport(:,1),w2022_tke_transport(:,2),'b*','LineWidth',2,'MarkerSize',8);
    plot(w2022_tke_production(:,1),w2022_tke_production(:,2),'g*','LineWidth',2,'MarkerSize',8);
    plot(w2022_tke_dissipation(:,1),w2022_tke_dissipation(:,2),'r*','LineWidth',2,'MarkerSize',8);
    
    legend([h1,h2,h3,h4,h5], {"$\mbox{Present}$", ...
            "$\mbox{Pantano \& Sarkar (2002)}$", ...
            "$\mbox{Rogers \& Moser (1994)}$", ...
            "$\mbox{Vaghefi (2014)}$", ...
            "$\mbox{Wang et al. (2022)}$"}, ...
            'interpreter','latex','location','northeast');

    saveas(f1, "results/tke_budget/tstep_"+string(timestep),"png");
    close(f1);
end

% Compute and plot growth rate
function plot_growth_rate(time, mth)

    % Growth rate
    dmth = (mth(2:end) - mth(1:end-1)) ./ (time(2:end) - time(1:end-1));

    % Normalization
    dmth = dmth / 2; % (dmth/dt) * (1/Delta U)

    % Integrated growth rate
    f_growth_rate = figure("DefaultAxesFontSize", 18); 

    h1 = plot(time(1:end-1),dmth,'-ko','LineWidth',2); hold on; grid on;
    h2 = plot([0 max(time)],[0.012 0.012], 'r--','LineWidth',1);
    h3 = plot([0 max(time)],[0.0135 0.0135], 'b--','LineWidth',1);
    h4 = plot([0 max(time)],[0.014 0.014], 'g--','LineWidth',1);
    h5 = plot([0 max(time)],[0.015 0.015], 'c--','LineWidth',1);
    h6 = plot([0 max(time)],[0.017 0.017], 'm--.','LineWidth',1);
    axis([0 max(time) 0 0.04]);
    xlabel('$t U_1 / \delta_{\omega h}^0$','interpreter','latex');
    ylabel('$\dot{\delta}_{\theta} / \Delta U$','interpreter','latex');
    legend([h1,h2,h3,h4,h5,h6],{"$\mbox{Present}$",...
            "$0.012^{[1]}$",...     % [1] Baltzer & Livescu (2020, JFM) 
            "$0.0135^{[2]}$",...    % [2] Vaghefi (2014, Thesis) 
            "$0.014^{[3,4]}$",...   % [3] Rogers & Moser (1993, PoF) [4] Blakeley et al. (2023, JFM)
            "$0.015^{[5]}$",...     % [5] Sharan, Matheou & Dimotakis (2019, JFM)
            "$0.017^{[6]}$"},...    % [6] Almagro et al. (2017, JFM)
            'Interpreter','latex','location','northeast');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f_growth_rate, "results/growth_rate.png"); 
    close(f_growth_rate); 

end
