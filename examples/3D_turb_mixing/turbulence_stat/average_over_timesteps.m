close all; clear all;

t_step_start = 0;
t_step_end = 0;



f_average_Reynolds_stress()
f_average_tke_budget()








%% Functions

%
function f_average_Reynolds_stress()


end

% f_average_tke_budget
function f_average_tke_budget()

    load variables/user_inputs.mat;
    load variables/grid.mat;

    ybeg = -5; yend = 5; ny = 101;
    y = linspace(ybeg,yend,ny);

    T1_averaged = zeros(ny,1);
    T2_averaged = zeros(ny,1);
    T3_averaged = zeros(ny,1);
    P1_averaged = zeros(ny,1);
    P2_averaged = zeros(ny,1);
    P3_averaged = zeros(ny,1);
    D_averaged = zeros(ny,1);
    
    for q = start_idx:Nfiles
        load("results/tke_budget_data/tke_budget_"+string(timesteps(q))+".mat");
        i_start = 1;
        for j = 1:ny
            for i = i_start:length(y_norm_mth) - 1
                if (y_norm_mth(i) <= y(j) && y_norm_mth(i+1) > y(j))
                    T1_averaged(j) = T1_averaged(j) + ((T1(i+1) - T1(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T1(i))/(Nfiles - start_idx + 1);
                    T2_averaged(j) = T2_averaged(j) + ((T2(i+1) - T2(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T2(i))/(Nfiles - start_idx + 1);
                    T3_averaged(j) = T3_averaged(j) + ((T3(i+1) - T3(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T3(i))/(Nfiles - start_idx + 1);
                    P1_averaged(j) = P1_averaged(j) + ((P1(i+1) - P1(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P1(i))/(Nfiles - start_idx + 1);
                    P2_averaged(j) = P2_averaged(j) + ((P2(i+1) - P2(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P2(i))/(Nfiles - start_idx + 1);
                    P3_averaged(j) = P3_averaged(j) + ((P3(i+1) - P3(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P3(i))/(Nfiles - start_idx + 1);
                    D_averaged(j) = D_averaged(j) + ((D(i+1) - D(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + D(i))/(Nfiles - start_idx + 1);
                    i_start = i;
                    break;
                end
            end
        end
    end 

    T0_averaged = T1_averaged + T2_averaged + T3_averaged;
    T_averaged = f_compute_derivative_1d(T0_averaged,y*mth); 
    P_averaged = P1_averaged + P2_averaged + P3_averaged;

    %
    f1 = figure("DefaultAxesFontSize",18);
    plot(y,T_averaged,'-bo','LineWidth',2); hold on; grid on;
    plot(y,P_averaged,'-go','LineWidth',2);
    plot(y,D_averaged,'-ro','LineWidth',2);
    % xlim([-5 5]); %ylim([-0.002 0.002]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    saveas(f1, "results/tke_budget/tke_budget_self_similar","png");
    close(f1);
end
