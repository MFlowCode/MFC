close all; clear all;


% Setup
disp("Start average_tke_over_self_similar ..."); tic;
set_user_inputs(); load variables/user_inputs.mat;

ybeg = -5; yend = 5; ny = 101;
y = linspace(ybeg,yend,ny);

% Array
T0_averaged = zeros(ny,1);
P_averaged = zeros(ny,1);
D_averaged = zeros(ny,1);

% Compute averaged TKE budget
for q = 1:Nfiles
    load("results/tke_budget_data/tstep_"+string(timesteps(q))+".mat");

    % Normalization
    T0 = T0 / (8/mth);  % T / (Delta U^3 / mth)
    P = P / (8/mth);    % P / (Delta U^3 / mth)
    D = D / (8/mth);    % D / (Delta U^3 / mth)

    % Interpolation
    i_start = 1;
    for j = 1:ny
        for i = i_start:length(y_norm_mth) - 1
            if (y_norm_mth(i) <= y(j) && y_norm_mth(i+1) > y(j))
                T0_averaged(j) = T0_averaged(j) + ((T0(i+1) - T0(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T0(i))/Nfiles;
                P_averaged(j) = P_averaged(j) + ((P(i+1) - P(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P(i))/Nfiles;
                D_averaged(j) = D_averaged(j) + ((D(i+1) - D(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + D(i))/Nfiles;
                i_start = i;
                break;
            end
        end
    end
end 

T_averaged = f_compute_derivative_1d(T0_averaged,y*mth); 

% Plot
plot_tke_budget(T_averaged, P_averaged, D_averaged, y, mth);

disp("End of program"); toc;

%% FUNCTIONS
% Compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = f_compute_derivative_1d(fun,s)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector

    % Compute one-sided derivative at the bottom boundary
    dfunds(1) = (fun(2) - fun(1)) / (s(2) - s(1));
    
    % Compute one-sided derivative at the top boundary
    dfunds(end) = (fun(end) - fun(end-1)) / (s(end) - s(end-1));
    
    % Compute two-sided derivatives for interior points
    for i = 2:length(s)-1
        dfunds(i) = (fun(i+1) - fun(i-1)) / (s(i+1) - s(i-1));
    end
end

% Plot TKE budget
function plot_tke_budget(T, P, D, y_norm_mth, mth)

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

    saveas(f1, "results/tke_budget/avg_self_similar","png");
    close(f1);
end
