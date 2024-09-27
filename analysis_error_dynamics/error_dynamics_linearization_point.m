%% Summary
% more detailed analytical analysis of the linearization error dyanmics of the VdP
% the linearization point is in this analysis the independet variable

clear;
close all;
clc;

syms x1_l_sym x2_l_sym u_l_sym x1_sym x2_sym u_sym
x_l_sym = [x1_l_sym; x2_l_sym];
x_sym = [x1_sym; x2_sym];
mu = 0.5;

% symbolic VdP
f_sym = [x2_sym;
    mu*(1 - x1_sym^2)*x2_sym - x1_sym + u_sym];

f_l_sym = [x2_l_sym;
    mu*(1 - x1_l_sym^2)*x2_l_sym - x1_l_sym + u_l_sym];

% Compute the Jacobians symbolically
grad_f_x_sym = jacobian(f_sym, x_sym);
grad_f_u_sym = jacobian(f_sym, u_sym);

f_linearized_sym = subs(grad_f_x_sym,x_sym,x_l_sym)*(x_sym-x_l_sym) + subs(grad_f_u_sym,u_sym,u_l_sym)*(u_sym-u_l_sym) + f_l_sym;
grad_f_linearized_x_l_sym = jacobian(f_linearized_sym, x_l_sym);

grad_e_x_l_sym = -grad_f_linearized_x_l_sym;

x = [-6;9];
u = 0;
grad_e_x_l_sym = subs(grad_e_x_l_sym,[x_sym;u_sym],[x;u]);

e_sym = f_sym - f_linearized_sym;
e_sym = subs(e_sym,[x_sym;u_sym;u_l_sym],[x;u;0]);


% Visualization

figure;
tl = tiledlayout(2,2);
title(tl,['Gradient of the linearization error evaluated at [',num2str(x(1)),',',num2str(x(2)),']'], ...
    'interpreter','latex');
nexttile;
    hold on;
    grid on;
    fsurf(grad_e_x_l_sym(1,1));
    colorbar;
    colormap('turbo');
    xlabel('$x_{l1}$','interpreter','latex','FontSize',14);
    ylabel('$x_{l2}$','interpreter','latex','FontSize',14);
    zlabel('$\nabla e_{1,1}$','interpreter','latex','FontSize',14);
nexttile;
    hold on;
    grid on;
    fsurf(grad_e_x_l_sym(1,2));
    colorbar;
    colormap('turbo');
    xlabel('$x_{l1}$','interpreter','latex','FontSize',14);
    ylabel('$x_{l2}$','interpreter','latex','FontSize',14);
    zlabel('$\nabla e_{1,2}$','interpreter','latex','FontSize',14);
nexttile;
    hold on;
    grid on;
    fsurf(grad_e_x_l_sym(2,1));
    colorbar;
    colormap('turbo');
    xlabel('$x_{l1}$','interpreter','latex','FontSize',14);
    ylabel('$x_{l2}$','interpreter','latex','FontSize',14);
    zlabel('$\nabla e_{2,1}$','interpreter','latex','FontSize',14);
nexttile;
    hold on;
    grid on;
    fsurf(grad_e_x_l_sym(2,2));
    colorbar;
    colormap('turbo');
    xlabel('$x_{l1}$','interpreter','latex','FontSize',14);
    ylabel('$x_{l2}$','interpreter','latex','FontSize',14);
    zlabel('$\nabla e_{2,2}$','interpreter','latex','FontSize',14);

figure;
hold on;
tl = tiledlayout(2,1);
title(tl,['$|e|$ evaluated at [',num2str(x(1)),',',num2str(x(2)),']'],'interpreter','latex');
nexttile;
hold on;
    grid on;
    fsurf(abs(e_sym(1)),[-10 10],'ShowContours','on');
    colorbar;
    colormap('turbo');
    xlabel('$x_{l1}$','interpreter','latex','FontSize',14);
    ylabel('$x_{l2}$','interpreter','latex','FontSize',14);
    zlabel('$|e_1|$','interpreter','latex','FontSize',14);
nexttile;
hold on;
    grid on;
    fsurf(abs(e_sym(2)),[-10 10],'ShowContours','on');
    colorbar;
    colormap('turbo');
    xlabel('$x_{l1}$','interpreter','latex','FontSize',14);
    ylabel('$x_{l2}$','interpreter','latex','FontSize',14);
    zlabel('$|e_2|$','interpreter','latex','FontSize',14);










