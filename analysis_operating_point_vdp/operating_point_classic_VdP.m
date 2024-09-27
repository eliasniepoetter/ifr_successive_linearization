%% Description
% Operating Point compution for the classic VdP
% Compensation of the dynamics at the operating point

close all;
clear;
clc;

% simulation time and steps
tend = 100;
dt = 0.001;

% initial state & parameters for the simulation
x0 = [6;1.5];
mu = 1;

% the equlibrium with a stady state input u_bar is exactly at x_star = [u_bar,0]
% disadvantage: operating point / setpoint constrained to x2 = 0
u_star = 5;

% controller design
% analytical solution for A and B are parametrized
A = [0,1;-1,mu*(1-u_star^2)];
B = [0;1];
K = lqr(A,B,eye(1)*1,1);

% check controllability via Kalman criteria
if rank([B, A*B]) == 2
    disp('linearized system is controllable!');
else
    disp('linearized system is uncontrollable!');
end

% simulation log
sim = struct;
sim.states = zeros(2,length(0:dt:tend)-1);
sim.states(:,1) = x0;
sim.input = zeros(1,length(0:dt:tend)-1);

% run the simulation
for n = 1 : length(0:dt:tend)-1
    x_i = sim.states(:,n);
    % input computation with steady state input added to shift equlibrium
    u = u_star - K*(x_i-[u_star;0]);
    sim.input = u;
    % call RK4 integration scheme
    x_ip1 = dynamics_step(x_i,mu,u,dt);
    sim.states(:,n+1) = x_ip1;
end

% state space trajectory
figure;
set(gca,'fontsize', 14);
hold on;
title(['Van der Pol Oscillator with $\mu=$',num2str(mu)],'Interpreter','latex');
scatter(x0(1),x0(2),70,'filled','DisplayName','Initial','Marker','square','MarkerFaceColor',[0,0,0]);
scatter(0,0,80,'filled','DisplayName','Origin','Marker','diamond','MarkerFaceColor',[0,0,0]);
plot(sim.states(1,1:n),sim.states(2,1:n),'DisplayName','Trajectory','LineStyle','-','Color',[0 0 0],'LineWidth',1.25);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
axis equal;
legend('Location','southwest','Interpreter','latex');
xlim([-10 10]);
ylim([-10 10]);
grid on;
hold off;




