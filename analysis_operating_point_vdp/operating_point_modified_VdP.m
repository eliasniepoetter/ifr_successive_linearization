%% Description
% Operating Point compution for the modified VdP
% Compensation of the dynamics at the operating point

close all;
clear;
clc;

% Define symbolic variables
syms x1_sym x2_sym u1_sym u2_sym
x_sym = [x1_sym; x2_sym];
u_sym = [u1_sym; u2_sym];
mu = 1;

% modified Van der Pol equations
f_sym = [x2_sym + u1_sym;
    mu*(1 - x1_sym^2)*x2_sym - x1_sym + u2_sym];

% equilirbium equations
eq1 = f_sym(1) == 0;
eq2 = f_sym(2) == 0;

% analytical symbolic  solutions of the equlibrium equations
equilibrium_x2 = solve(eq1,x2_sym);
equilibrium_x1 = solve(eq2,x1_sym);

% select setpoint x2 and compute steady state input u1_star
setpoint_x2 = -1;
u1_star = double(solve(equilibrium_x2 == setpoint_x2,u1_sym));

% substitue setpoint x2 into equlibrium eq. for x1
equilibrium_x1 = subs(equilibrium_x1,x2_sym,setpoint_x2);

% Plot the bifurcation diagram of x1 to select setpoint graphically
figure;
set(gca,'fontsize', 14);
hold on;
fplot(equilibrium_x1,'LineWidth',1.25);
title('Bifurcation Diagram of the modified Van der Pol with $x_2^*=-1$ fixed','Interpreter','latex');
ylabel('equlibrium state $x_1^*$','Interpreter','latex');
xlabel('steady state input $u_2^*$','Interpreter','latex');
grid on;

% examplary setpoint x1 and stedy state input u2
setpoint_x1 = 2;%-1
u2_star = -1;

%% Simulation of the modified Van der Pol Oscillator
% analogue to classic VdP with one more input

% simulation time and steps
tend = 100;
dt = 0.05;

% initial state & parameters
x0 = [-0.5;-1.5];
u_star = [u1_star;u2_star];
x_star = [setpoint_x1;setpoint_x2];
mu = 1;

% controller design
A = [0,1;-1,mu*(1-setpoint_x1^2)];
B = [1,0;0,1];
K = lqr(A,B,eye(1)*10,1);

if rank([B, A*B]) == 2
    disp('linearized system is controllable!');
else
    disp('linearized system is uncontrollable!');
end

% simulation log
sim = struct;
sim.states = zeros(2,length(0:dt:tend)-1);
sim.states(:,1) = x0;
sim.input = zeros(2,length(0:dt:tend)-1);

for n = 1 : length(0:dt:tend)-1
    x_i = sim.states(:,n);
    u = u_star - K*(x_i-x_star);
    sim.input(:,n) = u;
    x_ip1 = dynamics_modified_step(x_i,mu,u,dt);
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

