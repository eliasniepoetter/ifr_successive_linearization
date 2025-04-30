%% Description
% Script to test the estimation of a terminal region resp. terminal
% ingredients (for the duffing oscillator) proposed by Köhler et al.

close all;
clear;
clc;


% setup sim struct
sim = struct;
sim.c = 2;
sim.m = 2;
sim.k1 = 10;
sim.k3 = 5;
sim.F_duffing = @(x1,x2,u,dt,c,m,k1,k3) [
x1 + dt*x2;
x2 + dt*(-(c/m)*x2 - (k1/m)*x1 - (k3/m)*x1^3 + u);
];
sim.N = 10000;
sim.dt = 0.001;
sim.x_init = [0.5;2];
sim.x = zeros(2,sim.N);
sim.x_lin = zeros(2,sim.N);
sim.x(:,1) = sim.x_init;
sim.x_lin(:,1) = sim.x_init;
sim.u = zeros(1,sim.N-1);
sim.u_lin = zeros(1,sim.N-1);

sim.box_constraint = [-5, 5];

% linearization
sim.x_r = [0;0];
sim.u_r = 0;
[sim.Ad,sim.Bd,sim.Zd,sim.T] = duffing_linearization(sim.x_r,sim.u_r,sim.dt,sim.c,sim.m,sim.k1,sim.k3,sim.box_constraint);

% controller and terminal set synthesis
sim.Q = 1e0*eye(2);
sim.R = 1e1;
[sim.K,sim.S,~] = dlqr(sim.Ad, sim.Bd, sim.Q, sim.R);
sim.Adc = sim.Ad-sim.Bd*sim.K;
sim.eps = 1e4;

sim.Q_star = sim.Q + sim.K'*sim.R*sim.K + sim.eps*eye(2);
sim.P = dlyap(sim.Adc,sim.Q_star);

sim.cu = max(eig(sim.P));
sim.cl = min(eig(sim.P));
sim.ku = norm(sim.K);
sim.cu2 = max(eig(sim.P - sim.Q_star));

sim.L_star = sqrt((sim.cu2+sim.eps)/sim.cu) - sqrt(sim.cu2/sim.cu);

sim.alpha1 = sim.cl * (sim.L_star / (sim.T*(1+sim.ku^2)))^2;
disp(['alpha1: ',num2str(sim.alpha1)]);


%% simulation

for i = 1 : sim.N-1
    sim.x(:,i+1) = sim.F_duffing(sim.x(1,i),sim.x(2,i),sim.u(i),sim.dt,sim.c,sim.m,sim.k1,sim.k3);
    sim.x_lin(:,i+1) = sim.Adc*sim.x_lin(:,i);
    sim.u_lin(i) = -sim.K*sim.x_lin(:,i);
end


%% postprocessing

figure;
hold on;
grid minor;
set(gca,'FontSize',15);
set(gca,'TickLabelInterpreter','latex');
scatter(sim.x_init(1),sim.x_init(2),50,'black','filled','square','DisplayName','initial state')
stairs(sim.x(1,:),sim.x(2,:),'LineWidth',1.25,'Color',[0 0 0],'DisplayName','nonlinear unforced trajectory');
stairs(sim.x_lin(1,:),sim.x_lin(2,:),'LineWidth',1.25,'Color',[0 0 1],'DisplayName','linear dLQR trajectory');
plot([sim.box_constraint(1) sim.box_constraint(1)],[sim.box_constraint(1) sim.box_constraint(2)],'LineWidth',1.125,'Color',[1 0 0],'DisplayName','box constraint');
plot([sim.box_constraint(1) sim.box_constraint(2)],[sim.box_constraint(1) sim.box_constraint(1)],'LineWidth',1.125,'Color',[1 0 0],'HandleVisibility','off');
plot([sim.box_constraint(2) sim.box_constraint(2)],[sim.box_constraint(2) sim.box_constraint(1)],'LineWidth',1.125,'Color',[1 0 0],'HandleVisibility','off');
plot([sim.box_constraint(2) sim.box_constraint(1)],[sim.box_constraint(2) sim.box_constraint(2)],'LineWidth',1.125,'Color',[1 0 0],'HandleVisibility','off');
xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend('Interpreter','latex','Location','northwest');
axis equal;
hold off;
%%

figure;
tiledlayout(3,1);
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:sim.N-1,sim.u_lin,'LineWidth',1.25,'Color',[0 0 0]);
    xlabel('step',Interpreter='latex');
    ylabel('$u$',Interpreter='latex');
    hold off;
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:sim.N,sim.x(1,:),'LineWidth',1.25,'Color',[0 0 0]);
    xlabel('step',Interpreter='latex');
    ylabel('$x_1$',Interpreter='latex');
    hold off;
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:sim.N,sim.x(2,:),'LineWidth',1.25,'Color',[0 0 0]);
    xlabel('step',Interpreter='latex');
    ylabel('$x_2$',Interpreter='latex');
    hold off;


figure;
tiledlayout(1,2);
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(sim.x(1,1:end-1),sim.x(1,2:end),'LineWidth',1.25,'Color',[0 0 0]);
    plot([0.9*min(sim.x(1,:)), 1.1*max(sim.x(1,:))],[0.9*min(sim.x(1,:)), 1.1*max(sim.x(1,:))],'LineWidth',1.0,'Color',[0 0 0],'HandleVisibility','off');
    xlabel('$x_1^{i}$',Interpreter='latex');
    ylabel('$x_1^{i+1}$',Interpreter='latex');
    axis equal;
    hold off;
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(sim.x(2,1:end-1),sim.x(2,2:end),'LineWidth',1.25,'Color',[0 0 0]);
    plot([0.9*min(sim.x(2,:)), 1.1*max(sim.x(2,:))],[0.9*min(sim.x(2,:)), 1.1*max(sim.x(2,:))],'LineWidth',1.0,'Color',[0 0 0],'HandleVisibility','off');
    xlabel('$x_2^{i}$',Interpreter='latex');
    ylabel('$x_2^{i+1}$',Interpreter='latex');
    axis equal;
    hold off;








