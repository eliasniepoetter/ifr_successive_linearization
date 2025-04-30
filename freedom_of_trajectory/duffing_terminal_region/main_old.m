%% Description
% ...

close all;
clear;
clc;

% nonlinear duffing equation
f_duffing = @(x_1,x_2,u) [
x_2;
-(c/m)*x_2 - (k_1/mu)*x_1 - (k_3/m)*x_1^3 + u;
];

% setup sim struct
sim = struct;
sim.N = 100;
sim.Ad = [0.9507 -0.3609; 0.0097 0.9330];
sim.Bd = [0.3093; 0.3281];
sim.Q = 1*eye(2);
sim.R = 100;
[sim.K,sim.S,~] = dlqr(sim.Ad, sim.Bd, sim.Q, sim.R);
sim.Adc = sim.Ad-sim.Bd*sim.K;
sim.x_init = [0.2;0.5];
sim.x = zeros(2,sim.N);
sim.x(:,1) = sim.x_init;
sim.eps = 1e-2;

sim.Q_star = sim.Q + sim.K'*sim.R*sim.K + sim.eps*eye(2);
sim.P = dlyap(sim.Adc,sim.Q_star);

sim.cu = max(eig(sim.P));
sim.cl = min(eig(sim.P));
sim.ku = norm(sim.K);
sim.cu2 = max(eig(sim.P - sim.Q_star));

sim.L_star = sqrt((sim.cu2+sim.eps)/sim.cu) - sqrt(sim.cu2/sim.cu);
sim.T = 0.1;

sim.alpha1 = sim.cl * (sim.L_star / (sim.T*(1+sim.ku^2)))^2;


%% simulation

for i = 1 : sim.N-1
    sim.x(:,i+1) = sim.Adc*sim.x(:,i);
end


%% postprocessing

figure;
hold on;
grid minor;
set(gca,'FontSize',15);
set(gca,'TickLabelInterpreter','latex');
scatter(sim.x_init(1),sim.x_init(2),50,'black','filled','square','DisplayName','initial state')
stairs(sim.x(1,:),sim.x(2,:),'LineWidth',1.25,'Color',[0 0 0],'DisplayName','trajectory');
xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend('Interpreter','latex','Location','northwest');
hold off;


figure;
tiledlayout(3,1);
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:sim.N,-sim.K*sim.x,'LineWidth',1.25,'Color',[0 0 0]);
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
    plot([-1, 0.2],[-1, 0.2],'LineWidth',1.0,'Color',[0 0 0],'HandleVisibility','off');
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
    plot([-0.1, 0.6],[-0.1, 0.6],'LineWidth',1.0,'Color',[0 0 0],'HandleVisibility','off');
    xlabel('$x_2^{i}$',Interpreter='latex');
    ylabel('$x_2^{i+1}$',Interpreter='latex');
    axis equal;
    hold off;




