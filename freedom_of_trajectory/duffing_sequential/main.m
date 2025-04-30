%% Description
% Main simulation of the successive linearization framework for the
% discrete duffing oscillator. At each setpoit an equilibrium is guaranteed
% by enforcing x2_sp = 0, choosing x1_sp and solving for u_sp.

close all;
clear;
clc;


%% setup variables and parameters

% system parameters
c = 2;
m = 2;
k1 = 10;
k3 = -5;

% discrete dynamics
F_duffing = @(x1,x2,u,dt,c,m,k1,k3) [
x1 + dt*x2;
x2 + dt*(-(c/m)*x2 - (k1/m)*x1 - (k3/m)*x1^3 + u);
];

% simulation parameters
T = 10;
dt = 0.001;
N = T/dt;

% sequqnetial controller parameters
switchingDistance = 0.01;
stepLength = 0.1;

% linearization for pure linear controller
x_0 = [0;0];
u_0 = 0;
[Ad,Bd,Zd] = duffing_linearization(x_0,u_0,dt,c,m,k1,k3);

% initialization of variables
x = zeros(2,N);
x_lc = zeros(2,N);
dx_lc = zeros(2,N);
x_slc = zeros(2,N);
u = zeros(1,N-1);
u_lc = zeros(1,N-1);
u_slc = zeros(1,N-1);

% allocating first entries
x_init = [5;0];
x(:,1) = x_init;
x_lc(:,1) = x_init;
x_slc(:,1) = x_init;
dx_lc(:,1) = x_init-x_0;

% controller synthesis
Q = 5e1*eye(2);
R = 1e1;
[K,S,~] = dlqr(Ad, Bd, Q, R);
Adc = Ad-Bd*K;


%% simulation of discrete dynamics

for i = 1 : N-1
    % update unforced dynamics
    x(:,i+1) = F_duffing(x(1,i),x(2,i),u(i),dt,c,m,k1,k3);

    % update linearized controlled dynamics
    dx_lc(:,i+1) = Adc*dx_lc(:,i);
    u_lc(i) = -K*dx_lc(:,i) + u_0;
    x_lc(:,i+1) = F_duffing(x_lc(1,i),x_lc(2,i),u_lc(i),dt,c,m,k1,k3);

    % update sequentially linearized controller
    Q_seq = 5e1*eye(2);
    R_seq = 1e1;
    if i == 1
        x_sp(:,i) = duffing_setpoint_generation(x_slc(:,i),stepLength);
        u_sp(i) = (k1/m)*x_sp(1,i) + (k3/m)*x_sp(1,i)^3;
        [Ad_seq,Bd_seq,Zd_seq] = duffing_linearization(x_sp(:,i),u_sp(i),dt,c,m,k1,k3);
        [K_seq,S_seq,~] = dlqr(Ad_seq, Bd_seq, Q_seq, R_seq);
        Adc_seq = Ad_seq-Bd_seq*K_seq;
    else
        if (norm(x_sp(:,i-1) - x_slc(:,i-1)) < switchingDistance) && norm(x_sp(:,i-1)) > 0
            x_sp(:,i) = duffing_setpoint_generation(x_slc(:,i),stepLength);
            u_sp(i) = (k1/m)*x_sp(1,i) + (k3/m)*x_sp(1,i)^3;
            [Ad_seq,Bd_seq,Zd_seq] = duffing_linearization(x_sp(:,i),u_sp(i),dt,c,m,k1,k3);
            [K_seq,S_seq,~] = dlqr(Ad_seq, Bd_seq, Q_seq, R_seq);
            Adc_seq = Ad_seq-Bd_seq*K_seq;
        else
            x_sp(:,i) = x_sp(:,i-1);
            u_sp(i) = (k1/m)*x_sp(1,i) + (k3/m)*x_sp(1,i)^3;
        end
    end
    u_slc(i) = -K_seq*(x_slc(:,i)-x_sp(:,i)) + u_sp(i);
    x_slc(:,i+1) = F_duffing(x_slc(1,i),x_slc(2,i),u_slc(i),dt,c,m,k1,k3);

    % evaluate dissipativity
    e(:,i) = F_duffing(x_slc(1,i),x_slc(2,i),u_slc(i),dt,c,m,k1,k3) - x_slc(:,i) - Adc_seq*dx_lc(:,i) - Zd_seq;
    V_dot(i) = x_slc(:,i)'*(K_seq'*R_seq*K_seq+Q_seq)*x_slc(:,i) + 2*x_slc(:,i)'*S_seq*e(:,i);

    progressBar(i, N-1);

end


%% postprocessing

% trajectories
figure;
hold on;
grid minor;
set(gca,'FontSize',15);
set(gca,'TickLabelInterpreter','latex');
scatter(x_init(1),x_init(2),50,'black','filled','square','DisplayName','initial state');
stairs(x(1,:),x(2,:),'LineWidth',1.25,'Color',[0 0 0],'DisplayName','nonlinear unforced system');
stairs(x_lc(1,:),x_lc(2,:),'LineWidth',1.25,'Color',[0 0 1],'DisplayName','nonlinear dLQR system');
stairs(x_slc(1,:),x_slc(2,:),'LineWidth',1.25,'Color',[1 0 0],'DisplayName','nonlinear sequential dLQR system');
scatter(x_sp(1,:),x_sp(2,:),30,'red','filled','DisplayName','setpoint trajectory');
xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend('Interpreter','latex','Location','northwest');
axis equal;
xlim([-5 5]);
ylim([-5 5]);
hold off;

% traces
force_unforced_system = -(c/m)*x(2,:) - (k1/m)*x(1,:) - (k3/m)*x(1,:).^3;
figure;
tiledlayout(3,1);
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:N,100*(force_unforced_system)/max(abs(force_unforced_system)),'LineWidth',1.25,'Color',[0 0 0],'DisplayName',' force unforced system');
    stairs(1:N-1,100*(u_lc)/max(abs(force_unforced_system)),'LineWidth',1.25,'Color',[0 0 1],'DisplayName','dLQR');
    stairs(1:N-1,100*(u_slc)/max(abs(force_unforced_system)),'LineWidth',1.25,'Color',[1 0 0],'DisplayName','sequential dLQR');
    xlabel('step',Interpreter='latex');
    ylabel('$u$ $\left[\%_{|u_{max,uf}|}\right]$',Interpreter='latex');
    legend('Interpreter','latex','Location','northeast');
    hold off;
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:N,x_lc(1,:),'LineWidth',1.25,'Color',[0 0 1],'DisplayName','dLQR');
    stairs(1:N,x_slc(1,:),'LineWidth',1.25,'Color',[1 0 0],'DisplayName','sequential dLQR');
    xlabel('step',Interpreter='latex');
    ylabel('$x_1$',Interpreter='latex');
    hold off;
nexttile;
    hold on;
    grid minor;
    set(gca,'FontSize',15);
    set(gca,'TickLabelInterpreter','latex');
    stairs(1:N,x_lc(2,:),'LineWidth',1.25,'Color',[0 0 1],'DisplayName','dLQR');
    stairs(1:N,x_slc(2,:),'LineWidth',1.25,'Color',[1 0 0],'DisplayName','sequential dLQR');
    xlabel('step',Interpreter='latex');
    ylabel('$x_2$',Interpreter='latex');
    hold off;

% % cobweb diagrams
% figure;
% tiledlayout(1,2);
% nexttile;
%     hold on;
%     grid minor;
%     set(gca,'FontSize',15);
%     set(gca,'TickLabelInterpreter','latex');
%     stairs(x_lin(1,1:end-1),x_lin(1,2:end),'LineWidth',1.25,'Color',[0 0 0]);
%     plot([0.9*min(x(1,:)), 1.1*max(x(1,:))],[0.9*min(x(1,:)), 1.1*max(x(1,:))],'LineWidth',1.0,'Color',[0 0 0],'HandleVisibility','off');
%     xlabel('$x_1^{i}$',Interpreter='latex');
%     ylabel('$x_1^{i+1}$',Interpreter='latex');
%     axis equal;
%     hold off;
% nexttile;
%     hold on;
%     grid minor;
%     set(gca,'FontSize',15);
%     set(gca,'TickLabelInterpreter','latex');
%     stairs(x_lin(2,1:end-1),x_lin(2,2:end),'LineWidth',1.25,'Color',[0 0 0]);
%     plot([0.9*min(x(2,:)), 1.1*max(x(2,:))],[0.9*min(x(2,:)), 1.1*max(x(2,:))],'LineWidth',1.0,'Color',[0 0 0],'HandleVisibility','off');
%     xlabel('$x_2^{i}$',Interpreter='latex');
%     ylabel('$x_2^{i+1}$',Interpreter='latex');
%     axis equal;
%     hold off;


%% animation

% duffing_animate(x_slc,dt,c,m,k1,k3,100);
% create_gif('animation.gif','animation/',x_lin,dt,c,m,k1,k3,500);


%% energetic analysis

% KE = 0.5*m*x(2,:).^2;
% PE = 0.5*k1*x(1,:).^2 + 0.25*k3*x(1,:).^4;
% L = KE-PE;
% 
% xspace = linspace(-5,5,999);
% yspace = linspace(-5,5,999);
% [Xmesh,Ymesh] = meshgrid(xspace,yspace);
% lagragian = 0.5*m*Ymesh.^2 - (0.5*k1*Xmesh.^2 + 0.25*k3*Xmesh.^4);
% 
% figure;
% hold on;
% contourf(Xmesh,Ymesh,lagragian);
% cb = colorbar;
% colormap parula;
% stairs(x(1,:),x(2,:),'LineWidth',1.25,'Color',[0 0 0],'DisplayName','nonlinear unforced trajectory');
% hold off;


%% Helper functions

function progressBar(currentIteration, totalIterations, barWidth)
    % Create a single progress bar with in-place updates
    % Inputs:
    %   currentIteration  - Current iteration number
    %   totalIterations   - Total number of iterations
    %   barWidth          - Width of the progress bar (optional, default = 50)
    
    % Set default bar width if not provided
    if nargin < 3
        barWidth = 50;
    end
    
    % Persistent variable to track first call
    persistent firstCall
    if isempty(firstCall)
        firstCall = true;
    end
    
    % Calculate progress percentage
    percentComplete = currentIteration / totalIterations;
    
    % Calculate number of filled and empty segments
    filledLength = round(percentComplete * barWidth);
    emptyLength = barWidth - filledLength;
    
    % Create progress bar string
    progressBarString = [repmat('=', 1, filledLength), repmat(' ', 1, emptyLength)];
    
    % Calculate percentage as a string
    percentString = sprintf('%6.2f%%', percentComplete * 100);
    
    % Clear previous line if not first call
    if ~firstCall
        fprintf(repmat('\b', 1, 100));  % Adjust number of backspaces as needed
    else
        firstCall = false;
    end
    
    % Construct and print progress bar output
    outputString = sprintf('[%s] %s (%d/%d)', ...
        progressBarString, percentString, currentIteration, totalIterations);
    
    fprintf('%s', outputString);
    
    % Add final newline when complete
    if currentIteration == totalIterations
        fprintf('\n');
        % Reset persistent variable for potential reuse
        firstCall = true;
    end
end