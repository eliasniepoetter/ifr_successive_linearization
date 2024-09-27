%% Summary
% Main simulation of the successive linearization framework for the classic
% VdP with naive setpoint selection (no compensation of the setpoint dynamics)



close all;
clear;
clc;


% simulation time and steps
tend = 100;
dt = 0.01;


% allocation of arrays
valid_return_memory = cell(1,100);
system_memory = cell(4,100);
optimization_log_memory = cell(1,100);
u_memory = zeros(1,length(0:dt:tend));
x_memory = zeros(2,length(0:dt:tend));
x_memory_openloop = zeros(2,length(0:dt:tend));
setpoint_memory = zeros(2,length(0:dt:tend));


% initial state & parameters
x0 = [-8;-9];                                                                     % initial state
x_memory(:,1) = x0;
setpoint0 = [0;0];                                                              % initial setpoint
setpoint_memory(:,1) = setpoint0;
x_memory_openloop(:,1) = x0;
mu = 0.5;                                                                       % Van der Pol parameter
alpha = 0.2;                                                                    % threshold for origin as setpoint


% Simulation of the inverted Van der Pol
wb = waitbar(0,'simulation in progress');
counter = 1;
for n = 1 : length(0:dt:tend)-1
    x_i = x_memory(:,n);                                                        % allocate current state
    setpoint = setpoint_memory(:,n);                                            % allocate current setpoint
    x_i_openloop = x_memory_openloop(:,n);                                      % allocate openloop state

    %%% Start of Successive LQR - Input: x,setpoint %%%

    % Compute setpoint + corresponding linear systems + controller
    if norm(setpoint-x_i) < (0.01)
        [setpoint,optimization_log] = setpoint_generation_v4(x_i,alpha,mu);          % compute setpoint
        setpoint_input = 0;
        [A,B,Z] = linearization_scheme(setpoint,mu,setpoint_input);             % compute linearized system for setpoint
                % replacing setpoint with x_i???
        K = controller_synthesis(A,B);                                          % compute LQR controller
        system_memory{1,counter} = A;
        system_memory{2,counter} = B;
        system_memory{3,counter} = Z;
        system_memory{4,counter} = K;
        optimization_log_memory{counter} = optimization_log;
        %valid_return_memory{counter} = valid_return;
        counter = counter + 1;
    elseif n == 1                                                               % ensuring new setpoint in iter 1
        [setpoint,optimization_log] = setpoint_generation_v4(x_i,alpha,mu);
        setpoint_input = 0;
        [A,B,Z] = linearization_scheme(setpoint,mu,setpoint_input);
        K = controller_synthesis(A,B);
        system_memory{1,counter} = A;
        system_memory{2,counter} = B;
        system_memory{3,counter} = Z;
        system_memory{4,counter} = K;
        optimization_log_memory{counter} = optimization_log;
        %valid_return_memory{counter} = valid_return;
        counter = counter + 1;
    end

    % Compute control input -> state is tracking error
    %u0 = B\(-Z-A*setpoint);
    u = -K*(x_i-setpoint);                                                      % tracking formulation to enable LQR
    u_memory(n) = u;

    %%% End of Successive LQR - Output: u %%%
    
    % Update Dynamics -> true state = tracking error + setpoint
    x_ip1_star = dynamics_step((x_i-setpoint),mu,u,dt);                         % call integration scheme on tracking error
    x_ip1 = x_ip1_star + setpoint;                                              % update true state 
            % x_ip1_star = dynamics_step((x_i),mu,u,dt)???
            % x_ip1 = x_ip1_star???
    
    % Update openloop states and memories
    x_memory(:,n+1) = x_ip1;
    setpoint_memory(:,n+1) = setpoint;
    x_ip1_openloop = dynamics_step(x_i_openloop,mu,0,dt);
    x_memory_openloop(:,n+1) = x_ip1_openloop;

    % stopping criteria
    if norm(setpoint-x_i) > 1000
        break;
    elseif norm(x_i) < 0.01
        break
    end

    % update waitbar
    waitbar(n/(length(0:dt:tend)-1));
end
close(wb);


%% Visualization
 
[x1, x2] = meshgrid(-2:0.1:2, -2:0.1:2);
dx1 = x2;
dx2 = mu*(1 - x1.^2).*x2 - x1;

% state space trajectory
figure;
set(gca,'fontsize', 14);
hold on;
%title(['Van der Pol Oscillator with $\mu=$',num2str(mu)],'Interpreter','latex');
scatter(x0(1),x0(2),70,'filled','DisplayName','Initial','Marker','square','MarkerFaceColor',[0,0,0]);
scatter(setpoint_memory(1,1:n),setpoint_memory(2,1:n),'filled','DisplayName','Setpoints','MarkerFaceColor',[0,0,0]);
scatter(0,0,80,'filled','DisplayName','Target','Marker','diamond','MarkerFaceColor',[0,0,0]);
plot(x_memory(1,1:n),x_memory(2,1:n),'DisplayName','Trajectory','LineStyle','-','Color',[0 0 0],'LineWidth',1.25);
plot(x_memory_openloop(1,1:n),x_memory_openloop(2,1:n),'DisplayName','Openloop','LineStyle','-.','Color',[0 0 0],'LineWidth',1.25);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
axis equal;
legend('Location','southwest','Interpreter','latex');
xlim([-5 5]);
ylim([-5 5]);
grid on;
hold off;

% plots over time
figure;
tl = tiledlayout(3,1);
%title(tl,['Van der Pol Oscillator with $\mu=$',num2str(mu)],'Interpreter','latex','fontsize',16);
nexttile;
set(gca,'fontsize', 14);
    hold on;
    grid on;
    time_u = 0:dt:(length(u_memory)*dt)-dt;
    plot(time_u(1:n),u_memory(1:n),'Color','black','LineWidth',1.25);
    %title('Input $u$','Interpreter','latex');
    xlabel('Time [s]','Interpreter','latex');
    ylabel('$u$','Interpreter','latex');
    hold off;
nexttile;
set(gca,'fontsize', 14);
    hold on;
    grid on;
    time = 0:dt:(length(u_memory)*dt);
    plot(time(1:n),x_memory(1,1:n),'Color','black','LineWidth',1.25);
    %title('State $x_1$','Interpreter','latex');
    xlabel('Time [s]','Interpreter','latex');
    ylabel('$x_1$','Interpreter','latex');
    hold off;
nexttile;
set(gca,'fontsize', 14);
    hold on;
    grid on;
    plot(time(1:n),x_memory(2,1:n),'Color','black','LineWidth',1.25);
    %title('State $x_2$','Interpreter','latex');
    xlabel('Time [s]','Interpreter','latex');
    ylabel('$x_2$','Interpreter','latex');
    hold off;

%%

optimization_log_memory = optimization_log_memory(~cellfun('isempty',optimization_log_memory));

figure;
tiledlayout(2,ceil(length(optimization_log_memory)/2));
for i = 1 : length(optimization_log_memory)
    log = optimization_log_memory{i};
    cost_function_sym = log.cost_function_sym;
    nexttile;
    hold on;
    grid on;
    %fsurf(cost_function_sym,[-10 10],'ShowContours','on');
    fcontour(cost_function_sym,[-10 10],'Fill','off','LevelStep',5);
    colorbar;
    scatter3(log.x_eval(1),log.x_eval(2),log.cost_function(log.x_eval),50,'red','filled');
    scatter3(log.x_min(1),log.x_min(2),log.fval,50,'green','filled');
    legend('obective contour','initial','minimum');
end








