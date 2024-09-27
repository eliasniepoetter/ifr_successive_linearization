function [new_setpoint,valid_return,search_time] = setpoint_generation_v2(x_i,alpha,mu)

% Definition of symbolic variables and parameters
syms x1 x2 u
x = [x1;x2];

% Dynamics of the Van der Pol
f = @(x,u)[
    x(2);
    mu*(1-x(1)^2)*x(2)-x(1)+u;
];

% Jacobians of f(x,u)
grad_f_x = jacobian(f(x,u),x);
grad_f_u = jacobian(f(x,u),u);

% current state
x_c = x_i;

% Search Algorithm
n_samples = 100;
samples = zeros(2,n_samples);
sigma1 = 0.5;
sigma2 = 0.5;
rng(0,'twister');
samples(1,:) = normrnd(x_c(1),sigma1,[1,n_samples]);
samples(2,:) = normrnd(x_c(2),sigma2,[1,n_samples]);

errors = zeros(1,n_samples);

tic;
for i = 1 : n_samples
    x_l = samples(:,i);
    u_l = 0;
    e = @(x,u) f(x,u) - double(subs(grad_f_x,[x;u],[x_l;u_l]))*x - double(subs(grad_f_u,[x;u],[x_l;u_l]))*u - f(x_l,u_l);
    
    % Jacobian of e(x,u) w.r.t x
    grad_e_x = jacobian(e(x,u),x);
    
    % Compute induced 2-Norm of de/dx
    E = grad_e_x*grad_e_x';
    error_lipschitz = sqrt(max(eig(E)));
    errors(i) = double(subs(error_lipschitz,[x;u],[x_c;0]));
end
search_time = toc;

eps = 1;
valid_samples = samples(:,errors<eps);

dist = 1e6;
best_sample = [1000;1000];
x_t = [0;0];
for i = 1 : length(valid_samples)
    dist_new = norm(valid_samples(:,i)-x_t,2);
    if dist_new < dist
        dist = dist_new;
        best_sample = valid_samples(:,i);
    end
end

valid_e = errors(errors<eps);
valid_return = [valid_samples;valid_e];

if norm(x_i) > alpha
    new_setpoint = best_sample;
else
    new_setpoint = [0;0];
end

end

