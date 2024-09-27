function [x_r,valid_return] = setpoint_generation_v1(x_i,alpha)

% unitDirection = x_i / norm(x_i);
% if norm(x_i) > alpha
%     x_r = x_i-alpha*unitDirection;
% else
%     x_r = [0;0];
% end

% Definition of symbolic variables and parameters
syms x1 x2 u
x = [x1;x2];
mu = -0.75;

% Dynamics of the Van der Pol
f = @(x,u)[
    x(2);
    mu*(1-x(1)^2)*x(2)-x(1)+u;
];

% Jacobians of f(x,u)
grad_f_x = jacobian(f(x,u),x);
% grad_f_u = jacobian(f(x,u),u);

% linearization error dynamics
x_l = x_i;
u_l = 0;
e = @(x,u) f(x,u) - double(subs(grad_f_x,[x;u],[x_l;u_l]))*x - f(x_l,u_l);

% Jacobian of e(x,u) w.r.t x
grad_e_x = jacobian(e(x,u),x);

% Compute induced 2-Norm of de/dx
E = grad_e_x*grad_e_x';
error_lipschitz = sqrt(max(eig(E)));

% Search Algorithm
n_samples = 1000;
samples = zeros(2,n_samples);
sigma1 = 0.5;
sigma2 = 0.5;
samples(1,:) = normrnd(x_l(1),sigma1,[1,n_samples]);
samples(2,:) = normrnd(x_l(2),sigma2,[1,n_samples]);

e = zeros(1,n_samples);
for i = 1 : n_samples
    e(i) = double(subs(error_lipschitz,x,samples(:,i)));
end

eps = 10;
valid_samples = samples(:,e<eps);

dist = 1e6;
best_sample = [1000;1000];
x_t = [0;0];
for i = 1 : size(valid_samples,2)
    dist_new = norm(valid_samples(:,i)-x_t,2);
    if dist_new < dist
        dist = dist_new;
        best_sample = valid_samples(:,i);
    end
end

valid_e = e(e<eps);
valid_return = [valid_samples;valid_e];

if norm(x_i) > alpha
    x_r = best_sample;
else
    x_r = [0;0];
end


end

