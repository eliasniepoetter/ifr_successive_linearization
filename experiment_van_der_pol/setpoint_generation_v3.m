function [new_setpoint,valid_return,search_time] = setpoint_generation_v3(x_i,alpha,mu)

syms x1 x2 u real
x_sym = [x1; x2];  % Symbolic state vector
%mu = 0.5;

% Define the symbolic vector field
f_sym = [x2;
    mu*(1 - x1^2)*x2 - x1 + u^2];

% Compute the Jacobians symbolically
grad_f_x_sym = jacobian(f_sym, x_sym);
grad_f_u_sym = jacobian(f_sym, u);

% Now define the fully numeric version of f(x, u)
f = @(x,u) [
    x(2);
    mu*(1 - x(1)^2) * x(2) - x(1) + u;
];

% current state
x_c = x_i;

% Search Algorithm
n_samples = 500;
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

    % Substitute the nominal values into the Jacobians
    grad_f_x_numeric = double(subs(grad_f_x_sym, [x1, x2, u], [x_l' u_l]));
    grad_f_u_numeric = double(subs(grad_f_u_sym, [x1, x2, u], [x_l' u_l]));

    % Define the error function e(x,u)
    e = @(x,u) f(x,u) - grad_f_x_numeric * x - grad_f_u_numeric * u - f(x_l, u_l);
    errors(i) = norm(e(x_c-x_l,0),2);
end
search_time = toc;





%eps = mean(errors)*0.5;
valid_samples = samples(:,errors < (mean(errors)-var(errors)));

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

valid_e = errors(errors<eps);
valid_return = [valid_samples;valid_e];

if norm(x_i) > alpha
    new_setpoint = best_sample;
else
    new_setpoint = [0;0];
end

end

