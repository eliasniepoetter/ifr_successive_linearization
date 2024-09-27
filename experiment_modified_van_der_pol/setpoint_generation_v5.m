function [new_setpoint,optimization_log] = setpoint_generation_v4(x_i,alpha,mu)

syms x1_l_sym x2_l_sym u_l_sym x1_sym x2_sym u_sym
x_l_sym = [x1_l_sym; x2_l_sym];
x_sym = [x1_sym; x2_sym];

x_eval = x_i;
u_eval = 0;

f_sym = [x2_sym;
    mu*(1 - x1_sym^2)*x2_sym - x1_sym + u_sym];

f_l_sym = [x2_l_sym;
    mu*(1 - x1_l_sym^2)*x2_l_sym - x1_l_sym + u_l_sym];

grad_f_x_sym = jacobian(f_sym, x_sym);
grad_f_u_sym = jacobian(f_sym, u_sym);

f_linearized_sym = subs(grad_f_x_sym,x_sym,x_l_sym)*x_sym + subs(grad_f_u_sym,u_sym,u_l_sym)*u_sym + f_l_sym;

e_sym = f_sym - f_linearized_sym;
e_sym = subs(e_sym,[x_sym;u_sym;u_l_sym],[x_eval;u_eval;0]);
e_numeric_norm = @(x) norm(double(subs(e_sym,x_l_sym,x)),2);
cost_function = @(x) e_numeric_norm(x) + 1e-1*norm(x,2);
cost_function_sym = norm(e_sym,2) + norm(x_l_sym,2);

options = optimoptions('fmincon','Display','off','Algorithm','interior-point','OptimalityTolerance',1e-6);
search_range = 1;

%@(x)con(x,x_i)

tic;
[optimized_setpoint,fval,exitflag] = fmincon(cost_function,x_i,[],[],[],[], ...
    [], ...
    [],@(x)con(x,x_i,search_range),options);
t_search = toc;

optimization_log = struct;
optimization_log.fval = fval;
optimization_log.exitflag = exitflag;
optimization_log.t_search = t_search;
optimization_log.cost_function = cost_function;
optimization_log.cost_function_sym = cost_function_sym;
optimization_log.e_sym = e_sym;
optimization_log.x_min = optimized_setpoint;
optimization_log.x_eval = x_eval;


if norm(x_i) < alpha
    new_setpoint = [0;0];
elseif norm(optimized_setpoint,2) < alpha
    new_setpoint = [0;0];
else
    new_setpoint = optimized_setpoint;
end







end