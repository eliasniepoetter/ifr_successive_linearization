function [A,B,Z] = linearization_scheme(x_l,mu,u_l)

% Definition of symbolic variables
syms x1 x2 u
x = [x1;x2];

% Dynamics
f = @(x,u)[
    x(2);
    mu*(1-x(1)^2)*x(2)-x(1)+u;
];

grad_f_x = jacobian(f(x,u),x);
grad_f_u = jacobian(f(x,u),u);

A = double(subs(grad_f_x,[x;u],[x_l;u_l]));
B = double(subs(grad_f_u,[x;u],[x_l;u_l]));
Z = f(x_l,u_l);

end