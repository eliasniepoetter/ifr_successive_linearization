function [A,B,Z] = duffing_linearization(x_l,u_l,dt,c,m,k1,k3)

% Definition of symbolic variables
syms x1 x2 u
x = [x1;x2];

% Dynamics of the Duffing Equation
F_duffing = @(x,u,dt,c,m,k1,k3) [
x(1) + dt*x(2);
x(2) + dt*(-(c/m)*x(2) - (k1/m)*x(1) - (k3/m)*x(1)^3 + u);
];

% Gradients of the dynamics
grad_F_x = jacobian(F_duffing(x,u,dt,c,m,k1,k3),x);
grad_F_u = jacobian(F_duffing(x,u,dt,c,m,k1,k3),u);

% linear matrices
A = double(subs(grad_F_x,[x;u],[x_l;u_l]));
B = double(subs(grad_F_u,[x;u],[x_l;u_l]));
Z = F_duffing(x_l,u_l,dt,c,m,k1,k3);

% Compute Lipschitz constant on constrains set
Phi = @(x,u) F_duffing(x,u,dt,c,m,k1,k3) - double(subs(grad_F_x,[x;u],[x_l;u_l]))*(x-x_l) - double(subs(grad_F_u,[x;u],[x_l;u_l]))*(u-u_l) - F_duffing(x_l,u_l,dt,c,m,k1,k3);
grad_Phi_x = jacobian(Phi(x,u),x);
PHI = grad_Phi_x*grad_Phi_x';
Phi_lip_fun = sqrt(max(eig(PHI)));

end