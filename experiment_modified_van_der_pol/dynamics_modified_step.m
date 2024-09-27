function [x_ip1] = dynamics_step(x_i,mu,u,dt)

% Dynamics
f = @(x,u)[
    x(2)+u(1);
    mu*(1-x(1)^2)*x(2)-x(1)+u(2);
];

% Runge-Kutta 4 Integration Scheme
k1 = f(x_i,u);
k2 = f(x_i+dt/2*k1,u);
k3 = f(x_i+dt/2*k2,u);
k4 = f(x_i+dt*k3,u);
x_ip1 = x_i + dt/6*(k1+2*k2+2*k3+k4);

end