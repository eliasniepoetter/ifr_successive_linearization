%% Summary
% Examplary Optimization Routine for the modified VdP
% The cost function is the L2 norm of the linearization error dynamics
% Additionally a nonlinear constraint is added to enforce an equlibrium at the optmimized linearization point

clear;
close all;
clc;

syms x1_l x2_l u1_l u2_l x1 x2 u1 u2
x_l_sym = [x1_l;x2_l];
u_l_sym = [u1_l;u2_l];
x_sym = [x1;x2];
u_sym = [u1;u2];

x_eval = [2;-1];
u_eval = [1;-1];
mu = 1;

f_sym = [x2+u1;
    mu*(1 - x1^2)*x2 - x1 + u2];

grad_f_x_sym = jacobian(f_sym, x_sym);
grad_f_u_sym = jacobian(f_sym, u_sym);

f_linearized_sym = subs(grad_f_x_sym,x_sym,x_l_sym)*(x_sym-x_l_sym) + subs(grad_f_u_sym,u_sym,u_l_sym)*(u_sym-u_l_sym);

e_sym = f_sym - f_linearized_sym;
e_sym = subs(e_sym,x_sym,x_eval);
cost_function_sym = norm(e_sym,2) + norm(x_l_sym,2) + norm(u_l_sym,2);

% x is denoting the optmimzation variables --> [x1_l; x2_l; u1_l; u2_l]
cost_function = @(x) double(subs(cost_function_sym,[x1_l;x2_l;u1_l;u2_l],x));

options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point');
[xopt,fval,exitflag] = fmincon(cost_function,[x_eval-[1;1];u_eval-[0.5;0.5]],[],[],[],[],[],[],@(x)nonlcon(x,mu),options);


grid_val = 3;
grid_step = 0.25;

x1_grid = -grid_val:grid_step:grid_val;
x2_grid = -grid_val:grid_step:grid_val;
x3_grid = -grid_val:grid_step:grid_val;
x4_grid = -grid_val:grid_step:grid_val;

[X1, X2, X3, X4] = ndgrid(x1_grid, x2_grid, x3_grid, x4_grid);
n = length(x1_grid);
costs = zeros(n,n,n,n);
wb = waitbar(0,'grid creation in progress');
counter = 0;
for i = 1:n
    for j = 1:n
        for k = 1:n
            for l = 1:n
                costs(i,j,k,l) = cost_function([x1_grid(i);x2_grid(j);x3_grid(k);x4_grid(l)]);
                counter = counter + 1;
                waitbar(counter/n^4);
            end
        end
    end
end
close(wb);

%%



c = costs(:,:,5,12);
[X1_2D, X2_2D] = meshgrid(x1_grid, x2_grid);
[X1_res, X2_res] = meshgrid(-3:0.01:3, -3:0.01:3);
vq = griddata(X1_2D(:),X2_2D(:),c(:),X1_res, X2_res,'v4');

%%

u1_select = xopt(3);
u2_select = xopt(4);
figure;
hold on;
grid on;
% F_slice = costs(:,:,2, 1);
% surf(X1_2D, X2_2D, F_slice);
contour(X1_res,X2_res,vq,'LevelStep',1);
fcontour(@(xx,yy)yy+u1_select,[-3 3 -3 3],'LevelList',0,'LineWidth',1.25,'LineColor','black');
fcontour(@(xx,yy)u2_select - xx - yy.*(xx.^2 - 1),[-3 3 -3 3],'LevelList',0,'LineWidth',1.25,'LineColor','black');
scatter(xopt(1),xopt(2),80,'black','filled','diamond');
scatter(x_eval(1),x_eval(2),80,'black','filled','o');
axis equal;





















