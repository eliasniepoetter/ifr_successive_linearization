function [c,ceq] = nonlcon(x,mu)

    % Variable:
    % x = [x1_l; x2_l; u1_l; u2_l]

    % inequality constraints
    c = [];

    % equality constraints - equlibrium constraints
    ceq = [x(2)+x(3);
        mu*(1 - x(1)^2)*x(2) - x(1) + x(4)];

    % guarantee type double for fmincon
    ceq = double(ceq);
    
end
