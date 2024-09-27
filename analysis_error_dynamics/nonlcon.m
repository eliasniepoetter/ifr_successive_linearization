function [c,ceq] = nonlcon(x,mu)
    c = [];
    ceq = [x(2)+x(3);
        mu*(1 - x(1)^2)*x(2) - x(1) + x(4)];
    ceq = double(ceq);
end
