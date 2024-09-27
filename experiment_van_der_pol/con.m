function [c,ceq] = con(x,x_i,search_range)
% c = double(norm(x_i,2) - norm(x,2) > -1000);
% ceq = [];
c = double((x(1)-x_i(1))^2 + (x(2)-x_i(2))^2 - search_range^2);
ceq = [];
end