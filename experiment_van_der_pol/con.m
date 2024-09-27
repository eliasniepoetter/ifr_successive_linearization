function [c,ceq] = con(x,x_i,search_range)

% inequality constraints: searchable set if a disk with radius search_range
% around last setpoint x_i
c = double((x(1)-x_i(1))^2 + (x(2)-x_i(2))^2 - search_range^2);

% equality constraints
ceq = [];
end