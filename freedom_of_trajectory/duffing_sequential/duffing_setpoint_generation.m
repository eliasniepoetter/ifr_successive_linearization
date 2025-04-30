function [x_sp] = duffing_setpoint_generation(x,stepLength)

% compute unit vector to determine next setpoint direction
unitDirection = x / norm(x);

% compute next setpoint
if norm(x) > 1.5*stepLength
    x_sp = x-stepLength*unitDirection;
    x_sp(2) = 0;                        % zero is enforced to an equalibrium, otherwise a second input is needed
else
    x_sp = [0;0];                       % origin if close enough to it
end

end