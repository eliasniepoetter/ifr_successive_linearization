function duffing_animate(x,dt,c,m,k1,k3,skip)
    
% allocate data
positions = x(1,1:skip:end);
velocities = x(2,1:skip:end);
time = 0:dt:length(x);
t = time(1:skip:end);
zeta_lin = c / (2 * sqrt(m*k1));
zeta_nlin = c / (2 * sqrt(m*k3));

% Visualization setup
figure('Color', 'white', 'Position', [100, 100, 1000, 600]);
x_lim = [min(positions)-0.5, max(positions)+0.5];
y_lim = [-1, 1.5];
ground_y = 0;
pause(2);

% animation loop
for frame = 1:length(positions)

    % clear figure to start with blank figure
    clf;
    hold on;

    % Plot ground line
    plot(x_lim, [ground_y, ground_y], 'k-', 'LineWidth', 2);
    
    % Draw spring (sinusoidal representation)
    curr_pos = positions(frame);
    spring_x = linspace(0, curr_pos, 150);
    spring_y = 0.1 * sin(linspace(0, 4*pi, 150)) + ground_y;
    plot(spring_x, spring_y, 'b-', 'LineWidth', 1.5);
    
    % Draw mass
    plot(curr_pos, ground_y, 'ro', 'MarkerSize', 30, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    
    % figure settings and information
    xlim(x_lim);
    ylim(y_lim);
    title(sprintf('Duffing Oscillator\nm = %.1f kg, k1 = %.1f N/m, k3 = %.1f N/m, c = %.1f Ns/m\nTime: %.3f s', ...
        m, k1, k3, c, t(frame)), 'FontSize', 12);
    xlabel('Displacement (m)');
    ylabel('System Configuration');
    info_text = sprintf(...
        'Position: %.3f m\nVelocity: %.3f m/s\nLiner Damping Ratio: %.3f\nNonlinear Damping Ratio: %.3f', ...
        positions(frame), velocities(frame), zeta_lin, zeta_nlin);
    text(x_lim(1), y_lim(2), info_text, ...
        'VerticalAlignment', 'top', 'FontSize', 10);
    drawnow;
    pause(0);
end
end