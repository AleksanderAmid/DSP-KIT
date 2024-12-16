% MATLAB code to analyze the effect of different \mu values on stability for given formula
% Formula: \hat{h}(n) = \hat{h}(n-1) + 2 * \mu * y(n) * (x(n) - \hat{h}^T(n-1) * y(n))

% Parameters
n_max = 100; % Number of iterations
h_dim = 1;   % Dimension of h vector
my_values = [0.01, 0.1, 1]; % Different \mu values
x_n = 1;    % Constant x(n)
y_n = 1;    % Constant y(n)

% Initialize figure
figure;

% Iterate over each \mu value
for i = 1:length(my_values)
    my = my_values(i);
    h = zeros(h_dim, 1); % Initialize h vector
    h_values = zeros(h_dim, n_max);

    % Simulation loop
    for n = 1:n_max
        % Calculate error term e(n)
        e_n = x_n - h' * y_n;
        
        % Update h(n)
        h = h + 2 * my * y_n * e_n;
        h_values(:, n) = h;
        
        % Plot updated h(n)
        subplot(3, 1, i);
        plot(1:n, h_values(1, 1:n), 'bo-', 'LineWidth', 1.5);
        title(['\mu = ', num2str(my)]);
        xlabel('n');
        ylabel('h(n)');
        grid on;
        pause(0.05); % Refresh plot in real-time
    end
end
