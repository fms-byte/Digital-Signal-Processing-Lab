% Sinc Function

x = -5:0.1:5;

% Compute sinc(x) while handling division by zero at x = 0
y = sin(x);
z = ones(size(x));  % Initialize z with ones to handle division by zero
idx = find(x ~= 0); % Find indices where x is not equal to 0
z(idx) = y(idx) ./ x(idx); % Compute sinc(x) for non-zero x

% Plotting
subplot(2,1,1);
plot(x, y, 'b');
grid on;
xlabel("x");
ylabel("y = sin(x)");

subplot(2,1,2);
plot(x, z, 'r');
grid on;
xlabel("x");
ylabel("z = sinc(x)");