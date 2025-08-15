clear; clc;

lambda_arrival = 4;     % Arrival rate
lambda_service = 100;   % Service rate
num_events = 5000;      % Number of customers

% Generate event times (first customer arrives at time 0)
interarrival_times = [0; exprnd(1/lambda_arrival, num_events-1, 1)];
work_times = exprnd(1/lambda_service, num_events, 1);
arrival_times = cumsum(interarrival_times);

fprintf('First 10 Arrival Times: ');
fprintf('%.4f ', arrival_times(1:10));
fprintf('\n');

served_times = zeros(num_events, 1);
finish_times = zeros(num_events, 1);

busy_until = -1.0;
for i = 1:num_events
    served_time = max(arrival_times(i), busy_until);
    busy_until = served_time + work_times(i);
    served_times(i) = served_time;
    finish_times(i) = busy_until;
end

fprintf('First 10 Service Start Times: ');
fprintf('%.4f ', served_times(1:10));
fprintf('\n');

fprintf('First 10 Service Finish Times: ');
fprintf('%.4f ', finish_times(1:10));
fprintf('\n');

fprintf('System busy until time: %.4f\n', busy_until);

current_time = 0;
time_step = 0.1;
num_samples = floor(busy_until / time_step);

fprintf('Taking %d samples with time step %.1f\n', num_samples, time_step);

num_customer_serving = zeros(num_samples, 1);
num_customer_waiting = zeros(num_samples, 1);
num_customer_sys = zeros(num_samples, 1);
time_points = zeros(num_samples, 1);

for t = 1:num_samples
    current_time = (t-1) * time_step;
    time_points(t) = current_time;
    serving = sum((served_times <= current_time) & ...
                  (current_time < finish_times));

    waiting = sum((arrival_times <= current_time) & ...
                  (current_time < served_times));

    num_customer_serving(t) = serving;
    num_customer_waiting(t) = waiting;
    num_customer_sys(t) = serving + waiting;
end

avg_in_system = mean(num_customer_sys);
avg_waiting = mean(num_customer_waiting);
avg_serving = mean(num_customer_serving);

fprintf('\n--- System Statistics ---\n');
fprintf('Average customers in system: %.4f\n', avg_in_system);
fprintf('Average customers waiting: %.4f\n', avg_waiting);
fprintf('Average customers being served: %.4f\n', avg_serving);

rho = lambda_arrival / lambda_service;
theoretical_L = rho / (1 - rho);
theoretical_Lq = rho^2 / (1 - rho);
theoretical_server_util = rho;

fprintf('\n--- Theoretical M/M/1 Values ---\n');
fprintf('Utilization (ρ): %.4f\n', rho);
fprintf('Expected customers in system (L): %.4f\n', theoretical_L);
fprintf('Expected customers waiting (Lq): %.4f\n', theoretical_Lq);
fprintf('Expected server utilization: %.4f\n', theoretical_server_util);

figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(time_points, num_customer_sys, 'b-', 'LineWidth', 1);
title('Total Customers in System');
xlabel('Time');
ylabel('Number of Customers');
grid on;

subplot(2,2,2);
plot(time_points, num_customer_waiting, 'r-', 'LineWidth', 1);
title('Customers Waiting in Queue');
xlabel('Time');
ylabel('Number of Customers');
grid on;

subplot(2,2,3);
plot(time_points, num_customer_serving, 'g-', 'LineWidth', 1);
title('Customers Being Served');
xlabel('Time');
ylabel('Number of Customers');
grid on;

subplot(2,2,4);
plot(time_points, num_customer_sys, 'b-', 'LineWidth', 1); hold on;
plot(time_points, num_customer_waiting, 'r-', 'LineWidth', 1);
plot(time_points, num_customer_serving, 'g-', 'LineWidth', 1);
yline(theoretical_L, 'b--', 'LineWidth', 2, 'DisplayName', 'Theoretical L');
yline(theoretical_Lq, 'r--', 'LineWidth', 2, 'DisplayName', 'Theoretical Lq');
title('System Overview');
xlabel('Time');
ylabel('Number of Customers');
legend('Total in System', 'Waiting', 'Being Served', 'Theoretical L', 'Theoretical Lq');
grid on;

sgtitle('M/M/1 Queue Simulation Results');