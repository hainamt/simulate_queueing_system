% clear
% close all;
function average_waiting_time = run_system
duration_simulation = 50000;
% rng(9);

lambda = 1/10;
mu1 = 1/20;
mu2 = mu1;

global C
C = 2;
rho = lambda/C*mu1
global K
global Arrivals
global Losses
Arrivals = 0;
Losses = 0;
K = 80;
global Previous_queue_change_instant
global State_residence_time_histogram
Previous_queue_change_instant = 0;
State_residence_time_histogram = zeros(1,K+1);

%% Initialization

xSIGMA = 0; %# of users in the system

xSERVER1 = 0; %free = 0, busy = 1 
xSERVER2 = 0;

tab_event =[
    1E20 %arrival
    1E20 %server1 finish
    1E20 %server2 finish
    ];


%% First event
clk = 0;
count_events = 1;
tab_event(1) = clk + rand_exp(lambda);

%these should be used for plotting
abscissa_t(count_events) = clk;
vect_event(count_events) = 0;
elem_in_syst(count_events) = xSIGMA;

%% Simulation loop
while clk < duration_simulation
    [clk_next_event, type_next_event] = min(tab_event);
    count_events = count_events +1;
    switch type_next_event
        case 1
            [xSIGMA,xSERVER1,xSERVER2,tab_event] = manage_event_A(xSIGMA,xSERVER1,xSERVER2,tab_event,lambda,mu1,mu2);
        case 2
            [xSIGMA,xSERVER1,xSERVER2,tab_event] = manage_event_S1(xSIGMA,xSERVER1,xSERVER2,tab_event,lambda,mu1,mu2);
        case 3
            [xSIGMA,xSERVER1,xSERVER2,tab_event] = manage_event_S2(xSIGMA,xSERVER1,xSERVER2,tab_event,lambda,mu1,mu2);
    end
    clk = clk_next_event;
    abscissa_t(count_events) = clk;
    vect_event(count_events) = type_next_event;
    elem_in_syst(count_events) = xSIGMA;
end

Loss_prob = Losses/Arrivals
% 
% subplot(3,1,1)
% %better for discrete qnt rather than plot
% stairs(abscissa_t,vect_event)
% yticks([0 1 2 3 4]);
% xlim([0 duration_simulation]);
% title('Event')
% 
% % figure(2)
% subplot(3,1,2)
% stairs(abscissa_t,elem_in_syst)
% yticks([0 max(elem_in_syst)])
% xlim([0 duration_simulation])
% title('Elements in the system')

state_residence_time_pdf = State_residence_time_histogram/sum(State_residence_time_histogram);

% subplot(3,1,3)
% bar(0:K,state_residence_time_pdf);
% title("State residence time pdf");
% %this can be done also by using .* (Adamard product) but we need to sum all
% %the value returned by the product, so that
% %averageNumUsers = sum(0:K).*StateResidenceTimePDF;
average_number_of_users = sum((0:K).*state_residence_time_pdf);

average_waiting_time = average_number_of_users/lambda;

% sys_point2(average_waiting_time,state_residence_time_pdf)
end