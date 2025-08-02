function [xSIGMA,xSERVER1,xSERVER2,tab_event] = manage_event_A(xSIGMA,xSERVER1,xSERVER2,tab_event,lambda,mu1,mu2)
global K
global Losses
global Arrivals
global C
global Previous_queue_change_instant
global State_residence_time_histogram

clk = min(tab_event);
Arrivals = Arrivals + 1;
State_residence_time_histogram(xSIGMA+1) = ...
    State_residence_time_histogram(xSIGMA+1)+(clk-Previous_queue_change_instant);
%check if queue is full
if xSIGMA < K
    xSIGMA = xSIGMA +1;
else
    Losses = Losses +1;
end

if xSERVER1 == 0
    xSERVER1 = 1;
    %service time for first server
    tab_event(2) = clk+rand_exp(mu1);
elseif xSERVER2 == 0
    xSERVER2 = 1;
    %service time for second server
    tab_event(3) = clk+rand_exp(mu2);
end
%anyway next arrival
    tab_event(1) = clk + rand_exp(lambda);
    Previous_queue_change_instant = clk;

end
