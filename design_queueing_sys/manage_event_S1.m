function [xSIGMA,xSERVER1,xSERVER2,tab_event] = manage_event_S1(xSIGMA,xSERVER1,xSERVER2,tab_event,lambda,mu1,mu2)
global Previous_queue_change_instant
global State_residence_time_histogram

clk =min(tab_event);

State_residence_time_histogram(xSIGMA+1) = ...
    State_residence_time_histogram(xSIGMA+1)+(clk-Previous_queue_change_instant);
xSIGMA = xSIGMA-1;
if xSIGMA <= 1 
    xSERVER1 = 0;
    tab_event(2) = 1E20;
else 
    xSERVER1 = 1;
    tab_event(2) = clk +rand_exp(mu1);
end
Previous_queue_change_instant = clk;
end