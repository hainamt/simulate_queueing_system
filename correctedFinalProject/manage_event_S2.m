function [Tab_events,xSIGMA,xOUT] = manage_event_S2(Tab_events,xSIGMA,xOUT,mu)
    global C;
    global State_residence_time_histogram;
    global Previous_queue_change_instant;

    [clk, i] = min(Tab_events);
    State_residence_time_histogram(xSIGMA+1) = ... %because MATLAB counts from 1
        State_residence_time_histogram(xSIGMA+1)+(clk-Previous_queue_change_instant);

    xSIGMA = xSIGMA -1;

    if xSIGMA <=(C-1)
        xOUT((i/2)-1) = 0;
        Tab_events(i-1) = inf;
    else
        xOUT((i/2)-1) = 1;
        Tab_events(i-1) = clk + rand_exp(2*mu);
    end
    Previous_queue_change_instant = clk;
    Tab_events(i) = inf;
end