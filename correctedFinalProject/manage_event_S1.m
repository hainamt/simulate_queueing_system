function [Tab_events,xOUT] = manage_event_S1(Tab_events,xOUT,mu)
    [clk,i] = min(Tab_events);

    xOUT((i-1)/2) = 2;

    Tab_events(i) = inf;
    Tab_events(i+1) = clk +rand_exp(2*mu);
end