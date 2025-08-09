function [Tab_events,xIN] = manage_event_A1(Tab_events,xIN,lambda)
    clk = min(Tab_events);
    xIN = 2;
    
    Tab_events(1) = inf;
    Tab_events(2) = clk +rand_exp(2*lambda);
end