function [Tab_events,xIN,xSIGMA,xOUT] = manage_event_A2(Tab_events,xIN,lambda,mu,xSIGMA,xOUT)
    global K;
    global C;
    global Losses;
    global Arrivals;
    global Previous_queue_change_instant;
    global State_residence_time_histogram;

    clk = min(Tab_events);
    Arrivals = Arrivals +1;
    State_residence_time_histogram(xSIGMA+1) = ... %because MATLAB counts from 1
        State_residence_time_histogram(xSIGMA+1)+(clk-Previous_queue_change_instant);
    if xSIGMA < K
        xSIGMA = xSIGMA+1;
        for i = 1:C
            if (xOUT(i) == 0) %check if a server is idle
                xOUT(i) = 1;
                Tab_events((2*i)+1) = clk +rand_exp(2*mu); %generate the next S1_i
                break;
            end
        end
    Previous_queue_change_instant = clk;
    else
        Losses = Losses+1;
    end

    xIN = 1;


    Tab_events(1)=clk+rand_exp(2*lambda);
    Tab_events(2)=inf;

end
