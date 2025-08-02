clear ; close ;
tic;
global C;
global K;
global Losses;
global Arrivals
global Previous_queue_change_instant;
global State_residence_time_histogram;

simulation_duration = 50000;
rng(6);
lambda = 1/10;
rho = 1.05;
simulations_number = 100;
max_queue_length = 30;
c_values = [2,3];
k_values = c_values(1):max_queue_length;

loss_probabilities = zeros(length(c_values),length(k_values));
average_number_of_users = zeros(length(c_values),length(k_values));
point_2_matrix=zeros(length(k_values),simulations_number);
% with this loop we evaluate the possibilities for the # of servers
for c_idx = 1:length(c_values)
    %service rate changes as the # of servers
    C = c_values(c_idx);
    mu = lambda/(C*rho);

    for k_idx = 1:length(k_values)
        K = k_values(k_idx);
        if K < C
            average_number_of_users(c_idx,k_idx) = NaN;
            loss_probabilities(c_idx,k_idx) = NaN;
            continue; % K must be at least C
        end
        loss_probabilities_sim = zeros(1,simulations_number);
        average_number_of_users_simulation =zeros(1,simulations_number);


        for sim_idx = 1:simulations_number
            %initialization before each simulation
            Tab_events= inf*ones((C+1)*2,1);
            clk = 0;
            count_events = 1;
            State_residence_time_histogram = zeros(1,K+1);
            Previous_queue_change_instant = 0;
            Arrivals = 0;
            Losses = 0;
            xIN = 1; %input E2
            xOUT = zeros(1,C);
            xSIGMA = 0;

            %first arrival
            Tab_events(1) = clk + rand_exp(2*lambda);
            elem_in_syst(count_events)=xSIGMA;

            % simulation loop
            while clk < simulation_duration
                [t_min_event, type_event] = min(Tab_events);
                count_events = count_events+1;
                switch type_event
                    case 1
                        [Tab_events,xIN] = manage_event_A1(Tab_events,xIN,lambda);

                    case 2
                        [Tab_events,xIN,xSIGMA,xOUT] = manage_event_A2(Tab_events,xIN,lambda,mu,xSIGMA,xOUT);

                    otherwise

                    if(type_event >2 && mod(type_event,2)==1) % xout=1
                        [Tab_events,xOUT] = manage_event_S1(Tab_events,xOUT,mu);
                    else % xout=2
                        [Tab_events,xSIGMA,xOUT] = manage_event_S2(Tab_events,xSIGMA,xOUT,mu);
                    end
                end
                clk = t_min_event; %forwarding to next event
                elem_in_syst(count_events)=xSIGMA;
            end


            loss_probabilities_sim(sim_idx) = Losses / Arrivals; %vector of all simulation losses with same K
            state_residence_time_pdf = State_residence_time_histogram/sum(State_residence_time_histogram);
            average_number_of_users_simulation(sim_idx) = (0:K)*state_residence_time_pdf';
            if C == 2
                point_2_matrix(k_idx,sim_idx) = average_number_of_users_simulation(sim_idx);
            end
        end
        loss_probabilities(c_idx,k_idx) = sum(loss_probabilities_sim)/simulations_number; %vector with same C
        average_number_of_users(c_idx,k_idx) = mean(average_number_of_users_simulation);
        state_residence_time_pdf = [];
    end

end
%%
%subplot(2,1,1);
figure()
hold on;
title("Average number of users");
plot(k_values,average_number_of_users(1,:),LineWidth=2);
plot(k_values,average_number_of_users(2,:),LineWidth=2);
xlabel("K");
xticks(k_values);
ylabel("Number of users");
legend(["C=2", "C=3"])

grid on;
hold off;

%subplot(2,1,2);
figure()
hold on;
title("Loss probability");
plot(k_values,loss_probabilities(1,:),LineWidth=2);
grid on;
plot(k_values,loss_probabilities(2,:),LineWidth=2);
legend(["C=2", "C=3"])
xticks(k_values);
xlabel("K");
ylabel("Loss probability");
hold off;
%
%figure();
%boxchart(point_2_matrix');
%title("Average number of users CI");
%xlabel("K");
%xticks(1:length(k_values));
%xticklabels(string(k_values));
%ylabel("Users");



%% Point 2
p = 0.95;
alpha =(1-p)/2;
ni =simulations_number-1;

sample_mean = mean(point_2_matrix,2);
sample_var = var(point_2_matrix,0,2);
c = tinv(1-alpha,ni);

INF_Confidence_Interval = sample_mean-c*sqrt(sample_var)/(sqrt(simulations_number));
SUP_Confidence_Interval = sample_mean+c*sqrt(sample_var)/(sqrt(simulations_number));
Width_Confidence_Interval = SUP_Confidence_Interval-INF_Confidence_Interval;
Relative_width = Width_Confidence_Interval/sample_mean;

%%
figure();
errorbar(k_values,sample_mean,sample_mean-INF_Confidence_Interval,SUP_Confidence_Interval-sample_mean,".");
xticks(k_values);
yticks(2:2:26);
xlabel("K");
ylabel("Average number of users");
title("Average number of users (CI)");

toc
