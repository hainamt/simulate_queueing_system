function sys_point2(average_w_time,state_residence_time_pdf)
    variance_value = var(state_residence_time_pdf);
    average_value = average_w_time;
    p = 0.95;
    N_MAX = 1000;
    alpha = (1-p)/2;
    threshold = 0.1;
    std_value =sqrt(variance_value);
    A = average_value+rand*std_value;

    for n=2:N_MAX
        ni = n-1;
        A(n) = run_simulation
        sample_mean = mean(A);
        sample_variance = var(A);
        c= tinv(1-alpha,ni);
        INF_Confidence_Interval(ni)=sample_mean-c*sqrt(sample_variance)/sqrt(n)
        SUP_Confidence_Interval(ni)=sample_mean+c*sqrt(sample_variance)/sqrt(n)
        Width_Confidence_Interval(ni)=SUP_Confidence_Interval(ni)-INF_Confidence_Interval(ni);
        Relative_width(ni)=Width_Confidence_Interval(ni)/sample_mean;
        if Relative_width(ni)<threshold %we reached the goal
            End=1
            n
            Relative_width(ni)
        break
        end
    end
    figure
    plot(Width_Confidence_Interval)
    figure
    plot(Relative_width)
    figure
    plot(INF_Confidence_Interval)
    hold on
    plot(SUP_Confidence_Interval)
end
