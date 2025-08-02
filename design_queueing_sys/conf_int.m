clear
close all

p = 0.95;
N_MAX = 1000;
alpha = (1-p)/2;
threshold = 0.1;
rng(9);

A = run_system;

for n=2:N_MAX
    ni = n-1;
    A(n) = run_system;
    sample_mean = mean(A);
    sample_variance=var(A);
    c=tinv(1-alpha,ni);
    INF_Confidence_Interval(ni)=sample_mean-c*sqrt(sample_variance)/sqrt(n)
    SUP_Confidence_Interval(ni)=sample_mean+c*sqrt(sample_variance)/sqrt(n)
    Width_Confidence_Interval(ni)=SUP_Confidence_Interval(ni)-INF_Confidence_Interval(ni);
    Relative_width(ni)=Width_Confidence_Interval(ni)/sample_mean; 
    if Relative_width(ni)<threshold
        END = 1
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