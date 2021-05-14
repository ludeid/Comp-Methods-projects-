%%
close all

Y = importfile("mixture-observations.csv", [1, Inf]);

%histogram(Y)

theta_1 = 0.5;
n = length(Y);
N = 10^2;

theta= zeros(N,1);
theta(1) = theta_1;
for time = 2:N
    f_vec = zeros(n,1);
    for i = 1:n
        f_vec(i) = theta(time-1)*exp(-(Y(i)-1)^2/8)/(2*exp(-(Y(i)^2)/2)*(1-theta(time-1)) + theta(time-1)*exp(-(Y(i)-1)^2/8));
    end
    
    theta(time) = sum(f_vec)/n;
end
%%
plot(theta)
legend(['theta = ' num2str(theta(end))]) 