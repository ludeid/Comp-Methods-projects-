%% Problem 1.c
%investigate final t vector for different d
tau_obj = readtable('coal-mine.csv');
tau = table2array(tau_obj);

d_vec = [2 3 4 5 6];
m = 10^3;
rho = 0.1;
N = 1;
varTheta = 4;

T = cell(length(d_vec),1);
idx = 1;
for d = d_vec
    [t, lambda, theta] = prob_1_b(N,m,rho, d, varTheta,tau);
    T{idx} = t(:,1,end);
    idx = idx +1;
end
for i=1:length(T)
    disp(T{i})
end

%% 1.d


tau_obj = readtable('coal-mine.csv');
tau = table2array(tau_obj);

var_points = 3;
varT_vec = linspace(0.1,10,var_points);

burnIn = 10^2;
d = 2;
m = 10^4;
rho = 10;
N = 1;

l = zeros(length(varT_vec),d,m-burnIn+1);
breakPoint = zeros(length(varT_vec),m-burnIn+1);
T = zeros(m,var_points);
idx = 1;
for varTheta= varT_vec
    [t, lambda, theta] = prob_1_b(N,m,rho, d, varTheta,tau);
    for i = 1:d
        l(idx, i, :) = lambda(i,1,burnIn:m);
    end
    if d == 2
        breakPoint(idx,:) = t(2,1,burnIn:m);
    end
    %disp(num2str(idx))
    T(:,idx) = theta;
    idx = idx+1;
end


for i = 1:d
    figure(i)
    hold on
    histogram(l(1,i,:))
    histogram(l(round(var_points/2),i,:))
    histogram(l(var_points,i,:))
    %xlim([0 2.5*10^2]);
    legend(['varTheta = ' num2str(varT_vec(1))],['varTheta = ' num2str(varT_vec(round(var_points/2)))],['varTheta = ' num2str(varT_vec(var_points))])
    hold off
end


figure(5)
if d == 2
    hold on
    histogram(breakPoint(1,:),30)
    histogram(breakPoint(round(var_points/2),:),30)
    histogram(breakPoint(var_points,:),30)
    hold off
    legend(['varTheta = ' num2str(varT_vec(1))],['varTheta = ' num2str(varT_vec(round(var_points/2)))],['varTheta = ' num2str(varT_vec(var_points))])
end

figure(6)
hold on
histogram(T(burnIn:m,1))
histogram(T(burnIn:m,round(var_points/2)))
histogram(T(burnIn:m,var_points), 'FaceColor', 'black')
hold off
legend(['varTheta = ' num2str(varT_vec(1))],['varTheta = ' num2str(varT_vec(round(var_points/2)))],['varTheta = ' num2str(varT_vec(var_points))])


%% 1.e
clear all

tau_obj = readtable('coal-mine.csv');
tau = table2array(tau_obj);

rho_points = 3;
rho_vec = linspace(0.1,10,rho_points);

burnIn = 10^2;
d = 2;
m = 10^4;
N = 1;
varTheta= 5;

l = zeros(length(rho_vec),d,m-burnIn+1);
breakPoint = zeros(length(rho_vec),m-burnIn+1);
T = zeros(m,rho_points);
idx = 1;
for rho= rho_vec
    [t, lambda, theta] = prob_1_b(N,m,rho, d, varTheta,tau);
    for i = 1:d
        l(idx, i, :) = lambda(i,1,burnIn:m);
    end
    if d == 2
        breakPoint(idx,:) = t(2,1,burnIn:m);
    end
    disp(num2str(idx))
    T(:,idx) = theta;
    idx = idx+1;
end

for i = 1:d
    figure(i)
    hold on
    histogram(l(1,i,:))
    histogram(l(round(rho_points/2),i,:))
    histogram(l(rho_points,i,:))
    %xlim([0 2.5*10^2]);
    legend(['rho = ' num2str(rho_vec(1))],['rho = ' num2str(rho_vec(round(rho_points/2)))],['rho = ' num2str(rho_vec(rho_points))])
    hold off
end


figure(5)
if d == 2
    hold on
    histogram(breakPoint(1,:))
    histogram(breakPoint(round(rho_points/2),:))
    histogram(breakPoint(rho_points,:))
    hold off
    legend(['rho = ' num2str(rho_vec(1))],['rho = ' num2str(rho_vec(round(rho_points/2)))],['rho = ' num2str(rho_vec(rho_points))])
end

figure(6)
hold on
histogram(T(burnIn:m,1))
histogram(T(burnIn:m,round(rho_points/2)))
histogram(T(burnIn:m,rho_points)) %, 'FaceColor', 'black'
hold off
legend(['rho = ' num2str(rho_vec(1))],['rho = ' num2str(rho_vec(round(rho_points/2)))],['rho = ' num2str(rho_vec(rho_points))])

%mixing
figure(7)
plot(breakPoint(1,:))
legend(['rho = ' num2str(rho_vec(1))])
figure(8)
plot(breakPoint(round(rho_points/2),:))
legend(['rho = ' num2str(rho_vec(round(rho_points/2)))])
figure(9)
plot(breakPoint(rho_points,:))
legend(['rho = ' num2str(rho_vec(rho_points))])

figure(10)
autocorr(breakPoint(1,:));
legend(['rho = ' num2str(rho_vec(1))])
figure(11)
autocorr(breakPoint(round(rho_points/2),:))
legend(['rho = ' num2str(rho_vec(round(rho_points/2)))])
figure(12)
autocorr(breakPoint(rho_points,:))
legend(['rho = ' num2str(rho_vec(rho_points))])

%% PROBLEM 2 c)

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

plot(theta)
legend(['theta = ' num2str(theta(end))]) 
