%% investigate varTheta, rho fixed
clear all
close all

tau_obj = readtable('coal-mine.csv');
tau = table2array(tau_obj);

var_points = 3;
varT_vec = linspace(0.1,10,var_points);

burnIn = 10^2;
d = 2;
m = 10^5;
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
    disp(num2str(idx))
    T(:,idx) = theta;
    idx = idx+1;
end

%%
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
