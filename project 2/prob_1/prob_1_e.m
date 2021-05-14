%% investigate rho, rho fixed
clear all
close all

tau_obj = readtable('coal-mine.csv');
tau = table2array(tau_obj);

rho_points = 3;
rho_vec = linspace(0.01,0.15,rho_points);
rho_vec = [0.01 0.1 5];
burnIn = 10^2;
d = 2;
m = 10^5;
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

%%
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
autocorr(breakPoint(1,:))
figure(11)
autocorr(breakPoint(2,:))
figure(12)
autocorr(breakPoint(3,:))
% figure(10)
% [r,lags]= autocorr(breakPoint(1,:));
% plot(lags, r)
% xlim([0 10100])
% legend(['rho = ' num2str(rho_vec(1))])
% figure(11)
% [r,lags]= xcorr(breakPoint(1,:));
% plot(lags, r)
% xlim([0 10100])
% breakPoint(round(rho_points/2),:))
% legend(['rho = ' num2str(rho_vec(round(rho_points/2)))])
% figure(12)
% autocorr(breakPoint(rho_points,:))
% legend(['rho = ' num2str(rho_vec(rho_points))])