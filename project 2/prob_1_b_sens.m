% varTheta
var_points = 10;
varT_vec = linspace(0.0002,1000,var_points);
% rho
rho_points = 10;
rho_vec = 0.01;
rho_vec = linspace(0.00001,1000,rho_points);

N = 1;
m = 10^3;
d = 4;

%% investigate both
data = zeros(d+1,length(varT_vec)*length(rho_vec)); 

idx = 1;
for varTheta= varT_vec
    for rho = rho_vec
        [t, lambda, theta] = prob_1_b(N,m,rho, d, varTheta);
        data(:,idx) = t(:,1,m);
        idx = idx +1;
%         if idx == 13
%             disp(rho)
%         end
    end
    disp([num2str(idx) ' out of ' num2str(var_points*rho_points)])
end

%% investigate varTheta, rho fixed
burnIn = 10^2;
d = 2;
m = 10^4;
rho = 0.01;
N = 1;

l = zeros(length(varT_vec),d,m-burnIn+1);
breakPoint = zeros(length(varT_vec),m-burnIn+1);
idx = 1;
for varTheta= varT_vec
    [t, lambda, theta] = prob_1_b(N,m,rho, d, varTheta);
    for i = 1:d
        l(idx, i, :) = lambda(i,1,burnIn:m);
    end
    if d == 2
        breakPoint(idx,:) = t(2,1,burnIn:m);
    end
    disp(num2str(idx))
    idx = idx+1;
end

%%
for i = 1:d
    figure(i)
    hold on
    histogram(l(15,i,:))
    histogram(l(50,i,:))
    histogram(l(99,i,:))
%    xlim([0 40]);
    legend('15','50','99')
    hold off
end
figure(6)
histogram(l(5,1,:))

figure(5)
if d == 2
    hold on
    histogram(breakPoint(10,:),1000)
    histogram(breakPoint(50,:),1000)
    histogram(breakPoint(100,:),1000)
    hold off
    legend('10','50','100')
end
