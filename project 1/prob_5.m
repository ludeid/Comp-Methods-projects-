clear
N = 10^3;
m = 501;

zeta_vec = linspace(0,3,D+2);
zeta_vec = zeta_vec(2:D+1);

l = zeros(D,m);
for idx=1:D
    disp(['D= ' num2str(idx)])
    [temp, weights] = SISR(zeta_vec(idx),N, m);
    for time=1:m
        summed = log(cumsum(sum(weights(:,:),1).^(1/time)));
        l(idx,time) = (-time-1)/time*log(N) + sum(summed(1:time));
    end
end
%%
[m_val, m_idx] = max(l(:,500));

disp(['Running for maximized zeta=' num2str(zeta_vec(m_idx))])
[weights, X] = SISR(zeta_vec(m_idx), N,m);


tau_1 = zeros(1,m);
tau_2 = zeros(1,m);

for time = 1:m
    big_omega = sum(weights(:,time));
   
    tau_1(time) = sum(weights(:,time).*X(1,:,time)')/big_omega; 
    tau_2(time) = sum(weights(:,time).*X(4,:,time)')/big_omega;
end

stations = matfile('stations.mat');
pos_vec = stations.pos_vec;
hold on
plot(tau_1,tau_2)
plot(pos_vec(1,:),pos_vec(2,:),'*')
title(['Estimated path for zeta= ' num2str(zeta_vec(m_idx))])
hold off