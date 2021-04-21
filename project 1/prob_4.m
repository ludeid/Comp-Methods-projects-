%% Constants
clear
m = 501;
N = 10^3;
dt = 0.5;
alpha = 0.6;
sigma = 0.5;
v = 90;
eta = 3;
zeta = 1.5;
resamplerate = 10;



P = [
16 1 1 1 1;
1 16 1 1 1;
1 1 16 1 1;
1 1 1 16 1;
1 1 1 1 16].*(1/20);

phi = [
1 dt dt^2/2 0 0 0;
0 1 dt 0 0 0;
0 0 alpha 0 0 0;
0 0 0 1 dt dt^2/2;
0 0 0 0 1 dt;
0 0 0 0 0 alpha];

psi_z = [dt^2/2 0;
    dt 0;
    0 0;
    0 dt^2/2;
    0 dt;
    0 0];
psi_w = [dt^2/2 0;
    dt 0;
    1 0;
    0 dt^2/2;
    0 dt;
    0 1];

Z = [
    0 3.5 0 0 -3.5;
    0 0 3.5 -3.5 0];

stations = matfile('stations.mat');
pos_vec = stations.pos_vec;
RSSI_obj = matfile('RSSI-measurements.mat');
Y = RSSI_obj.Y;
%%
X = zeros(6,N,m);
X(:,:,1) = mvnrnd(zeros(6,1), [500 0 0 0 0 0;0 5 0 0 0 0;0 0 5 0 0 0;0 0 0 200 0 0;0 0 0 0 5 0;0 0 0 0 0 5], N)'; %X0
weights = zeros(N,m); %weight0
for part=1:N
    weights(part,1) = mvnpdf(X(:,part,1),zeros(6,1), [500 0 0 0 0 0;0 5 0 0 0 0;0 0 5 0 0 0;0 0 0 200 0 0;0 0 0 0 5 0;0 0 0 0 0 5]);  %xhi(0)*p0
    mu = zeros(6,1);
    for i=1:6
        mu(i) = v -10*eta*log10(norm( [X(1,part,1); X(4,part,1)] - pos_vec(:,i)));
    end
    weights(part,1) = weights(part,1)*mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2));
    weights(part,1) = weights(part,1)/( mvnpdf(X(:,part,1), zeros(6,1), [500 0 0 0 0 0;0 5 0 0 0 0;0 0 5 0 0 0;0 0 0 200 0 0;0 0 0 0 5 0;0 0 0 0 0 5]));
end

mw = max(weights(:,1));
weights(:,1) =weights(:,1)/mw;
%%
driver = round(rand()*4)+1;  %Z0 index


count = 0;
for time=2:m
    for part= 1:N        
        W = randn(2,1).*sigma; %Wn+1
        W = psi_w*W;
        X(:,part,time) = phi*X(:,part,time-1) + psi_z*Z(:,driver) + W;
        %Driver n+1
        driver = randsample(5,1,true,P(:,driver));
        
        %weights
        for i=1:6
            mu(i) = v -10*eta*log10(norm( [X(1,part,time); X(4,part,time)] - pos_vec(:,i)));
        end
        weights(part,time) = weights(part,time-1)*mvnpdf(Y(:,time), mu, diag(ones(6,1)*zeta^2));
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
    end
    if count >= resamplerate
        count=0;
        weights(:,time) = randsample(N,N,true, weights(:,time));
    else
        count = count+1;
    end
end

tau_1 = zeros(1,m);
tau_2 = zeros(1,m);

for time = 1:m
    big_omega = sum(weights(:,time));

    tau_1(time) = sum(weights(:,time).*X(1,:,time)')/big_omega; 
    tau_2(time) = sum(weights(:,time).*X(4,:,time)')/big_omega;
end

hold on
plot([1:1:501],tau_1)
plot([1:1:501],tau_2)
hold off
% %%
% plot(tau_1,tau_2)