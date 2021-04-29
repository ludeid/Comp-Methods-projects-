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
eff_vec = [2 5 10 round(m/10) round(3*m/10) round(m/2) round(8*m/10) m];
CV = zeros(length(eff_vec), 1);
ess = zeros(length(eff_vec), 1);
resamplerate = 1;



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
diagg = [500 0 0 0 0 0;0 5 0 0 0 0;0 0 5 0 0 0;0 0 0 200 0 0;0 0 0 0 5 0;0 0 0 0 0 5];


stations = matfile('stations.mat');
pos_vec = stations.pos_vec;
RSSI_obj = matfile('RSSI-measurements.mat');
Y = RSSI_obj.Y;
%% Generate X_0,w_0
X = zeros(6,N,m);
X(:,:,1) = mvnrnd(zeros(6,1), diagg, N)'; %X0
weights = zeros(N,m); %weight0
for part=1:N
    weights(part,1) = mvnpdf(X(:,part,1),zeros(6,1), diagg);  %xhi(0)*p0
    mu = zeros(6,1);
    for i=1:6
        mu(i) = v -10*eta*log10(norm( [X(1,part,1); X(4,part,1)] - pos_vec(:,i)));
    end
    weights(part,1) = weights(part,1)*mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2));
    weights(part,1) = weights(part,1)/( mvnpdf(X(:,part,1), zeros(6,1), diagg));
end

mw = max(weights(:,1));
weights(:,1) =weights(:,1)/mw;
%%
driver_hist = zeros(N,m);
driver = randsample(5,N, true);  %Z0 index
driver_hist(:,1) = driver;
driver_plot = zeros(m,1);
driver_plot(1) = mode(driver);

idx = 1;
count = 0;
countrate = round(N/10);

for time=2:m
    if mod(time, resamplerate) == 0 
        %resample
        indices = randsample(N,N,true,weights(:,time-1));
        X(:,:,:) = X(:,indices,:);
        
        %draw as usual
        for part= 1:N
            W = randn(2,1).*sigma; %Wn+1
            X(:,part,time) = phi*X(:,part,time-1) + psi_z*Z(:,driver(part)) + psi_w*W;
            %Driver n+1
            driver(part) = randsample(5,1,true,P(:,driver(part)));
            
            %weights
            for i=1:6
                mu(i) = v -10*eta*log10(norm( [X(1,part,time); X(4,part,time)] - pos_vec(:,i)));
            end
            weights(part,time) = mvnpdf(Y(:,time), mu, diag(ones(6,1)*zeta^2));
        end
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
        
    else
        for part= 1:N        
            W = randn(2,1).*sigma; %Wn+1
            W = psi_w*W;
            X(:,part,time) = phi*X(:,part,time-1) + psi_z*Z(:,driver(part)) + psi_w*W;
            %Driver n+1
            driver(part) = randsample(5,1,true,P(:,driver(part)));
            
            %weights
            for i=1:6
                mu(i) = v -10*eta*log10(norm( [X(1,part,time); X(4,part,time)] - pos_vec(:,i)));
            end
            weights(part,time) = weights(part,time-1)*mvnpdf(Y(:,time), mu, diag(ones(6,1)*zeta^2));
            
        end
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
        
        %driver_hist(:,time) = driver;
        driver_plot(time) =mode(driver);
    end
    

    
    
    if time == eff_vec(idx)
        big_omega = sum(weights(:,time));

        CV(idx) = sqrt(N)*norm(weights(:,time)./big_omega - 1/N);
        %sqrt(N)*sqrt(sum((weights(:,time)./big_omega - 1/N).^2));
        ess(idx) = N/(1+CV(idx)^2);
        
        idx = idx +1;
    end
    
    if count >= countrate
        disp(time);
        count = 0;
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

%% Create semilog histograms
count = 1;
for i =[1 5 10]
    figure(count)
    count = count + 1;
    [~,edges] = histcounts(log10(weights(:,i)));
    histogram(weights(:,i),10.^edges)
    set(gca, 'xscale','log')
end

%figure(count)
% hold on
% plot([1:1:501],tau_1)
% plot([1:1:501],tau_2)
% hold off
%% PLOT average path
figure(count+1)
hold on
plot(tau_1,tau_2)
plot(pos_vec(1,:),pos_vec(2,:),'*')
hold off
%% Effective sample size
figure(count +2)
plot(eff_vec, ess)
%% Driver
