function [tau_1, tau_2] = prob_3(N, generated_data, eff_vec, hist_vec)
%% Constants
m = 501;
dt = 0.5;
alpha = 0.6;
sigma = 0.5;
v = 90;
eta = 3;
zeta = 1.5;
CV = zeros(length(eff_vec), 1);
ess = zeros(length(eff_vec), 1);

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
%% Self 
if generated_data
    [X_true, Y] = generate_data();
end
    %%
X = zeros(6,N,m);
X(:,:,1) = mvnrnd(zeros(6,1), diagg, N)'; %X0
weights = zeros(N,m); %weight0

%% Vectorized. Do not need chi(X_0) since g_0=xhi_0, z_0=xhi_0*p_0 --> z_0=p(X_0,Y_0)
mu = zeros(6,N);
for station=1:6
    mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,1); X(4,:,1)] - pos_vec(:,station),2,1));
end
%weights(:,1) = cellfun(@mvnpdf,num2cell(X(:,:,1),1),
%repmat(num2cell(zeros(6,1),[1 2]),1,N), repmat(num2cell(diagg,[1
%2]),1,N));  = chi(X_0)
zeta_matrices = repmat(num2cell(diag(ones(6,1)*zeta^2),[1 2]),1,N);
weights(:,1) = cellfun(@mvnpdf,num2cell(repmat(Y(:,1),1,N),1), num2cell(mu,1) , zeta_matrices); %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
mw = max(weights(:,1));
weights(:,1) =weights(:,1)/mw;
%%
driver = randsample(5,N, true);  %Z0 index

idx = 1;
count = 0;
countrate = round(m/10);
for time=2:m
    W = psi_w*randn(2,N)*sigma; %Wn+1 N times
    X(:,:,time) = phi*X(:,:,time-1) + psi_z*Z(:,driver) + W;    
    driver = cellfun(@randsample,num2cell(ones(1,N)*5,1),num2cell(ones(1,N),1),num2cell(true(1,N)),num2cell(P(:,driver),1)); %driver n+1
    
    for station=1:6
        mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,time); X(4,:,time)] - pos_vec(:,station),2,1));
    end
    weights(:,time) = weights(:,time-1).*cellfun(@mvnpdf,num2cell(repmat(Y(:,time),1,N),1), num2cell(mu,1) ,zeta_matrices)'; %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
    mw = max(weights(:,time));
    weights(:,time) =weights(:,time)/mw;
    
    if time == eff_vec(idx)
        big_omega = sum(weights(:,time));

        CV(idx) = sqrt(N)*norm(weights(:,time)./big_omega - 1/N);
        %sqrt(N)*sqrt(sum((weights(:,time)./big_omega - 1/N).^2));
        ess(idx) = N/(1+CV(idx)^2);
        
        idx = idx +1;
    end
    
    if count >= countrate
        disp(['Time = ' num2str(time) ' out of ' num2str(m)]);
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
for idx = 1:length(hist_vec)
    i = hist_vec(idx);
    figure(count)
    count = count + 1;
    [~,edges] = histcounts(log10(weights(:,i)));
    histogram(weights(:,i),10.^edges)
    set(gca, 'xscale','log')
    title(['Time = ' num2str(i)])
end
%% PLOT average path
figure(count+1)
hold on
plot(tau_1,tau_2)
plot(pos_vec(1,:),pos_vec(2,:),'*')
%%
if generated_data
    plot(X_true(1,:), X_true(4,:))
end
%%
legend('Approx','Stations', 'True');
hold off
%%
figure(count +2)
plot(eff_vec, ess)
title('ESS')
end