function [weights, X] = SISR(zeta, N, m)
%% Constants
dt = 0.5;
alpha = 0.6;
sigma = 0.5;
v = 90;
eta = 3;
%zeta = 1.5;
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
RSSI_obj = matfile('RSSI-measurements-unknown-sigma.mat');
Y = RSSI_obj.Y;

%% Generate X_0,w_0
zeta_matrices = repmat(num2cell(diag(ones(6,1)*zeta^2),[1 2]),1,N);

X = zeros(6,N,m);
X(:,:,1) = mvnrnd(zeros(6,1), diagg, N)'; %X0
weights = zeros(N,m); %weight0

mu = zeros(6,N);
for station=1:6
    mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,1); X(4,:,1)] - pos_vec(:,station),2,1));
end

weights(:,1) = cellfun(@mvnpdf,num2cell(repmat(Y(:,1),1,N),1), num2cell(mu,1) , zeta_matrices); %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
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
countrate = round(m/5);

for time=2:m
    if mod(time, resamplerate) == 0 
        %resample
        indices = randsample(N,N,true,weights(:,time-1));
        X(:,:,:) = X(:,indices,:);
        
        %draw as usual but without mult. weights
        W = psi_w*randn(2,N)*sigma; %Wn+1 N times
        X(:,:,time) = phi*X(:,:,time-1) + psi_z*Z(:,driver) + W;    
        driver = cellfun(@randsample,num2cell(ones(1,N)*5,1),num2cell(ones(1,N),1),num2cell(true(1,N)),num2cell(P(:,driver),1)); %driver n+1

        for station=1:6
            mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,time); X(4,:,time)] - pos_vec(:,station),2,1));
        end
        weights(:,time) = cellfun(@mvnpdf,num2cell(repmat(Y(:,time),1,N),1), num2cell(mu,1) ,zeta_matrices)'; %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
        driver_plot(time) =mode(driver);
        
        
        
    else
        W = psi_w*randn(2,N)*sigma; %Wn+1 N times
        X(:,:,time) = phi*X(:,:,time-1) + psi_z*Z(:,driver) + W;    
        driver = cellfun(@randsample,num2cell(ones(1,N)*5,1),num2cell(ones(1,N),1),num2cell(true(1,N)),num2cell(P(:,driver),1)); %driver n+1

        for station=1:6
            mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,time); X(4,:,time)] - pos_vec(:,station),2,1));
        end
        weights(:,time) = weights(:,time-1).*cellfun(@mvnpdf,num2cell(repmat(Y(:,time),1,N),1), num2cell(mu,1) ,zeta_matrices)'; %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
        
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
        
        %driver_hist(:,time) = driver;
        driver_plot(time) =mode(driver);
    end
    

    
    
%     if time == eff_vec(idx)
%         big_omega = sum(weights(:,time));
% 
%         CV(idx) = sqrt(N)*norm(weights(:,time)./big_omega - 1/N);
%         %sqrt(N)*sqrt(sum((weights(:,time)./big_omega - 1/N).^2));
%         ess(idx) = N/(1+CV(idx)^2);
%         
%         idx = idx +1;
%     end
    
    if count >= countrate
        disp(['Time= ' num2str(time)]);
        count = 0;
    else
        count = count+1;
    end
end
end