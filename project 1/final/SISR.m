function [weights, X] = SISR(zeta, N, m)
%% Constants
dt = 0.5;
alpha = 0.6;
sigma = 0.5;
v = 90;
eta = 3;
%zeta = 1.5;
%resamplerate = 1;



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
%% Driver
disp('Driver')
driver_hist = zeros(N,m);
driver = randsample(5,N, true);
driver_hist(:,1) = driver;

for time=2:m
    for part=1:N
        driver(part) = randsample(5,1,true,P(:,driver(part)));
    end
    driver_hist(:,time)=driver;
end
% W
W = randn(2,N,m)*sigma;

% X
X=zeros(6,N,m);
X(:,:,1) = mvnrnd(zeros(6,1), diagg, N)'; %X0

%Weights
zeta_matrices = repmat(num2cell(diag(ones(6,1)*zeta^2),[1 2]),1,N);
weights = zeros(N,m); %weight0
mu = zeros(6,N);
for station=1:6
    mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,1); X(4,:,1)] - pos_vec(:,station),2,1));
end
weights(:,1) = cellfun(@mvnpdf,num2cell(repmat(Y(:,1),1,N),1), num2cell(mu,1) , zeta_matrices); %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
mw = max(weights(:,1));
weights(:,1) =weights(:,1)/mw;

dz = det(zeta^2*eye(6))^(-1/2);
invZ = inv(zeta^2*eye(6));

disp('Particles')
for time = 2:m
    indices = randsample(N,N,true,weights(:,time-1));
    X(:,:,time-1) = X(:,indices,time-1);
    driver = driver(indices);
    
    X(:,:,time) = phi*X(:,:,time-1) + psi_z*Z(:,driver_hist(:,time)) + psi_w*W(:,:,time);
    
    mu = zeros(6,N);
    for station=1:6
        mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,time); X(4,:,time)] - pos_vec(:,station),2,1));
    end
    
    weights(:,time) = (2*pi)^(-6/2)*dz*exp(-1/2*( dot(Y(:,time)- mu, invZ*(Y(:,time)- mu))));
    mw = max(weights(:,time));
    weights(:,1) =weights(:,time)/mw;
end
end