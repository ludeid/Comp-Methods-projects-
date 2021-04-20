%% Constants
clear
m = 10000;
N = 10^4;
dt = 0.01;
alpha = 0.6;
sigma = 0.5;


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

%% 
stations = matfile('stations.mat');
pos_vec = stations.pos_vec;
RSSI_obj = matfile('RSSI-measurements.mat');
RSSI_meas = RSSI_obj.Y;

X = zeros(6,N,m);
X(:,:,1) = mvnrnd(zeros(6,1), [500 0 0 0 0 0;0 5 0 0 0 0;0 0 5 0 0 0;0 0 0 200 0 0;0 0 0 0 5 0;0 0 0 0 0 5], N)'; %X0
weights = zeros(N,m); %weight0

for time=1:m
    for part= 1:N
        mu = -phi*X(:,part,time) - psi_z*
        s2 =
        X(:,part,time) = mvnrnd(mu,s2,1);    %draw from q
        
    end
end