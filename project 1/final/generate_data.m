function [X Y] = generate_data()
m = 501;
dt = 0.5;
alpha = 0.6;
sigma = 0.5;
v = 90;
eta = 3;
zeta = 1.5;

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

X = zeros(6,m);
X(:,1) = mvnrnd(zeros(6,1), diagg, 1)'; %X0
driver = randsample(5,1, true);  %Z0 index


for time  = 2:m
    W = randn(2,1).*sigma; %Wn+1
    X(:,time) = phi*X(:,time-1) + psi_z*Z(:,driver) + psi_w*W;
    %Driver n+1
    driver = randsample(5,1,true,P(:,driver));
end

Y = zeros(6,m);
for time=1:m
    for station =1:6
        Y(station,time) = v - 10*eta*log10(norm([X(1,time); X(4,time)] - pos_vec(:,station))) + randn()*zeta;
    end
end

end