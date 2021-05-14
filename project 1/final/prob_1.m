function prob_1(m)
%% Constants
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

%% 1.0: Trajectory
X = randn(6,1).*sqrt([500;5;5;200;5;5]); %X0
driver = randsample(5,1);  %Z0 index

X_hist = zeros(2, m);

for i=1:m
    W = randn(2,1).*sigma; %Wn+1
    W = psi_w*W;
    
    X = phi*X + psi_z*Z(:,driver) + W;
    
    %Driver n+1
    driver = randsample(5,1,true,P(:,driver));
    
    X_hist(:,i) = [X(1);X(4)];
end
%% PLot trajectory
plot(X_hist(1,:),X_hist(2,:))

end

