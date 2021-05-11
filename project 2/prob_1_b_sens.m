clear all

% varTheta
var_points = 50;
varT_vec = linspace(0.002,100,var_points);
% rho
rho_points = 50;
rho_vec = linspace(0.0005,1,rho_points);

N = 1;
m = 10^3;
d = 4;

data = zeros(d+1,length(varT_vec)*length(rho_vec)); 

idx = 1;
for varTheta= varT_vec
    for rho = rho_vec
        t= prob_1_b(N,m,rho, d, varTheta);
        data(:,idx) = t(:,1,m);
        idx = idx +1;
    end
    disp([num2str(idx) ' out of ' num2str(var_points*rho_points)])
end
