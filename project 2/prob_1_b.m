function t = prob_1_b(N,m,rho, d, varTheta)
%% Prob 1
% N = 1;
% m = 10^5;
% d =5;
% rho = 0.01;
% varTheta = 5;

t_1 = 1851;
t_end = 1963;

%%
theta = zeros(N,m);
lambda = zeros(d,N,m);
t = zeros(d+1,N,m);
t(1,:,:) = ones(N,m)*t_1;
t(d+1,:,:) = ones(N,m)*t_end;

%initial States
for part = 1:N
    t(:,part,1) = linspace(t_1,t_end,d+1);
    theta(part,1) = gamrnd(2,1/varTheta);
    lambda(:,part,1) = gamrnd(2,1/theta(part,1),d,1);
end
%%
for time = 2:m
    for part = 1:N
        % theta
        theta(part, time) =gamrnd(2,1/varTheta);
        % lambda
        lambda(:,part,time) = gamrnd(2,1/theta(part,time),d,1);
        % MH-t using RV
        cand_vec = zeros(d+1,1);
        cand_vec(1) = t_1;
        cand_vec(d+1) = t_end;
        for i = 2:d
            R = rho*(t(i+1,part,time-1) -t(i-1,part,time-1));
            t(i,part,time) =t(i,part,time-1) + 2*R*rand() -R;
        end
        alpha = min(1,f(cand_vec)/f(t(:,part,time)));
        if rand()<=alpha
            t(:,part,time) = cand_vec;
        else
            t(:,part,time) = t(:,part,time-1);
        end
    end
end

t = round(t);
%disp(t(:,1,m))
end

