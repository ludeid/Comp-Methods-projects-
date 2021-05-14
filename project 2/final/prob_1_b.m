function [t,lambda, theta] = prob_1_b(N,m,rho, d, varTheta,tau)
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
    theta(part,1) = gamrnd(2*d +2,1/(varTheta));
    for i = 1:d
        lambda(i,part,1) = gamrnd(n_func(i,t(:,part,1),tau) +2 , 1/theta(part,1));
    end
end
%%
for time = 2:m
    for part = 1:N
        % theta
        theta(part, time) =gamrnd(2*d +2,1/(varTheta +sum(lambda(:,part,time-1))));
        % lambda
        for i = 1:d
            lambda(i,part,time) = gamrnd(n_func(i,t(:,part,time-1),tau) +2 , 1/(theta(part,time)+ t(i+1, part, time-1)- t(i,part, time-1)));
        end
        % MH-t using RV
        cand_vec = zeros(d+1,1);
        cand_vec(1) = t_1;
        cand_vec(d+1) = t_end;
        for i = 2:d
            R = rho*(t(i+1,part,time-1) -t(i-1,part,time-1));
            cand_vec(i) =t(i,part,time-1) + 2*R*rand() -R;
        end
        temp = f_t_cond(cand_vec, lambda(:,part,time), tau);
        alpha = min(1,temp/f_t_cond(t(:,part,time-1),  lambda(:,part,time), tau));
        %disp(['temp ' num2str(temp)])
        %disp(['old ' num2str(f_t_cond(t(:,part,time-1),  lambda(:,part,time), tau))])
        %disp(['alpha ' num2str(alpha)])
        if (rand()<=alpha && temp ~= 0 )
            t(:,part,time) = cand_vec;
%             if (round(cand_vec(2)) <= 1851 || round(cand_vec(2))>= 1963)
%                 disp('__')
%                 disp(num2str(round(cand_vec(2))))
%                 disp(['temp ' num2str(temp)])
%                 disp(['old ' num2str(f_t_cond(t(:,part,time-1),  lambda(:,part,time), tau))])
%                 disp(['alpha ' num2str(alpha)])
%             end
        else
            t(:,part,time) = t(:,part,time-1);
        end
    end
end

t = round(t);
%disp('t')
%disp(t(:,1,m))
end

