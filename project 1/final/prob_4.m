function prob_4(N,generated_data, plot_histograms)
%% Constants
m = 501;
dt = 0.5;
alpha = 0.6;
sigma = 0.5;
v = 90;
eta = 3;
zeta = 1.5;
eff_vec = 1:1:m; % times at which we compute Efficient sample size
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
%%
if generated_data
    [X_true, Y] = generate_data();
end
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
driver_hist = zeros(5,m);
driver = randsample(5,N, true);  %Z0 index
driver_hist2= zeros(N,501);
driver_hist2(:,1)= driver;

driver_plot = zeros(m,1);
driver_plot(1) = mode(driver);
[temp1, temp2, driver_hist(:,1)]= groupcounts(driver);

idx = 1;
count = 0;
countrate = round(m/5);

for time=2:m
    if mod(time, resamplerate) == 0 
        %resample
        indices = randsample(N,N,true,weights(:,time-1));
        X(:,:,time-1) = X(:,indices,time-1);
        driver = driver(indices);
        
        %draw as usual but without mult. weights
        W = psi_w*randn(2,N)*sigma; %Wn+1 N times
        X(:,:,time) = phi*X(:,:,time-1) + psi_z*Z(:,driver) + W;    
        driver = cellfun(@randsample,num2cell(ones(1,N)*5,1),num2cell(ones(1,N),1),num2cell(true(1,N)),num2cell(P(:,driver),1)); %driver n+1
        driver_hist2(:,time)= driver;

        for station=1:6
            mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,time); X(4,:,time)] - pos_vec(:,station),2,1));
        end
        weights(:,time) = cellfun(@mvnpdf,num2cell(repmat(Y(:,time),1,N),1), num2cell(mu,1) ,zeta_matrices)'; %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
        driver_plot(time) =mode(driver);
        [temp1, temp2, driver_hist(:,time)]= groupcounts(driver');        
        
    else
        W = psi_w*randn(2,N)*sigma; %Wn+1 N times
        X(:,:,time) = phi*X(:,:,time-1) + psi_z*Z(:,driver) + W;    
        driver = cellfun(@randsample,num2cell(ones(1,N)*5,1),num2cell(ones(1,N),1),num2cell(true(1,N)),num2cell(P(:,driver),1)); %driver n+1
        driver_hist2(:,time)= driver;

        for station=1:6
            mu(station,:) = v -10*eta*log10(vecnorm( [X(1,:,time); X(4,:,time)] - pos_vec(:,station),2,1));
        end
        weights(:,time) = weights(:,time-1).*cellfun(@mvnpdf,num2cell(repmat(Y(:,time),1,N),1), num2cell(mu,1) ,zeta_matrices)'; %mvnpdf(Y(:,1), mu, diag(ones(6,1)*zeta^2)) (for each particle where mu=mu(part))
        
        mw = max(weights(:,time));
        weights(:,time) =weights(:,time)/mw;
        
        driver_plot(time) =mode(driver);
        [temp1, temp2, driver_hist(:,time)]= groupcounts(driver');        
    end
    

    
    
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
if plot_histograms
for i =[1 10 50 100 400]
    figure(count)
    count = count + 1;
    [~,edges] = histcounts(log10(weights(:,i)));
    histogram(weights(:,i), 10.^edges)
    set(gca, 'xscale', 'log')
    title(['Time =' num2str(i)])
end
end
%% PLOT average path
figure(count+1)
hold on
plot(tau_1,tau_2)
plot(pos_vec(1,:),pos_vec(2,:),'*')
title('Estimated path')
%%
if generated_data
    plot(X_true(1,:), X_true(4,:))
end
legend('Approx','Stations', 'True');

hold off

max_driver = zeros(1,m);
for time = 1:m
    max_driver(time) = max(driver_hist(:,time));
end
%plot(max_driver)
%driver_hist

 %%
% comm = zeros(5,m);
% for time = 1:m
%     for part = 1:N
%         comm(driver_hist2(part, time),time) = comm(driver_hist2(part, time),time) + weights(part, time);
%     end
%     comm(:,time)=comm(:,time)/sum(comm(:,time));
% end
% 
% 
% prob_comm = zeros(1,m);
% for time = 1:m
%     [val, I]= max(comm(:,time));
%     prob_comm(time) = val;
% end
% figure(count +3)
% area(1:1:m,comm')

mp_comm = zeros(m,1);
for time = 1:m
    comm_vec = zeros(5,1);
    for part = 1:N
        curr_dc = driver_hist2(part,time);
        comm_vec(curr_dc) = comm_vec(curr_dc) + weights(part,time);
    end
    [M,I] = max(comm_vec);
    mp_comm(time) = I;
end

plot(mp_comm,'*')
[GC,GR] = groupcounts(mp_comm)
title(['Count for Stat, East, Norh, South, West = '  num2str(GC')])
yticks([1 2 3 4 5])
yticklabels({'Stationary','East','North', 'South', 'West'})
xlabel('Time')
end