clear
close all

%% Problem 1
m = 1000;

prob_1(m)  %Simulating trajectory with m timesteps



%% Problem 3
N = 10^4; %number of particles
m = 501; 
generated_data = false; %true if one want to check correctedness of algorithm against data personally generated 
eff_vec = 1:1:m; % times at which we compute Efficient sample size
hist_vec = [1 5 100 400]; %times to show histograms of weights

[tau_1, tau_2] = prob_3(N,generated_data, eff_vec, hist_vec);



%% Problem 4
N = 10^2; %number of particles
plot_histograms = false; %plot histograms for (already) selected timepoints
generated_data = true; %same as above


prob_4(N,generated_data, plot_histograms)