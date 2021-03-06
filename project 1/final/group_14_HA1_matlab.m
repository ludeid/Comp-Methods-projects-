clear
close all

%% Problem 1
m = 1000;

prob_1(m)  %Simulating trajectory with m timesteps



%% Problem 3
N = 10^2; %number of particles
generated_data = false; %true if one want to check correctedness of algorithm against data personally generated 
hist_vec = [1 5 100 400]; %times to show (semilog-)histograms of weights

[tau_1, tau_2] = prob_3(N,generated_data, hist_vec);



%% Problem 4
N = 10^4; %number of particles
plot_histograms = true; %plot histograms for (already) selected timepoints
generated_data = false; %same as above

prob_4(N,generated_data, plot_histograms)



%% Problem 5
N = 10^4; %%number of particles in each run of zeta
D = 0;
%choose to either run algotirhm with zeta over a evenly spaced mesh with D
%points or specified zeta_vec. If zeta_vec is passed as an empty array then
%D-mesh will be used.

zeta_vec = [2.25 2.3 2.35 2.4 2.45 2.5];

%or

%D = 10;    %number of points in mesh evenly spaced inside (0,3). D = 3 yields mesh [0.75    1.50    2.25]
%zeta_vec = [];

prob_5(N,D, zeta_vec)