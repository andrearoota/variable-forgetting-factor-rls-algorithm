%% MAIN.M
% This script simulates and compares different Recursive Least Squares (RLS)
% algorithms for system identification.
%
% Algorithms compared:
% 1. Classic RLS with fixed forgetting factor.
% 2. Classic RLS with no forgetting (lambda = 1).
% 3. Variable Forgetting Factor RLS (VFF-RLS).
%
% The system simulates a time-varying channel with abrupt changes and
% continuous variations to test the tracking capabilities of the algorithms.

clc
clearvars

% Set random seed for reproducibility
rng(1);

%% Simulation Parameters
P_delta = 10^6;    % Regularization parameter for initializing the inverse correlation matrix P
N = 1000;          % Number of samples in the simulation
L = 2;             % Number of system coefficients (filter taps)
lambda_rls = 0.95; % Fixed forgetting factor for the standard RLS comparison
start_sin = 200;   % Sample index where sinusoidal variation starts

%% System Generation
% Define noise power
sigma_v = 1;
% Generate Gaussian white noise
noise = randn(N,1) * sqrt(sigma_v);

% Generate true system parameters (time-varying weights)
true_weights = generate_true_weights(N, L, start_sin);

% Generate input signal (regressor)
input_matrix = generate_input(N, L);

% Generate observed output signal (with noise)
observed_signal = generate_output(true_weights, input_matrix, noise);

% Generate clean output signal (noiseless) for reference/VFF calculation
clean_signal = generate_output(true_weights, input_matrix, zeros(N,1));

% Visualize the noiseless desired signal
% figure;
% plot(clean_signal);
% title('Desired Signal (Noiseless)');
% grid on;

%% Algorithm Execution

% 1. RLS with fixed forgetting factor (lambda < 1)
disp('Running RLS (lambda=0.95)...');
[weights_classic] = classic_rls(input_matrix, observed_signal, P_delta, lambda_rls);
mis_classic = misalignment(weights_classic, true_weights);

% 2. Classic RLS with infinite memory (lambda = 1)
disp('Running RLS (lambda=1)...');
[weights_plain] = classic_rls(input_matrix, observed_signal, P_delta, 1);
mis_plain = misalignment(weights_plain, true_weights);

% 3. Variable Forgetting Factor RLS (VFF-RLS)
disp('Running VFF-RLS...');
% Note: VFF-RLS here uses the clean signal for optimal parameter tuning/oracle comparison
[weights_vff, lambda_vff] = vff_rls(input_matrix, observed_signal, clean_signal, P_delta, 2, 10);
mis_vff = misalignment(weights_vff, true_weights);

%% Plotting Results
disp('Generating plots...');
generate_plots(true_weights, weights_classic, weights_vff, mis_classic, mis_vff, ...
               ones(1,N)*lambda_rls, lambda_vff, weights_plain, mis_plain, ones(1,N));

