%% MAIN.M
% This script simulates and compares different Recursive Least Squares (RLS)
% algorithms for system identification.
%
% Algorithms compared:
% 1. Classic RLS with fixed forgetting factor.
% 2. Classic RLS with no forgetting (lambda = 1).
% 3. Variable Forgetting Factor RLS (VFF-RLS).
% 4. Gradient Variable Forgetting Factor RLS (GVFF-RLS).
%
% The system simulates a time-varying channel with abrupt changes and
% continuous variations to test the tracking capabilities of the algorithms.

clc
clearvars

% Set random seed for reproducibility
rng(1);

%% Simulation Parameters
% SNR = 20;       % Signal-to-Noise Ratio (dB) (Unused in current setup, noise variance is fixed)
P_delta = 10^6;   % Regularization parameter for initializing the inverse correlation matrix P
N = 1000;         % Number of samples in the simulation
L = 2;            % Number of system coefficients (filter taps)
lambda_rls = 0.95; % Fixed forgetting factor for the standard RLS comparison

with_gvff = false; % Flag to enable/disable GVFF-RLS plotting if needed
output_dir = "main"; % Name of the experiment for output file naming

% Ensure output directory exists for saving plots
if ~exist("output/"+output_dir, 'dir')
    mkdir("output/"+output_dir)
end

%% System Generation
% Define noise power
sigma_v = 1;
% Generate Gaussian white noise
noise = randn(N,1) * sqrt(sigma_v);

% Generate true system parameters (time-varying weights)
true_weights = course_theta(N, L);

% Generate input signal (regressor)
input_matrix = course_input(N, L);

% Generate observed output signal (with noise)
observed_signal = ma_realizations(true_weights, input_matrix, noise);

% Generate clean output signal (noiseless) for reference/VFF calculation
clean_signal = ma_realizations(true_weights, input_matrix, zeros(N,1));

% Visualize the noiseless desired signal
figure;
plot(clean_signal);
title('Desired Signal (Noiseless)');
grid on;

%% Algorithm Execution

% 1. RLS with fixed forgetting factor (lambda < 1)
disp('Running RLS (lambda=0.95)...');
[h_classic] = classic_rls(input_matrix, observed_signal, P_delta, lambda_rls);
mis_classic = misalignment(h_classic, true_weights);

% 2. Classic RLS with infinite memory (lambda = 1)
disp('Running RLS (lambda=1)...');
[h_plain] = classic_rls(input_matrix, observed_signal, P_delta, 1);
mis_plain = misalignment(h_plain, true_weights);

% 3. Variable Forgetting Factor RLS (VFF-RLS)
disp('Running VFF-RLS...');
% Note: VFF-RLS here uses the clean signal for optimal parameter tuning/oracle comparison
[h_vff, lambda_vff] = vff_rls(input_matrix, observed_signal, clean_signal, P_delta, 2, 10);
mis_vff = misalignment(h_vff, true_weights);

% 4. Gradient VFF RLS (GVFF-RLS)
disp('Running GVFF-RLS...');
[h_gvff, lambda_gvff] = gvff_rls(input_matrix, observed_signal, P_delta, 2, 10);
mis_gvff = misalignment(h_gvff, true_weights);

%% Plotting Results
disp('Generating plots...');
plot_true_parameters(true_weights, output_dir);

for i = 1:L
    plot_parameter(i, true_weights, h_classic, h_vff, h_gvff, with_gvff, output_dir, h_plain, lambda_rls)
end
plot_mis_and_lambda(mis_classic, mis_vff, mis_gvff, ones(1,N)*lambda_rls, lambda_vff, lambda_gvff, with_gvff, output_dir, mis_plain, ones(1,N))

%% Helper Functions

function x = course_input(N, L)
% COURSE_INPUT Generates the input signal matrix.
%
% Inputs:
%   N - Number of samples.
%   L - Number of input channels (taps).
%
% Outputs:
%   x - Input matrix (N x L) with independent Gaussian entries.

% Variance of the input signal
sigma_u = sqrt(100);
x = randn(N, L) * sigma_u;
end

function h = course_theta(N, L)
% COURSE_THETA Generates the true time-varying system parameters.
%
% Inputs:
%   N - Number of samples.
%   L - Number of coefficients.
%
% Outputs:
%   h - Matrix of true weights (N x L).

h = zeros(N, L);
t = (1:N)';

% Coefficient 1: STEP CHANGE
% Tests the algorithm's ability to react to abrupt changes.
% Value changes from 1 to -0.5 at N/2.
h(:,1) = 1;
h(floor(N/2):end, 1) = -0.5;

% Coefficient 2: SINE WAVE
% Tests the algorithm's ability to track continuous variations.
% Sinusoidal variation starts at sample 200.
h(:,2) = 0.5;
start_sin = 200;
h(start_sin:end, 2) = 0.5 + 0.3 * sin(2*pi*0.005*t(start_sin:end));

% Coefficient 3: LINEAR RAMP (Commented out)
% Tests tracking of a constant drift.
% h(:,3) = linspace(-0.5, 0.5, N)';

% Coefficient 4: CONSTANT (Commented out)
% Tests stability.
% h(:,4) = -0.8;
end
