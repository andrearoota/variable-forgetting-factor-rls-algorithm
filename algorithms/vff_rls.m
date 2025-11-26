function [weights_history, lambda_history] = vff_rls(input_matrix, desired_signal, desired_signal_clean, delta, K_alpha, K_beta, gamma_threshold, epsilon_small, lambda_max)
% VFF_RLS Implements Variable Forgetting Factor RLS.
%
% This algorithm adaptively adjusts the forgetting factor based on the
% prediction error power, input signal correlation, and noise variance.
% It balances tracking speed and noise sensitivity for time-varying systems.
%
% Inputs:
%   input_matrix         - Input signal matrix (N x L), where N = samples, L = taps.
%   desired_signal       - Desired signal vector (N x 1) with measurement noise.
%   desired_signal_clean - Clean desired signal vector (N x 1) for initialization.
%                          NOTE: This is oracle knowledge not available in practice.
%   delta                - Regularization parameter for initializing P (e.g., 10^6).
%   K_alpha              - Window length parameter for error power estimation.
%   K_beta               - Window length parameter for noise variance estimation.
%   gamma_threshold      - Threshold for λ adaptation (optional, default: 1.5).
%   epsilon_small        - Small constant to avoid division by zero (optional, default: 1e-8).
%   lambda_max           - Maximum forgetting factor (optional, default: 0.999999).
%
% Outputs:
%   weights_history - Estimated filter weights over time (N x L).
%   lambda_history  - Forgetting factor evolution over time (N x 1).
%
% Algorithm:
%   The forgetting factor λ is adapted based on:
%     - σ_e: Prediction error power (exponentially weighted)
%     - σ_q: Input signal correlation measure
%     - σ_v: Noise variance estimate
%   When σ_e > γ·σ_v (high tracking error), λ decreases for faster adaptation.
%   Otherwise, λ = λ_max for stability.
%
% Reference:
%   Paleologu, C., Benesty, J., & Ciochină, S. (2008). A Robust Variable
%   Forgetting Factor Recursive Least-Squares Algorithm for System Identification.
%   IEEE Signal Processing Letters, 15, 597-600.

% Set default values for optional parameters
if nargin < 7 || isempty(gamma_threshold)
    gamma_threshold = 1.5;
end
if nargin < 8 || isempty(epsilon_small)
    epsilon_small = 1e-8;
end
if nargin < 9 || isempty(lambda_max)
    lambda_max = 0.999999;
end

[N, L] = size(input_matrix);

% Initialize adaptive filter weights
weights = zeros(L, 1);

% Initialize inverse correlation matrix P = δ·I
% P tracks the inverse of the input signal autocorrelation
inv_corr_matrix = eye(L) * delta;

% Initialize forgetting factor at maximum (most conservative)
lambda = lambda_max;

% Compute smoothing coefficients for exponential averaging
% Controls how quickly estimates adapt to new observations
alpha_error = 1 - 1/(K_alpha * L);  % Error power smoothing
beta_noise = 1 - 1/(K_beta * L);    % Noise variance smoothing

% Initialize power estimates for the adaptation mechanism
% correlation_init: Initial correlation metric from second sample
correlation_init = input_matrix(2,:) * inv_corr_matrix * input_matrix(2,:)';
sigma_correlation = sqrt(var(correlation_init));  % Std dev of input correlation

% σ_e: Prediction error power (initialized using clean signal energy)
sigma_error = sqrt(mean(desired_signal_clean.^2));

% σ_v: Noise variance estimate (initialized to known value)
sigma_noise = 1;

% Pre-allocate storage arrays
weights_history = zeros(N, L);  % Weight evolution over time
lambda_history = zeros(N, 1);   % Forgetting factor evolution
output_signal = zeros(N, 1);    % Filter output for MSE calculation

% RLS with Variable Forgetting Factor - Main Loop
for n = 1:N
    % Extract current input vector (regressor)
    input_vector = input_matrix(n, :)';

    % Compute a priori filter output
    output_signal(n) = weights' * input_vector;

    % Compute a priori prediction error
    error_signal = desired_signal(n) - output_signal(n);

    % === RLS Update ===
    % Compute intermediate vector: P·x
    p_times_x = inv_corr_matrix * input_vector;

    % Compute gain denominator: λ + x'·P·x
    gain_denominator = lambda + input_vector' * p_times_x;

    % Compute Kalman gain: k = (P·x) / (λ + x'·P·x)
    gain_vector = p_times_x / gain_denominator;

    % Update filter weights: w = w + k·e
    weights = weights + gain_vector * error_signal;

    % Compute correlation metric: x'·P·x (used for λ adaptation)
    correlation_metric = input_vector' * p_times_x;

    % Update inverse correlation matrix using matrix inversion lemma:
    % P = (1/λ)·(P - k·x'·P)
    inv_corr_matrix = (inv_corr_matrix - (gain_vector * input_vector' * inv_corr_matrix)) / lambda;

    % === Forgetting Factor Adaptation ===
    % Update prediction error power using exponential smoothing
    sigma_error_squared = alpha_error * sigma_error^2 + (1 - alpha_error) * error_signal^2;
    sigma_error = sqrt(sigma_error_squared);

    % Update input correlation metric using exponential smoothing
    sigma_correlation_squared = alpha_error * sigma_correlation^2 + (1 - alpha_error) * correlation_metric^2;
    sigma_correlation = sqrt(sigma_correlation_squared);

    % Update noise variance estimate using exponential smoothing
    sigma_noise_squared = beta_noise * sigma_noise^2 + (1 - beta_noise) * error_signal^2;
    sigma_noise = sqrt(sigma_noise_squared);

    % Adapt forgetting factor based on prediction error
    if sigma_error > gamma_threshold * sigma_noise
        % High error detected: reduce λ for faster tracking
        % Formula: λ = min(σ_correlation·σ_noise / |σ_error - σ_noise|, λ_max)
        lambda = min((sigma_correlation * sigma_noise) / (epsilon_small + abs(sigma_error - sigma_noise)), lambda_max);
    else
        % Low error: use maximum λ for stability
        lambda = lambda_max;
    end

    % Store results for this iteration
    weights_history(n, :) = weights;
    lambda_history(n) = lambda;
end

% Display final Mean Squared Error
disp("MSE:")
disp(mean((output_signal - desired_signal).^2))
end