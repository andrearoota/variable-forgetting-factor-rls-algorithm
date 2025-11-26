function [weights_history] = classic_rls(input_matrix, desired_signal, delta, lambda)
% CLASSIC_RLS Implements the Classic Recursive Least Squares (RLS) algorithm.
%
% This is the standard RLS algorithm with a fixed forgetting factor.
% It provides exponentially weighted least squares estimation, giving
% more weight to recent data while gradually forgetting older measurements.
%
% Inputs:
%   input_matrix   - Input signal matrix (N x L), where N = samples, L = taps.
%   desired_signal - Desired signal vector (N x 1) with measurement noise.
%   delta          - Regularization parameter for initializing P (e.g., 10^6).
%                    Large values give more initial uncertainty.
%   lambda         - Forgetting factor (0 < lambda <= 1).
%                    lambda = 1: infinite memory (standard least squares)
%                    lambda < 1: exponential forgetting (tracks time-varying systems)
%
% Outputs:
%   weights_history - Estimated filter weights over time (N x L).
%
% Algorithm:
%   Uses the matrix inversion lemma for efficient recursive updates of
%   the inverse correlation matrix P. The forgetting factor λ controls
%   the effective memory length: smaller λ = faster tracking but more noise.
%
% Reference:
%   Haykin, S. (2002). Adaptive Filter Theory (4th ed.). Prentice Hall.

[N, L] = size(input_matrix);

% Initialize inverse correlation matrix P = δ·I
% P represents the inverse of the weighted input autocorrelation matrix
inv_corr_matrix = delta * eye(L);

% Initialize filter weights to zero
weights = zeros(L, 1);

% Pre-allocate storage arrays for efficiency
weights_history = zeros(N, L);  % Weight evolution over time
output_signal = zeros(N, 1);    % Filter output for MSE calculation

% RLS Main Loop
for n = 1:N
    % Extract current input vector (regressor)
    input_vector = input_matrix(n, :)';
    
    % === RLS Update ===
    % Compute intermediate vector: P·x
    p_times_x = inv_corr_matrix * input_vector;
    
    % Compute normalization factor: λ + x'·P·x
    normalization = lambda + input_vector' * p_times_x;
    
    % Update inverse correlation matrix using matrix inversion lemma:
    % P(n) = (1/λ)·(P(n-1) - (P·x·x'·P) / (λ + x'·P·x))
    inv_corr_matrix = (1/lambda) * (inv_corr_matrix - (p_times_x * p_times_x') / normalization);
    
    % Compute Kalman gain: k = P·x
    gain_vector = inv_corr_matrix * input_vector;
    
    % Compute a priori filter output
    output_signal(n) = input_vector' * weights;
    
    % Compute a priori prediction error
    error_signal = desired_signal(n) - output_signal(n);
    
    % Update filter weights: w = w + k·e
    weights = weights + gain_vector * error_signal;
    
    % Store weights for this iteration
    weights_history(n, :) = weights.';
end

% Display final Mean Squared Error
disp("MSE:")
disp(mean((output_signal - desired_signal).^2))
end
