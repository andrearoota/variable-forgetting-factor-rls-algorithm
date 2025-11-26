function [weights_history, lambda_history] = vff_rls(input_matrix, desired_signal, desired_signal_clean, delta, K_alpha, K_beta, gamma_threshold, epsilon_small, lambda_max)
% VFF_RLS Implements Variable Forgetting Factor RLS.
%
% Inputs:
%   input_matrix         - Input signal matrix (N x L).
%   desired_signal       - Desired signal vector (N x 1) (noisy).
%   desired_signal_clean - Clean desired signal vector (N x 1) (for initialization).
%   delta                - Regularization parameter.
%   K_alpha, K_beta      - Parameters for window lengths.
%   gamma_threshold      - Threshold for forgetting factor adaptation (optional, default: 1.5).
%   epsilon_small        - Small constant to avoid division by zero (optional, default: 1e-8).
%   lambda_max           - Maximum forgetting factor (optional, default: 0.999999).
%
% Outputs:
%   weights_history - History of estimated weights.
%   lambda_history  - History of forgetting factors.

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

[N, num_taps] = size(input_matrix);
weights = zeros(num_taps, 1);                   % Adaptive filter weights
inv_corr_matrix = eye(num_taps) * delta;        % Inverse of the input auto-correlation matrix
lambda = lambda_max;                            % Variable Forgetting Factor

alpha = 1 - 1/(K_alpha * num_taps);
beta = 1 - 1/(K_beta * num_taps);

% Initialization of power estimates
% Note: q calculation here seems to be a placeholder or specific initialization
q_init = input_matrix(2,:) * inv_corr_matrix * input_matrix(2,:)';
sigma_q = sqrt(var(q_init));                % prediction error associated with the input (starts at 0 if scalar)
sigma_e = sqrt(mean(desired_signal_clean.^2));  % power of the a priori error signal
sigma_v = 1; % Initial noise power estimate

% Storage for tracking
weights_history = zeros(N, num_taps);
lambda_history = zeros(N, 1);
output_signal = zeros(N, 1);

% RLS with Variable Forgetting Factor
for n = 1:N
    input_vector = input_matrix(n, :)';  % Input vector

    output_signal(n) = weights' * input_vector;
    error_signal = desired_signal(n) - output_signal(n); % A priori error

    % Kalman gain vector
    % k_n = (P * x_n) / (lambda + x_n' * P * x_n)
    pi_vector = inv_corr_matrix * input_vector;
    denom = lambda + input_vector' * pi_vector;
    gain_vector = pi_vector / denom;

    weights = weights + gain_vector * error_signal;

    q_n = input_vector' * pi_vector; % x_n' * P * x_n

    % Update inverse correlation matrix
    % P = (P - k_n * x_n' * P) / lambda
    inv_corr_matrix = (inv_corr_matrix - (gain_vector * input_vector' * inv_corr_matrix)) / lambda;

    % Update forgetting factor
    sigma2_e = alpha * sigma_e^2 + (1 - alpha) * error_signal^2;
    sigma_e = sqrt(sigma2_e);

    sigma2_q = alpha * sigma_q^2 + (1 - alpha) * q_n^2;
    sigma_q = sqrt(sigma2_q);

    sigma2_v = beta * sigma_v^2 + (1 - beta) * error_signal^2;
    sigma_v = sqrt(sigma2_v);

    if(sigma_e > gamma_threshold * sigma_v)
        lambda = min((sigma_q * sigma_v) / (epsilon_small + abs(sigma_e - sigma_v)), lambda_max);
    else
        lambda = lambda_max;
    end

    % Store metrics
    weights_history(n, :) = weights;
    lambda_history(n) = lambda;
end
disp("MSE")
disp(mean((output_signal - desired_signal).^2))
end