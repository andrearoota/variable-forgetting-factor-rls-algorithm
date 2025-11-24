function [weights_history] = classic_rls(input_matrix, desired_signal, delta, lambda)
% CLASSIC_RLS Implements the Classic Recursive Least Squares (RLS) algorithm.
%
% Inputs:
%   input_matrix   - Input signal matrix (N x L), where N is the number of samples
%                    and L is the number of filter coefficients (taps).
%   desired_signal - Desired signal vector (N x 1).
%   delta          - Regularization parameter for initializing the inverse correlation matrix P.
%   lambda         - Forgetting factor (0 < lambda <= 1).
%
% Outputs:
%   weights_history - History of estimated filter weights (N x L).

[N, num_taps] = size(input_matrix);

% Initialize inverse correlation matrix P
% P is initialized as delta * I
inv_corr_matrix = delta * eye(num_taps);

% Initialize filter weights
weights = zeros(num_taps, 1);

% Pre-allocate memory for history and output
weights_history = zeros(N, num_taps);
output_signal = zeros(N, 1);

for n = 1:N
    % Extract current input vector (regressor)
    input_vector = input_matrix(n, :)';

    % Calculate intermediate term for P update
    % Using the original implementation logic where P is updated before k
    % beta = lambda + phi'*P*phi
    pi_vector = inv_corr_matrix * input_vector;
    beta = lambda + input_vector' * pi_vector;

    % Update inverse correlation matrix P
    % P(n) = (1/lambda) * (P(n-1) - (P(n-1)*u(n)*u(n)'*P(n-1)) / beta)
    inv_corr_matrix = (1/lambda) * (inv_corr_matrix - (pi_vector * pi_vector') / beta);

    % Calculate Kalman gain
    gain_vector = inv_corr_matrix * input_vector;

    % Calculate filter output (a priori estimation)
    output_signal(n) = input_vector' * weights;

    % Calculate a priori error
    error_signal = desired_signal(n) - output_signal(n);

    % Update filter weights
    weights = weights + gain_vector * error_signal;

    % Store weights history
    weights_history(n, :) = weights.';
end

disp("MSE")
disp(mean((output_signal - desired_signal).^2))
end

