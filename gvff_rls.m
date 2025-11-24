function [weights_history, lambda_history] = gvff_rls(input_matrix, desired_signal, delta, K_alpha, K_beta)
% GVFF_RLS Implements Gradient Variable Forgetting Factor RLS.
%
% Inputs:
%   input_matrix   - Input signal matrix (N x L).
%   desired_signal - Desired signal vector (N x 1).
%   delta          - Regularization parameter.
%   K_alpha, K_beta - Parameters for window lengths.
%
% Outputs:
%   weights_history - History of estimated weights.
%   lambda_history  - History of forgetting factors.

[N, num_taps] = size(input_matrix);
weights = zeros(num_taps, 1);              % Initial weights
inv_corr_matrix = delta * eye(num_taps);   % Inverse correlation matrix
mu = 0.1;                                  % Gradient step size
lambda_max = 0.999999;
lambda = lambda_max;
rho = 1;
rho_tilde = 1;
sigma_e = 0;  % Power of the a priori error signal
sigma_v = 0;
sigma2_e_partial = 0;

alpha = 1 - (1/(K_alpha * num_taps));
beta = 1 - (1/(K_beta * num_taps));

weights_history = zeros(N, num_taps);
lambda_history = zeros(N, 1);
% history = zeros(N, 1); % Unused output? It was storing zeta.

for n = 1:N
    input_vector = input_matrix(n, :)';          % Input vector (L x 1)
    desired_sample = desired_signal(n);          % Desired signal sample

    % Error
    error_signal = desired_sample - weights' * input_vector;

    % Gradient-based Lambda Update Logic
    % Calculate discriminant D
    D = (num_taps - 2)^2 * rho^2  - 8 * (num_taps + 2) * ((num_taps + 1) * rho_tilde + rho^2);

    if(D > 0)
        lambda_min = ((num_taps - 2) * rho + sqrt(D)) / (4 * ((num_taps + 1) * rho_tilde + rho^2));
    else
        lambda_min = 0.000001;
    end

    rho = 1 + lambda * rho;
    rho_partial = rho; % This seems to be just rho? Wait, rho is updated before this line.

    rho_tilde = 1 + lambda^2 * rho_tilde;
    rho_tilde_partial = 2 * lambda * rho_tilde;

    sigma2_e = alpha * sigma_e^2 + (1 - alpha) * error_signal^2;
    sigma_e = sqrt(sigma2_e);
    sigma2_v = beta * sigma_v^2 + (1 - beta) * error_signal^2;
    sigma_v = sqrt(sigma2_v);

    zeta = 1 - ((2 * ((num_taps + 1) * rho_tilde + rho^2) - (num_taps + 2) * rho) / (rho * ((num_taps + 1) * rho_tilde + rho^2)));

    temp1 = 2 * rho_partial / rho^2;
    temp2 = ((num_taps + 1) * rho_tilde + rho^2)^2;
    temp3 = ((num_taps + 1) * rho_tilde_partial + 2 * rho * rho_partial);

    zeta_partial = temp1 - ((num_taps + 2) / temp2 * temp3);
    h_partial = temp1 - (2 / temp2 * temp3);

    sigma2_e_partial = zeta * sigma2_e_partial + zeta_partial * sigma_e + h_partial * sigma_v;

    % Update lambda
    lambda = lambda - (mu * sigma2_e_partial / (1 - lambda));
    lambda = max(2 * lambda_min, min(lambda_max, lambda));


    % Gain
    pi_vector = inv_corr_matrix * input_vector;
    denom = lambda + input_vector' * pi_vector;
    gain_vector = pi_vector / denom;

    % Update weights
    weights = weights + gain_vector * error_signal;

    % Update inverse correlation matrix P
    inv_corr_matrix = (inv_corr_matrix - gain_vector * input_vector' * inv_corr_matrix) / lambda;

    % Save history
    % history(n) = zeta;
    weights_history(n, :) = weights;
    lambda_history(n) = lambda;
end
end