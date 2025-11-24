function misalignment_metric = misalignment(estimated_weights, true_weights)
% MISALIGNMENT Computes the normalized misalignment in dB.
%
% Inputs:
%   estimated_weights - Estimated weights matrix (N x L).
%   true_weights      - True system weights matrix (N x L).
%
% Outputs:
%   misalignment_metric - Misalignment vector (N x 1) in dB.

arguments
    estimated_weights (:,:) double
    true_weights (:,:) double
end

N = size(estimated_weights, 1);
misalignment_metric = zeros(N, 1);

for n = 1:N
    norm_error = norm(true_weights(n, :) - estimated_weights(n, :));
    norm_true = norm(true_weights(n, :));

    misalignment_metric(n) = 20 * log10(norm_error / norm_true);
end

disp("Mean Misalignment [dB]:")
disp(mean(misalignment_metric))
end