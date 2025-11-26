function output_signal = generate_output(weights, input_matrix, noise)
% GENERATE_OUTPUT Generates the output signal of a time-varying system.
%
% Inputs:
%   weights      - Time-varying filter weights (N x L).
%   input_matrix - Input signal matrix (N x L).
%   noise        - Additive noise vector (N x 1).
%
% Outputs:
%   output_signal - Generated output signal (N x 1).

N = size(noise, 1);
output_signal = zeros(N, 1);

for n = 1:N
    output_signal(n) = input_matrix(n, :) * weights(n, :)' + noise(n);
end
end