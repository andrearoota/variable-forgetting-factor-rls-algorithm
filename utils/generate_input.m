function input_matrix = generate_input(N, L, sigma_u)
% GENERATE_INPUT Generates the input signal matrix.
%
% Inputs:
%   N       - Number of samples.
%   L       - Number of input channels (taps).
%   sigma_u - Standard deviation of the input signal (optional, default: sqrt(100)).
%
% Outputs:
%   input_matrix - Input matrix (N x L) with independent Gaussian entries.

if nargin < 3
    sigma_u = sqrt(100); % Default variance
end

input_matrix = randn(N, L) * sigma_u;
end
