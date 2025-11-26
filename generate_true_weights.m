function weights = generate_true_weights(N, L, start_sin)
% GENERATE_TRUE_WEIGHTS Generates the true time-varying system parameters.
%
% Inputs:
%   N         - Number of samples.
%   L         - Number of coefficients.
%   start_sin - Sample index where sinusoidal variation starts (optional, default: 200).
%
% Outputs:
%   weights - Matrix of true weights (N x L).

if nargin < 3
    start_sin = 200; % Default starting point for sine wave
end

weights = zeros(N, L);
t = (1:N)';

% Coefficient 1: STEP CHANGE
% Tests the algorithm's ability to react to abrupt changes.
% Value changes from 1 to -0.5 at N/2.
weights(:,1) = 1;
weights(floor(N/2):end, 1) = -0.5;

% Coefficient 2: SINE WAVE
% Tests the algorithm's ability to track continuous variations.
% Sinusoidal variation starts at the specified sample.
weights(:,2) = 0.5;
weights(start_sin:end, 2) = 0.5 + 0.3 * sin(2*pi*0.005*t(start_sin:end));

end
