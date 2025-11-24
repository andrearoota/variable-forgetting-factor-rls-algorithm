# Variable Forgetting Factor RLS Algorithm

*Adaptive Learning, Estimation and Supervision of Dynamical Systems (ALES) for the Master's Degree Course in Computer Engineering at the University of Bergamo.*

A MATLAB implementation and comparison of Recursive Least Squares (RLS) algorithms with variable forgetting factors for adaptive system identification.

## Overview

This project implements and compares four different RLS algorithms for tracking time-varying systems:

1. **Classic RLS** - Standard RLS with fixed forgetting factor (λ = 0.95)
2. **RLS with Infinite Memory** - RLS without forgetting (λ = 1)
3. **VFF-RLS** - Variable Forgetting Factor RLS
4. **GVFF-RLS** - Gradient-based Variable Forgetting Factor RLS

The algorithms are tested on a time-varying system with both abrupt changes (step) and continuous variations (sinusoidal) to evaluate their tracking capabilities.

## Requirements

- MATLAB (tested on R2025b)

## Usage

### Basic Simulation

Simply run the main script:

```matlab
main
```

This will:
1. Generate a time-varying system with N=1000 samples
2. Run all four RLS algorithms
3. Generate comparison plots in the `output/{output_dir}/` directory
4. Display MSE (Mean Squared Error) for each algorithm

### Configuration

You can modify simulation parameters in `main.m`:

```matlab
P_delta = 10^6;      % Regularization parameter
N = 1000;            % Number of samples
L = 2;               % Number of filter taps
lambda_rls = 0.95;   % Fixed forgetting factor
sigma_v = 1;         % Noise variance

with_gvff = false;    % Flag to enable/disable GVFF-RLS plotting if needed
output_dir = "main"; % Name of the experiment for output file naming
```

### System Definition

The time-varying system is defined in `course_theta()` with:
- **Coefficient 1**: Step change from 1 to -0.5 at N/2
- **Coefficient 2**: Sinusoidal variation: 0.5 + 0.3·sin(2π·0.005·t) starting at sample 200

## Output

The simulation generates plots showing:
- True parameter evolution
- Estimated parameters from each algorithm
- Misalignment (tracking error) over time
- Forgetting factor λ evolution for adaptive methods

All plots are saved to `output/<output_dir>/` directory.

## References

This implementation is based on the following paper:

**Paleologu, C., Benesty, J., & Ciochină, S.** (2008). A Robust Variable Forgetting Factor Recursive Least-Squares Algorithm for System Identification. *IEEE Signal Processing Letters*, 15, 597-600.  
https://doi.org/10.1109/LSP.2008.2001559

The work applies adaptive filtering theory and variable forgetting factor techniques for non-stationary signal processing.