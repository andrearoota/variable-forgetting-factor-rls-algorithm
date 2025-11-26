function generate_plots(true_weights, weights_classic, weights_vff, mis_classic, mis_vff, lambda_classic, lambda_vff, weights_plain, mis_plain, lambda_plain)
% GENERATE_PLOTS Creates all visualization plots for RLS algorithm comparison.
%
% This function generates comprehensive plots including:
%   - True parameter evolution for each coefficient
%   - Estimated vs true parameters for each coefficient
%   - Misalignment comparison across algorithms
%   - Forgetting factor evolution over time
%
% Inputs:
%   true_weights   - True parameter matrix (N x L).
%   weights_classic - Estimated parameters from classic RLS (N x L).
%   weights_vff    - Estimated parameters from VFF-RLS (N x L).
%   mis_classic    - Misalignment vector for classic RLS (N x 1) in dB.
%   mis_vff        - Misalignment vector for VFF-RLS (N x 1) in dB.
%   lambda_classic - Forgetting factor for classic RLS (N x 1).
%   lambda_vff     - Forgetting factor for VFF-RLS (N x 1).
%   weights_plain  - (Optional) Estimated parameters from RLS with λ=1 (N x L).
%   mis_plain      - (Optional) Misalignment for RLS with λ=1 (N x 1) in dB.
%   lambda_plain   - (Optional) Forgetting factor for RLS with λ=1 (N x 1).
%
% Output:
%   Multiple PDF figures saved to 'figures/' directory:
%     - true_theta_<i>.pdf: True parameter evolution
%     - theta_<i>.pdf: Parameter tracking comparison
%     - misalignment.pdf: Tracking error comparison
%     - lambda.pdf: Forgetting factor evolution

[N, L] = size(true_weights);

% Extract lambda value for legend
lambda_val = lambda_classic(1);
legend_label = sprintf('RLS $\\lambda = %.3g$', lambda_val);

% Check for optional plain RLS arguments
has_plain = nargin >= 8 && ~isempty(weights_plain);

%% Plot 1: True Parameters
for i = 1:L
    figure('Position', [100, 100, 800, 400]);
    plot(true_weights(:,i), 'Color', 'black', 'LineWidth', 2);
    xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
    grid on;
    set(gca, 'FontSize', 16);
    ylim([-1, 1.5]);
    xlim([0, N]);

    exportgraphics(gcf, 'figures/true_theta_'+string(i)+'.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');
    close(gcf);
end

%% Plot 2: Parameter Tracking
for i = 1:L
    figure('Position', [100, 100, 800, 400]);
    plot(true_weights(:,i), 'Color','black', 'DisplayName', ['True $\theta_' num2str(i) '$'], 'LineWidth', 2);
    hold on;
    plot(weights_classic(:,i), 'b', 'DisplayName', ['$\theta_' num2str(i) '$ ' legend_label], 'LineWidth', 1.5);
    plot(weights_vff(:,i), 'r', 'DisplayName', '$\theta_'+string(i)+'$ VFF', 'LineWidth', 1.5);

    if has_plain
        plot(weights_plain(:,i), 'k--', 'DisplayName', 'RLS $\lambda = 1$', 'LineWidth', 1.5);
    end

    legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 14);
    xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
    grid on;
    set(gca, 'FontSize', 16);
    ylim([-1, 1.5]);
    xlim([0, N]);

    exportgraphics(gcf, 'figures/theta_'+string(i)+'.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');
    close(gcf);
end

%% Plot 3: Misalignment
figure('Position', [100, 100, 800, 400]);
plot(mis_classic, 'b', 'DisplayName', legend_label, 'LineWidth', 1.5);
hold on;
plot(mis_vff, 'r', 'DisplayName', 'VFF-RLS', 'LineWidth', 1.5);

if has_plain && ~isempty(mis_plain)
    plot(mis_plain, 'k--', 'DisplayName', 'RLS $\lambda = 1$', 'LineWidth', 1.5);
end

xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Misalignment [dB]', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);
ylim([-60, 10]);
xlim([0, N]);

exportgraphics(gcf, 'figures/misalignment.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');
close(gcf);

%% Plot 4: Lambda Evolution
figure('Position', [100, 100, 800, 400]);
plot(lambda_classic, 'b', 'DisplayName', legend_label, 'LineWidth', 1.5);
hold on;
plot(lambda_vff, 'r', 'DisplayName', 'VFF-RLS', 'LineWidth', 1.5);

if has_plain && ~isempty(lambda_plain)
    plot(lambda_plain, 'k--', 'DisplayName', 'RLS $\lambda = 1$', 'LineWidth', 1.5);
end

xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);
ylim([-0.05, 1.05]);
xlim([0, N]);

exportgraphics(gcf, 'figures/lambda.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector');
close(gcf);

end
