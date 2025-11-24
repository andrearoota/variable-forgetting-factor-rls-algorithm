function plot_mis_and_lambda(mis_classic,mis_vff,mis_gvff,lambda_classic,lambda_vff,lambda_gvff, enable_gvff, output_dir, mis_plain, lambda_plain)
is_gvff = "without_gvff";
% Extract lambda value for legend
lambda_val = lambda_classic(1);
legend_label = sprintf('RLS $\\lambda = %.3g$', lambda_val);

%% Plot 1: Misalignment
figure('Position', [100, 100, 800, 400]);
plot(mis_classic, 'b', 'DisplayName', legend_label, 'LineWidth', 1.5);
hold on;
plot(mis_vff, 'r', 'DisplayName', 'VFF-RLS', 'LineWidth', 1.5);
if enable_gvff
    plot(mis_gvff, 'g', 'DisplayName', 'GVFF-RLS', 'LineWidth', 1.5);
    is_gvff = "with_gvff";
end
if nargin > 8 && ~isempty(mis_plain)
    plot(mis_plain, 'k--', 'DisplayName', 'RLS $\lambda = 1$', 'LineWidth', 1.5);
end
xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Misalignment [dB]', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);
ylim([-60, 10])
xlim([0, length(mis_classic)])

exportgraphics(gcf, 'output/'+output_dir+'/misalignment_'+is_gvff+'.pdf', 'Resolution',300, 'BackgroundColor', 'none', 'ContentType', 'vector');

%% Plot 2: Lambda
figure('Position', [100, 100, 800, 400]);
plot(lambda_classic,'b', 'DisplayName', legend_label, 'LineWidth', 1.5);
hold on
plot(lambda_vff,'r', 'DisplayName', 'VFF-RLS', 'LineWidth', 1.5);
if enable_gvff
    plot(lambda_gvff,'g', 'DisplayName', 'GVFF-RLS', 'LineWidth', 1.5);
end
if nargin > 9 && ~isempty(lambda_plain)
    plot(lambda_plain, 'k--', 'DisplayName', 'RLS $\lambda = 1$', 'LineWidth', 1.5);
end
xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 16);
ylim([-0.05, 1.05])
xlim([0, length(lambda_classic)])

exportgraphics(gcf, 'output/'+output_dir+'/lambda_'+is_gvff+'.pdf', 'Resolution',300, 'BackgroundColor', 'none', 'ContentType', 'vector');
end

