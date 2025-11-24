function plot_parameter(idx,real,est_classic,est_vff,est_gvff,enable_gvff, output_dir, est_plain, lambda_val)
is_gvff = "without_gvff";

if nargin < 9 || isempty(lambda_val)
    lambda_val = 0.95; % Default fallback
end
legend_label = sprintf('RLS $\\lambda = %.3g$', lambda_val);

figure('Position', [100, 100, 800, 400]);
plot(real(:,idx), 'Color','black', 'DisplayName', ['True $\theta_' num2str(idx) '$'], 'LineWidth', 2);
hold on;
plot(est_classic(:,idx), 'b', 'DisplayName', ['Theta ' num2str(idx) ' ' legend_label], 'LineWidth', 1.5);
plot(est_vff(:,idx), 'r', 'DisplayName', 'Theta '+string(idx)+' VFF', 'LineWidth', 1.5);
if enable_gvff
    plot(est_gvff(:,idx), 'g', 'DisplayName', 'GVFF-RLS', 'LineWidth', 1.5);
    is_gvff = "with_gvff";
end
if nargin > 7 && ~isempty(est_plain)
    plot(est_plain(:,idx), 'k--', 'DisplayName', 'RLS $\lambda = 1$', 'LineWidth', 1.5);
end
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 14);
% title removed
xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 18);
grid on
set(gca, 'FontSize', 16);
ylim([-1, 1.5])
xlim([0, size(real,1)])

exportgraphics(gcf, 'output/'+output_dir+'/theta_'+string(idx)+'_'+is_gvff+'.pdf', 'Resolution',300, 'BackgroundColor', 'none', 'ContentType', 'vector');
end