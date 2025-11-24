function plot_true_parameters(h, output_dir)
[N, L] = size(h);
for i = 1:L
    figure('Position', [100, 100, 800, 400]);
    plot(h(:,i), 'Color', 'black', 'LineWidth', 2);
    % title removed for presentation
    xlabel('Iterations', 'Interpreter', 'latex', 'FontSize', 18);
    grid on;
    set(gca, 'FontSize', 16);
    ylim([-1, 1.5]);
    xlim([0, N]);

    % Save with transparent background
    exportgraphics(gcf, 'output/'+output_dir+'/true_theta_'+string(i)+'.pdf', 'Resolution', 300, 'BackgroundColor', 'none', 'ContentType', 'vector');
end
end
