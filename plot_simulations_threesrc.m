function plot_simulations_threesrc(r, var, detect)

    c = lines(8);
    mc = size(r, 3);

    figure
    subplot(4, 4, 1)
    histogram(r(1, 1, :), 'BinWidth', 0.002, 'FaceColor', c(1, :), ...
        'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
    hold on
    grid on
    plot([median(r(1, 1, :)) median(r(1, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.06])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    ylabel('Number of simulations')
    set(gca, 'FontSize', 18)
    title('Estimation bias')

    subplot(4, 4, 2)
    histogram(var(1, 1, :), 'BinWidth', 0.001, 'FaceColor', c(1, :), ...
        'FaceAlpha', 0.8, 'LineWidth', 2)
    hold on
    grid on
    plot([median(var(1, 1, :)) median(var(1, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.02])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    set(gca, 'FontSize', 18)
    title('Point spreading')

    subplot(4, 4, 5)
    histogram(r(4, 1, :), 'BinWidth', 0.002, 'FaceColor', c(6, :), ...
        'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
    hold on
    grid on
    plot([median(r(4, 1, :)) median(r(4, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.06])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    ylabel('Number of simulations')
    set(gca, 'FontSize', 18)

    subplot(4, 4, 6)
    histogram(var(4, 1, :), 'BinWidth', 0.001, 'FaceColor', c(6, :), ...
        'FaceAlpha', 0.8, 'LineWidth', 2)
    hold on
    grid on
    plot([median(var(4, 1, :)) median(var(4, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.02])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    set(gca, 'FontSize', 18)

    subplot(4, 4, 9)
    histogram(r(2, 1, :), 'BinWidth', 0.002, 'FaceColor', c(2, :), ...
        'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
    hold on
    grid on
    plot([median(r(2, 1, :)) median(r(2, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.06])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    ylabel('Number of simulations')
    set(gca, 'FontSize', 18)

    subplot(4, 4, 10)
    histogram(var(2, 1, :), 'BinWidth', 0.001, 'FaceColor', c(2, :), ...
        'FaceAlpha', 0.8, 'LineWidth', 2)
    hold on
    grid on
    plot([median(var(2, 1, :)) median(var(2, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.02])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    set(gca, 'FontSize', 18)

    subplot(4, 4, 13)
    histogram(r(3, 1, :), 'BinWidth', 0.002, 'FaceColor', c(4, :), ...
        'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
    hold on
    grid on
    plot([median(r(3, 1, :)) median(r(3, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.06])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    ylabel('Number of simulations')
    set(gca, 'FontSize', 18)

    subplot(4, 4, 14)
    histogram(var(3, 1, :), 'BinWidth', 0.001, 'FaceColor', c(4, :), ...
        'FaceAlpha', 0.8, 'LineWidth', 2)
    hold on
    grid on
    plot([median(var(3, 1, :)) median(var(3, 1, :))], [0 mc / 4], '--k', 'LineWidth', 2)
    xlim([0 0.02])
    ylim([0 mc / 4])
    xlabel('Distance, m')
    set(gca, 'FontSize', 18)

    for i = 1:4
        for j = 0:3
            ratio_detected(i, j+1) = sum((detect(i, 1, :) == j)) / mc * 100;
        end
        ratio_detected(i, :) = fliplr(ratio_detected(i, :));
        s(i, :) = cumsum(ratio_detected(i, :));
    end

    subplot(4, 4, [11, 12, 15, 16])
    plot(s(1, 1:3), '-o', 'Color', c(1, :), 'LineWidth', 3, 'MarkerSize', 10, ...
        'MarkerEdgeColor', c(1, :), 'MarkerFaceColor', c(1, :))
    hold on
    plot(s(4, 1:3), '-s', 'Color', c(6, :), 'LineWidth', 3, 'MarkerSize', 10, ...
        'MarkerEdgeColor', c(6, :),'MarkerFaceColor', c(6, :))
    plot(s(2, 1:3), '-p', 'Color', c(2, :), 'LineWidth', 3, 'MarkerSize', 10, ...
        'MarkerEdgeColor', c(2, :), 'MarkerFaceColor', c(2, :))
    plot(s(3, 1:3), '-*', 'Color', c(4, :), 'LineWidth', 3, 'MarkerSize', 10, ...
        'MarkerEdgeColor', c(4, :), 'MarkerFaceColor', c(4, :))
    ylabel('Ratio from all simulations, %')
    grid on
    ylim([0 100])
    xticks([1 2 3])
    xticklabels({'3', '\geq 2', '\geq 1'})
    xlabel('Number of detected sources')
    set(gca, 'FontSize', 18)
    title('Detection ratio')
    legend('ReciPSIICOS', 'Whitened ReciPSIICOS', 'LCMV', 'MNE', ...
        'Location', 'southeast')

end