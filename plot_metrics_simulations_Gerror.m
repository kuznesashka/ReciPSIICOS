function [r, var] = plot_metrics_simulations_Gerror(G3, G3_red, ...
    Z_total, picked_src)

%--------------------------------------------------------------------------
% Calculate point spreading and bias for precomputed `Z_total_G.mat` and
% `picked_src_G.mat` and plot them.
%
% Parameters
% ----------
% G3 : Forward operator for dense cortical model
% G3_red : Forward operator for reduced cortical model
% Z_total : (method, synch, length(rank), Nmc, Nsrc), where methods: 
%       [ReciPSIICOS, WReciPSIICOS, LCMV, MNE]
% picked_src : (Nmc, 2) pairs of simulated sources
%
% Returns
% -------
% r : bias
% var : point spreading
%--------------------------------------------------------------------------

    R = G3.GridLoc; % source location, dense matrix
    R_red = G3_red.GridLoc; % source location reduced matrix

    Nmc = size(picked_src, 1);
    % 2. SIMULATIONS 

    src_left_red = find(R_red(:, 2) > 0);
    src_right_red = find(R_red(:, 2) < 0);

    fm_error = [0, 0.05, 0.1, 0.2];
    Range_frac = [0.65; 0.25];
    thresh = [0.02, 0.02];

    bad_detection = zeros(4, 2, length(fm_error), size(Range_frac, 2));

    for fm_error_i = 1:length(fm_error) 
        for mc = 1:Nmc
            ind_generated = picked_src(mc, :);
            for synch = 1:2
                for method = 1:4
                    for frac = 1:size(Range_frac, 2)
                        Z = squeeze(Z_total(method, synch, fm_error_i, mc, :))';
                        [max_val, ~] = max(Z);
                        if method == 4
                        Z(Z < 0.65 * max_val) = 0;
                        else
                        Z(Z < Range_frac(synch, frac) * max_val) = 0;
                        end

                        clear dist_l dist_r
                        Z_left = Z(src_left_red); % values from left hemisphere
                        Z_right = Z(src_right_red); % values from right hemisphere
                        [~, max_l] = max(Z_left); % maximum in the left hemisphere
                        max_ind_l = src_left_red(max_l);
                        [~, max_r] = max(Z_right);
                        max_ind_r = src_right_red(max_r);

                        if (length(Z_left(Z_left > 0)) == 0) | ...
                                (length(Z_right(Z_right > 0)) == 0)
                            r(method, synch, fm_error_i, frac, mc) = NaN;
                            var(method, synch, fm_error_i, frac, mc) = NaN;
                            bad_detection(method, synch, fm_error_i, frac) = ...
                                bad_detection(method, synch, fm_error_i, frac) + 1; 
                        else
                            r(method, synch, fm_error_i, frac, mc) = ...
                                (norm(R_red(max_ind_l, :) - R(ind_generated(1), :))+ ...
                                norm(R_red(max_ind_r, :) - R(ind_generated(2), :))) / 2;
                            if r(method, synch, fm_error_i, frac, mc) > thresh(synch)
                               bad_detection(method, synch, fm_error_i, frac) = ...
                                   bad_detection(method, synch, fm_error_i, frac) + 1; 
                            end

                            Z_left_norm = Z_left ./ sum(Z_left);
                            Z_right_norm = Z_right ./ sum(Z_right);

                            coord_active_l = R_red(src_left_red(Z_left_norm > 0), :);
                            for i = 1:size(coord_active_l, 1)
                                dist_l(i) = norm(coord_active_l(i, :) - R_red(max_ind_l, :));
                            end

                            coord_active_r = R_red(src_right_red(Z_right_norm > 0), :);
                            for i = 1:size(coord_active_r, 1)
                                dist_r(i) = norm(coord_active_r(i, :) - R_red(max_ind_r, :));
                            end

                            var(method, synch, fm_error_i, frac, mc) = ...
                                (sum(Z_left_norm(Z_left_norm > 0) .* dist_l) + ...
                                sum(Z_right_norm(Z_right_norm > 0) .* dist_r)) / 2;
                        end
                    end
                end        
            end
            mc
        end
    end

    mean_r = median(r, 5, 'omitnan');
    mean_var = median(var, 5, 'omitnan');


    c = lines(7);
    col = c([1, 6, 2, 4], :);
    figure
    subplot(3, 2, 1)
    hold on
    grid on
    for i = 1:4
        plot(fm_error, squeeze(mean_r(i, 1, :)), 'Color', col(i, :), 'LineWidth', 3)
    end
    scatter(fm_error, squeeze(mean_r(1, 1, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_r(2, 1, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_r(3, 1, :)), ...
        150, 'p', 'MarkerFaceColor', col(3, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_r(4, 1, :)), ...
        150, 'h', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Error rate, \delta')
    ylabel('Distance, m')
    ylim([0 0.07])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = 2;

    subplot(3, 2, 3)
    hold on
    grid on
    for i = 1:4
        plot(fm_error, squeeze(mean_var(i, 1, :)), 'Color', col(i, :), 'LineWidth', 3)
    end
    scatter(fm_error, squeeze(mean_var(1, 1, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_var(2, 1, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_var(3, 1, :)), ...
        150, 'p', 'MarkerFaceColor', col(3, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_var(4, 1, :)), ...
        150, 'h', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Error rate, \delta')
    ylabel('Distance, m')
    ylim([0 0.03])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = 2;

    subplot(3, 2, 5)
    hold on
    grid on
    for i = 1:4
        plot(fm_error, (1 - (squeeze(bad_detection(i, 1, :, 1)) ./ Nmc)) * 100, ...
            'Color', col(i, :), 'LineWidth', 3)
    end
    scatter(fm_error, (1 - (squeeze(bad_detection(1, 1, :, 1)) ./ Nmc)) * 100, ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, (1 - (squeeze(bad_detection(2, 1, :, 1)) ./ Nmc)) * 100, ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, (1 - (squeeze(bad_detection(3, 1, :, 1)) ./ Nmc)) * 100, ...
        150, 'p', 'MarkerFaceColor', col(3, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, (1 - (squeeze(bad_detection(4, 1, :, 1)) ./ Nmc)) * 100, ...
        150, 'h', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Error rate, \delta')
    ylabel('Percents, %')
    ylim([0 100])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = 2;

    subplot(3, 2, 2)
    hold on
    grid on
    for i = 1:3
        plot(fm_error, squeeze(mean_r(i, 2, :)), 'Color', col(i, :), 'LineWidth', 3)
    end
    scatter(fm_error, squeeze(mean_r(1, 2, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_r(2, 2, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_r(3, 2, :)), ...
        150, 'p', 'MarkerFaceColor', col(3, :), 'MarkerEdgeColor', 'k')
    plot(fm_error, squeeze(mean_r(4, 2, :)), 'Color', col(4, :), 'LineWidth', 3)
    scatter(fm_error, squeeze(mean_r(4, 2, :)), ...
        150, 'h', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Error rate, \delta')
    ylabel('Distance, m')
    ylim([0.01 0.015])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = 2;

    subplot(3, 2, 4)
    hold on
    grid on
    for i = 1:3
        plot(fm_error, squeeze(mean_var(i, 2, :)), 'Color', col(i, :), 'LineWidth', 3)
    end
    scatter(fm_error, squeeze(mean_var(1, 2, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_var(2, 2, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, squeeze(mean_var(3, 2, :)), ...
        150, 'p', 'MarkerFaceColor', col(3, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    plot(fm_error, squeeze(mean_var(4, 2, :)), 'Color', col(4, :), 'LineWidth', 3)
    scatter(fm_error, squeeze(mean_var(4, 2, :)), ...
        150, 'h', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
    xlabel('Error rate, \delta')
    ylabel('Distance, m')
    ylim([0 0.015])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = 2;

    subplot(3, 2, 6)
    hold on
    grid on
    for i = 1:3
        plot(fm_error, (1 - (squeeze(bad_detection(i, 2, :, 1)) ./ Nmc)) * 100, ...
            'Color', col(i, :), 'LineWidth', 3)
    end
    scatter(fm_error, (1 - (squeeze(bad_detection(1, 2, :, 1)) ./ Nmc)) * 100, ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, (1 - (squeeze(bad_detection(2, 2, :, 1)) ./ Nmc)) * 100, ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    scatter(fm_error, (1 - (squeeze(bad_detection(3, 2, :, 1)) ./ Nmc)) * 100, ...
        150, 'p', 'MarkerFaceColor', col(3, :), 'MarkerEdgeColor', 'k')
    plot(fm_error, (1 - (squeeze(bad_detection(4, 2, :, 1)) ./ Nmc)) * 100, ...
        'Color', col(4, :), 'LineWidth', 3)
    scatter(fm_error, (1 - (squeeze(bad_detection(4, 2, :, 1)) ./ Nmc)) * 100, ...
        150, 'h', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Error rate, \delta')
    ylabel('Percents, %')
    ylim([0 100])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = 2;
end
