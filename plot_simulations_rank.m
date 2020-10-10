function [r, var] = plot_simulations_rank(G3, G3_red, ...
    Z_total, picked_src, rank)

%--------------------------------------------------------------------------
% Calculate point spreading and bias for precomputed `Z_total_rank.mat` and
% `picked_src_rank.mat` and plot them.
%
% Parameters
% ----------
% G3 : Forward operator for dense cortical model
% G3_red : Forward operator for reduced cortical model
% Z_total : (method, synch, length(rank), Nmc, Nsrc), where methods: 
%       [ReciPSIICOS, WReciPSIICOS, LCMV, MNE]
% picked_src : (Nmc, 2) pairs of simulated sources
% rank : list of projection rank values
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
    % Z_total = (method, synch, length(rank), Nmc, Nsrc)
    % methods: ReciPSIICOS, WReciPSIICOS, LCMV, MNE

    src_left_red = find(R_red(:, 2) > 0);
    src_right_red = find(R_red(:, 2) < 0);

    Range_frac = [0.65; 0.25];
    thresh = [0.02, 0.02];

    bad_detection = zeros(2, 2, length(rank), size(Range_frac, 2));
    clear r var
    for rank_i = 1:length(rank) 
        for mc = 1:Nmc
            ind_generated = picked_src(mc, :);
            for synch = 1:2
                for method = 1:2
                    for frac = 1:size(Range_frac, 2)
                        Z = squeeze(Z_total(method, synch, rank_i, mc, :))';
                        [max_val, ~] = max(Z);
                        Z(Z < Range_frac(synch, frac) * max_val) = 0;

                        clear dist_l dist_r
                        Z_left = Z(src_left_red); % values from left hemisphere
                        Z_right = Z(src_right_red); % values from right hemisphere
                        [~, max_l] = max(Z_left); % maximum in the left hemisphere
                        max_ind_l = src_left_red(max_l);
                        [~, max_r] = max(Z_right);
                        max_ind_r = src_right_red(max_r);

                        if (length(Z_left(Z_left > 0)) == 0) | ...
                                (length(Z_right(Z_right > 0)) == 0)
                            r(method, synch, rank_i, frac, mc) = NaN;
                            var(method, synch, rank_i, frac, mc) = NaN;
                            bad_detection(method, synch, rank_i, frac) = ...
                                bad_detection(method, synch, rank_i, frac) + 1; 
                        else
                            r(method, synch, rank_i, frac, mc) = ...
                                (norm(R_red(max_ind_l, :) - R(ind_generated(1), :))+ ...
                                norm(R_red(max_ind_r, :) - R(ind_generated(2), :))) / 2;
                            if r(method, synch, rank_i, frac, mc) > thresh(synch)
                               bad_detection(method, synch, rank_i, frac) = ...
                                   bad_detection(method, synch, rank_i, frac) + 1; 
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

                            var(method, synch, rank_i, frac, mc) = ...
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
    grid on
    hold on
    for i = 1:2
        plot(rank, squeeze(mean_r(i, 1, :)), 'Color', col(i, :), 'LineWidth', 3)
        ylim([0 0.03])
        xlim([0 1100])
    end
    scatter(rank, squeeze(mean_r(1, 1, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(rank, squeeze(mean_r(2, 1, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    set(gca,'FontSize', 24)
    xlabel('Projection rank')
    ylabel('Distance, m')
    plot([500 500], [0 0.03], '--', 'Color', c(1, :), 'LineWidth', 2)
    plot([1000 1000], [0 0.03], '--', 'Color', c(6, :), 'LineWidth', 2)

    subplot(3, 2, 3)
    grid on
    hold on
    for i = 1:2
        plot(rank, squeeze(mean_var(i, 1, :)), 'Color', col(i, :), 'LineWidth', 3)
         ylim([0 0.02])
        xlim([0 1100])
    end
    scatter(rank, squeeze(mean_var(1, 1, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(rank, squeeze(mean_var(2,1,:)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Projection rank')
    ylabel('Distance, m')
    plot([500 500], [0 0.03], '--', 'Color', c(1, :), 'LineWidth', 2)
    plot([1000 1000], [0 0.03], '--', 'Color', c(6, :), 'LineWidth', 2)

    subplot(3, 2, 5)
    grid on
    hold on
    for i = 1:2
        plot(rank, (1 - (squeeze(bad_detection(i, 1, :, 1)) ./ Nmc)) * 100, ...
            'Color', col(i, :), 'LineWidth', 3)
        ylim([30 100])
        xlim([0 1100])
    end
    scatter(rank, (1 - (squeeze(bad_detection(1, 1, :, 1)) ./ Nmc)) * 100, ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(rank, (1 - (squeeze(bad_detection(2, 1, :, 1)) ./ Nmc)) * 100, ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Projection rank')
    ylabel('Percents, %')
    plot([500 500], [0 100], '--', 'Color', c(1, :), 'LineWidth', 2)
    plot([1000 1000], [0 100], '--', 'Color', c(6, :), 'LineWidth', 2)

    subplot(3, 2, 2)
    grid on
    hold on
    for i = 1:2
        plot(rank, squeeze(mean_r(i, 2, :)), 'Color', col(i,:), 'LineWidth', 3)
        ylim([0 0.03])
        xlim([0 1100])
    end
    scatter(rank, squeeze(mean_r(1, 2, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(rank, squeeze(mean_r(2, 2, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Projection rank')
    ylabel('Distance, m')
    plot([500 500], [0 0.03], '--', 'Color', c(1, :), 'LineWidth', 2)
    plot([1000 1000], [0 0.03], '--', 'Color', c(6, :), 'LineWidth', 2)

    subplot(3, 2, 4)
    grid on
    hold on
    for i = 1:2
        plot(rank, squeeze(mean_var(i, 2, :)), 'Color', col(i, :), 'LineWidth', 3)
        ylim([0 0.02])
        xlim([0 1100])
    end
    scatter(rank, squeeze(mean_var(1, 2, :)), ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(rank, squeeze(mean_var(2, 2, :)), ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Projection rank')
    ylabel('Distance, m')
    plot([500 500], [0 0.03], '--', 'Color', c(1, :), 'LineWidth', 2)
    plot([1000 1000], [0 0.03], '--', 'Color', c(6, :), 'LineWidth', 2)

    subplot(3, 2, 6)
    grid on
    hold on
    for i = 1:2
        plot(rank, (1 - (squeeze(bad_detection(i, 2, :, 1)) ./ Nmc)) * 100, ...
            'Color', col(i, :), 'LineWidth', 3)
        ylim([30 100])
        xlim([0 1100])
    end
    scatter(rank, (1 - (squeeze(bad_detection(1, 2, :, 1)) ./ Nmc)) * 100, ...
        150, 'o', 'MarkerFaceColor', col(1, :), 'MarkerEdgeColor', 'k')
    scatter(rank, (1 - (squeeze(bad_detection(2, 2, :, 1)) ./ Nmc)) * 100, ...
        150, 's', 'MarkerFaceColor', col(2, :), 'MarkerEdgeColor', 'k')
    set(gca, 'FontSize', 24)
    xlabel('Projection rank')
    ylabel('Percents, %')
    plot([500 500], [0 100], '--', 'Color', c(1, :), 'LineWidth', 2)
    plot([1000 1000], [0 100], '--', 'Color', c(6, :), 'LineWidth', 2)
    
end




       
