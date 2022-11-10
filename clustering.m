function [cluster, center, maxind, clustnum] = clustering(Z, R_red, thr_dist)
    
    for frac = 1:size(Z, 2)
        src_idx = find(Z(:, frac) > 0)';

        for i = 1:size(src_idx, 2) % distances between each vertex
            for j = 1:size(src_idx, 2)
                dist(i, j) = norm(R_red(src_idx(i), :) - R_red(src_idx(j), :));
            end
        end

        Nmin = 1; % minimal number of entries in cluster
        fl = 1; 
        k = 1;
        while fl == 1
            dst = sum(dist < thr_dist, 2); 
            [val, ind] = max(dst); % vertex with the highest number of close neighbours
            if val >= Nmin
                ind_nbh = find(dist(ind, :) < thr_dist); % neighbours
                cluster{frac, k} = src_idx(:, ind_nbh);
                center{frac, k} = mean(R_red(cluster{frac, k}, :));
                [val, indm] = max(Z(cluster{frac, k}));
                maxind{frac, k} = cluster{frac, k}(indm);
                src_idx = src_idx(:, setdiff(1:size(dist, 1), ind_nbh));
                dist = dist(setdiff(1:size(dist, 1), ind_nbh), setdiff(1:size(dist, 1), ind_nbh));
                k = k + 1;
                fl = 1;
            else
                fl = 0;
                clustnum(frac) = k - 1;
            end
        end
    end
end
