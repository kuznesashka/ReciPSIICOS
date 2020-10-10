function [r, var, detect] = ...
    save_simulations_three_src_asynch(G3, G3_red, chans, Nmc)

%--------------------------------------------------------------------------
% Generate simulations with fixed SNR for three asynchronous sources. Then
% run source recounstrunction and calculate bias, point spreading and
% detection ratio for four methods.
%
% Parameters
% ----------
% G3 : Forward operator for dense cortical model
% G3_red : Forward operator for reduced cortical model
% chans : channel locations
% Nmc : number of Monte Carlo simulations
%
% Returns
% -------
% r : (method, length(rank), Nmc), where methods: 
%       [ReciPSIICOS, WReciPSIICOS, LCMV, MNE], bias
% var : (method, length(rank), Nmc), point spreading
% detect : (method, length(rank), Nmc), detection ratio
%--------------------------------------------------------------------------

    R = G3.GridLoc; % source location, dense matrix
    R_red = G3_red.GridLoc; % source location reduced matrix
    % set to use gradiometers only
    ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD'));
    Nch = length(ChUsed);

    [~, G2d0, Nsites] = G3toG2(G3, ChUsed);
    [G2d_red, G2d0_red, Nsites_red] = G3toG2(G3_red, ChUsed);

    % 2. REDUCING SENSOR SPACE for sparse matrix
    % 0.05 results into 47 eigensensors and makes it run faster 
    % but produces less contrasting subcorr scans for a more reliable 
    % performance use 0.01 to get all the sensor on board but be ready to wait
    
    GainSVDTh = 0.001; 
    [ug, ~, ~] = spm_svd(G2d_red * G2d_red', GainSVDTh);
    UP = ug'; % direction of dimention reduction
    G2dU_red = UP * G2d_red;
    G2d0U_red = UP * G2d0_red;

    % 3. PSIICOS projection for sparse matrix
    [Upwr, ~, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, 1500, 0);
    Rnk_rp = 500;

    % 4. Simulations
    snr = [5*3/2]; % snr level in the data
    c = lines(50); % colors
    Fs = 500; % sampling frequency
    Ntr = 100; % number of simulated trials
    T = Fs; % number of time points in one trial
    t = linspace(1, 2 * pi, T);

    RankG = size(G2d0U_red, 1);
    Wks =  load('C_re.mat', 'C_re');
    Ccorr = Wks.C_re;
    Cpwr = Apwr * Apwr';
    Wpwr = sqrtm(inv(Cpwr + 0.001 * trace(Cpwr) / (RankG^2) * eye(size(Cpwr))));
    WCcorrW = Wpwr * double(Ccorr) * Wpwr';
    ProjRnkMax = 1000;
    [u, ~] = eigs(WCcorrW, ProjRnkMax);
    UcorrW_rnk = u(:, 1:ProjRnkMax);

    PrFromCorr_W = inv(Wpwr + 0.001 * trace(Wpwr) / size(Wpwr, 1) * eye(size(Wpwr, 1)))* ...
        (eye(size(UcorrW_rnk, 1)) - UcorrW_rnk(:, 1:1000) * UcorrW_rnk(:, 1:1000)') * Wpwr;

    r = zeros(4, length(snr), Nmc);
    var = zeros(4, length(snr), Nmc);

    for noise_i = 1:length(snr) 
        for mc = 1:Nmc
            % generate brain noise with dense matrix
            range = 1:T;
            for tr = 1:Ntr
                Noise(:, range) = GenerateBrainNoise(G2d0, T, 200, 500, Fs);
                range = range + T;
            end
            Noise_av = mean(reshape(Noise, [Nch, T, Ntr]), 3); % average by trial
            Noise_av_0 = Noise_av / norm(Noise_av); % normalized average noise

            % generate signal 
            Farsites = find(R(:, 2) > 0.05); % sources from left hemisphere > 2 cm far from midline (6542/10001)
            ii = randperm(length(Farsites));
            src(1) = Farsites(ii(1));  % pick the first source from the left hemisphere > 2 cm far from midline 

            Farsites = find(R(:, 2) < -0.05); % sources from left hemisphere > 2 cm far from midline (6542/10001)
            ii = randperm(length(Farsites));
            src(2) = Farsites(ii(1));  % pick the first source from the left hemisphere > 2 cm far from midline 

            clear dd
            for i = 1:Nsites
                for j = 1:2
                    dd(i, j) = norm(R(src(j), :) - R(i, :));
                end
            end
            range = 1:Nsites;
            far_third = range((dd(:, 1) > 0.06) & (dd(:, 2) > 0.06));
            ii = randperm(length(far_third));
            src(3) = far_third(ii(1));
            I = src * 2;

%             figure
%             scatter3(R_red(:,1),R_red(:,2),R_red(:,3))
%             hold on
%             scatter3(R(src,1),R(src,2),R(src,3), 'r', 'filled')

            fi = randn;
            S(1, :) = sin(fi + t);
            S(2, :) = sin(fi + t + pi/3);
            S(3, :) = sin(fi + t + 2*pi/3);

%             figure
%             plot(1:2:1000, S(1, :), 'Color', c(7, :), 'LineWidth', 3)
%             hold on
%             plot(1:2:1000, S(2, :), 'Color', c(1, :), 'LineWidth', 3)
%             plot(1:2:1000, S(3, :), 'Color', c(5, :), 'LineWidth', 3)
%             xlabel('Time, ms')
%             ylabel('Activation amplitude, a.u.')
%             set(gca, 'FontSize', 18)

            X = G2d0(:, I(1)) * S(1, :) + G2d0(:, I(2)) * S(2, :) + ...
                G2d0(:, I(3)) * S(3, :); % signal, activations only by y-axis of G
            X_av_0 = X / norm(X);

            Data = snr(noise_i) * X_av_0 + Noise_av_0; % add noise to the data
            Ca = UP * Data * Data' * UP'; % compute covariance in the virtual sensor space

            % AntiPSIICOS
            Zp = Apsiicos(G2dU_red, Ca, Upwr, Rnk_rp, 1, 'power');
            % LCMV
            Zbf = lcmv(G2dU_red, Ca);
            % MNE
            Zmne = mne(G2dU_red, Ca, 0.1);

            % Whitened reciPsiicos
            Cap = reshape(PrFromCorr_W * Ca(:), size(Ca));
            [e, a] = eig(Cap);
            Cap = e * abs(a) * e';
            iCap = tihinv(Cap, 0.01);
            range2d = 1:2;
            for i = 1:Nsites_red
                g = G2dU_red(:, range2d);
                m = inv(g' * iCap * g);
                [~, ss, ~] = svd(m);
                Zpw(i) = ss(1, 1);
                range2d = range2d + 2;
            end

            % closest sources
            clear dd
            for i = 1:3
                for j = 1:size(R_red, 1)
                    dd(i, j) = norm(R(src(i), :) - R_red(j, :));
                end
            end

            % GOODNESS OF DETECTION METRICS
            InvSol = [Zp; Zbf; Zmne; Zpw];
            frac_all = {[0.025:0.05:0.5], [0.1:0.1:0.7], [0.1:0.1:0.7], ...
                [0.025:0.05:0.5]};        

            for method = 1:size(InvSol, 1)
                frac = frac_all{method};
                Zinit = InvSol(method, :);
                max_val = max(Zinit);
                Z = repmat(Zinit, length(frac), 1)';
                Z(Z < [frac * max_val]) = 0;
                [cluster, ~, maxind, ~] = clustering(Z, R_red, 0.025);

                clear val_detect_all src_detect_all detect_num
                for frac_num = 1:length(frac)
                    clear dd
                    % number of clusters found for this frac
                    for i = 1:length(find(~cellfun(@isempty, cluster(frac_num, :))))
                        for j = 1:length(src)
                            dd(i, j) = norm(R_red(maxind{frac_num, i}, :) - R(src(j), :));
                        end
                    end
                    [val_detect_all{frac_num}, src_detect_all{frac_num}] = min(dd');
                    detect_num(frac_num) = length(unique(src_detect_all{frac_num}));
                end

                best_ind = max(find(detect_num == 3));
                if isempty(best_ind)
                    best_ind = max(find(detect_num == 2));
                end
                if isempty(best_ind)
                    best_ind = max(find(detect_num == 1));
                end

                Z = Zinit;
                Z(Z < frac(best_ind) * max_val) = 0;

                best_cluster = cluster(best_ind, :);
                best_cluster = best_cluster(1, find(~cellfun(@isempty, best_cluster)));
                best_maxind = maxind(best_ind, :);
                src_detect = src_detect_all{best_ind};
                val_detect = val_detect_all{best_ind};

%                 Inv.ImageGridAmp = Z';
%                 Inv.ImageGridAmp = [];
%                 Inv.ImageGridAmp = zeros(5001, 1);
%                 Inv.ImageGridAmp(src_detect) = 100;
%                             
%                 figure
%                 scatter3(R_red(:, 1), R_red(:, 2), R_red(:, 3), 'MarkerEdgeAlpha', 0.2)
%                 hold on
%                 for i = 1:size(best_cluster, 2)
%                     scatter3(R_red(best_cluster{i}, 1), R_red(best_cluster{i}, 2), ...
%                         R_red(best_cluster{i}, 3), 'filled', 'MarkerFaceColor', c(i, :))
%                 end
%                 scatter3(R(src, 1), R(src, 2), R(src, 3), ...
%                     'filled', 'MarkerFaceColor', 'r')
%                 grid off
%                 view(270, 90)

                if size(best_cluster, 1) == 0
                    r(method, noise_i, mc) = NaN;
                    var(method, noise_i, mc) = NaN;
                else
                    detect(method, noise_i, mc) = length(unique(src_detect(val_detect < 0.02)));
                    for i = 1:size(best_cluster, 2)
                        r(method, noise_i, mc) = r(method, noise_i, mc) + ...
                            norm(R_red(best_maxind{i}, :) - R(src(src_detect(i)), :));
                    end
                    r(method, noise_i, mc) = r(method, noise_i, mc) / size(best_cluster, 2);

                    clear dist pointspr
                    for i = 1:size(best_cluster, 2)
                        Z_norm{i} = Z(best_cluster{i}) ./ sum(Z(best_cluster{i}));
                        coord_active = R_red(best_cluster{i}, :);
                        for j = 1:size(coord_active, 1)
                            dist{i}(j) = norm(coord_active(j, :) - R_red(best_maxind{i}, :));
                        end
                        pointspr{i} = sum(Z_norm{i} .* dist{i});
                        var(method, noise_i, mc) = var(method, noise_i, mc) + pointspr{i};
                    end
                        var(method, noise_i, mc) = var(method, noise_i, mc) / size(best_cluster, 2);
                end
                method
            end
            mc
        end
    end

    figure 
    plot(1:2:1000, S(1, 1:T), 'Color', c(7, :), 'LineWidth', 3)
    hold on
    plot(1:2:1000, S(2, 1:T), 'Color', c(1, :), 'LineWidth', 3)
    plot(1:2:1000, S(3, 1:T), 'Color', c(5, :), 'LineWidth', 3)
    xlabel('Time, ms')
    ylabel('Activation amplitude, a.u.')
    set(gca, 'FontSize', 18)

end

