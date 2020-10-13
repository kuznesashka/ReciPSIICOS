function [Z_total, picked_src] = save_simulations_snr(G3, G3_red, chans, ...
    snr, Nmc, load_corr, load_src)

%--------------------------------------------------------------------------
% Generate simulations with different SNRs, then apply source reconstruction 
% with PSIICOS, ReciPSIICOS, LCMV, MNE techniques.
%
% Parameters
% ----------
% G3 : Forward operator for dense cortical model
% G3_red : Forward operator for reduced cortical model
% chans : channel locations
% snr : list of SNR values
% Nmc : number of Monte Carlo simulations
% load_corr : bool
%       If true, load the precomputed correlation matrix
% load_src : bool
%       If true, generate simulations for precomputed symmetrical sources
%
% Returns
% -------
% Z_total : (method, synch, length(rank), Nmc, Nsrc), where methods: 
%       [ReciPSIICOS, WReciPSIICOS, LCMV, MNE]
% picked_src : (Nmc, 2) pairs of simulated sources
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
    % but produces less contrasting subcorr scans. For a more reliable 
    % preformance use 0.01 to get all the sensor on board but be ready to wait

    GainSVDTh = 0.001; 
    [ug, ~, ~] = spm_svd(G2d_red * G2d_red', GainSVDTh);
    UP = ug'; % direction of dimention reduction
    G2dU_red = UP * G2d_red;
    G2d0U_red = UP * G2d0_red;

    % 3. PROJECTION MATRICES
    % now find span of the target space
    
    Swp = zeros(4);
    Swp(2, 3) = 1;
    Swp(3, 2) = 1;
    Swp(1, 1) = 1; 
    Swp(4, 4) = 1;
    RankG = size(G2d0U_red, 1);

   if (~load_corr)
        C_re = zeros(RankG^2, RankG^2);
        parfor i = 1:Nsites
             range_i = i * 2 - 1:i * 2;
             ai = G2d0U_red(:, range_i);
             rng = 1:4;
             X = zeros(RankG^2, 4 * (Nsites - i), 'single');
             for j = i + 1:Nsites
                range_j = j * 2 - 1:j * 2;
                aj = G2d0U_red(:, range_j);
                X(:, rng) = kron(ai, aj) + kron(aj, ai) * Swp;
                rng = rng + 4;
             end
             C_re = C_re + X*X';
             i
        end
        save C_re C_re;
    end

    % 3. PSIICOS projection for sparse matrix
    [Upwr, ~, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, 1500, 0);

    % power subspace correlation matrix
    Cpwr = Apwr * Apwr';
    Wks = load('C_re.mat','C_re');
    
    % correlation subspace correlation matrix
    Ccorr = Wks.C_re;

    % whitener wrt to power subspace 
    Wpwr = sqrtm(inv(Cpwr + 0.001 * trace(Cpwr) / (RankG^2) * eye(size(Cpwr))));
    WCcorrW = Wpwr * double(Ccorr) * Wpwr';

    ProjRnkMax = 1280;
    [u, ~] = eigs(WCcorrW, ProjRnkMax);
    UcorrW_rnk = u(:, 1:ProjRnkMax);

    PrFromCorr_W = inv(Wpwr + 0.001 * trace(Wpwr)/size(Wpwr, 1) * eye(size(Wpwr, 1)))* ...
        (eye(size(UcorrW_rnk, 1)) - UcorrW_rnk(:, 1:500) * UcorrW_rnk(:, 1:500)') * Wpwr;
    PrPwr = Upwr(:, 1:500) * Upwr(:, 1:500)';

    % 5. SIMULATIONS
    Fs = 500; % sampling frequency
    Ntr = 100; % number of simulated trials
    T = Fs; % number of time points in one trial
    t = 1:T;

    Zp = zeros(2, length(snr), Nmc, Nsites_red);
    Zp_two = zeros(2, length(snr), Nmc, Nsites_red);
    Zpw = zeros(2, length(snr), Nmc, Nsites_red);
    Zbf = zeros(2, length(snr), Nmc, Nsites_red);
    Zmne = zeros(2, length(snr), Nmc, Nsites_red);

    if (load_src)
        picked_src = load('picked_src_snr');
        picked_src = picked_src.picked_src;
    end

    for mc = 1:Nmc
        % generate brain noise with dense matrix
        range = 1:T;
        for tr = 1:Ntr
            Noise(:, range) = GenerateBrainNoise(G2d0, T, 200, 500, Fs);
            range = range + T;
        end
        Noise_av = mean(reshape(Noise, [Nch, T, Ntr]), 3); % average by trial
        Noise_av_0 = Noise_av / norm(Noise_av); % normalized average noise

        % generate signal for dense forward model matrix

        if (~load_src)
            % sources from left hemisphere > 2 cm far from midline (6542/10001)
            far_sites = find(R(:, 2) > 0.02);
            rand_idx = randperm(length(far_sites));
            % pick the first source from the left hemisphere > 2 cm far from midline
            picked_src(mc, 1) = far_sites(rand_idx(1)); 

            for i = 1:Nsites % find the symmetrical source
                dd_sym(i) = norm([R(picked_src(mc, 1), 1), ...
                    -R(picked_src(mc, 1), 2), R(picked_src(mc, 1), 3)] - R(i, :));
            end
            [~, picked_src(mc, 2)] = min(dd_sym);
        end

        picked_src_oriented = picked_src(mc, :) * 2;

        range = 1:T; % activation functions
        for tr = 1:Ntr
            S(1, range) = sin(2 * pi * t / Fs + randn * 0.1);
            S(2, range) = sin(2 * pi * t / Fs);
            S(3, range) = cos(2 * pi * t / Fs);
            range = range + T;
        end

        for noise_i = 1:length(snr) 
            for synch = 1:2
                X = G2d0(:, picked_src_oriented(1)) * S(synch, :) + ...
                    G2d0(:, picked_src_oriented(2)) * S(synch + 1, :);
                X_av = mean(reshape(X, [Nch, T, Ntr]), 3);

                X_av_0 = X_av / norm(X_av);
                % add noise to the data
                Data  = snr(noise_i) * X_av_0 + Noise_av_0;
                % compute covariance in the virtual sensor space
                Ca = UP * Data * Data' * UP';

                % LCMV BF
                Zbf(synch, noise_i, mc, :) = lcmv(G2dU_red, Ca);

                % Whitened ReciPSIICOS beamformer
                Cap =reshape(PrFromCorr_W * Ca(:), size(Ca));
                [e, a] = eig(Cap);
                Cap = e * abs(a) * e';
                iCap = tihinv(Cap, 0.01);

                range2d = 1:2;
                for i = 1:Nsites_red
                    g = G2dU_red(:, range2d);
                    m = inv(g' * iCap * g);
                    [~, ss, ~] = svd(m);
                    Zpw(synch, noise_i, mc, i) = ss(1, 1);
                    range2d = range2d + 2;
                end

                % ReciPSIICOS beamformer
                Cap = reshape(PrPwr * Ca(:), size(Ca));
                [e, a] = eig(Cap);
                Cap = e * abs(a) * e';
                iCap = tihinv(Cap, 0.01);

                range2d = 1:2;
                for i = 1:Nsites_red
                    g = G2dU_red(:, range2d);
                    m = inv(g' * iCap * g);
                    [~, ss, ~] = svd(m);
                    Zp(synch, noise_i, mc, i) = ss(1, 1);
                    range2d = range2d + 2;
                end

                % Whitened ReciPSIICOS beamformer
                Cap = reshape(PrPwr * PrFromCorr_W * Ca(:), size(Ca));
                [e, a] = eig(Cap);
                Cap = e * abs(a) * e';
                iCap = tihinv(Cap, 0.01);

                range2d = 1:2;
                for i = 1:Nsites_red
                    g = G2dU_red(:, range2d);
                    m = inv(g' * iCap * g);
                    [~, ss, ~] = svd(m);
                    Zp_two(synch, noise_i, mc, i) = ss(1, 1);
                    range2d = range2d + 2;
                end

                % Minimum norm estimate
                lambda = 0.1;
                Zmne(synch, noise_i, mc, :) = mne(G2dU_red, Ca, lambda);
            end
        end
        mc
    end

    Z_total = zeros(5, 2, length(snr), Nmc, 5001);
    Z_total(1, :, :, :, :) = Zp;
    Z_total(2, :, :, :, :) = Zpw;
    Z_total(3, :, :, :, :) = Zbf;
    Z_total(4, :, :, :, :) = Zmne;
    Z_total(5, :, :, :, :) = Zp_two;

end