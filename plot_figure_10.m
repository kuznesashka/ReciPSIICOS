G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channel locations

Nmc = 1;
% 1.1 Generate simulations for different methods (MNE, LCMV beamformer, 
% ReciPSIICOS, WReciPSIICOS) and different SNRs

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix
% use gradiometers only
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
snr = [5 * 3 / 2]; % snr level in the data
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

PrFromCorr_W = inv(Wpwr + 0.001 * trace(Wpwr)/size(Wpwr, 1) * eye(size(Wpwr, 1)))* ...
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
        % sources from left hemisphere > 2 cm far from midline (6542/10001)
        Farsites = find(R(:, 2) > 0.05);
        ii = randperm(length(Farsites));
        % pick the first source from the left hemisphere > 2 cm far from midline
        src(1) = Farsites(ii(1));  

        % sources from left hemisphere > 2 cm far from midline (6542/10001)
        Farsites = find(R(:, 2) < -0.05);
        ii = randperm(length(Farsites));
        % pick the first source from the left hemisphere > 2 cm far from midline 
        src(2) = Farsites(ii(1));

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
%             scatter3(R_red(:, 1), R_red(:, 2), R_red(:, 3))
%             hold on
%             scatter3(R(src, 1), R(src, 2), R(src, 3), 'r', 'filled')

        fi = randn;
        S(1, :) = sin(0.1 * fi + t);
        S(2, :) = sin(0.05 * fi + t);
        S(3, :) = sin(t);

%             figure
%             plot(1:2:1000, S(1, :), 'Color', c(7, :), 'LineWidth', 3)
%             hold on
%             plot(1:2:1000, S(2, :), 'Color', c(1, :), 'LineWidth', 3)
%             plot(1:2:1000, S(3, :), 'Color', c(5, :), 'LineWidth', 3)
%             xlabel('Time, ms')
%             ylabel('Activation amplitude, a.u.')
%             set(gca, 'FontSize', 18)

        % signal, activations only by y-axis of G
        X = G2d0(:, I(1)) * S(1, :) + G2d0(:, I(2)) * S(2, :) + ...
            G2d0(:, I(3)) * S(3, :);
        X_av_0 = X / norm(X);

        % add noise to the data
        Data = snr(noise_i) * X_av_0 + Noise_av_0;
        % compute covariance in the virtual sensor space
        Ca = UP * Data * Data' * UP';

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

    end
        mc
end


R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix

coord1 = R(src(1), :);
coord2 = R(src(2), :);
coord3 = R(src(3), :);

dist = zeros(2, size(R_red, 1));
for i = 1:size(R_red, 1)
    dist(1, i) = norm(R_red(i, :) - coord1);
    dist(2, i) = norm(R_red(i, :) - coord2);
    dist(3, i) = norm(R_red(i, :) - coord3);
end
[~, ind1] = min(dist(1, :));
[~, ind2] = min(dist(2, :));
[~, ind3] = min(dist(3, :));

Z = Zp;
Inv.ImageGridAmp = Z';
%             Inv.ImageGridAmp = [];

Inv.ImageGridAmp = zeros(5001, 1);
Inv.ImageGridAmp([ind1, ind2, ind3]) = 100;

