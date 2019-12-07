% 1. FORWARD MODEL
% for dense and sparse matrices

G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channels

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix
ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % set to use gradiometers only
Nch = length(ChUsed);

Noise_forw_ratio = [0.05, 0.1, 0.2];
Nsites_red = size(G3_red.Gain, 2)/3;
noise_forw = rand(Nch, 3*Nsites_red) - 0.5;
Gain_in = G3_red.Gain(ChUsed,:);

for fw_error = 1:length(Noise_forward_ratio)
    G3_red.Gain(ChUsed,:) = G3_red.Gain(ChUsed,:) + Noise_forw_ratio(fw_error)*noise_forw*...
        norm(G3_red.Gain(ChUsed,:),'fro')/norm(noise_forw,'fro');

    % norm(Noise_forw_ratio*noise_forw*norm(G3_red.Gain(ChUsed,:),'fro')/norm(noise_forw,'fro'), 'fro')/norm(Gain_in, 'fro')

    [G2d, G2d0, Nsites] = G3toG2(G3, ChUsed);
    [G2d_red, G2d0_red, Nsites_red] = G3toG2(G3_red, ChUsed);

    % 2. REDUCING SENSOR SPACE for sparse matrix
    GainSVDTh = 0.001; % 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
                   % for a more reliable preformance use 0.01 to get all the sensor on board but be ready to wait;
    [ug sg vg] = spm_svd(G2d_red*G2d_red',GainSVDTh);
    UP = ug'; % direction of dimention reduction
    G2dU_red = UP*G2d_red;
    G2d0U_red = UP*G2d0_red;

    % 3. PROJECTION MATRICES
    RankG = size(G2d0U_red,1);

    % now find span of the target space
    bLoad = true;
    Nsrc = size(G2d0U_red, 2)/2;
    Swp = zeros(4);
    Swp(2, 3) = 1;
    Swp(3, 2) = 1;
    Swp(1, 1) = 1; 
    Swp(4, 4) = 1;
    RankG = size(G2d0U_red, 1);
    NSites = fix(size(G2d0U_red, 2)/2);

    errname = string(Noise_forw_ratio(fw_error))
    errname = erase(errname, '.')
    C_re = zeros(RankG^2, RankG^2);
    if(bLoad)
      gname = strcat('C_re_', errname, '.mat')
      Wks =  load(strcat('C_re_', errname, '.mat'),'C_re');
        % this is correlation subspace correlation matrix
       Ccorr = Wks.C_re;
    else
        parfor i = 1:NSites
             range_i = i*2-1:i*2;
             ai = G2d0U_red(:, range_i);
             rng = 1:4;
             X = zeros(RankG^2, 4*(NSites-i), 'single');
             for j = i+1:NSites
                range_j = j*2-1:j*2;
                aj = G2d0U_red(:, range_j);
                X(:, rng) = kron(ai,aj) + kron(aj,ai)*Swp;
                rng = rng + 4;
             end
            C_re = C_re + X*X';
            i
        end
        save C_re_005 С_re
    end

    % 3. PSIICOS projection for sparse matrix
    [Upwr, ds, Apwr]  = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, 1500, 0);
    Rnk = 500;

    % this is power subspace correlation matrix
    Cpwr = Apwr*Apwr';

    % Wks =  load('C_re_G_02.mat','C_re');
    % this is correlation subspace correlation matrix
    % Ccorr = Wks.C_re;

    %this is a whitener wrt to power subspace 
    Wpwr = sqrtm(inv(Cpwr+0.001*trace(Cpwr)/(RankG^2)*eye(size(Cpwr))));
    WCcorrW = Wpwr*double(Ccorr)*Wpwr';

    ProjRnkMax = 1280;
    [u s] = eigs(WCcorrW,ProjRnkMax);
    UcorrW_rnk = u(:,1:ProjRnkMax);

    [u s] = eigs(double(Ccorr),ProjRnkMax);
    Ucorr_rnk = u(:,1:ProjRnkMax);

    PrFromCorr_W = inv(Wpwr+0.05*trace(Wpwr)/size(Wpwr,1))*(eye(size(UcorrW_rnk,1))-UcorrW_rnk(:,1:500)*UcorrW_rnk(:,1:500)')*Wpwr;
    PrPwr = Upwr(:,1:500)*Upwr(:,1:500)';

    % 5. SIMULATIONS
    snr = [5]; % snr level in the data
    c = lines(7); % colors
    Fs = 500; % sampling frequency
    Ntr = 100; % number of simulated trials
    T = Fs; % number of time points in one trial
    t = 1:T;
    Nmc = 500;

    Zp = zeros(2, length(snr), Nmc, Nsites_red);
    Zpw = zeros(2, length(snr), Nmc, Nsites_red);
    Zbf = zeros(2, length(snr), Nmc, Nsites_red);
    Zmne = zeros(2, length(snr), Nmc, Nsites_red);

    load picked_src

    for mc = 1:Nmc
        % generate brain noise with dense matrix
        range = 1:T;
        for tr = 1:Ntr
            Noise(:, range) = GenerateBrainNoise(G2d0, T, 200, 500, Fs);
            range = range + T;
        end
        Noise_av = mean(reshape(Noise, [Nch, T, Ntr]), 3); % average by trial
        Noise_av_0 = Noise_av/norm(Noise_av); % normalized average noise

        % generate signal for dense forward model matrix

    %     far_sites = find(R(:,2)>0.02); % sources from left hemisphere >2 cm far from midline (6542/10001)
    %     rand_idx = randperm(length(far_sites));
    %     picked_src(mc,1) = far_sites(rand_idx(1));  % pick the first source from the left hemisphere >2 cm far from midline 
    % 
    %     clear dd_sym
    %     for i = 1:Nsites % find the symmetrical source
    %         dd_sym(i) = norm([R(picked_src(mc,1),1), ...
    %             -R(picked_src(mc,1),2), R(picked_src(mc,1),3)]-R(i,:));
    %     end
    %     [val, picked_src(mc,2)] = min(dd_sym);
        picked_src_oriented = picked_src(mc,:)*2;

        %         figure
        %         scatter3(R_red(:,1),R_red(:,2),R_red(:,3))
        %         hold on
        %         scatter3(R(picked_src,1), R(picked_src,2), R(picked_src,3), 'filled', 'r')

        range = 1:T; % activation functions
        for tr = 1:Ntr
            S(1,range) = sin(2*pi*t/Fs + randn*0.1);
            S(2,range) = sin(2*pi*t/Fs + randn*0);
            S(3,range) = cos(2*pi*t/Fs + randn*0);
            range = range+T;
        end

        for noise_i = 1:length(snr) 
            for synch = 1:2
                X = G2d0(:,picked_src_oriented(1))*S(synch,:) + G2d0(:,picked_src_oriented(2))*S(synch+1,:);
                X_av = mean(reshape(X, [Nch, T, Ntr]), 3);

                X_av_0 = X_av/norm(X_av);
                Data  = snr(noise_i)*X_av_0 + Noise_av_0; % add noise to the data
                Ca = UP*Data*Data'*UP'; % compute covariance in the virtual sensor space

                % LCMV BF
                Zbf(synch, noise_i, mc, :) = lcmv(G2dU_red, Ca);

                % Whitened ReciPSIICOS beamformer
                Cap =reshape(PrPwr*PrFromCorr_W*Ca(:), size(Ca));
                [e a] = eig(Cap);
                Cap = e*abs(a)*e';
                iCap = tihinv(Cap, 0.01);

                range2d = 1:2;
                for i=1:Nsites_red
                    g = G2dU_red(:,range2d);
                    m = inv(g'*iCap*g);
                    [u ss v] = svd(m);
                    Zpw(synch, noise_i, mc, i) = ss(1,1);
                    range2d = range2d+2;
                end

                % ReciPSIICOS beamformer
                Cap = reshape(PrPwr*Ca(:), size(Ca));
                [e a] = eig(Cap);
                Cap = e*abs(a)*e';
                iCap = tihinv(Cap, 0.01);

                range2d = 1:2;
                for i=1:Nsites_red
                    g = G2dU_red(:,range2d);
                    m = inv(g'*iCap*g);
                    [u ss v] = svd(m);
                    Zp(synch, noise_i, mc, i) = ss(1,1);
                    range2d = range2d+2;
                end

                % Minimum norm estimate
                lambda = 0.1;
                Zmne(synch, noise_i, mc, :) = mne(G2dU_red, Ca, lambda);
            end
        end
        mc
    end

    Z_total_02 = zeros(4, 2, 1, 500, 5001);
    Z_total_02(1,:,:,:,:) = Zp;
    Z_total_02(2,:,:,:,:) = Zpw;
    Z_total_02(3,:,:,:,:) = Zbf;
    Z_total_02(4,:,:,:,:) = Zmne;

%     save(strcat('Z_total_', errname, '.mat'), Z_total_02)
end

