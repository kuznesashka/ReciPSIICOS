% 1. FORWARD MODEL

G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channels

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix
ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % set to use gradiometers only

[Nch, Nsites] = size(G3.Gain(ChUsed,1:3:end));
[Nch, Nsites_red] = size(G3_red.Gain(ChUsed,1:3:end));

G_pure = G3.Gain(ChUsed,:); % 2D dense forward matrix 
G2d = zeros(Nch,Nsites*2);
G2d0 = zeros(Nch,Nsites*2);
range = 1:2;
for i = 1:Nsites
    g = [G_pure(:,1+3*(i-1)) G_pure(:,2+3*(i-1)) G_pure(:,3+3*(i-1))];
    [u sv v] = svd(g);
    gt = g*v(:,1:2);
    G2d(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
    G2d0(:,range) = gt;
    range = range + 2;
end

G = G3_red.Gain(ChUsed,:); % 2D sparse forward matrix
G2d_red = zeros(Nch,Nsites_red*2);
G2d0_red = zeros(Nch,Nsites_red*2);
range = 1:2;
for i = 1:Nsites_red
    g = [G(:,1+3*(i-1)) G(:,2+3*(i-1)) G(:,3+3*(i-1))];
    [u sv v] = svd(g);
    gt = g*v(:,1:2);
    G2d_red(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
    G2d0_red(:,range) = gt;
    range = range + 2;
end


% 2. REDUCING SENSOR SPACE for sparse matrix
GainSVDTh = 0.001; % 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
                   % for a more reliable preformance use 0.01 to get all the sensor on board but be ready to wait;
[ug sg vg] = spm_svd(G2d_red*G2d_red',GainSVDTh);
UP = ug'; % direction of dimention reduction
G2dU_red = UP*G2d_red;
G2d0U_red = UP*G2d0_red;

RankG = size(G2d0U_red,1);

% now find span of the target space
bLoad = true;
Nsrc = size(G2d0U_red,2)/2;
Swp = zeros(4);
Swp(2,3) = 1;
Swp(3,2) = 1;
Swp(1,1) = 1; 
Swp(4,4) = 1;
RankG = size(G2d0U_red,1);
NSites = fix( size(G2d0U_red,2)/2);

if(bLoad)
  Wks = load('C_re.mat','C_re');
% this is correlation subspace correlation matrix
   Ccorr = Wks.C_re;
else
    parfor i = 1:NSites
         range_i = i*2-1:i*2;
         ai = G2d0U_red(:,range_i);
         rng = 1:4;
         X = zeros(RankG^2,4*(NSites-i),'single');
         for j = i+1:NSites
            range_j = j*2-1:j*2;
            aj = G2d0U_red(:,range_j);
            X(:,rng) = kron(ai,aj)+kron(aj,ai)*Swp;
            rng = rng+4;
         end
        C_re = C_re + X*X';
        i
    end
    save C_re;
end

% 3. PSIICOS projection for sparse matrix
[Upwr, ds, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, 1500, 0);
Rnk = 500;

% this is power subspace correlation matrix
Cpwr = Apwr*Apwr';

Wks =  load('C_re.mat','C_re');
% this is correlation subspace correlation matrix
Ccorr = Wks.C_re;

%this is a whitener wrt to power subspace 
Wpwr = sqrtm(inv(Cpwr+0.001*trace(Cpwr)/(RankG^2)*eye(size(Cpwr))));

WCcorrW = Wpwr*double(Ccorr)*Wpwr';

ProjRnkMax = 1280;
[u s] = eigs(WCcorrW,ProjRnkMax);
UcorrW_rnk = u(:,1:ProjRnkMax);

[u s] = eigs(double(Ccorr),ProjRnkMax);
Ucorr_rnk = u(:,1:ProjRnkMax);

%this is a projector away from corr. subspace with whitening w.r.t to power subspace  
PrFromCorr_W = inv(Wpwr+0.05*trace(Wpwr)/size(Wpwr,1))*(eye(size(UcorrW_rnk,1))-UcorrW_rnk(:,1:1000)*UcorrW_rnk(:,1:1000)')*Wpwr;
PrPwr = Upwr(:, 1:500)*Upwr(:, 1:500)';


GG = zeros(size(G2dU_red,1)^2, Nsites_red*4);

range2d = 1:2;
range4d = 1:4;
for i = 1:Nsites_red
    g = G2dU_red(:,range2d);
    GG(:,range4d) = [kron(g(:,1),g(:,1)),kron(g(:,1),g(:,2)),kron(g(:,2),g(:,1)) kron(g(:,2),g(:,2))];
    range2d = range2d+2;
    range4d = range4d+4;
    i
end
PrGG = PrFromCorr_W*GG;
PrPwrGG = PrPwr*GG;

normGG = sqrt(sum(GG.^2,1));
normPrGG = sqrt(sum(PrGG.^2,1));
normPrPwrGG = sqrt(sum(PrPwrGG.^2,1));
figure
plot(normPrPwrGG./normGG)
hold on
plot(normPrGG./normGG)

% 4. Simulations

% left and right hemispheres index
src_left = find(R(:,2)>0);
src_right = find(R(:,2)<0);
src_left_red = find(R_red(:,2)>0);
src_right_red = find(R_red(:,2)<0);

snr = [5]; % snr level in the data
c = lines(7); % colors
Fs = 500; % sampling frequency
Ntr = 50; % number of simulated trials
T = Fs; % number of time points in one trial
t = 1:T;
frac = 0.65; % threshold for amplitude
plotfigs = false;
symmetrical = 1;

for noise_i = 1:length(snr) 
    for mc = 1:100
      
        % generate brain noise with dense matrix
        range = 1:T;
        for tr = 1:Ntr
            Noise(:,range) = GenerateBrainNoise(G2d0,T,200,500,Fs);
            range = range+T;
        end
        Noise_av = mean(reshape(Noise,[Nch, T, Ntr]), 3); % average by trial
        Noise_av_0 = Noise_av/norm(Noise_av); % normalized average noise

        % generate signal from dense forward model matrix
        
        Farsites = find(R(:,2)>0.02); % sources from left hemisphere >2 cm far from midline (6542/10001)
                
        ii = randperm(length(Farsites));
        src = Farsites(ii(1));  % pick the first source from the left hemisphere >2 cm far from midline 

        if symmetrical == 1
            for i = 1:Nsites % find the symmetrical source in the right hemisphere
                dd_sym(i) = norm([R(src,1),-R(src,2), R(src,3)]-R(i,:));
            end
            [val, ind_sym] = min(dd_sym);
        else
            ii = randperm(length(src_right));
            ind_sym = ii(1);
        end
            
        ind_generated = [src, ind_sym]; % active symmetrical sources
        I = ind_generated*2;
                
        range = 1:T; % activation functions
        for tr = 1:Ntr
            S(1,range) = sin(2*pi*t/Fs + randn*0.1);
            S(2,range) = sin(2*pi*t/Fs + randn*0);
            range = range+T;
        end
        X = G2d0(:,I(1))*S(1,:) + G2d0(:,I(2))*S(2,:); % signal, activations only by y-axis of G
        X_av = mean(reshape(X, [Nch, T, Ntr]), 3);
        
        X_av_0 = X_av/norm(X_av);
        Data  = snr(noise_i)*X_av_0 + Noise_av_0; % add noise to the data
        
        Ca = UP*Data*Data'*UP'; % compute covariance in the virtual sensor space

        % LCMV BF
        iCa = tihinv(Ca, 0.01);
        range2d = 1:2;
        for i = 1:Nsites_red
            g = G2dU_red(:,range2d);
            m = inv(g'*iCa*g);
            [u ss v] = svd(m);
            Z(i) = ss(1,1);
            range2d = range2d+2;
        end

        % PSIICOS modified beamformer
        % create projection matrix
        
        Cap =reshape(PrFromCorr_W*Ca(:), size(Ca));
        [e a] = eig(Cap);
        Cap = e*abs(a)*e';
        iCap = tihinv(Cap, 0.01);
        
        range2d = 1:2;
        for i=1:Nsites_red
            g = G2dU_red(:,range2d);
            m = inv(g'*iCap*g);
            [u ss v] = svd(m);
            Zpw(i) = ss(1,1);
            range2d = range2d+2;
        end
        
        % PSIICOS modified beamformer
       
        Cap = reshape(PrPwr*Ca(:), size(Ca));
        [e a] = eig(Cap);
        Cap = e*abs(a)*e';
        iCap = tihinv(Cap, 0.01);
        
        range2d = 1:2;
        for i=1:Nsites_red
            g = G2dU_red(:,range2d);
            m = inv(g'*iCap*g);
            [u ss v] = svd(m);
            Zp(i) = ss(1,1);
            range2d = range2d+2;
        end
       
 
%       Minimum norm estimate

        lambd = 0.1;
        Cs = eye(Nsites_red*2);
        Cn = eye(size(G2dU_red,1));
        GGt = G2dU_red*Cs*G2dU_red';
        Wmne = G2dU_red'*inv(G2dU_red*Cs*G2dU_red'+lambd*trace(GGt)/size(GGt,1)*Cn);

        clear Zmne;
        range2d = 1:2;
        for i = 1:Nsites_red
            w = Wmne(range2d,:);
            m = w*Ca*w';
            [u ss v] = svd(m);
            Zmne(i) = ss(1,1);
            range2d = range2d+2;
        end
        
        
        % find the closest sources in the reduced model
                
        for j = 1:length(src_left_red)
            dd_left(j) = norm(R_red(src_left_red(j),:)-R(ind_generated(1),:));
        end        
        [val, ind_new_left] = min(dd_left);
        ind_left_red = src_left_red(ind_new_left);
       
        for j = 1:length(src_right_red)            
            dd_right(j) = norm(R_red(src_right_red(j),:)-R(ind_generated(2),:));
        end
        [val, ind_new_right] = min(dd_right);
        ind_right_red = src_right_red(ind_new_right);
        
        ind_generated_red = [ind_left_red, ind_right_red];
        XYZGen_red = R_red(ind_generated_red,:);
        
%         figure
%         scatter3(R(:,1),R(:,2),R(:,3))
%         hold on
%         scatter3(R(ind_generated,1),R(ind_generated,2),R(ind_generated,3), 'r', 'filled')
%         scatter3(R_red(ind_generated_red,1),R_red(ind_generated_red,2),R_red(ind_generated_red,3), 'g', 'filled')
        %range = 1:Nsites_red;
        if(plotfigs)
            figure
            subplot(2,2,1)
            stem(Zmne, 'LineWidth', 2)
            hold on
            plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(2,:))
            plot([ind_generated_red(2),ind_generated_red(2)], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(2,:))
            %stem(range(Zmne>frac*max(Zmne)),Zmne(Zmne>frac*max(Zmne)))
            ylim([0 max(Zmne)])
            title('mne')

            subplot(2,2,2)
            stem(squeeze(Z(1,:,:)), 'LineWidth', 2)
            hold on
            plot([ind_generated_red(1),ind_generated_red(1)], [0,max(squeeze(Z(1,:,:)))], 'LineWidth', 2, 'Color', c(2,:))
            plot([ind_generated_red(2),ind_generated_red(2)], [0,max(squeeze(Z(1,:,:)))], 'LineWidth', 2, 'Color', c(2,:))
            ylim([0 max(Z(1,:,:))])
            title('LCMV')

            subplot(2,2,3)
            stem(squeeze(Zp(1,:,:)), 'LineWidth', 2)
            hold on
            plot([ind_generated_red(1),ind_generated_red(1)], [0,max(squeeze(Zp(1,:,:)))], 'LineWidth', 2, 'Color', c(2,:))
            plot([ind_generated_red(2),ind_generated_red(2)], [0,max(squeeze(Zp(1,:,:)))], 'LineWidth', 2, 'Color', c(2,:))
            ylim([0 max(Zp(1,:,:))])
            title('PSIICOS')
            
            subplot(2,2,4)
            stem(squeeze(Zpw(1,:,:)), 'LineWidth', 2)
            hold on
            plot([ind_generated_red(1),ind_generated_red(1)], [0,max(squeeze(Zpw(1,:,:)))], 'LineWidth', 2, 'Color', c(2,:))
            plot([ind_generated_red(2),ind_generated_red(2)], [0,max(squeeze(Zpw(1,:,:)))], 'LineWidth', 2, 'Color', c(2,:))
            ylim([0 max(Zpw(1,:,:))])
            title('Whitened PSIICOS')
        end
        
       
        % for MNE
        clear dist_l dist_r
        Zmne_left = Zmne(src_left_red); % values from left hemisphere
        Zmne_right = Zmne(src_right_red); % values from right hemisphere
        [max_val_l max_l] = max(Zmne_left); % maximum in the left hemisphere
        max_ind_l = src_left_red(max_l);
        [max_val_r max_r] = max(Zmne_right);
        max_ind_r = src_right_red(max_r);
        
        
%         figure
%         subplot(1,2,1)
%         stem(src_left_red, Zmne_left, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(2,:))
%         plot([max_ind_l,max_ind_l], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(3,:))
%         ylim([0 max(Zmne)])
%         title('Left hemisphere')
%         subplot(1,2,2)
%         stem(src_right_red, Zmne_right, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(2),ind_generated_red(2)], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(2,:))
%         plot([max_ind_r,max_ind_r], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(3,:))
%         title('Right hemisphere')
%         ylim([0 max(Zmne)])
      
        max_val = max(max_val_l, max_val_r);
        Zmne_left(Zmne_left<frac*max_val) = 0;
        Zmne_right(Zmne_right<frac*max_val) = 0;

%         figure
%         subplot(1,2,1)
%         stem(Zmne_left, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(2,:))
%         ylim([0 max(Zmne)])
%         subplot(1,2,2)
%         stem(Zmne_right, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Zmne)], 'LineWidth', 2, 'Color', c(2,:))
%         ylim([0 max(Zmne)])
        
        
        if (length(Zmne_left(Zmne_left>frac*max_val)) == 0)|(length(Zmne_right(Zmne_right>frac*max_val)) == 0)
            r_mne(noise_i, mc) = NaN;
            var_mne(noise_i, mc) = NaN;
        else
 
            r_mne(noise_i, mc) = (norm(R_red(max_ind_l,:)-R_red(ind_generated_red(1),:))+ ...
                norm(R_red(max_ind_r,:)-R_red(ind_generated_red(2),:)))/2;

            Zmne_left_norm = Zmne_left./sum(Zmne_left);
            Zmne_right_norm = Zmne_right./sum(Zmne_right);

            coord_active_l = R_red(src_left_red(Zmne_left_norm>0),:);
            for i = 1:size(coord_active_l,1)
                dist_l(i) = norm(coord_active_l(i,:)-R_red(max_ind_l,:));
            end

            coord_active_r = R_red(src_right_red(Zmne_right_norm>0),:);
            for i = 1:size(coord_active_r,1)
                dist_r(i) = norm(coord_active_r(i,:)-R_red(max_ind_r,:));
            end

            var_mne(noise_i, mc) = (sum(Zmne_left_norm(Zmne_left_norm>0).*dist_l)+ ...
                sum(Zmne_right_norm(Zmne_right_norm>0).*dist_r))/2;
        end

      
        
        % for LCMV
        clear dist_l dist_r
        Z_left = Z(src_left_red); % values from left hemisphere
        Z_right = Z(src_right_red); % values from right hemisphere
        [max_val_l max_l] = max(Z_left); % maximum in the left hemisphere
        max_ind_l = src_left_red(max_l);
        [max_val_r max_r] = max(Z_right);
        max_ind_r = src_right_red(max_r);
        
         
%         figure
%         subplot(1,2,1)
%         stem(src_left_red, Z_left, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Z)], 'LineWidth', 2, 'Color', c(2,:))
%         plot([max_ind_l,max_ind_l], [0,max(Z)], 'LineWidth', 2, 'Color', c(3,:))
%         ylim([0 max(Z)])
%         title('Left hemisphere')
%         subplot(1,2,2)
%         stem(src_right_red, Z_right, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(2),ind_generated_red(2)], [0,max(Z)], 'LineWidth', 2, 'Color', c(2,:))
%         plot([max_ind_r,max_ind_r], [0,max(Z)], 'LineWidth', 2, 'Color', c(3,:))
%         title('Right hemisphere')
%         ylim([0 max(Z)])
%       
        max_val = max(max_val_l, max_val_r);
        Z_left(Z_left<frac*max_val) = 0;
        Z_right(Z_right<frac*max_val) = 0;

%         figure
%         stem(Z_left, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Z)], 'LineWidth', 2, 'Color', c(2,:))
%         ylim([0 max(Zmne)])
        if (length(Z_left(Z_left>frac*max_val)) == 0)|(length(Z_right(Z_right>frac*max_val)) == 0)
            r_bf(noise_i, mc) = NaN;
            var_bf(noise_i,mc) = NaN;
        else
 
            r_bf(noise_i,mc) = (norm(R_red(max_ind_l,:)-R_red(ind_generated_red(1),:))+ ...
                norm(R_red(max_ind_r,:)-R_red(ind_generated_red(2),:)))/2;

        
            Z_left_norm = Z_left./sum(Z_left);
            Z_right_norm = Z_right./sum(Z_right);

            coord_active_l = R_red(src_left_red(Z_left_norm>0),:);
            for i = 1:size(coord_active_l,1)
                dist_l(i) = norm(coord_active_l(i,:)-R_red(max_ind_l,:));
            end

            coord_active_r = R_red(src_right_red(Z_right_norm>0),:);
            for i = 1:size(coord_active_r,1)
                dist_r(i) = norm(coord_active_r(i,:)-R_red(max_ind_r,:));
            end

            var_bf(noise_i,mc) = (sum(Z_left_norm(Z_left_norm>0).*dist_l)+ ...
                sum(Z_right_norm(Z_right_norm>0).*dist_r))/2;
        end

        
         % for Antipsicos
        clear dist_l dist_r
        Zp_left = Zp(src_left_red); % values from left hemisphere
        Zp_right = Zp(src_right_red); % values from right hemisphere
        [max_val_l max_l] = max(Zp_left); % maximum in the left hemisphere
        max_ind_l = src_left_red(max_l);
        [max_val_r max_r] = max(Zp_right);
        max_ind_r = src_right_red(max_r);
        
      
%         figure
%         subplot(1,2,1)
%         stem(src_left_red, Zp_left, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Zp)], 'LineWidth', 2, 'Color', c(2,:))
%         plot([max_ind_l,max_ind_l], [0,max(Zp)], 'LineWidth', 2, 'Color', c(3,:))
%         ylim([0 max(Zp)])
%         title('Left hemisphere')
%         subplot(1,2,2)
%         stem(src_right_red, Zp_right, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(2),ind_generated_red(2)], [0,max(Zp)], 'LineWidth', 2, 'Color', c(2,:))
%         plot([max_ind_r,max_ind_r], [0,max(Zp)], 'LineWidth', 2, 'Color', c(3,:))
%         title('Right hemisphere')
%         ylim([0 max(Zp)])
%       
        max_val = max(max_val_l, max_val_r);
        Zp_left(Zp_left<frac*max_val) = 0;
        Zp_right(Zp_right<frac*max_val) = 0;

%         figure
%         stem(Z_left, 'LineWidth', 2)
%         hold on
%         plot([ind_generated_red(1),ind_generated_red(1)], [0,max(Z)], 'LineWidth', 2, 'Color', c(2,:))
%         ylim([0 max(Zmne)])
        if (length(Zp_left(Zp_left>frac*max_val)) == 0)|(length(Zp_right(Zp_right>frac*max_val)) == 0)
            r_ps(noise_i,mc) = NaN;
            var_ps(noise_i,mc) = NaN;
        else
 
            r_ps(noise_i,mc) = (norm(R_red(max_ind_l,:)-R_red(ind_generated_red(1),:))+ ...
                norm(R_red(max_ind_r,:)-R_red(ind_generated_red(2),:)))/2;
            Zp_left_norm = Zp_left./sum(Zp_left);
            Zp_right_norm = Zp_right./sum(Zp_right);

            coord_active_l = R_red(src_left_red(Zp_left_norm>0),:);
            for i = 1:size(coord_active_l,1)
                dist_l(i) = norm(coord_active_l(i,:)-R_red(max_ind_l,:));
            end

            coord_active_r = R_red(src_right_red(Zp_right_norm>0),:);
            for i = 1:size(coord_active_r,1)
                dist_r(i) = norm(coord_active_r(i,:)-R_red(max_ind_r,:));
            end

            var_ps(noise_i,mc) = (sum(Zp_left_norm(Zp_left_norm>0).*dist_l)+ ...
                sum(Zp_right_norm(Zp_right_norm>0).*dist_r))/2;
        end

         % for Antipsicos
        clear dist_l dist_r
        Zpw_left = Zpw(src_left_red); % values from left hemisphere
        Zpw_right = Zpw(src_right_red); % values from right hemisphere
        [max_val_l max_l] = max(Zpw_left); % maximum in the left hemisphere
        max_ind_l = src_left_red(max_l);
        [max_val_r max_r] = max(Zpw_right);
        max_ind_r = src_right_red(max_r);
        
        max_val = max(max_val_l, max_val_r);
        Zpw_left(Zpw_left<frac*max_val) = 0;
        Zpw_right(Zpw_right<frac*max_val) = 0;

        if (length(Zpw_left(Zpw_left>frac*max_val)) == 0)|(length(Zpw_right(Zpw_right>frac*max_val)) == 0)
            r_psw(noise_i,mc) = NaN;
            var_psw(noise_i,mc) = NaN;
        else
 
            r_psw(noise_i,mc) = (norm(R_red(max_ind_l,:)-R_red(ind_generated_red(1),:))+ ...
                norm(R_red(max_ind_r,:)-R_red(ind_generated_red(2),:)))/2;
            Zpw_left_norm = Zpw_left./sum(Zpw_left);
            Zpw_right_norm = Zpw_right./sum(Zpw_right);

            coord_active_l = R_red(src_left_red(Zpw_left_norm>0),:);
            for i = 1:size(coord_active_l,1)
                dist_l(i) = norm(coord_active_l(i,:)-R_red(max_ind_l,:));
            end

            coord_active_r = R_red(src_right_red(Zpw_right_norm>0),:);
            for i = 1:size(coord_active_r,1)
                dist_r(i) = norm(coord_active_r(i,:)-R_red(max_ind_r,:));
            end

            var_psw(noise_i,mc) = (sum(Zpw_left_norm(Zpw_left_norm>0).*dist_l)+ ...
                sum(Zpw_right_norm(Zpw_right_norm>0).*dist_r))/2;
        end
    
    mc
    
    end

    
end


minimum_r = min([min(r_mne), min(r_bf), min(r_ps)]);
maximum_r = max([max(r_mne), max(r_bf), max(r_ps)]);
minimum_v = min([min(var_mne), min(var_bf), min(var_ps)]);
maximum_v = max([max(var_mne), max(var_bf), max(var_ps)]);


figure
subplot(4,2,1)
histogram(r_mne, 20)
xlim([minimum_r, maximum_r])
ylim([0 25])
title('Bias, MNE')
subplot(4,2,3)
histogram(r_bf,20)
xlim([minimum_r, maximum_r])
ylim([0 25])
title('Bias, LCMV')
subplot(4,2,5)
histogram(r_ps,20)
xlim([minimum_r, maximum_r])
ylim([0 25])
title('Bias, APSIICOS')
subplot(4,2,7)
histogram(r_psw,20)
xlim([minimum_r, maximum_r])
ylim([0 25])
title('Bias, WAPSIICOS')
subplot(4,2,2)
histogram(var_mne,20)
xlim([minimum_v, maximum_v])
ylim([0 25])
title('PointSpread, MNE')
subplot(4,2,4)
histogram(var_bf,20)
xlim([minimum_v, maximum_v])
ylim([0 25])
title('PointSpread, LCMV')
subplot(4,2,6)
histogram(var_ps,20)
xlim([minimum_v, maximum_v])
ylim([0 25])
title('PointSpread, APSIICOS')
subplot(4,2,8)
histogram(var_psw,20)
xlim([minimum_v, maximum_v])
ylim([0 25])
title('PointSpread, WAPSIICOS')

sum(isnan(r_psw))
