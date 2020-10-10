function fname_results = save_subspace_suppression_curves(G3_red, chans, SensorsType, ... 
    bNormalize, bLoad, GainSVDTh, bPlot)
%--------------------------------------------------------------------------
% Parameters
% ----------
% GainSVDTh : 0.001 % Threshold of sensor count reduction step; 0.01 results
                      % into 50 eigensensors and makes it run faster but produces
                      % less contrasting subcorr scans
                      % for a more reliable preformance use 0.001 to get all the 
                      % sensor on board but be ready to wait;
% bLoad : 0 % load source covariance space correlation matrix
                      % set to false, if need to recompute. takes time!    
% SensorsType : 0     % 0 - MAG, 1 - GRAD, 2 - MAG & GRAD                      
% bNormalize : 0      % normalize topographies
%--------------------------------------------------------------------------

    if (nargin < 1)
        bPlot = true;
    end

    % 0. SETTINGS
    if (SensorsType == 0)
        RPRank = 300;  % reciPSIICOSprojection rank
        wRPRank = 500;  % whitened reciPSIICOS projection rank
        PwrRankMax = 400;
    elseif(SensorsType == 1) 
        RPRank = 300;  % reciPSIICOSprojection rank
        wRPRank = 1020;  % whitened reciPSIICOS projection rank
        PwrRankMax = 1500;
    elseif(SensorsType == 2)
        RPRank = 300;  % reciPSIICOSprojection rank
        wRPRank = 1020;  % whitened reciPSIICOS projection rank
        PwrRankMax = 1500;
    end

    % 1. FORWARD MODEL & channels selection
    % for sparse matrices
    if (SensorsType == 0)
        ChUsed = find(strcmp({chans.Channel.Type}, 'MEG MAG')); % use magnitometers only
    elseif (SensorsType == 1) 
        ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % use gradiometers only
    elseif (SensorsType == 2)
        ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD') | ...
                      strcmp({chans.Channel.Type}, 'MEG MAG')); % use gradiometers + magnitometers
        ChUsedMAG = find(strcmp({chans.Channel.Type}, 'MEG MAG')); % use magnitometers only
        ChUsedGRAD = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % use gradiometers only
        G3_red.Gain(ChUsedMAG, :) = G3_red.Gain(ChUsedMAG, :) * 20;
    end

    % 3D -> 2D tangent space, xxx0 stands for not normalized topographies
    [~, G2d0_red, ~] = G3toG2(G3_red, ChUsed);

    % 2. REDUCING SENSOR SPACE for 2 versions of the FM
    % non-normalized reduced src. space
    [ug, ~, ~] = spm_svd(G2d0_red * G2d0_red', GainSVDTh);
    UP0_red = ug'; % dimension reduction transform matrix

    % compute low dmensional sversion of the forward model
    G2d0U_red  = UP0_red * G2d0_red; 

    % 3. SET SOURCE COV. SPACE CORR. MATRIX
    if (bLoad == 1)
      fname = sprintf('Ccorr_nrm_%d_type_%d_th_%2.3f.mat', bNormalize, ...
          SensorsType, GainSVDTh);
      Tmp = load(fname, 'Ccorr');
      Ccorr = Tmp.Ccorr;
      clear Tmp;
    else % compute source covariance space correlation matrix (takes time!)   
      Ccorr = ComputeSynchSpaceCorrelation(G2d0U_red, bNormalize, true); 
      fname = sprintf('Ccorr_nrm_%d_type_%d_th_%2.3f.mat', bNormalize, ...
          SensorsType, GainSVDTh);
      save(fname, 'Ccorr');
    end

    % 4. COMPUTE reciPSIICOS PROJECTION MATRIX
    RankG = size(G2d0U_red, 1);
    
    % use non-normalized FM matrix to form PSIICOS projection operator 
    % normalization happens inside the function 
    [Upwr, ~, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U_red, ...
        PwrRankMax, 0, bNormalize);
    
    % Apwr - matrix of normalized source power space topographies, need for
    % wreciPSIICOS
    % Form reciPSIICOS projection matrix
    Pr_rp = Upwr(:, 1:RPRank) * Upwr(:, 1:RPRank)';

    % 7. COMPUTE WHITENED reciPSIICOS PROJECTION MATRIX  
    % form  power subspace correlation matrix
    Cpwr = 1 / size(Apwr, 2) * (Apwr * Apwr');
    % compute whitening wrt to power subspace operator 
    Wpwr = sqrtm(inv(Cpwr + 0.001 * trace(Cpwr) / (RankG^2) * eye(size(Cpwr))));
    % take source space cov. matrip into the whitened space
    WCcorrW = Wpwr * double(Ccorr) * Wpwr';
    % compute principal subspace of the whitened source corr. covariance matrix  
    [u, ~] = eigs(0.5 * WCcorrW + 0.5 * WCcorrW', size(WCcorrW, 1));
    UcorrW = u;
    
    % compute whitened reciPSCIICOS projection in the orgiganal space
    Pr_wrp = inv(Wpwr) * ...
        (eye(size(UcorrW, 1)) - UcorrW(:, 1:wRPRank) * UcorrW(:, 1:wRPRank)') * ...
        Wpwr;

    % 8. DO DIAGNOSTICS I
    % explore corr. suppression vs pwr. depletion phenomenon 
    clear PcorrRat_wrp PpwrRat_wrp PcorrRat_rp PpwrRat_rp Pr_wrp_rnk Pr_rp_rnk
    % reciPSIICOS
    clear  PcorrRat_rp PpwrRat_rp PrRnk_rp
    
    i_rnk = 1;
    for prrnk = 1:2:PwrRankMax
        Pr_rp_rnk = Upwr(:, 1:prrnk) * Upwr(:, 1:prrnk)';
        PcorrRat_rp(i_rnk) = trace(Pr_rp_rnk * double(Ccorr) * Pr_rp_rnk) / ...
                             trace(double(Ccorr));
        PpwrRat_rp(i_rnk) = trace(Pr_rp_rnk * Cpwr * Pr_rp_rnk) / trace(Cpwr);
        PrRnk_rp(i_rnk) = prrnk;
        i_rnk = i_rnk + 1;
    end

    % whitened reciPSIICOS
    i_rnk = 1;
    for prrnk = 1:2:size(UcorrW, 2)
        UcorrW_rnk = UcorrW(:, 1:prrnk);
        Pr_wrp_rnk = inv(Wpwr) * ...
                     (eye(size(UcorrW, 1)) - UcorrW(:, 1:prrnk) * UcorrW(:, 1:prrnk)') * ...
                     Wpwr;
        PcorrRat_wrp(i_rnk) = trace(Pr_wrp_rnk * double(Ccorr) * Pr_wrp_rnk) / ...
            trace(double(Ccorr));
        PpwrRat_wrp(i_rnk) = trace(Pr_wrp_rnk * Cpwr * Pr_wrp_rnk) / trace(Cpwr);
        i_rnk = i_rnk + 1;
    end

    % calculate the operating point corresponding to the selected projection ranks
    % whitened reciPSIICOS
    PcorrRat_wrp_cur = trace(Pr_wrp * double(Ccorr) * Pr_wrp) / trace(double(Ccorr));
    PpwrRat_wrp_cur = trace(Pr_wrp * Cpwr * Pr_wrp) / trace(Cpwr);
    % reciPSIICOS
    PcorrRat_rp_cur = trace(Pr_rp * double(Ccorr) * Pr_rp) / trace(double(Ccorr));
    PpwrRat_rp_cur = trace(Pr_rp * Cpwr * Pr_rp) / trace(Cpwr);
    fname_results = sprintf('SS_Suppression_nrm_%d_type_%d_th_%2.3f.mat', ...
        bNormalize, SensorsType, GainSVDTh);
    save(fname_results, 'PpwrRat_rp','PcorrRat_rp', 'PpwrRat_wrp','PcorrRat_wrp',...
        'PcorrRat_wrp_cur','PpwrRat_wrp_cur','PcorrRat_rp_cur','PpwrRat_rp_cur');

    % plot the suppression vs depletion curves and mark the operating points
    if (bPlot == 1)
        % plot results
        figure
        col_rp = [0 0.4470 0.7410];
        col_wrp = 1.3*[0 0.4470 0.7410];
        h1 = plot(PpwrRat_rp([1:15, 20:10:length(PpwrRat_rp), length(PpwrRat_rp)]), ...
            PcorrRat_rp([1:15, 20:10:length(PpwrRat_rp), length(PpwrRat_rp)]), 'o');
        h1.MarkerFaceColor = col_rp;
        h1.Color = col_rp;
        h1.Marker = 'o';
        h1.MarkerSize = 8;
        hold on
        h2 = plot(PpwrRat_wrp([1:15, 20:10:length(PpwrRat_wrp), length(PpwrRat_wrp)]), ...
            PcorrRat_wrp([1:15, 20:10:length(PpwrRat_wrp), length(PpwrRat_wrp)]), 'o');
        h2.MarkerFaceColor = col_wrp;
        h2.Color = col_wrp;
        h2.Marker = 's';
        h2.MarkerSize = 8;
        h3 = plot(PpwrRat_rp, PcorrRat_rp);
        h3.LineWidth = 2;
        h3.Color = col_rp;
        h4 = plot(PpwrRat_wrp, PcorrRat_wrp, '.-');
        h4.LineWidth = 2;
        h4.Color = col_wrp;
        xlabel('Fraction of power subspace energy left');
        ylabel('Fraction of corr. subspace energy left');
        grid on
        legend('Prj. into power subspace', 'Whitened prj. away from corr. subspace ');

        close all
        dPwr  = diff(PpwrRat_rp);
        dCorr = diff(PcorrRat_rp);
        ind = find(dPwr ./ dCorr <= 0);
        rp_rat = log(dPwr) - log(dCorr);
        rp_rat(ind) = NaN;
        h5 = plot(rp_rat);
        hold on
        dPwr = diff(PpwrRat_wrp);
        dCorr = diff(PcorrRat_wrp);
        ind = find(dPwr ./ dCorr <= 0);
        wrp_rat = log(dPwr) - log(dCorr);
        wrp_rat(ind) = NaN;
        h5 = plot(wrp_rat);
    end
end
 
   
  
