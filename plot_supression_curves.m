
SensorsType = [1 2];
bLoad = 0;
GainSVDTh = 3.5e-3;
bPlot = 1;

col_rp = [7 110 187] / 255;
col_wrp = 0.9 * [74 198 255] / 255;

bNormalize = false;

for st = 1:length(SensorsType)
    fname_results = sprintf('SS_Suppression_nrm_%d_type_%d_th_%2.3f.mat', ...
        bNormalize, SensorsType(st), GainSVDTh);
    load(fname_results)
    
    figure
    xyi_rp = interparc(90, PpwrRat_rp, PcorrRat_rp, 'lin');
    h1 = plot(xyi_rp(1:3:end, 1), xyi_rp(1:3:end, 2), 'o');
    h1.MarkerFaceColor = col_rp;
    h1.Color = col_rp;
    h1.Marker = 'o';
    h1.MarkerSize = 10;
    h1.Parent.FontSize = 16;
    hold on
    xyi_wrp = interparc(90, PpwrRat_wrp, PcorrRat_wrp, 'lin');
    h2 = plot(xyi_wrp(1:3:end, 1), xyi_wrp(1:3:end, 2), 's');
    h2.MarkerFaceColor = col_wrp;
    h2.Color = col_wrp;
    h2.Marker = 's';
    h2.MarkerSize = 10;
    h2.LineWidth = 2;
    h3 = plot(xyi_rp(:, 1), xyi_rp(:, 2));
    h3.LineWidth = 2;
    h3.Color = col_rp;
    h3.Parent.FontSize = 16;
    h4 = plot(xyi_wrp(:, 1), xyi_wrp(:, 2));
    h4.LineWidth = 2;
    h4.Color = col_wrp;
    xlabel('Power subspace fraction','FontSize', 1);
    ylabel('Correlation subspace fraction','FontSize', 16);
    
    axis([0 1 0 1]);
    dPwr = diff(PpwrRat_rp);
    dCorr = diff(PcorrRat_rp);
    ind = find(dPwr ./ dCorr <= 0);
    rp_rat = log(dPwr(1:ind(1) - 1) ./ dCorr(1:ind(1) - 1));
    rnk_step = 2;
    xyi_rp = interparc(50, 1:2:rnk_step * length(rp_rat), rp_rat, 'lin');
    
    figure
    h5 = plot(xyi_rp(:, 1), xyi_rp(:, 2));
    h5.Color = col_rp;
    h5.LineWidth = 2;
    h3.Parent.FontSize = 16;
    hold on
    h6 = plot(xyi_rp(1:10:end, 1), xyi_rp(1:10:end, 2), 'o');
    h6.Color = col_rp;
    h6.LineWidth = 2;
    h6.Marker = 'o';
    h6.MarkerSize = 10;
    h6.MarkerFaceColor = col_rp;
    
    dPwr = diff(PpwrRat_wrp);
    dCorr = diff(PcorrRat_wrp);
    ind = find(dPwr ./ dCorr <= 0);
    rp_rat = log(dPwr(1:ind(1) - 1) ./ dCorr(1:ind(1) - 1));
    rnk_step = 2;
    xyi_wrp = interparc(50, 1:2:rnk_step * length(rp_rat), rp_rat, 'lin');
    hold on
    h7 = plot(xyi_wrp(:, 1), -xyi_wrp(:, 2));
    h7.Color = col_wrp;
    h7.LineWidth = 2;
    h3.Parent.FontSize = 18;
    hold on
    h8 = plot(xyi_wrp(1:10:end, 1), -xyi_wrp(1:10:end, 2), 's');
    h8.Color = col_wrp;
    h8.LineWidth = 2;
    h8.Marker = 's';
    h8.MarkerSize = 10;
    h8.MarkerFaceColor = col_wrp;

    xlabel('Projection rank, K','FontSize',18)
    ylabel('Logarithm of marginal gain','FontSize',18);
    axis tight
end
   
