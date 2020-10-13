G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channel locations

% 1.1 Generate simulations for different methods (MNE, LCMV beamformer, 
% ReciPSIICOS, WReciPSIICOS) and different SNRs

snr = [5];
Nmc = 1;
load_corr = true;
load_src = false;

[Z_total, picked_src] = save_simulations_snr(G3, G3_red, chans, ...
    snr, Nmc, load_corr, load_src);

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix

coord1 = R(picked_src(1), :);
coord2 = R(picked_src(2), :);

dist = zeros(2, size(R_red, 1));
for i = 1:size(R_red, 1)
    dist(1, i) = norm(R_red(i, :) - coord1);
    dist(2, i) = norm(R_red(i, :) - coord2);
end
[~, ind1] = min(dist(1, :));
[~, ind2] = min(dist(2, :));

Z = squeeze(Z_total(3, 2, :, :, :));
Inv.ImageGridAmp = Z;
%             Inv.ImageGridAmp = [];

Inv.ImageGridAmp = zeros(5001, 1);
Inv.ImageGridAmp([ind1, ind2]) = 100;

c = lines(7);
figure
subplot(2, 2, 1)
Z = squeeze(Z_total(1, 1, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(1, :), 'MarkerEdgeColor', 'k')
ylim([0 5])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

subplot(2, 2, 2)
Z = squeeze(Z_total(2, 1, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(6, :), 'MarkerEdgeColor', 'k')
ylim([0 5])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

subplot(2, 2, 3)
Z = squeeze(Z_total(3, 1, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(2, :), 'MarkerEdgeColor', 'k')
ylim([0 0.45])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

subplot(2, 2, 4)
Z = squeeze(Z_total(4, 1, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(4, :), 'MarkerEdgeColor', 'k')
ylim([0 8*10^(-4)])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')




c = lines(7);
figure
subplot(2, 2, 1)
Z = squeeze(Z_total(1, 2, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(1, :), 'MarkerEdgeColor', 'k')
ylim([0 5])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

subplot(2, 2, 2)
Z = squeeze(Z_total(2, 2, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(6, :), 'MarkerEdgeColor', 'k')
ylim([0 5])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

subplot(2, 2, 3)
Z = squeeze(Z_total(3, 2, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(2, :), 'MarkerEdgeColor', 'k')
ylim([0 5])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

subplot(2, 2, 4)
Z = squeeze(Z_total(4, 2, :, :, :));
scatter(1:5001, Z, 'MarkerFaceColor', c(4, :), 'MarkerEdgeColor', 'k')
ylim([0 8*10^(-4)])
xlim([0 5001])
xlabel('Source index')
ylabel('Activation amplitude, a.u.')

