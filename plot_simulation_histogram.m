% Input:
% r -- array with all biases of size (number_of_methods x 2 x
%   number_of_snrs x 1 x number_of_simulations), where:
%   number_of_methods = 4 (ReciPSIICOS, wReciPSIICOS, LCMV, MNE),
%   2 corresponds to synchronous/asynchronous
%   among all snrs used you need to choose the desired index as noise_i
% var -- array with all spreading values, the same atructure as r

% BIG FIGURE
c = lines(8);
noise_i = 4;
mc = size(r, 5);

% SYNCHRONOUS SOURCES
figure
subplot(4,1,1)
histogram(r(1,1,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(1,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(1,1,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(1,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)

subplot(4,1,2)
histogram(r(2,1,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(6,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(2,1,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(6,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)

subplot(4,1,3)
histogram(r(3,1,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(2,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(3,1,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(2,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)

subplot(4,1,4)
histogram(r(4,1,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(4,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(4,1,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(4,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)


% ASYNCHRONOUS SOURCES

figure
subplot(4,1,1)
histogram(r(1,2,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(1,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(1,2,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(1,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)

subplot(4,1,2)
histogram(r(2,2,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(6,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(2,2,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(6,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)

subplot(4,1,3)
histogram(r(3,2,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(2,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(3,2,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(2,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)

subplot(4,1,4)
histogram(r(4,1,noise_i,1,:), 'BinWidth', 0.002, 'FaceColor', c(4,:), ...
    'EdgeColor', 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3, 'LineWidth', 2)
hold on
histogram(var(4,1,noise_i,1,:), 'BinWidth', 0.001, 'FaceColor', c(4,:), ...
    'FaceAlpha', 0.8, 'LineWidth', 2)
xlim([0, 0.08])
ylim([0 mc/4])
xlabel('Meters, m')
ylabel('Number of simulations')
set(gca,'FontSize', 18)
legend({'Estimation bias', 'Point spreading'}, 'FontSize', 24)
