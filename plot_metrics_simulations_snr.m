% 1. FORWARD MODEL

G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channels

R = G3.GridLoc; % source location, dense matrix
R_red = G3_red.GridLoc; % source location reduced matrix
ChUsed = find(strcmp({chans.Channel.Type}, 'MEG GRAD')); % set to use gradiometers only

[Nch, Nsites] = size(G3.Gain(ChUsed,1:3:end));
[Nch, Nsites_red] = size(G3_red.Gain(ChUsed,1:3:end));

% 2. SIMULATIONS 

% Z_total = (method, synch, noise_i, mc, src)
% methods: ReciPSIICOS, WReciPSIICOS, LCMV, MNE
load Z_total_full
load picked_src

src_left_red = find(R_red(:,2)>0);
src_right_red = find(R_red(:,2)<0);

snr = [0, 0.5, 0.7, 1.5, 1.7, 2, 4, 5, 7];
Nmc = 500;
Range_frac = [0.65; 0.25];

clear r var
for noise_i = 1:length(snr) 
    for mc = 1:Nmc
        ind_generated = picked_src(mc,:);
        for synch = 1:2
            for method = 1:4
                for frac = 1:size(Range_frac, 2)
                    Z = squeeze(Z_total(method,synch,noise_i,mc,:))';
                    [max_val max_ind] = max(Z);
                    Z(Z < Range_frac(synch, frac)*max_val) = 0;
          
                    clear dist_l dist_r
                    Z_left = Z(src_left_red); % values from left hemisphere
                    Z_right = Z(src_right_red); % values from right hemisphere
                    [max_val_l max_l] = max(Z_left); % maximum in the left hemisphere
                    max_ind_l = src_left_red(max_l);
                    [max_val_r max_r] = max(Z_right);
                    max_ind_r = src_right_red(max_r);
        

                    if (length(Z_left(Z_left>0)) == 0)|(length(Z_right(Z_right>0)) == 0)
                        r(method, synch, noise_i, frac, mc) = NaN;
                        var(method, synch, noise_i, frac, mc) = NaN;
%                         coord_active_l = R_red(src_left_red(Z_left>0),:);
%                         coord_active_r = R_red(src_right_red(Z_right>0),:);
%                         figure
%                         scatter3(R_red(:,1), R_red(:,2), R_red(:,3))
%                         hold on
%                         scatter3(coord_active_r(:,1), coord_active_r(:,2), ...
%                             coord_active_r(:,3), 'filled', 'r')
%                         scatter3(coord_active_l(:,1), coord_active_l(:,2), ...
%                             coord_active_l(:,3), 'filled', 'g')
%                         scatter3(R(ind_generated,1), R(ind_generated,2), ...
%                             R(ind_generated,3), 'filled', 'k')
                    else

                        r(method, synch, noise_i, frac, mc) = (norm(R_red(max_ind_l,:)-R(ind_generated(1),:))+ ...
                            norm(R_red(max_ind_r,:)-R(ind_generated(2),:)))/2;

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
                        
%                         figure
%                         scatter3(R_red(:,1), R_red(:,2), R_red(:,3))
%                         hold on
%                         scatter3(coord_active_r(:,1), coord_active_r(:,2), ...
%                             coord_active_r(:,3), 'filled', 'r')
%                         scatter3(coord_active_l(:,1), coord_active_l(:,2), ...
%                             coord_active_l(:,3), 'filled', 'g')
%                         scatter3(R(ind_generated,1), R(ind_generated,2), ...
%                             R(ind_generated,3), 'filled', 'k')
% 

                        var(method, synch, noise_i, frac, mc) = (sum(Z_left_norm(Z_left_norm>0).*dist_l)+ ...
                            sum(Z_right_norm(Z_right_norm>0).*dist_r))/2;
                    end
                end
            end        
        end
        mc
    end
    
%     minimum_r = min(min(r(:,1,noise_i,:)));
%     maximum_r = max(max(r(:,1,noise_i,:)));
%     minimum_v = min(min(var(:,1,noise_i,:)));
%     maximum_v = max(max(var(:,1,noise_i,:)));
% 
% 
%     figure
%     subplot(4,2,1)
%     histogram(squeeze(r(4,1,noise_i,1,:)), 20)
%     xlim([minimum_r, maximum_r])
%     ylim([0 125])
%     title('Bias, MNE')
%     subplot(4,2,3)
%     histogram(squeeze(r(3,1,noise_i,1,:)),20)
%     xlim([minimum_r, maximum_r])
%     ylim([0 125])
%     title('Bias, LCMV')
%     subplot(4,2,5)
%     histogram(squeeze(r(1,1,noise_i,1,:)),20)
%     xlim([minimum_r, maximum_r])
%     ylim([0 125])
%     title('Bias, APSIICOS')
%     subplot(4,2,7)
%     histogram(squeeze(r(2,1,noise_i,1,:)),20)
%     xlim([minimum_r, maximum_r])
%     ylim([0 125])
%     title('Bias, WAPSIICOS')
%     subplot(4,2,2)
%     histogram(squeeze(var(4,1,noise_i,1,:)),20)
%     xlim([minimum_v, maximum_v])
%     ylim([0 125])
%     title('PointSpread, MNE')
%     subplot(4,2,4)
%     histogram(squeeze(var(3,1,noise_i,1,:)),20)
%     xlim([minimum_v, maximum_v])
%     ylim([0 125])
%     title('PointSpread, LCMV')
%     subplot(4,2,6)
%     histogram(squeeze(var(1,1,noise_i,1,:)),20)
%     xlim([minimum_v, maximum_v])
%     ylim([0 125])
%     title('PointSpread, APSIICOS')
%     subplot(4,2,8)
%     histogram(squeeze(r(2,1,noise_i,1,:)),20)
%     xlim([minimum_v, maximum_v])
%     ylim([0 125])
%     title('PointSpread, WAPSIICOS')

end

detect = sum(isnan(r), 5);
mean_r = mean(r, 5, 'omitnan');
mean_var = mean(var, 5, 'omitnan');

c = lines(7);
col = c([1, 6, 2, 4], :);
figure
subplot(3, 2, 1)
hold on
for i = 1:4
    plot(snr, squeeze(mean_r(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(snr, squeeze(mean_r(i,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim([0 0.07])
    xlim([1 10])
end
plot([4, 4], [0 0.07], '--', 'Color', 'k')
set(gca,'FontSize', 24)
xlabel('Signal-to-noise ratio')
ylabel('Meters, m')
subplot(3,2,3)
hold on
for i = 1:4
    plot(snr, squeeze(mean_var(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(snr, squeeze(mean_var(i,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim([0 0.03])
    xlim([1 10])
end
set(gca,'FontSize', 24)
plot([4, 4], [0 0.03], '--', 'Color', 'k')
xlabel('Signal-to-noise ratio')
ylabel('Meters, m')
subplot(3, 2, 5)
hold on
for i = 1:4
    plot(snr,(1-(squeeze(detect(i,1,:))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
    scatter(snr, (1-(squeeze(detect(i,1,:))./Nmc))*100, ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim([60 100])
    xlim([1 10])
end
plot([4, 4], [60 100], '--', 'Color', 'k')
set(gca,'FontSize', 24)
xlabel('Signal-to-noise ratio')
ylabel('Percents, %')


subplot(3, 2, 2)
hold on
for i = 1:3
    plot(snr, squeeze(mean_r(i,2,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(snr, squeeze(mean_r(i,2,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim([0 0.05])
    xlim([1 10])
end
plot([4, 4], [0 0.05], '--', 'Color', 'k')
plot(snr, squeeze(mean_r(4,1,:)), 'Color', col(4,:), 'LineWidth', 3)
scatter(snr, squeeze(mean_r(4,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(4,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Signal-to-noise ratio')
ylabel('Meters, m')

subplot(3, 2, 4)
hold on
for i = 1:3
    plot(snr, squeeze(mean_var(i,2,:)), 'Color', col(i,:), 'LineWidth', 3)
    scatter(snr, squeeze(mean_var(i,2,:)), ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim([0 0.03])
    xlim([1 10])
end
set(gca,'FontSize', 24)
plot([4, 4], [0 0.03], '--', 'Color', 'k')
plot(snr, squeeze(mean_var(4, 1, :)), 'Color', col(4, :), 'LineWidth', 3)
scatter(snr, squeeze(mean_var(4, 1, :)), ...
    150, 'o', 'MarkerFaceColor', col(4, :), 'MarkerEdgeColor', 'k')
xlabel('Signal-to-noise ratio')
ylabel('Meters, m')

subplot(3, 2, 6)
hold on
for i = 1:3
    plot(snr,(1-(squeeze(detect(i,2,:))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
    scatter(snr, (1-(squeeze(detect(i,2,:))./Nmc))*100, ...
    150, 'o', 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k')
    ylim([60 100])
    xlim([1 10])
end
plot([4, 4], [60 100], '--', 'Color', 'k')
plot(snr,(1-(squeeze(detect(4,1,:))./Nmc))*100, 'Color', col(4,:), 'LineWidth', 3)
    scatter(snr, (1-(squeeze(detect(4,1,:))./Nmc))*100, ...
    150, 'o', 'MarkerFaceColor', col(4,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Signal-to-noise ratio')
ylabel('Percents, %')


       
