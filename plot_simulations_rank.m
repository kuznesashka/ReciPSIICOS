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
load Z_total_rank
load picked_src_rank

src_left_red = find(R_red(:,2)>0);
src_right_red = find(R_red(:,2)<0);

snr = [100, 300, 500, 700, 1000];
Nmc = 500;
Range_frac = [0.65; 0.25];
thresh = [0.02, 0.02];

bad_detection = zeros(2, 2, length(snr), size(Range_frac, 2));
clear r var
for noise_i = 1:length(snr) 
    for mc = 1:Nmc
        ind_generated = picked_src(mc,:);
        for synch = 1:2
            for method = 1:2
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
                        bad_detection(method, synch, noise_i, frac) = bad_detection(method, synch, noise_i, frac) + 1; 
                    else

                        r(method, synch, noise_i, frac, mc) = (norm(R_red(max_ind_l,:)-R(ind_generated(1),:))+ ...
                            norm(R_red(max_ind_r,:)-R(ind_generated(2),:)))/2;
                        if r(method, synch, noise_i, frac, mc) > thresh(synch)
                           bad_detection(method, synch, noise_i, frac) = bad_detection(method, synch, noise_i, frac) + 1; 
                        end

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

mean_r = median(r, 5, 'omitnan');
mean_var = median(var, 5, 'omitnan');

c = lines(7);
col = c([1,6,2,4],:);
figure
subplot(3,2,1)
hold on
for i = 1:2
    plot(snr, squeeze(mean_r(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
    ylim([0 0.03])
    xlim([0 1100])
end
scatter(snr, squeeze(mean_r(1,1,:)), ...
    150, 'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'k')
scatter(snr, squeeze(mean_r(2,1,:)), ...
    150, 's', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'k')
%title('Mean bias')
set(gca,'FontSize', 24)
xlabel('Projection rank')
ylabel('Meters, m')
plot([500 500], [0 0.03], '--', 'Color', c(1,:), 'LineWidth', 2)
plot([1000 1000], [0 0.03], '--', 'Color', c(6,:), 'LineWidth', 2)

subplot(3,2,3)
hold on
for i = 1:2
    plot(snr, squeeze(mean_var(i,1,:)), 'Color', col(i,:), 'LineWidth', 3)
     ylim([0 0.02])
    xlim([0 1100])
end
scatter(snr, squeeze(mean_var(1,1,:)), ...
150, 'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'k')
scatter(snr, squeeze(mean_var(2,1,:)), ...
150, 's', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Projection rank')
ylabel('Meters, m')
plot([500 500], [0 0.03], '--', 'Color', c(1,:), 'LineWidth', 2)
plot([1000 1000], [0 0.03], '--', 'Color', c(6,:), 'LineWidth', 2)

subplot(3,2,5)
hold on
for i = 1:2
    plot(snr,(1-(squeeze(bad_detection(i,1,:,1))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
     ylim([30 100])
    xlim([0 1100])
end
scatter(snr, (1-(squeeze(bad_detection(1,1,:,1))./Nmc))*100, ...
150, 'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'k')
scatter(snr, (1-(squeeze(bad_detection(2,1,:,1))./Nmc))*100, ...
150, 's', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Projection rank')
ylabel('Percents, %')
plot([500 500], [0 100], '--', 'Color', c(1,:), 'LineWidth', 2)
plot([1000 1000], [0 100], '--', 'Color', c(6,:), 'LineWidth', 2)


subplot(3,2,2)
hold on
for i = 1:2
    plot(snr, squeeze(mean_r(i,2,:)), 'Color', col(i,:), 'LineWidth', 3)
    ylim([0 0.03])
    xlim([0 1100])
end
scatter(snr, squeeze(mean_r(1,2,:)), ...
150, 'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'k')
scatter(snr, squeeze(mean_r(2,2,:)), ...
150, 's', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Projection rank')
ylabel('Meters, m')
plot([500 500], [0 0.03], '--', 'Color', c(1,:), 'LineWidth', 2)
plot([1000 1000], [0 0.03], '--', 'Color', c(6,:), 'LineWidth', 2)


subplot(3,2,4)
hold on
for i = 1:2
    plot(snr, squeeze(mean_var(i,2,:)), 'Color', col(i,:), 'LineWidth', 3)
    ylim([0 0.02])
    xlim([0 1100])
end
scatter(snr, squeeze(mean_var(1,2,:)), ...
150, 'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'k')
scatter(snr, squeeze(mean_var(2,2,:)), ...
150, 's', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Projection rank')
ylabel('Meters, m')
plot([500 500], [0 0.03], '--', 'Color', c(1,:), 'LineWidth', 2)
plot([1000 1000], [0 0.03], '--', 'Color', c(6,:), 'LineWidth', 2)


subplot(3,2,6)
hold on
for i = 1:2
    plot(snr,(1-(squeeze(bad_detection(i,2,:,1))./Nmc))*100, 'Color', col(i,:), 'LineWidth', 3)
    ylim([30 100])
    xlim([0 1100])
end
scatter(snr, (1-(squeeze(bad_detection(1,2,:,1))./Nmc))*100, ...
150, 'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'k')
scatter(snr, (1-(squeeze(bad_detection(2,2,:,1))./Nmc))*100, ...
150, 's', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'k')
set(gca,'FontSize', 24)
xlabel('Projection rank')
ylabel('Percents, %')
plot([500 500], [0 100], '--', 'Color', c(1,:), 'LineWidth', 2)
plot([1000 1000], [0 100], '--', 'Color', c(6,:), 'LineWidth', 2)





       
