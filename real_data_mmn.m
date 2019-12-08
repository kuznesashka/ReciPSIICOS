load G3_mmn
load mmn

% 1. FORWARD MODEL
R = G3.GridLoc; 
ChUsed = 31:304;

% calculate tangential plane dipoles
[Nch, Nsites] = size(G3.Gain(ChUsed,1:3:end));
G2d = zeros(Nch,Nsites*2);
G2d0 = zeros(Nch,Nsites*2);
range = 1:2;
for i = 1:Nsites
    g = [G3.Gain(ChUsed,1+3*(i-1)) G3.Gain(ChUsed,2+3*(i-1)) G3.Gain(ChUsed,3+3*(i-1))];
    [u sv v] = svd(g);
    gt = g*v(:,1:2);
    G2d(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
    G2d0(:,range) = gt;
    range = range + 2;
end

% 2. REDUCING SENSOR SPACE
GainSVDTh = 0.0012; 
[ug sg vg] = spm_svd(G2d*G2d', GainSVDTh);
UP = ug'; % direction of dimention reduction
G2dU = UP*G2d;
G2d0U = UP*G2d0;

% 3. PSIICOS projector
[Upwr, ds, Apwr] = ProjectorOnlyAwayFromPowerComplete(G2d0U,1000);
Rnk = 500;

% 4. whitened projector
% Wks =  load('C_MMN.mat','C_re');
% Ccorr = Wks.C_re;
% Cpwr = Apwr*Apwr';
% RankG = size(G2dU,1);
% Wpwr = sqrtm(inv(Cpwr+0.001*trace(Cpwr)/(RankG^2)*eye(size(Cpwr))));
% WCcorrW = Wpwr*double(Ccorr)*Wpwr';
% 
% ProjRnkMax = 1000;
% [u s] = eigs(WCcorrW,ProjRnkMax);
% UcorrW_rnk = u(:,1:ProjRnkMax);
% 
% PrFromCorr_W = inv(Wpwr+0.05*trace(Wpwr)/size(Wpwr,1))*(eye(size(UcorrW_rnk,1))-UcorrW_rnk(:,1:1000)*UcorrW_rnk(:,1:1000)')*Wpwr;
% 
% 4. DATA loading
% protocolPath = '/media/ksasha/F21EF9851EF94361/For_new_comp/BST_protocols/AntiPSIICOS/data/Subject01/S01_AEF_20131218_01_600Hz_band_022';
% cd(protocolPath)
% std_cond = ['*standard_fix_trial*.mat'];
% Std = dir(std_cond);
%     
% for i = 1:length(Std)
%     Data = load([Std(i).folder, '/', Std(i).name]);
%     Data_std(:,:,i) = Data.F(ChUsed,:); 
% end
% aStd = mean(Data_std, 3);
% 
% dev_cond = ['*deviant_fix_trial*.mat'];
% Dev = dir(dev_cond);
%     
% for i = 1:length(Dev)
%     Data = load([Dev(i).folder, '/', Dev(i).name]);
%     Data_dev(:,:,i) = Data.F(ChUsed,:); 
% end
% aDev = mean(Data_dev, 3);

% c = lines(8);
% figure
% plot(-200:(1/2400*1000):800,aStd'*10^13, 'k')
% hold on
% plot([150,150], [-2.5, 2.5],'Color',c(2,:), 'LineWidth', 2)
% plot([0,0], [-2.5, 2.5], '--k')
% xlim([-200,800])
% ylim([-2.5 2.5])
% xlabel('Time [ms]')

% 5. COMPUTE COVARIANCE in the virtual sensors space
% mmn = aDev-aStd;
mmn = mmn(ChUsed,:);
Ca = UP*mmn*mmn'*UP';
Nsrc = size(G2dU,2);

cd('/home/ksasha/Projects/AntiPSIICOS_beamformer/For_paper/Real_data_MMN')
% 6. BEAMFORMER
iCa = tihinv(Ca, 0.01); % regularized correlation matrix
range2d = 1:2;
clear w;
for i = 1:Nsites
    g = G2dU(:,range2d);
    m = inv(g'*iCa*g);
    w(:,:,i) = m*g'*iCa;
    [u ss v] = svd(m);
    Z(i) = ss(1,1);
    range2d = range2d+2;
end

% 7. PSIICOS MODIFIED BEAMFORMER
% project the covariance matrix away from zero-phase synchrony subspace
Cap = reshape(Upwr(:,1:Rnk)*Upwr(:,1:Rnk)'*Ca(:),size(Ca));
% fix negative eigenvalues (after projection) issue
[e a] = eig(Cap);

Cap = e*abs(a)*e';
iCap = tihinv(Cap, 0.01);
clear wp;
range2d = 1:2;
for i = 1:Nsites
    g = G2dU(:,range2d);
    m = inv(g'*iCap*g);
    [u ss v] = svd(m);
    Zp(i) = ss(1,1);
    wp(:,:,i) = m*g'*iCap;
    range2d = range2d+2;
end


% 8. Whitened PSIICOS
% project the covariance matrix away from zero-phase synchrony subspace
% Cap = reshape(Upwr(:,1:Rnk)*Upwr(:,1:Rnk)'*PrFromCorr_W*Ca(:),size(Ca));
% fix negative eigenvalues (after projection) issue
% [e a] = eig(Cap);
% 
% Cap = e*abs(a)*e';
% iCap = tihinv(Cap, 0.01);
% clear wps;
% range2d = 1:2;
% for i = 1:Nsites
%     g = G2dU(:,range2d);
%     m = inv(g'*iCap*g);
%     [u ss v] = svd(m);
%     Zps(i) = ss(1,1);
%     wps(:,:,i) = m*g'*iCap;
%     range2d = range2d+2;
% end

% 8. MNE
lambda = 0.1;
Cs = eye(Nsites*2);
Cn = eye(size(G2dU,1));
GGt = G2dU*Cs*G2dU';
Wmne = G2dU'*inv(G2dU*Cs*G2dU'+lambda*trace(GGt)/size(GGt,1)*Cn);
clear Zmne;
range2d = 1:2;
for i = 1:Nsites
    wm(:,:,i) = Wmne(range2d,:);
    m = wm(:,:,i)*Ca*wm(:,:,i)';
    [u ss v] = svd(m);
    Zmne(i) = ss(1,1);
    range2d = range2d+2;
end

% 9. APPLY TO THE DATA
mmn_low = UP*mmn; % new sensor number

% sum of squares

data_bf = ((squeeze(w(1,:,:))'*mmn_low).^2 + ...
    (squeeze(w(2,:,:))'*mmn_low).^2)';
data_ps = ((squeeze(wp(1,:,:))'*mmn_low).^2 + ...
    (squeeze(wp(2,:,:))'*mmn_low).^2)';
data_mne = ((squeeze(wm(1,:,:))'*mmn_low).^2 + ...
    (squeeze(wm(2,:,:))'*mmn_low).^2)'; 
%     data_wps(i,:) = (squeeze(wps(1,:,:))'*mmn_low(:,i)).^2 + ...
%         (squeeze(wps(2,:,:))'*mmn_low(:,i)).^2; 


% InvSol.ImageGridAmp = [];
% InvSol.ImageGridAmp = data_ps';
% 
% InvSol.ImageGridAmp = zeros(Nsites,1);
% InvSol.ImageGridAmp(src_max) = 100;

% figure
% hold on
% stem(Zp)
% hold on
% stem(Z)
% 
% stem(src_max, Z(src_max), 'r')
% stem(src_max, Zp(src_max), 'r')

%3963
src_max = 3961;
time_max = 159;

% ind = 10000:12000;
% ind2 = intersect(find(data_ps(259,10000:12000)>5*10^(-26)),find(data_bf(259,10000:12000)>1.5*10^(-26)));
% [val, ind_m] = max(data_ps(259,ind(ind2)))

% figure
% subplot(1,2,1)
% hold on
% stem(data_ps(259,:))
% stem(ind(ind2), data_ps(259, ind(ind2)), 'r')
% subplot(1,2,2)
% hold on
% stem(data_bf(259,:))
% stem(ind(ind2), data_bf(259, ind(ind2)), 'r')

% src_max = ind(ind2(ind_m));

% figure
% scatter3(R(:,1), R(:,2), R(:,3))
% hold on
% scatter3(R(src_max,1),R(src_max,2),R(src_max,3), 'r', 'filled')

for i = 1:Nsites
    norm_ps(i) = (norm(wp(1,:,i))+norm(wp(2,:,i)))/2;
    norm_bf(i) = (norm(w(1,:,i))+norm(w(2,:,i)))/2;
    norm_mne(i) = (norm(wm(1,:,i))+norm(wm(2,:,i)))/2;
end

% for i = 1:Nsites
%     dd(i) = norm(R(i,:)-R(src_max,:));
% end
% 
% InvSol.ImageGridAmp(dd<0.005) = 100;

% title('ERP on standard auditory stimulus, 1-45Hz')
% subplot(4,1,3)
% stem(data_ps(252,:))
% hold on
% stem(data_bf(252,:))
% plot([src_max, src_max], [0, 2*10^(-25)], 'k', 'LineWidth', 2)
% legend({'AntiPSIICOS variance', 'LCMV variance', 'Picked source'})
% xlabel('Source index')
% title('Picked source for analysis')
% xlim([1,Nsites])


c = lines(7);

figure
subplot(2,1,1)
plot(-100:500,mmn*10^15, 'k')
hold on
plot([0,0], [-300, 300], '--k')
plot([time_max, time_max], [-300, 300], 'k', 'LineWidth', 2)
xlim([-100,500])
ylim([-300 300])
ylabel('Signal amplitude (fT)')
xlabel('Time [ms]')
% ylabel('\muA')
set(gca,'FontSize',20)
subplot(2,1,2)
plot(-100:500, data_ps(:,src_max), 'LineWidth', 3, 'Color', c(1,:))
hold on
plot(-100:500, data_bf(:,src_max), 'LineWidth', 3, 'Color', c(2,:))
% area(-100:500, data_ps(257,src_max)/data_mne(257,src_max)*data_mne(:,src_max), ...
%     'FaceAlpha', 0.1, 'FaceColor', c(4,:), 'EdgeColor', c(4,:),'LineStyle','--')
area(-100:500, (data_ps(257,src_max)/data_bf(257,src_max)).*data_bf(:,src_max),...
    'FaceAlpha', 0.2, 'FaceColor', c(2,:), 'EdgeColor', c(2,:),'LineStyle','--')
plot([0,0], [0, 7*10^(-25)], '--k')
% plot([159,159], [0, 7*10^(-25)], 'k', 'LineWidth', 2)
ylim([0, 7*10^(-25)])
legend({'ReciPSIICOS', 'LCMV', 'Scaled LCMV',  'Stimulus onset'}, 'Location',...
    'west', 'FontSize', 28)
legend('boxoff')
xlabel('Time [ms]')
ylabel('Estimated source amplitude (a.u.)')
set(gca,'FontSize',20)

norm_ratio = (norm_bf./norm_ps)';
figure
boxplot((norm_bf./norm_ps)', 'Colors', c(1,:), 'Symbol', 'ko')
set(gca,'FontSize',20)


clear m
for i = 1:15002
    m(i) = norm(R(i,:)-R(src_max,:));
end
ind_close = find(m<0.002);

InvSol.ImageGridAmp = zeros(1,15002)';
InvSol.ImageGridAmp(ind_close) = 100;




ind = 10000:12000;
ind2 = find(data_ps(259,10000:12000)>5*10^(-26));
[val, max_src] = max(data_ps(259,10000:12000));

src_max = ind(max_src);

c = lines(7);
src_max = ind(ind2(i));
figure
subplot(2,1,1)
plot(-100:500,mmn*10^13, 'k')
hold on
plot([0,0], [-2.5, 2.5], '--k')
plot([159,159], [-2.5, 2.5], 'k', 'LineWidth', 2)
xlim([-100,500])
ylim([-2.5 2.5])
xlabel('Time [ms]')
set(gca,'FontSize',20)
subplot(2,1,2)
plot(-100:500, data_ps(:,src_max), 'LineWidth', 3, 'Color', c(1,:))
% plot([0,0], [0, 1.5*10^(-25)], '--k')
plot([159,159], [0, 1.5*10^(-25)], 'k', 'LineWidth', 2)
ylim([0, 1.5*10^(-25)])
% legend({'ReciPSIICOS', 'LCMV', 'Scaled MNE', 'Scaled LCMV',  'Stimulus onset', 'MMN peak latency'})
% legend('boxoff')
xlabel('Time [ms]')
ylabel('a.u.')
set(gca,'FontSize',20)

       