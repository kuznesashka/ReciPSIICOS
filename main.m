
% 1. Upload forward model for dense and sparse matrices
G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channel locations

% 1.1 Generate simulations for different methods (MNE, LCMV beamformer, 
% ReciPSIICOS, WReciPSIICOS) and different SNRs

snr = [0, 0.5, 0.7, 1.5, 1.7, 2, 4, 5, 7];
Nmc = 500;
load_corr = true;
load_src = true;

[Z_total, picked_src] = save_simulations_snr(G3, G3_red, chans, ...
snr, Nmc, load_corr, load_src);

save Z_total_snr Z_total
if (~load_src)
    save picked_src_snr picked_src
end
 
% 1.2 Calculate point spreading and bias for precomputed `Z_total_full.mat` 
% and `picked_src.mat` and plot them
clear Z_total picked_src
load Z_total_snr
load picked_src_snr
snr = [0, 0.5, 0.7, 1.5, 1.7, 2, 4, 5, 7];
Nmc = 500;

[r, var] = plot_metrics_simulations_snr(G3, G3_red, ...
    Z_total, picked_src, snr, Nmc);
save r_snr r
save var_snr var


% 2.1 Generate simulations for different forward model inaccuracies
clear Z_total
load Z_total_snr

load_corr = true;
load_src = true;
Nmc = 500;
[Z_total_G, picked_src] = save_simulations_Gerror(G3, G3_red, chans, ...
    Nmc, load_corr, load_src, Z_total);

save Z_total_G Z_total_G
if (~load_src)
    save picked_src_G picked_src
end

% 2.2 Calculate bias, pointspreading, detection ratio and plot

clear Z_total picked_src
load Z_total_G
load picked_src_G

plot_metrics_simulations_Gerror(G3, G3_red, Z_total_G, picked_src);


% 3. Plot histograms for bias and pointspreading distribution for fixed SNR

clear r var 
load r_snr
load var_snr
plot_simulation_histogram(r, var)


% 4.1 Source reconstruction with PSIICOS and ReciPSIICOS for different
% projection ranks

rank = [100, 300, 500, 700, 1000];
Nmc = 500;
load_corr = true;
load_src = true;
[Z_total, picked_src] = save_simulations_rank(G3, G3_red, chans, ...
    rank, Nmc, load_corr, load_src);

save Z_total_rank Z_total 
if (~load_src)
    save picked_src_rank picked_src
end

% 4.2 Calculate point spreading and bias for precomputed 
% `Z_total_rank.mat` and `picked_src_rank.mat` and plot them

load Z_total_rank
load picked_src_rank
plot_simulations_rank(G3, G3_red, Z_total, picked_src, rank, Nmc);


% 5.1 Generate simulations for three asynchronous sources, calculate bias,
% point spreading and detection ratio

Nmc = 500;
[r, var, detect] = ...
    save_simulations_three_src_asynch(G3, G3_red, chans, Nmc);
save r_three_asynch r
save var_three_asynch var
save detect_three_asynch detect

% 5.2 Generate simulations for three synchronous sources, calculate bias,
% point spreading and detection ratio

Nmc = 500;
[r, var, detect] = ...
    save_simulations_three_src_synch(G3, G3_red, chans, Nmc);
save r_three_synch r
save var_three_synch var
save detect_three_synch detect

% 5.3 Plot simulations
% synchronous sources
load r_three_synch
load var_three_synch
load detect_three_synch

plot_simulations_threesrc(r, var, detect)

% asynchronous sources
load r_three_asynch
load var_three_asynch
load detect_three_asynch

plot_simulations_threesrc(r, var, detect)
