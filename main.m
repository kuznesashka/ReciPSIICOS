
% 1. UPLOAD FORWARD MODEL for dense and sparse matrices
G3 = load('R0015_G.mat'); % forward model matrix 20002 sources
G3_red = load('RR015_G_5k.mat'); % reduced forward model matrix 5001 sources
chans = load('R0015_chans.mat'); % channel locations


% Source reconstruction with PSIICOS and ReciPSIICOS for different
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

% Calculate point spreading and bias for precomputed 
% `Z_total_rank.mat` and `picked_src_rank.mat` and plot them 
load Z_total_rank
load picked_src_rank
plot_simulations_rank(G3, G3_red, chans, Z_total, picked_src, rank, Nmc);