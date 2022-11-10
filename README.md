# ReciPSIICOS

**Data source**: https://disk.yandex.ru/d/aTsr6XO4ZdMhLg

**Simulations**:

Run `main.m` script to perform the following steps:

1. Figure 7.A.
	* 1.1 `save_simulations_snr.m` — Generate simulations for different methods (MNE, LCMV beamformer, ReciPSIICOS, WReciPSIICOS) and different SNRs ([0, 0.5, 0.7, 1.5, 1.7, 2, 4, 5, 7]). The resulting activation maps saved at `Z_total_full.mat` and source locations in `picked_src.mat` (available in data folder).
	* 1.2 `plot_metrics_simulations_snr.m` — calculate point spreading and bias for precomputed `Z_total_full.mat` and `picked_src.mat` and plot them. Save collected bias and point spreading values as `r_snr.mat` and `var_snr.mat`.

2. Figure 7.B.
	* 2.1 `save_simulations_Gerror.m` — Generate simualtions for different forward model inaccuracies, using matrices `C_re_005`, `C_re_01`, `C_re_02` (data folder), save the results into `Z_total_G.mat` (data folder).
	* 2.2 `plot_metrics_simulations_Gerror.m` — calculate point spreading and bias for precomputed `Z_total_G.mat` and `picked_src_G.mat` and plot them.

3. Figure 8.
	* 3.1 `plot_simulation_histogram.m` — plot histograms for precomputed `r_snr.mat` and `var_snr.mat`.

4. Figure 9.
	* 4.1 `save_simulations_rank.m` — Generate simulations with fixed SNR, then apply source reconstruction with PSIICOS and ReciPSIICOS techniques with different projection ranks. The resulting activation maps saved at `Z_total_rank.mat` and source locations at `picked_src_rank.mat` (uploaded in data folder).
	* 4.2 `plot_simulations_rank.m` — calculate point spreading and bias for precomputed `Z_total_rank.mat` and `picked_src_rank.mat` and plot them.

5. Figures 10, 11, 12.
	* 5.1 `save_simulations_three_src_asynch.m` — generate simulations with three active moderatly related sources. Calculate bias, pointspreading variance, detection rate.
	* 5.2 `save_simulations_three_src_synch.m` — generate simulations with three synchronous sources.
	* 5.3 `plot_simulations_threesrc.m` — visualize simulations with three sources (input data uploaded in data folder) — for Figures 11 and 12.

**Real data**:

1. MMN analysis:
	* real_data_mmn.m — analysis of auditory evoked MMN, using `mmn.mat` and `G3_mmn.mat` (data folder).
