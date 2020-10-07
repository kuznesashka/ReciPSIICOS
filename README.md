# ReciPSIICOS

Data source: https://drive.google.com/drive/folders/1sOtgn66MiLxC1eer25z9Y1Vrzcc8Ga6F?usp=sharing

Simulations:

1. Figure 7.A.
	* 1.1 save_simulations_snr.m — Generate simulations for different methods (MNE, LCMV beamformer, ReciPSIICOS, WReciPSIICOS) and different SNRs ([0, 0.5, 0.7, 1.5, 1.7, 2, 4, 5, 7]). The resulting activation maps saved at `Z_total_full.mat` and source locations in `picked_src.mat` (available in GDrive folder).
	* 1.2 plot_metrics_simulations_snr.mat — calculate point spreading and bias for precomputed `Z_total_full.mat` and `picked_src.mat` and draw them on a graph.

2. Figure 7.B.
	* 2.1 save_simulations_Gerror.m — Generate simualtions for different forward model inaccuracies, using matrices `C_re_005`, `C_re_01`, `C_re_02` (GDrive), save the results into `Z_total_G.mat` (GDrive).
	* 2.2 plot_metrics_simulations_Gerror.mat — calculate point spreading and bias for precomputed `Z_total_G.mat` and `picked_src.mat` and draw them on a graph.

3. Figure 8.
	* 3.1 plot_simulation_histogram.m — plot histograms for precomputed `Z_total_full.mat`.

4. Figure 9.
	* 3.1 `save_simulations_rank.m` — Generate simulations with fixed SNR, then apply source reconstruction with PSIICOS and ReciPSIICOS techniques with different projection ranks. The resulting activation maps saved at `Z_total_rank.mat` and source locations at `picked_src_rank.mat` (uploaded in GoogleDrive folder).
	* 3.2 `plot_simulations_rank.mat` — calculate point spreading and bias for precomputed `Z_total_rank.mat` and `picked_src_rank.mat` and plot them.

5. Figures 10, 11, 12.
	* 4.1 save_simulations_three_src_asynch.m — generate simulations with three active moderatly related sources.
	* 4.2 save_simulations_three_src_synch.m — generate simulations with three synchronous sources.
	* 4.2 plot_simulations_threesrc.m — visualize simulations with three sources (input data on GDrive) — for Figures 11 and 12.

Real data:

1. MMN analysis:
	* real_data_mmn.m — analysis of auditory evoked MMN, using `mmn.mat` and `G3_mmn.mat` (GDrive).
