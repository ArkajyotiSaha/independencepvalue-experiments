# independencepvalue-experiments
Code to reproduce simulation results and gene expression data analysis results from the article.

# Package
independencepavlue_0.0.2.tar.gz contains the R package `independencepavlue` used to implement the proposed selective inference approach in the article. 

# Simulation
Simulation codes and results are stored in the folders Simulation_codes and Simulation_results, respectively. 

# Gene expression data analysis
We use the data from [DREAM5 network inference challenge](https://www.synapse.org/#!Synapse:syn2787209/wiki/70354). For user convenience we provide the relevant data in the subfolder DREAM5_data, under Real_data_codes. The subsubfolders DREAM5_data/Training and DREAM5_data/gold_standard_edges_only contain data pertaining to our article, and are downloaded from [training data](https://www.synapse.org/#!Synapse:syn2787212) and [Evaluation scripts](https://www.synapse.org/#!Synapse:syn2787219), respectively.
Running Real_data_analysis.R in Real_data_codes produces the gene expression data analysis results, that are stored in folder Real_data_results. 

# Figures
Figure codes and figures are stored in the folders Figures_codes and Figures respectively.

1. Figure 1(a)-(b): Figure_1(a)-(b).R produces this plot.
2. Figure 1(c): Running Simulation_1(c).R produces the simulation results. Figure_1(c).R produces Figure 1(c).
3. Figure 2: Running Simulation_2.R produces the simulation results. Figure_2.R produces Figure 2.
4. Figures 3-4: Running Simulation_3&4.R produces the simulation results. Figure_3.R and Figure_4.R produces Figures 3 and 4.
5. Figure 5: Figure_5.R produces Figure 5 using the data in folder Real_data_results. 
6. Figure 6: Figure_6.R produces Figure 6 using the data in folder Real_data_results.
7. Figure S1: Running Simulation_S1_variance_filtering.R and Simulation_S1_mean_filtering.R produces the simulation results and Figure_S1.R produces Figure S1. 
8. Figure S2: Running Simulation_S2.R produces the simulation results. Figure_S2.R produces Figure S2.
9. Figure S3: Running Simulation_S3.R produces the simulation results. Figure_S3.R produces Figure S3.
10. Figure S4: Running Simulation_S4.R produces the simulation results. Figure_S4.R produces Figure S4.
11. Figure S5: Running Simulation_S5.R produces the simulation results. Figure_S5.R produces Figure S5.
12. Figure S6: Running Simulation_S6.R produces the simulation results. Figure_S6.R produces Figure S6.
