# independencepvalue-experiments

Code to reproduce simulation results and gene expression data analysis results from the paper "Inferring independent sets of Gaussian variables after thresholding correlations" by Arkajyoti Saha, Daniela Witten, and Jacob Bien. 

# Figures for simulation results

1. Figure 1(a)-(b): Figure_1(a)-(b).R produces this plot.
2. Figure 1(c): Running Simulation_1(c).R with indica between 1-100 produces the simulation results, that are stored in folder Simulation_results/histogram. Using the data in Simulation_results/histogram, Figure_1(c).R produces Figure 1(c).
3. Figure 2: Running Simulation_2.R with indica between 1-300 produces the simulation results, that are stored in folder Simulation_results/global_null. Using the data in Simulation_results/global_null, Figure_2.R produces Figure 2.
4. Figure 3: Running Simulation_3.R with indica between 1-300 produces the simulation results, that are stored in folder Simulation_results/alternative. Using the data in Simulation_results/alternative Figure_3.R produces Figure 3.

# Figure for gene expression data analysis results
We use the data from [DREAM5 network inference challenge](https://www.synapse.org/#!Synapse:syn2787209/wiki/70354). For user convenience we provide the relevant data in folder DREAM5_data. The subfolders DREAM5_data/Training and DREAM5_data/gold_standard_edges_only contains data partainig to our article, and are downloaded from [training data](https://www.synapse.org/#!Synapse:syn2787212) and [Evaluation scrips](https://www.synapse.org/#!Synapse:syn2787219), respectively.

Running Real_data_analysis.R with indica between 1-240 produces the gene expression data analysis results, that are stored in folder Real_data_results. 

1. Figure 4: Figure_4.R produces Figure 4 using the data in folder Real_data_results. 
2. Figure 5: Figure_5.R produces Figure 5 using the data in folder Real_data_results.

The generated figures are stored in folder Figures. 
