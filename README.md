# Logistic-tree normal models for microbiome compositions --- reproducible workflow

## Data

The data of the T1D cohort of the DIABIMMUNE project (Kostic et al., 2015) is used in this paper. 

The folder "./res" contains the raw data. The contents of this folder are as follows:

- mmc2.xlsx: information on the 33 subjects in the study (publicly available [here](https://www.cell.com/cms/10.1016/j.chom.2015.01.001/attachment/1f0883f8-1df7-447d-a47b-c1aa2bb2bbaf/mmc2.xlsx) )
- diabimmune_t1d_16s_metadata.rdata: the covariate information on the 777 samples (publicly available [here](https://diabimmune.broadinstitute.org/diabimmune/t1d-cohort) )
- diabimmune_t1d_16s_otu_table.txt: the OTU table and taxonomy table on all OTUs (publicly available [here](https://diabimmune.broadinstitute.org/diabimmune/t1d-cohort) )



## Code

The folder "./src" contains all scripts to reproduce the results in the paper. The contents of this folder are as follows:

- ./src/data/process_data_otu50.R: pre-processes raw data
- ./src/data/process_data_simulation.R: creates the data that is used throughout the simulation studies  
- ./src/simulation/covariance_estimation/dtm_data_clrcov.R: estimates the true clr covariance for the DTM examples

- ./src/simulation/covariance_estimation/dtm_sim.R: simulates data from DTM model
- ./src/simulation/covariance_estimation/ln_sim.R: simulates data from LN model
- ./src/simulation/covariance_estimation/ln_graph.R: functions for generating sparse networks (adapted from [this file](https://github.com/yuanpeicao/COAT/blob/master/simulation.R) )
- ./src/simulation/covariance_estimation/dtm_fit_parse.R: fits the LTN model to DTM examples
- ./src/simulation/covariance_estimation/ln_fit_parse.R: fits the LTN model to LN examples
- ./src/simulation/covariance_estimation/COAT.R: functions for fitting COAT (Cao et al., 2019) (adapted from [this file](https://github.com/yuanpeicao/COAT/blob/master/coat.R) )  
- ./src/simulation/covariance_estimation/fit_COAT.R: fits COAT to the simulated data
- ./src/simulation/covariance_estimation/collect_results.R: collects the results and produces the tables shown in section 3.2 in the paper
- ./src/simulation/cross_group_comparison/single_otu_sim.R: generates simulated datasets where cross-group differences exist at a single OTU
- ./src/simulation/cross_group_comparison/multi_otu_sim.R: generates simulated datasets where cross-group differences exist at multiple OTUs
- ./src/simulation/cross_group_comparison/fit.R: fits LTN to the simulated data
- ./src/simulation/cross_group_comparison/roc.R: calculates the ROC curve of LTN
- ./src/simulation/cross_group_comparison/dirfactor: scripts for fitting DirFactor (Ren et al., 2020) to the simulated data (adapted from [this Github repository](https://github.com/boyuren158/DirFactor-fix) )
- ./src/simulation/cross_group_comparison/plot.R: produces figures shown in section 3.2 of the paper
-  ./src/application/application.R: fits LTN-based mixed-effects model to the T1D cohort of DIABIMMUNE data
- ./src/applicaation/figures_tables.R: produces all figures and tables shown in the case study section of the paper



## Workflow

The `LTN` package is used throughout the analysis and can be installed using `devtools`: 
```R
devtools::install_github('MaStatLab/LTN')
```

The simulation study and case study is run on computing clusters. Arguments like random seeds are passed to R scripts from command lines for easier use. The workflow is documented below. 

### Set the working directory

```bash
export WORK_DIR=.../LTN_reproducible/
```

### Data processing

```bash
$WORK_DIR/src/data/process_data_otu50.R $WORK_DIR
$WORK_DIR/src/data/process_data_simulation.R $WORK_DIR
```

### Simulation studies

#### Covariance estimation

Calculate DTM parameters and the true clr covariance

```bash
$WORK_DIR/src/simulation/covariance_estimation/dtm_data_clrcov.R $WORK_DIR 1000000
```

Run the following script for a single simulation with random seed 1 for each simulation scenario. To produce the results in the paper, set i=1,...,100​. Change the value of lambda for sensitivity analysis.

```bash
i=1
$WORK_DIR/src/simulation/covariance_estimation/dtm_sim.R $WORK_DIR $i
for j in `seq 1 3`;do
    $WORK_DIR/src/simulation/covariance_estimation/ln_sim.R $WORK_DIR $i $j
done
# LTN
lambda=10
$WORK_DIR/src/simulation/covariance_estimation/dtm_fit_parse.R --SEED $i --lambda $lambda --WORK_DIR $WORK_DIR
for j in `seq 1 3`;do
    $WORK_DIR/src/simulation/covariance_estimation/ln_fit_parse.R --SEED $i --modelCov $j --lambda $lambda --WORK_DIR $WORK_DIR
done
# COAT
$WORK_DIR/src/simulation/covariance_estimation/fit_COAT.R $WORK_DIR $i
```

Summarize results and make tables

```bash
lambda=10
$WORK_DIR/src/simulation/covariance_estimation/collect_results.R $WORK_DIR $lambda 100
```

The results shown in Table. 1 in the paper are saved in the folder "./results/covariance_estimation". 

#### Cross-group comparison with LTN-based mixed-effects model

Run the following bash script for a single simulation with random seed 1 for the two scenarios respectively. To reproduce the results shown in the paper, set i=1,...,1000. Change the value of lambda for sensitivity analysis. 

```bash
lambda=10
i=1
declare -a rvec=(2 5)
for h in `seq 0 1`;do
$WORK_DIR/src/simulation/cross_group_comparison/single_otu_sim.R $WORK_DIR $i $h
$WORK_DIR/src/simulation/cross_group_comparison/multi_otu_sim.R $WORK_DIR $i $h 
# LTN
for reffcov in `seq 1 2`;do
$WORK_DIR/src/simulation/cross_group_comparison/fit.R --reffcov $reffcov --h $h --lambda $lambda --i $i --scenario single_otu --niter 10000 --WORK_DIR $WORK_DIR
$WORK_DIR/src/simulation/cross_group_comparison/fit.R --reffcov $reffcov --h $h --lambda $lambda --i $i --scenario multi_otu --niter 10000 --WORK_DIR $WORK_DIR
done
# DirFactor
for r in "${rvec[@]}";do
$WORK_DIR/src/simulation/cross_group_comparison/dirfactor/R/fit_single_otu.R $WORK_DIR $i $r $h 100000
$WORK_DIR/src/simulation/cross_group_comparison/dirfactor/R/fit_multi_otu.R $WORK_DIR $i $r $h 100000
done
done
```

Summarize results and make figures

```bash
lambda=10
$WORK_DIR/src/simulation/cross_group_comparison/roc.R $WORK_DIR $lambda
$WORK_DIR/src/simulation/cross_group_comparison/dirfactor/roc_dirfactor.R $WORK_DIR
$WORK_DIR/src/simulation/cross_group_comparison/plot.R $WORK_DIR
```

The figures shown in section 3.2 in the paper are saved in the folder "./results/cross_group_comparison".

### Case Study

Run the following to reproduce the case study in the paper. The figures and tables are saved in the folders "./results/application/figures" and "./results/application/tables" respectively. Change the value of lambda for sensitivity analysis.

```bash
lambda=10
for i in `seq 1 10`; do
$WORK_DIR/src/application/application.R $WORK_DIR 10000 $i $lambda
done
$WORK_DIR/src/application/figures_tables.R $WORK_DIR
```



## References

A. D. Kostic, D. Gevers, H. Siljander, T. Vatanen, T. Hyo ̈tyl ̈ainen, A.-M. H ̈ama ̈la ̈inen, A. Peet, V. Tillmann, P. Po ̈h ̈o, I. Mattila, H. La ̈hdesm ̈aki, E. A. Franzosa, O. Vaar- ala, M. de Goffau, H. Harmsen, J. Ilonen, S. M. Virtanen, C. B. Clish, M. Oreˇsiˇc, C. Huttenhower, M. Knip, and R. J. Xavier. The dynamics of the human infant gut microbiome in development and in progression toward type 1 diabetes. Cell Host & Microbe, 17(2):260–273, 2020/12/18 2015. doi: 10.1016/j.chom.2015.01.001. URL https://doi.org/10.1016/j.chom.2015.01.001.

Y. Cao, W. Lin, and H. Li. Large covariance estimation for compositional data via composition-adjusted thresholding. Journal of the American Statistical Association, 114 (526):759–772, 2019. doi: 10.1080/01621459.2018.1442340. URL https://doi.org/10. 1080/01621459.2018.1442340.

B. Ren, S. Bacallado, S. Favaro, T. Vatanen, C. Huttenhower, and L. Trippa. Bayesian mixed effects models for zero-inflated compositions in microbiome data analysis. Ann. Appl. Stat., 14(1):494–517, 03 2020. doi: 10.1214/19-AOAS1295. URL https://doi. org/10.1214/19-AOAS1295.







