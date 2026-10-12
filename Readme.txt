Simulation Instructions

Install the required R packages: MASS, grpreg, fda, survcomp, locpol, survival, ggplot2, plyr, gridExtra, dplyr, purrr, readr, and tidyr.

Run the scripts in R or RStudio. Before each step, set the specified folder as the working directory. In RStudio, select Session > Set Working Directory > Choose Directory. Enter the source() commands below in the R console.

1. Run Scenarios 1-4

In each scenario folder, run the following files in order:

Working directory: Scenario1
Execution order: main.R -> saveRdata.R -> saveEPV_Scenario1.R

Working directory: Scenario2
Execution order: main.R -> saveRdata.R -> saveEPV_Scenario2.R

Working directory: Scenario3
Execution order: main.R -> saveRdata.R -> saveEPV_Scenario3.R

Working directory: Scenario4
Execution order: main.R -> saveRdata.R -> saveEPV_Scenario4.R

For example, with Scenario1 as the working directory:

source("main.R")
source("saveRdata.R")
source("saveEPV_Scenario1.R")

2. Run Comparisons without the Variability Component

Set each of the following folders as the working directory in turn:

- Scenario1/without_variability
- Scenario2/without_variability
- Scenario3/without_variability
- Scenario4/without_variability

Run these files in each folder:

source("main.R")
source("saveRdata_nodev.R")

Copy all generated files ending in nodev.RData to the corresponding parent scenario folder. For example, copy the *nodev.RData files from Scenario1/without_variability to Scenario1.

3. Run Scenarios 5-8

In each scenario folder, run the following files in order:

Working directory: Scenario5
Execution order: main.R -> saveRdata_Scenario5.R -> saveEPV_Scenario5.R

Working directory: Scenario6
Execution order: main.R -> saveRdata_Scenario6.R -> saveEPV_Scenario6.R

Working directory: Scenario7
Execution order: main.R -> saveRdata_Scenario7.R -> saveEPV_Scenario7.R

Working directory: Scenario8
Execution order: main.R -> saveRdata_Scenario8.R -> saveEPV_Scenario8.R

For example, with Scenario5 as the working directory:

source("main.R")
source("saveRdata_Scenario5.R")
source("saveEPV_Scenario5.R")

4. Generate Figures

After completing the preceding steps, run each plotting script in its corresponding scenario folder:

Working directory: Scenario1
Plotting script: Figure S1-w.R
Manuscript figure: Figures 1-2

Working directory: Scenario1
Plotting script: plot_curves.R
Manuscript figure: Figure A1

Working directory: Scenario2
Plotting script: Figure S2-w.R
Manuscript figure: Figures A2-A3

Working directory: Scenario3
Plotting script: Figure S3-w.R
Manuscript figure: Figures A4-A5

Working directory: Scenario4
Plotting script: Figure S4-w.R
Manuscript figure: Figures A6-A7

Working directory: Scenario5
Plotting script: FigureS5.R
Manuscript figure: Figure 3

Working directory: Scenario5
Plotting script: stabilityThreshold_r05c20varyn.R
Manuscript figure: Figure A8

Working directory: Scenario6
Plotting script: FigureS6.R
Manuscript figure: Figure A9

Working directory: Scenario7
Plotting script: FigureS7.R
Manuscript figure: Figure A10

Working directory: Scenario8
Plotting script: FigureS8.R
Manuscript figure: Figure A11

For example, with Scenario1 as the working directory:

source("Figure S1-w.R")
source("plot_curves.R")

Keep plot_curves.RData in Scenario1 when running plot_curves.R.

Keep the raw CSV files and generated RData files in their respective working directories. Do not combine CSV files from different scenario folders.