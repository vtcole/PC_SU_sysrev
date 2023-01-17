# Contents of data and code file
We have **data**, **code**, and a **codebook** that summarizes all of the variables in the dataset. We also have a **code walkthrough** that explains a lot of the different decisions made in each of the code files. 

## Data
In the `data` folder, we have:
- The main file `effectsheet.rds`. It is in RDS format because that's easier to read in. In order to read it into R, write: `readRDS(here("data/effectsheet.rds"))`
- Model output written out by the multiple imputation procedure (which have `table_imp<imputation>.csv` in the filename). It is not recommended to analyze these directly!
- Pooled model outputs, which correspond to the tables for the meta-regressions in the paper. These generally have `pooled` in the filename. 

## Codebook
In the `doc` folder, we have `codebook.xlsx`. Each variable in the `effectsheet.rds` is named in the codebook.

## Code

There are three main R scripts:
- `regressions.R` - This runs all the meta-analytic regression models, the one shown in Tables 1-4. It is the biggest file, because it contains all of the multiple imputation steps.
- `figure2.R` - This creates Figure 2, as well as the extension thereof in the supplemental material. 
- `supplemental.R` - This runs all of the supplemental analyses shown in Section I. 

We also have a walkthrough, an RMarkdown file called `allanalyses.Rmd` which renders all of the above code and gives some justification for some of the longer code snippets. 

## Contact

If you have any questions, ideas, feedback, or concerns, please email Veronica Cole at <colev@wfu.edu>. 
