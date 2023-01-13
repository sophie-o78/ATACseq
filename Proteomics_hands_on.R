library(protti)
library(magrittr)
library(dplyr)

#### Qualtiy Control (QC) Workflow

set.seed(123) 

# by setting the seed we are making sure that the random object generation can be reproduced

data <- create_synthetic_data(n_proteins = 100,
                              frac_change = 0.05,
                              n_replicates = 3,
                              n_conditions = 2,
                              method = "effect_random",
                              additional_metadata = TRUE)

# This function creates a synthetic limited proteolysis proteomics dataset that can be used to test functions while knowing the ground truth.

## Coeficient of variation:

input <- data %>%
  # as the data is log2 transformed, we need to transform it back before calculating the CVs
  mutate(raw_intensity = 2^peptide_intensity_missing) 

qc_cvs(data = input,
       grouping = peptide,
       condition = condition,
       intensity = raw_intensity, 
       plot = FALSE)
#> # A tibble: 2 × 3
#>   condition   median_cv median_cv_combined
#>   <chr>           <dbl>              <dbl>
#> 1 condition_2      6.06               7.49
#> 2 condition_1      6.07               7.49

qc_cvs(data = input,
       grouping = peptide,
       condition = condition,
       intensity = raw_intensity, 
       plot = TRUE, 
       plot_style = "violin")

## Number of identifications (IDs):

qc_ids(data = input,
       sample = sample,
       grouping = protein,
       intensity = peptide_intensity_missing,
       condition = condition, 
       plot = FALSE)
#> # A tibble: 6 × 3
#>   condition   sample   count
#>   <chr>       <chr>    <int>
#> 1 condition_1 sample_3    98
#> 2 condition_2 sample_4    97
#> 3 condition_2 sample_5    99
#> 4 condition_2 sample_6    98
#> 5 condition_1 sample_2    99
#> 6 condition_1 sample_1    98

qc_ids(data = input,
       sample = sample,
       grouping = protein,
       intensity = peptide_intensity_missing,
       condition = condition, 
       title = "Protein identifications per sample",
       plot = TRUE)

## Peptide types:

qc_peptide_type(data = input,
                sample = sample, 
                peptide = peptide, 
                pep_type = pep_type, 
                method = "intensity", 
                intensity = raw_intensity, 
                plot = TRUE, 
                interactive = FALSE)

## Run intensities:

qc_intensity_distribution(data = input,
                          sample = sample,
                          grouping = peptide,
                          intensity_log2 = peptide_intensity_missing,
                          plot_style = "boxplot")

## Charge states:

qc_charge_states(data = input, 
sample = sample, 
grouping = peptide, 
charge_states = charge, 
method = "intensity",
intensity = raw_intensity, 
plot = TRUE)

## Missed cleavages:

qc_missed_cleavages(data = input, 
                    sample = sample, 
                    grouping = peptide, 
                    missed_cleavages = n_missed_cleavage, 
                    method = "intensity",
                    intensity = raw_intensity, 
                    plot = TRUE)

## Sequence coverage:

qc_sequence_coverage(data = input, 
                     protein_identifier = protein, 
                     coverage = coverage)

## Sample correlation:

qc_sample_correlation(data = input,
                      sample = sample, 
                      grouping = peptide, 
                      intensity_log2 = peptide_intensity_missing, 
                      condition = condition, 
                      interactive = FALSE)

## Principal component analysis (PCA):

qc_pca(
  data = data,
  sample = sample,
  grouping = peptide,
  intensity = peptide_intensity_missing,
  condition = condition,
  digestion = NULL,
  plot_style = "scree"
)

qc_pca(
  data = data,
  sample = sample,
  grouping = peptide,
  intensity = peptide_intensity_missing,
  condition = condition,
  components = c("PC1", "PC2"), 
  plot_style = "pca"
)

#### Data analysis:

diff_abundance_data <- data %>%
  assign_missingness(
    sample = sample,
    condition = condition,
    grouping = peptide,
    intensity = peptide_intensity_missing,
    ref_condition = "condition_1",
    completeness_MAR = 0.7,
    completeness_MNAR = 0.25,
    
  ) %>%
  calculate_diff_abundance(
    sample = sample,
    condition = condition,
    grouping = peptide,
    intensity_log2 =peptide_intensity_missing ,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    
  ) 


## p-value distribution:

pval_distribution_plot(data = diff_abundance_data,
                       grouping = peptide,
                       pval = pval
)

## Volcano plot:

volcano_plot(
  data = diff_abundance_data,
  grouping = peptide,
  log2FC = diff,
  significance = pval,
  method = "target",
  target_column = peptide,
  target = "peptide_4_6",
  x_axis_label = "log2(fold change) condition_1 vs condition_2",
  significance_cutoff = c(0.05, "adj_pval") 
)