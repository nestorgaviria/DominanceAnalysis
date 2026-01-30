Overview
This repository contains R functions for performing category-level dominance analysis on environmental controls of denudation rates across the Chilean Andes. The approach combines ensemble model generation with LMG dominance analysis to robustly quantify the relative importance of five geomorphic process categories (Topography, Climate, Land Cover, Surface Materials, and Seismotectonics) while accounting for multicollinearity and limited sample sizes.

📂 Repository Structure

├── generate_models_ensemble_no_duplicates.R  # Model ensemble generator
├── Mdls_dominance_Analysis_OLS_AIC.R         # Dominance analysis workflow
├── run_analysis.R                             # Example script to run both functions
├── config.R                                   # Configuration parameters
└── README.md                                  # This file

🎯 Purpose
Problem: Environmental predictors in Earth surface systems are inherently correlated (e.g., precipitation, vegetation, and elevation co-vary). Standard regression approaches produce unstable coefficients and ambiguous interpretations when multicollinearity is present.
Solution: This workflow:

Generates a diverse ensemble of linear models with controlled complexity (5-6 predictors)
Selects top-performing models using Akaike Information Criterion (AIC)
Performs Lindeman-Merenda-Gold (LMG) dominance analysis on selected models
Aggregates results at the category level with uncertainty quantification


🔧 Functions

1. generate_models_ensemble_no_duplicates()
File: generate_models_ensemble_no_duplicates.R
Purpose: Generates a diverse ensemble of candidate models while preventing duplicates and enforcing scientific constraints.
Key Features:

Ensures minimum 1 predictor per category (all processes represented)
Limits maximum 3 predictors per category (prevents single-category dominance)
Constrains model size to 5-6 predictors (appropriate for small sample sizes)
Tracks predictor usage to balance selection across variables
Generates reference models (full and subset end-members per category)
Creates unique signatures to prevent duplicate models

Arguments:
generate_models_ensemble_no_duplicates(
  pool,              # Named list of predictors grouped by category
  n_random_models,   # Number of random interaction models to generate
  n_subsets = 2,     # Number of subset reference models per category
  seed               # Random seed for reproducibility
)
Returns:

Named list where each element is a character vector of predictor names
Model names indicate type: Ref_Full_[Category], Ref_Sub[N]_[Category], Rand_[N]

