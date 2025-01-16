
# FunRed

R package for calculating functional redundancy measures of a single trait based on relative entropy

## Installation

You can install this package from GitHub using the `devtools` package:

```R
# Install the devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the package
devtools::install_github("faesslerd/FunRed")
```

## Functions

### 1. fredundancy(abundance, functions, n_reference)

Compute measures of functional redundancy and functional interdependency.

#### Arguments
- `abundance`: A vector of non-zero abundances, corresponding to each species in the community.
- `functions`: The paired vector of functions of the abundance vector. Needs to be the same length as the abundance vector.
- `n_reference`: Integer value, corresponding to of the number of species in the reference that can perform the function (Optional).


#### Returns
- `sample_based`: A coefficient corresponding to taxon sample-based redundancy.
- `reference_based`: A coefficient corresponding to taxon reference-based redundancy (optional, if n_reference is set).
- `abundance_based`: A coefficient corresponding to abundance-based redundancy.
- `interdependency`: A coefficient corresponding to functional interdependency.

#### Example

```R
functions <- c(0.8, 0.1, 0.05, 0.05, 0)
abundance <- c(0.2, 0.1, 0.05, 0.05, 0.6)
n_reference <- 7
> fredundancy(functions, abundance, n_reference)
$sample_based
[1] -0.9010913

$reference_based
[1] -1.237564

$abundance_based
[1] -1.109035

$interdependency
[1] 0.1927448

```

