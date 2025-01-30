# FunRed

**R package** for computing **functional redundancy** and **functional interdependency** measures for a single trait using relative entropy.  
  
Relative Entropy or the Kullback-Leibler divergence is calculated using the [`philentropy`](https://cran.r-project.org/package=philentropy) package.


## Installation

You can install this package from GitHub using the `devtools` package:

```R
# Install required dependencies first
install.packages("philentropy")

# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install FunRed package from GitHub
devtools::install_github("SysPsyHertel/FunRed")

```

## Functions

### 1. fredundancy(functions, abundance, n_reference)

Compute measures of functional redundancy and functional interdependency.

#### **Arguments**
- `abundance`: A vector of non-zero abundances, corresponding to each species in the community.
- `functions`: The paired vector of functions of the abundance vector. Needs to be the same length as the abundance vector.
- `n_reference`: Integer value, corresponding to the number of species in the reference that can perform the function (optional).



#### **Returns**
- `sample_based`: A coefficient corresponding to taxon sample-based functional redundancy.
- `reference_based`: A coefficient corresponding to taxon reference-based functional (optional, if n_reference is set).
- `abundance_based`: A coefficient corresponding to abundance-based functional redundancy.
- `interdependency`: A coefficient corresponding to functional interdependency.

#### **Example Usage**

```R
# Assuming 'functions' and 'abundance' are defined:
result <- fredundancy(functions, abundance, n_reference = 7)
```
### Toy Example and description

Below is a toy example to demonstrate how the `fredundancy` function works.

```R
functions <- c(0.8, 0.1, 0.05, 0.05, 0)
abundance <- c(0.2, 0.1, 0.05, 0.05, 0.6)
n_reference <- 7

result <- fredundancy(functions, abundance, n_reference)
print(result)
```
**Output:**

```R
$sample_based
[1] -0.9010913

$reference_based
[1] -1.237564

$abundance_based
[1] -1.109035

$interdependency
[1] 0.1927448

```
#### **Description**
Note that, except when computing functional interdependency, the Kullback-Leibler divergence is negated in all calculations.  

- To compute sample taxon-based functional redundancy, we calculate the Kullback-Leibler divergence between the function vector
and a uniformly distributed vector of the same length as the number of species in a sample, i.e. the number of nonzero
elements in the abundance vector.  
In the toy example, we compute the Kullback-Leibler divergence for sample taxon-based
redundancy between  `(0.8, 0.1, 0.05, 0.05, 0)` and `(1/5, 1/5, 1/5, 1/5, 1/5)`. 
- For the calculation of reference taxon-based functional redundancy, the parameter `n_reference` corresponds to the total number of species
in the reference that can perform the function. The nonzero values in the function vector remain unchanged, while zeros are
added until its total length equals `n_reference`. Then, the Kullback-Leibler divergence is computed between this padded vector
and a uniform distribution of the same length as `n_reference`.  
In the toy example, we compute the Kullback-Leibler divergence for sample taxon-based redundancy between
`(0.8, 0.1, 0.05, 0.05, 0, 0 ,0)` and `(1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7)`.
- To compute abundance-based functional redundancy, we calculate the Kullback-Leibler divergence between the function vector and the abundance vector.  
In the toy example, we compute the Kullback-Leibler divergence for abundance-based redundancy 
between `(0.8, 0.1, 0.05, 0.05, 0)` and `(0.2, 0.1, 0.05, 0.05, 0.6)`.
- For the computation of functional interdependency only non-zero values of the function vector are considered.
For the abundance vector, only the abundances of species that have non-zero function are considered. The abundance
vector is then divided by the sum of the restricted vector, such that its components sum up to 1.
In the toy example, we compute the Kullback-Leibler divergence for sample functional interdependency between
`(0.8, 0.1, 0.05, 0.05)` and `(0.2, 0.1, 0.05, 0.05)/(0.2+0.1+0.05+0.05)`.
