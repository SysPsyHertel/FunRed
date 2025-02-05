#' Compute measures of functional redundancy and functional interdependency based on relative entropy.
#' Relative Entropy is calculated using the *KL* function of the [`philentropy`](https://cran.r-project.org/package=philentropy) package
#' with the following parameters: `test.na = TRUE` (to handle missing values), `unit = "log"` (to compute divergence using the natural logarithm), `est.prob = NULL` (indicating that empirical probabilities were used), and `epsilon = 1e-05` (a small constant used to avoid undefined log calculations caused by zero probabilities)
#' @param functions The paired vector of functions for the abundance vector. Needs to be the same size as the abundance vector.
#'                  The i-th component of \code{functions} corresponds to the i-th abundance of the species, encoded in the abundance vector.
#'                  Gets transformed into relative frequencies, such that the components sum up to 1 and its components correspond to relative frequencies. Might contain zeros.
#' @param abundance A vector of non-zero abundances representing each species in the community. The vector is checked to ensure the abundances sum to 1.
#'                    If they do not, the vector is normalized, provided the sum is within an acceptable tolerance (*|1 - sum(abundances)| < 1e-5*). If the sum falls outside this tolerance, an error message is returned.
#' @param n_reference Optional value to compute reference-based redundancy. It is an integer value of the number of species in the reference
#'                    that can perform the function.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sample_based}: A coefficient corresponding to taxon sample-based redundancy.
#'   \item \code{reference_based}: A coefficient corresponding to taxon reference-based redundancy (optional, if \code{n_reference} is set).
#'   \item \code{abundance_based}: A coefficient corresponding to abundance-based redundancy.
#'   \item \code{interdependency}: Interdependency measure.
#' }
#'
#' @examples
#' functions <- c(0.8, 0.1, 0.05, 0.05, 0)
#' abundance <- c(0.2, 0.1, 0.05, 0.05, 0.6)
#' n_reference <- 7
#' fredundancy(functions, abundance, n_reference)
#'
#' @importFrom philentropy KL
#' @export

fredundancy <- function(functions, abundance, n_reference = NULL) {  # n_reference is optional

  # Normalize function vector if the sum is not 1
  normalize_vector <- function(x) {
    if (sum(x) != 1) {
      return(x / sum(x))
    }
    return(x)
  }


  # Validate that abundance sums to 1
  validate_abundance <- function(abundance, tolerance = 1e-5) {
    total_abundance <- sum(abundance)

    # Check if the sum is exactly 1
    if (abs(total_abundance - 1) == 0) {
      return(abundance)
    }

    # Check if the sum is within the tolerance range
    if (abs(total_abundance - 1) <= tolerance) {
      abundance <- abundance / total_abundance
      return(abundance)
    }
    # If sum is outside the tolerance range, stop with an error
    stop("Abundances do not sum up to 1 within the allowed tolerance!")
  }

  # Helper function to calculate KL divergence and return the numeric result
  calculate_kl_divergence <- function(x1, negative = TRUE) {
    kl_value <- suppressMessages(KL(x1, test.na = TRUE, unit = "log", est.prob = NULL, epsilon = 1e-05))
    # Optionally apply the negative sign
    if (negative) {
      return(-as.numeric(kl_value))  # Negative by default
    } else {
      return(as.numeric(kl_value))  # No negative sign for interdependency
    }
  }

  # Normalize functions and check abundance
  functions <- normalize_vector(functions)
  abundance <- validate_abundance(abundance)
  functions_gz <- functions[functions > 0]  # Only functions > 0


  # Sample-based redundancy
  a_uniform_sample <- rep(1 / sum(abundance>0), sum(abundance>0))
  zero_sample <- rep(0, (sum(abundance>0))-length(functions_gz))  # Pad with zeros
  functions_sample <- c(functions_gz, zero_sample)
  x_sample <- rbind(functions_sample, a_uniform_sample)
  sample_based <- calculate_kl_divergence(x_sample)

  # Initialize reference-based redundancy as NULL (in case it's not calculated)
  reference_based <- NULL

  # Only calculate reference-based redundancy if n_reference is provided
  if (!is.null(n_reference)) {
    zero_reference <- rep(0, n_reference - length(functions_gz))  # Pad with zeros
    functions_reference <- c(functions_gz, zero_reference)
    a_uniform_reference <- rep(1 / n_reference, n_reference)
    x_reference <- rbind(functions_reference, a_uniform_reference)
    reference_based <- calculate_kl_divergence(x_reference)
  }

  # Abundance-based redundancy
  x_abundance <- rbind(functions, abundance)
  abundance_based <- calculate_kl_divergence(x_abundance)

  # Interdependency
  # Consider only abundances that have a function and normalize it that they sum up to 1
  abundance_I <- abundance[functions > 0]
  if(sum(abundance_I)==0){
    abundance_I <- 0
  } else {
      abundance_I <- abundance_I/sum(abundance_I)

  }
  x_interdependency <- rbind(functions_gz, abundance_I)
  interdependency <- calculate_kl_divergence(x_interdependency, negative = FALSE)

  # Return results in the specified order
  result <- list(
    sample_based = sample_based,
    reference_based = reference_based,  # Will be NULL if not calculated
    abundance_based = abundance_based,
    interdependency = interdependency
  )

  # Remove reference_based if it was not calculated (NULL)
  result <- result[!sapply(result, is.null)]

  return(result)
}
