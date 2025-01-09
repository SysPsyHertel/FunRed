# Function


fredundancy <- function(functions, abundance, n_reference = NULL) {  # n_reference is optional

  # Normalize function vector if the sum is not 1
  normalize_vector <- function(x) {
    if (sum(x) != 1) {
      return(x / sum(x))
    }
    return(x)
  }

  # Validate that abundance sums to 1
  validate_abundance <- function(abundance) {
    if (sum(abundance) != 1) {
      stop("Error: Abundances do not sum up to 1.")
    }
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
  validate_abundance(abundance)

  # Sample-based redundancy
  a_uniform_sample <- rep(1 / length(abundance), length(abundance))
  x_sample <- rbind(functions, a_uniform_sample)
  sample_based <- calculate_kl_divergence(x_sample)

  # Initialize reference-based redundancy as NULL (in case it's not calculated)
  reference_based <- NULL

  # Only calculate reference-based redundancy if n_reference is provided
  if (!is.null(n_reference)) {
    functions_gz <- functions[functions > 0]  # Only functions > 0
    zero_f <- rep(0, n_reference - length(functions_gz))  # Pad with zeros
    functions_reference <- c(functions_gz, zero_f)
    a_uniform_reference <- rep(1 / n_reference, n_reference)
    x_reference <- rbind(functions_reference, a_uniform_reference)
    reference_based <- calculate_kl_divergence(x_reference)
  }

  # Abundance-based redundancy
  x_abundance <- rbind(functions, abundance)
  abundance_based <- calculate_kl_divergence(x_abundance)

  # Interdependency
  abundance_I <- abundance[functions_gz > 0]
  x_interdependency <- rbind(functions_gz, abundance_I)
  interdependency <- calculate_kl_divergence(x_interdependency, negative = FALSE)

  # Return results as a list
  result <- list(
    sample_based = sample_based,
    abundance_based = abundance_based,
    interdependency = interdependency
  )

  # Add reference_based to the results if it was calculated
  if (!is.null(reference_based)) {
    result$reference_based <- reference_based
  }

  return(result)
}
