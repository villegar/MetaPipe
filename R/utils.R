#' Check variable types
#'
#' Check variable types in a data frame. Useful whenever the data must be 
#' filter prior some numerical analysis or procwessing
#' 
#' @param raw_data data frame containing the raw data 
#' @param excluded_columns list of excluded columns (won't be check)
#' @param numeric boolean to indicate the operation mode. If \code{TRUE} will
#'     look up for non-numeric variables and exclude them, otherwise it will
#'     filter numeric variables, including those with the following classes:
#'     - \code{complex} e.g. \code{1+i}
#'     - \code{integer} e.g. \code{1L}
#'     - \code{numeric} e.g. \code{1}, \code{1.01}
#' @param quiet boolean to hide warning messages
#'     
#' @return vector with column indices not meeting the checking mode, if 
#' \code{excluded_columns != NULL}, returns a vector with both originally 
#' excluded and newly found column indices
#' @export
#'
#' @examples
#' example_data <- data.frame(ID = c(1,2,3,4,5), 
#'                            P1 = c("one", "two", "three", "four", "five"), 
#'                            T1 = rnorm(5))
#' MetaPipe::check_types(example_data)
#' MetaPipe::check_types(example_data, numeric = FALSE)
#' MetaPipe::check_types(example_data, excluded_columns = 1, numeric = FALSE)
#' MetaPipe::check_types(example_data, quiet = FALSE)
check_types <- function(raw_data,
                        excluded_columns = NULL,
                        numeric = TRUE,
                        quiet = TRUE) {
  # Extract data type (class) of each variable (column)
  var_types <- sapply(raw_data, class)
  if (!is.null(excluded_columns)) {
    var_types <- var_types[-excluded_columns]
  }
  if (numeric) { # if numeric = TRUE, find non-numeri columns
    idx <- !(var_types %in% c("numeric", "integer", "complex"))
  } else { # if numeric = FALSE, find numeric columns
    idx <- var_types %in% c("numeric", "integer", "complex")
  }
  # If any variables were found to meet the search criteria
  if (any(idx)) {
    if (!quiet) {
      warning(paste0("\nThe following variable",
                     ifelse(sum(idx) > 1, "s ", " "), # Plural
                     "will be excluded, as ",
                     ifelse(sum(idx) > 1, "they are ", "it is "), # Plural
                     ifelse(numeric, "non-", ""), # Check mode
                     "numeric: ",
                     paste0("\n - ", names(var_types)[idx], collapse = "")))
    }
    # Append to the original list
    excluded_columns <- unique(c(excluded_columns, which(idx)))
  }
  return(excluded_columns)
}
