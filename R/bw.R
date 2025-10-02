#' Estimating Rotemberg Weights for a Bartik "Instrument"
#'
#' `bw()` estimates the Rotemberg weights for a Bartik "instrument" outlined in
#' Goldsmith-Pinkham, Sorkin and Swift (2018) and the just-identified IV
#' estimates using local industry shares, which are the actual instruments.
#'
#' The typical Bartik "instrument" has two components: the local industry
#' shares (usually at the level of location-industry or location-year-industry),
#' and the overall industry growth rates (usually at the level of industry or
#' industry-year). The outcome variable and the endogenous variable are at the
#' level of location or location-year.
#'
#' Because the key variables are at three different levels, `bw()` proceeds by
#' taking three different datasets: (1) a main `master` data frame containing
#' the dependent variable (`y`), the causal variable of interest (`x`), a set of
#' control variables (`controls`), and the weighted variables (`weight`); (2) a
#' `local` data frame containing the set of local industry shares (`Z`); and (3)
#' a `global` data frame containing the overall industry growth rates (`G`). At
#' the moment, it is necessary to transform the "local" dataset from long format
#' to wide format.
#'
#' @md
#' @param master The master data frame.
#' @param y A string for outcome variable. It should be a variable in `master`.
#' @param x A string for the endogenous variable. It should be a variable in
#'   `master`.
#' @param controls A string or character vector for the control variables.
#'   `controls` are optional and should be in `master`.
#' @param weight A string for the weighted variable. `weight` is optional and
#'   should be in `master`.
#' @param local The local data frame. It should be in wide format.
#' @param Z A string or character vector for the the local industry shares.
#' @param global The global data frame.
#' @param G A string for the the overall industry growth rates.
#'
#' @return A data frame built from `global` with two additional variables:
#'   `alpha` (the Rotemberg weights), and `beta` (the just-identified IV
#'   estimates).
#'
#' @importFrom tibble as_tibble
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartik.weight, .registration = TRUE
#'
#' @examples
#' library(bartik.weight)
#'
#' index = c("czone", "year")
#' y = "d_sh_empl_mfg"
#' x = "d_tradeusch_pw"
#' controls = c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
#'   "reg_escen", "reg_wscen", "reg_mount", "reg_pacif", "l_sh_popedu_c",
#'   "l_sh_popfborn", "l_sh_empl_f", "l_sh_routine33", "l_task_outsource",
#'   "t2", "l_shind_manuf_cbp")
#' weight = "timepwt48"
#' Z = setdiff(names(ADH_local_wide), index)
#' G = "trade_"
#'
#' bw(ADH_master, y, x, controls, weight, ADH_local_wide, Z, ADH_global, G)
#'
#' @export
bw = function(master, y, x, controls = NULL, weight = NULL, local, Z, global, G) {

    # Parsing the master file
    y = master[[y]]
    x = master[[x]]
    n = length(x)

    if (is.null(weight)) {
        weight_vec = rep(1.0, n)
    } else {
        weight_vec = as.numeric(master[[weight]])
    }

    if (is.null(controls)) {
        WW <- matrix(1, n, 1)
    } else {
        W <- as.matrix(master[controls])
        # If controls is an empty character vector, W will have 0 columns.
        WW <- cbind(W, rep(1.0, n))
    }
    # Force to plain numeric matrix (drop potential tibble classes)
    storage.mode(WW) <- "double"
    WW <- as.matrix(WW)
    if (!is.matrix(WW)) stop("Internal error: WW not a matrix after coercion.")

    # Parsing the local file
    Z <- as.matrix(local[Z])
    storage.mode(Z) <- "double"
    if (!is.matrix(Z)) stop("Z must subset to a matrix after coercion.")
    if (nrow(Z) != n) stop("Row mismatch: nrow(Z) != length(x).")

    # Parsing the global file
    G <- global[[G]]
    # Validate G length: should match number of instruments (columns of Z)
    if (length(G) != ncol(Z)) {
        stop(sprintf("Length of G (%d) must equal number of instruments / columns of Z (%d).",
                     length(G), ncol(Z)))
    }
    G <- as.numeric(G)

    # Compute the Rotemberg weights (alpha) and the just-identified coefficients (beta)
    alpha_beta <- ComputeAlphaBeta(y, x, WW, weight_vec, Z, G)

    # Return a tibble
    tibble::as_tibble(cbind(global, alpha = alpha_beta[[1]], beta = alpha_beta[[2]]))
}
