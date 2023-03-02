library(dplyr)
library(nimble)

# create object to store the output
brier_out <- list(
  brier_df = data.frame(
    model = c(
      "null", "covars", "fall","covar_fall","covar_inxs",  "covar_inxs_fall"
    ),
    score = rep(NA,6)
  ),
  brier_mcmc = list(
    null = NA,
    covars = NA,
    fall = NA,
    covar_fall = NA,
    covar_inxs = NA,
    covar_inxs_fall = NA
  ),
  script = c(
    "./R/intercept_brier_score.R",
    "./R/spatial_covariates_brier_score.R",
    "./R/fall_term_brier_score.R",
    "./R/spatial_covariates_fall_brier_score.R",
    "./R/spatial_covariates_urbless_brier_score.R",
    "./R/spatial_covariates_fall_urbless_brier_score.R"
  )
)

# iterate through each script, and then save the output in 
#  the correct location.

for(bs in 1:nrow(brier_out$brier_df)){
    cat(
      paste0(
        paste0(
          rep(
            "=",
            nchar(brier_out$script[bs])
          ),
          collapse = ""
        ),"\n",
        brier_out$script[bs],"\n",
        paste0(
          rep(
            "=",
            nchar(brier_out$script[bs])
          ),
          collapse = ""
        ),"\n"
      )
    )
    tryCatch(
      suppressWarnings(
        source(
          brier_out$script[bs]
        )
      ),
      error = function(e) cat("\nScript did not run :_( \n"),
      finally = cat("\nScript ran :-) \n")
    )
    if(
      exists(
        'brier_post'
      )
    ){
      brier_out$brier_df$score[bs] <- mean(brier_post, na.rm = TRUE)
      brier_out$brier_mcmc[[bs]] <- brier_post
    }
    objects_to_drop <- ls()[-grep('brier_out', ls())]
    rm(list = c(objects_to_drop, "objects_to_drop"))
    gc()
}
