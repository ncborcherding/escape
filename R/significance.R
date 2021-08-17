#' Perform significance testing between groups and enrichement scores.
#' 
#' This functions takes the enrichment scores and performs statistical 
#' testing to evaluate the difference by group selected. The function 
#' can perform 3 tests: 1) linear model based on the limma package, 
#' 2) Welch's T test, and 3) one-way ANOVA. The output includes 
#' adjusted p-values based on the Benjamini Hochberg method. 
#' 
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param group The parameter to group for the comparison, should a column of 
#' the enriched input
#' @param fit The test used for significance, either ANOVA or T.test
#'
#' @import limma
#' @importFrom dplyr select_if
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' output <- getSignificance(ES2, group = "Type", fit = "T.test")
#' 
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return Data frame of test statistics
getSignificance <- function(enriched, group = NULL, 
                        fit = NULL) {
    fit <- match.arg(fit,  choices = c("T.test", "ANOVA"))
    group2 <- enriched[,group]
    gr_names <- unique(group2)
    input <- select_if(enriched, is.numeric)
    output <- NULL
    if (fit == "T.test") {
        if (length(unique(group2)) != 2) {
            message("Ensure the group selection has only two levels for T.test 
                fit") 
        } else {
        out <- lapply(input, function(x) t.test(x ~ group2))
        for (i in seq_along(out)) {
            df <- out[[i]]
            mat <- c(df$statistic, df$conf.int[[1]], df$conf.int[[2]],
                        df$estimate[1], df$estimate[2], df$p.value)
            output <- rbind(output,mat)
        }
        output <- as.data.frame(output)
        colnames(output) <- c("T_statistic", "CI_low", "CI_high", 
            paste0("mean_", gr_names[1]), paste0("mean_", gr_names[2]), 
            "p.value")
        rownames(output) <- colnames(input)
        output$FDR <- p.adjust(output$p.value) }
    } else if (fit == "ANOVA") {
        if (length(unique(group2)) <= 2) {
            message("Ensure the group selection has more than two levels 
                for ANOVA fit") }
        out <- lapply(input, function(x) aov(x ~ group2))
        for (i in seq_along(out)) {
            df <- out[[i]]
            fval <- summary(df)[[1]]$'F value'[[1]]
            pval <- summary(df)[[1]]$'Pr(>F)'[[1]]
            output <- rbind(output, c(fval, pval))
        }
        output <- as.data.frame(output)
        colnames(output) <- c("f.value", "p.value")
        rownames(output) <- colnames(input)
        output$FDR <- p.adjust(output$p.value)
    }
    return(output)
}
