#' Perform significance testing between groups and enrichement scores.
#' 
#' This functions takes the enrichment scores and performs statistical 
#' testing to evaluate the difference by group selected. The function 
#' can perform 5 tests: 1) Welch's T test (T.test), 2) Logistic 
#' Regression (LR), 3) Wilcoxon Rank Sum Test (Wilcoxon), 
#' 4) one-way ANOVA (ANOVA), and 5) Kruskal-Wallis (KW). The latter 
#' two output will include the individual comparisons between groups 
#' using TukeyHSD for ANOVA and pairwise Wilcoxon Rank Sum Test 
#' for KW. The output includes adjusted p-values based on the 
#' Benjamini Hochberg method. 
#' 
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param group The parameter to group for the comparison, should a column of 
#' the enriched input
#' @param gene.sets Names of gene sets to compare
#' @param fit The test used for significance, 2 group: Wilcoxon, LR, T.test.
#' Multigroup: ANOVA or KW.
#' @importFrom dplyr select_if
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom stats TukeyHSD median glm wilcox.test pairwise.wilcox.test kruskal.test
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
                            gene.sets = NULL,
                            fit = NULL) {
    fit <- match.arg(fit,  choices = c("T.test", "ANOVA", "Wilcoxon", "LR", "KW"))
    group2 <- enriched[,group]
    gr_names <- unique(group2)
    if (!is.null(gene.sets)) {
        input <- enriched[,colnames(enriched) %in% gene.sets]
    } else {
        input <- select_if(enriched, is.numeric)
    }
    medians <- get.medians(input, group2)
    output <- NULL
    if (fit == "T.test" || fit == "Wilcoxon" || fit == "LR") {
        if (length(unique(group2)) != 2) {
            message("Ensure the group selection has only two levels for T.test 
                fit") 
        } else {
            if (fit == "T.test") {
                out <- lapply(input, function(x) t.test(x ~ group2))
                stat <- "T"
            } else if (fit == "Wilcoxon") {
                out <- lapply(input, function(x) wilcox.test(x ~ group2))
                stat <- "W"
            }  else if (fit == "LR") {
                group2 <- ifelse(group2 == gr_names[1], 0,1)
                out <- lapply(input, function(x) glm(group3 ~ x, family = "binomial"))
                out <- lapply(out, function(x) tidy(x)[2,])
                stat <- "L"
            }
        for (i in seq_along(out)) {
            df <- out[[i]]
            mat <- c(df$statistic, df$p.value)
            output <- rbind(output,mat)
        }
        output <- as.data.frame(output)
        colnames(output) <- c(paste0(stat, ".statistic"), "p.value")
        }
    } else if (fit == "ANOVA") {
        if (length(unique(group2)) <= 2) {
            message("Ensure the group selection has more than two levels 
                for ANOVA fit") }
        out <- lapply(input, function(x) aov(x ~ group2))
        for (i in seq_along(out)) {
            df <- out[[i]]
            tukey <- TukeyHSD(df)
            ind.p.values <- melt(tukey$group2[,4])
            names.ind.p.values <- gsub("-", "v", rownames(ind.p.values))
            names.ind.p.values <- paste0(names.ind.p.values,".p.value")
            fval <- summary(df)[[1]]$'F value'[[1]]
            pval <- summary(df)[[1]]$'Pr(>F)'[[1]]
            output <- rbind(output, c(fval, pval, t(ind.p.values)))
        }
        output <- as.data.frame(output)
        colnames(output) <- c("f.value", "p.value", names.ind.p.values)
    } else if (fit == "KW") {
        if (length(unique(group2)) <= 2) {
            message("Ensure the group selection has more than two levels 
                for Kruskal-Wallis test")}
        out <- lapply(input, function(x) kruskal.test(x ~ group2))
        out.ind <- lapply(input, function(x) pairwise.wilcox.test(x, group2, p.adjust.method = "BH"))
        for (i in seq_along(out)) {
            ind.p.values <- na.omit(melt(out.ind[[i]]$p.value))
            names.ind.p.values <- paste0(ind.p.values$Var1, "v", ind.p.values$Var2)
            names.ind.p.values <- paste0(names.ind.p.values,".p.value")
            ind.p.values <- ind.p.values[,3]
            Chi.squared <- out[[i]]$statistic 
            pval <- out[[i]]$p.value
            output <- rbind(output, c(fval, pval, t(ind.p.values)))
        }
        output <- as.data.frame(output)
        colnames(output) <- c("Chi.square", "p.value", names.ind.p.values)
    }
    rownames(output) <- colnames(input)
    output$FDR <- p.adjust(output$p.value) 
    output <- cbind.data.frame(output, medians)
    return(output)
}

get.medians<- function(input, group2) {
    input <- cbind.data.frame(group2, input)
    num <- ncol(input)-1
    med <- input %>%
        group_by(group2) %>%
        summarise(across(seq_len(all_of(num)), median))
    med <- as.data.frame(med[,seq_len(num) + 1])
    rownames(med) <- paste0("median.", unique(group2))
    med <- as.data.frame(t(med))
    return(med)
}
