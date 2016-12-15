#' Calculate the enrichment in the co-expression of interacting proteins 
#' relative to non-interacting protein pairs from the same dataset. 
#' 
#' @param interactions A list of interactions, where the first two columns
#' correspond to the accessions of the interacting proteins. 
#' @param expression A dataset of protein expression, used to calculate 
#' coexpression (see \code{\link{ReadExpression}})
#' @param threshold Pearson correlation coefficient to use when calculating fold
#' enrichment and hypergeometric P 
#' 
#' @return A data table containing the Pearson correlation coefficients of all 
#' interacting and non-interacting protein pairs
#' 
#' @examples 
#' interactions <- read.csv("interactions.csv")
#' expression <- ReadExpression("expression.csv")
#' coexpression <- CalculateCoexpression(interactions, expression)
#' enrichment <- CalculateCoexpressionEnrichment(coexpression, interactions, 
#' expression)
#' 
#' @export
CalculateCoexpressionEnrichment <-
  function(coexpression, interactions, expression, threshold = 0.75) {
  # Calculate coexpression for all random non-interacting pairs 
  new.dt <- function() return(
    data.frame(protein.A = character(), protein.B = character(), 
               pcc = numeric()))
  dts <- list()
  dt <- new.dt() 
  
  proteins <- unique(c(interactions[,1], interactions[,2]))
  n <- (length(proteins)^2 - length(proteins)) / 2  - nrow(coexpression)
  message("Calculating coexpression (PCC) for ", n, 
          " non-interacting protein pairs...")
  for (i in 1:length(proteins)) {
    for (j in 1:length(proteins)) {
      if (j <= i) next
      protein.A <- proteins[i]
      protein.B <- proteins[j]
      # Ignore interacting proteins
      if (nrow(interactions[
        (interactions[,1] == protein.A & interactions[,2] == protein.B) | 
        (interactions[,2] == protein.A & interactions[,1] == protein.B),]) > 0)
        next
      correlation <- CalculateCoexpressionPcc(protein.A, protein.B, expression)
      if (!is.na(correlation) & length(correlation) > 0)
        dt <- rbind(dt, list(protein.A, protein.B, correlation))
      if (nrow(dt) > 1000) {
        colnames(dt) <- c("protein.A", "protein.B", "pcc")
        dts[[length(dts) + 1]] <- dt
        dt <- new.dt()
      }
    }
  }
  
  # Merge data tables
  colnames(dt) <- c("protein.A", "protein.B", "pcc")
  dts[[length(dts) + 1]] <- dt
  non.interacting <- ldply(dts, data.frame)
  non.interacting$interact <- F

  # Merge interacting and non-interacting proteins 
  pcc.idx <- which(colnames(coexpression) == "pcc")
  interacting <- coexpression[,c(1:2, pcc.idx)]
  colnames(interacting) <- c("protein.A", "protein.B", "pcc")
  interacting$interact <- T
  compare <- rbind(interacting, non.interacting)
  
  # Calculate hypergeometric P at threshold
  sample.successes <- sum(compare$interact == T & compare$pcc >= threshold)
  pop.successes <- sum(compare$interact == F & compare$pcc >= threshold)
  sample.trials <- sum(compare$interact == T)
  pop.trials <- sum(compare$interact == F)
  p <- phyper(sample.successes, pop.successes, pop.trials - pop.successes, 
              sample.trials, lower.tail = F)
  
  # Print fold enrichment and hypergeometric P
  enrichment <- (sample.successes/sample.trials) / (pop.successes/pop.trials)
  message("Interacting protein pairs are ", format(enrichment, digits=2), 
          "-fold enriched for PCC >= ", threshold)
  message("  Interacting protein pairs: ", 
          format(sample.successes/sample.trials, digits = 4), "%")
  message("  Non-interacting protein pairs: ", 
          format(pop.successes/pop.trials, digits = 4), "%")
  message("  Hypergeometric P = ", format(p, digits=2))
  
  # Return the enrichment data 
  class(compare) <- append(class(compare), "FlavinCoexpressionEnrichment")
  return(compare)
}

# ==============================================================================
#' Plot enrichment in coexpression for interacting proteins 
#' 
#' Plot the fraction of proteins with Pearson correlation coefficient greater 
#' than a threshold relative to non-interacting protein pairs, for all 
#' thresholds in the range [0,1]. 
#'
#' @param enrichment A data table containing the Pearson correlation 
#' coefficients of all interacting and non-interacting protein pairs, as output 
#' by \code{\link{CalculateCoexpressionEnrichment}}
#'
#' @return A \code{ggplot2} plot
#' 
#' @examples
#' interactions <- read.csv("interactions.csv")
#' expression <- ReadExpression("expression.csv")
#' enrichment <- CalculateCoexpressionEnrichment(interactions, expression)
#' plot <- PlotCoexpressionEnrichment(enrichment)
#' plot(plot)
#'
#' @export
PlotCoexpressionEnrichment <- function(enrichment) {
  # Check input data 
  if (!("FlavinCoexpressionEnrichment" %in% class(enrichment)))
    stop("Must provide a data table created by CalculateCoexpressionEnrichment")
  
  # Calculate the fraction of protein pairs at each threshold T with PCC > T
  threshold.df <- data.frame(threshold = numeric(), interact = logical(), 
                             fraction = numeric())  
  for (threshold in seq(0,1,0.01)) {
    subset <- subset(enrichment, pcc >= threshold)
    interacting.fraction <- sum(subset$interact) / sum(enrichment$interact)
    not.interacting.fraction <- sum(!subset$interact) / 
      sum(!enrichment$interact)
    threshold.df[nrow(threshold.df) + 1,] <-
      list(threshold, T, interacting.fraction)
    threshold.df[nrow(threshold.df) + 1,] <- 
      list(threshold, F, not.interacting.fraction)
  }

  # Create plot
  plot <- ggplot(threshold.df, aes(x = threshold, y = fraction, 
                                   colour = interact, group = interact)) +
    geom_line() + 
    labs(x = "Threshold", y = "Percentage of proteins") + 
    scale_colour_manual(values = flavin_palette, name = "Interacting") +
    flavin_theme
  return(plot)
}

# ==============================================================================
#' Plot coexpression of interacting proteins and random non-interacting pairs
#' 
#' Plot the distribution of the Pearson correlation coefficient for the 
#' coexpression of interacting proteins, and a random sample of non-interacting
#' protein pairs of equal size. 
#' 
#' @param enrichment A data table containing the Pearson correlation 
#' coefficients of all interacting and non-interacting protein pairs, as output 
#' by \code{\link{CalculateCoexpressionEnrichment}}
#'
#' @return A \code{ggplot2} plot
#' 
#' @examples
#' interactions <- read.csv("interactions.csv")
#' expression <- ReadExpression("expression.csv")
#' enrichment <- CalculateCoexpressionEnrichment(interactions, expression)
#' plot <- PlotRandomCoexpression(enrichment)
#' plot(plot)
#'
#' @export 
PlotRandomCoexpression <- function(enrichment) {
  # Check input data 
  if (!("FlavinCoexpressionEnrichment" %in% class(enrichment)))
    stop("Must provide a data table created by CalculateCoexpressionEnrichment")
  
  # Get interacting proteins
  interacting <- subset(enrichment, interact == T)
  
  # Get a random sample of non-interacting proteins
  not.interacting <- subset(enrichment, interact == F)
  sample <- not.interacting[sample(1:nrow(not.interacting), nrow(interacting)),]
  
  # Merge 
  enrichment <- rbind(interacting, not.interacting)
  
  # Create plot 
  plot <- ggplot(enrichment, aes(x = pcc, colour = interact, fill = interact)) + 
    geom_density(alpha = 0.5) + 
    labs(x = "Tissue proteome abundance correlation (PCC)", y = "Density") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values = flavin_palette, name = "Interacting") +
    scale_colour_manual(values = flavin_palette, name = "Interacting") +
    flavin_theme
  return(plot)
}

# ==============================================================================
#' Calculate the correlation between coexpression and interaction score
#' 
#' Many software packages designed to predict protein-protein interactions 
#' assign a score to predicted interactions. This method calculates the 
#' correlation (Spearman's rho) between interaction score and coexpression as a
#' means of validating the score assigned to the interactions (under the 
#' assumption that higher-confidence interactions should be more likely to be 
#' coexpressed). 
#' 
#' @param coexpression A data frame output by 
#' \code{\link{CalculateCoexpression}}, with an additional column "score" noting
#' the score assigned to each interaction
#' 
#' @return none
#' 
#' @examples 
#' interactions <- read.csv("interactions.csv")
#' expression <- ReadExpression("expression.csv")
#' coexpression <- CalculateCoexpression(interactions, expression)
#' CalculateCoexpressionScoreCorrelation(coexpression)
#' 
#' @export
CalculateCoexpressionScoreCorrelation <- function(coexpression) {
  # Check input 
  if (!("FlavinCoexpression" %in% class(coexpression)))
    stop("Must provide a data table created by CalculateCoexpression")
  
  # Calculate correlation between interaction score and PCC
  cor <- suppressWarnings(cor.test(coexpression$score, coexpression$pcc, 
                                   method="spearman"))
  rho <- formatC(signif(cor$estimate, digits = 2), digits = 2,
                 format = "g", flag = "#")
  p <- formatC(signif(cor$p.value, digits = 2), digits = 2, 
               format = "g", flag = "#")
  message("Correlation between interaction score and coexpression ", 
          "(Pearson correlation coefficient): Spearman's rho = ", rho, ", P = ",
          p)
}

# ==============================================================================
#' Plot the correlation between coexpression and binned interaction score
#' 
#' Many software packages designed to predict protein-protein interactions 
#' assign a score to predicted interactions. This method creates a plot of
#' coexpression (Pearson correlation coefficient) vs. binned interaction score,
#' to visually assess the relationsihp between interaction score and 
#' coexpression. 
#' 
#' @param coexpression A data frame output by 
#' \code{\link{CalculateCoexpression}}, with an additional column "score" noting
#' the score assigned to each interaction
#' @param bins Then umber of bins to use on the X-axis 
#' 
#' @return A \code{ggplot2} plot
#' 
#' @examples 
#' interactions <- read.csv("interactions.csv")
#' expression <- ReadExpression("expression.csv")
#' coexpression <- CalculateCoexpression(interactions, expression)
#' plot <- PlotCoexpressionScoreCorrelation(coexpression)
#' plot(plot)
#' 
#' @export
PlotCoexpressionScoreCorrelation <- function(coexpression, bins = 10) {
  # Check input 
  if (!("FlavinCoexpression" %in% class(coexpression)))
    stop("Must provide a data table created by CalculateCoexpression")
  
  # Create plot 
  coexpression$score.cut <- cut(coexpression$score, bins)
  summary <- ddply(coexpression, .(score.cut), summarize, mean = mean(pcc), 
                   sd = sd(pcc), sem = sd(pcc)/sqrt(length(pcc)))
  plot <- ggplot(summary, aes(x = score.cut, y = mean)) +
    geom_line(group = 1) + 
    geom_point() + geom_linerange(aes(ymin = mean - sem, ymax = mean + sem)) +
    labs(x = "Interaction score", 
         y = "Coexpression (Pearson correlation coefficient)") + 
    flavin_theme
  coexpression <- subset(coexpression, select = -c(score.cut))
  return(plot)
}

# ==============================================================================
#' Calculate the coexpression of interacting proteins. 
#' 
#' Calculate the Pearson correlation coefficient for the coexpression of all
#' interacting protein pairs in an interaction dataset. 
#' 
#' @param interactions A list of interactions, where the first two columns
#' correspond to the accessions of the interacting proteins. 
#' @param expression A dataset of protein expression, used to calculate 
#' coexpression (see \code{\link{ReadExpression}})
#' 
#' @return The input table of interactions, with the PCC for each 
#' interacting pair in a new column
#' 
#' @examples 
#' interactions <- read.csv("interactions.csv")
#' expression <- ReadExpression("expression.csv")
#' coexpression <- CalculateCoexpression(interactions, expression)
#' 
#' @export 
CalculateCoexpression <- function(interactions, expression) {
  # Calculate Pearson correlation coefficient for interacting protein pairs
  pcc.df <- interactions
  pcc.df$pcc <- unlist(apply(interactions, 1, function(row) 
    CalculateCoexpressionPcc(row[1], row[2], expression)))
  
  # Ignore interactions for which PCC could not be calculated 
  message("Calculated coexpression (PCC) for ", sum(!is.na(pcc.df$pcc)), 
          " of ", nrow(pcc.df)," interactions")
  pcc.df <- subset(pcc.df, !is.na(pcc)) 
  
  # Return coexpression data 
  class(pcc.df) <- append(class(pcc.df), "FlavinCoexpression")
  return(pcc.df)
}

# ==============================================================================
#' Calculate protein expression Pearson correlation coefficient of two proteins
#' 
#' Calculates the Pearson correlation coefficient for the coexpression of two 
#' proteins based on an input dataset of gene or protein expression.
#' 
#' @param protein.A The accession of the first interacting protein 
#' @param protein.B The accession of the second interacting protein
#' @param expression A dataset of protein expression, used to calculate 
#' coexpression (see \code{\link{ReadExpression}})
#' 
#' @return the Pearson correlation coefficient for the coexpression of the two
#' proteins. 
CalculateCoexpressionPcc <- function(protein.A, protein.B, expression) {
  idx.A <- which(expression[,1] == protein.A)
  idx.B <- which(expression[,1] == protein.B)
  # Ignore absent mappings 
  if (length(idx.A) == 0 || length(idx.A) > 1 ||
      length(idx.B) == 0 || length(idx.B) > 1) return(NA)
  # Expression data is in all columns after 1
  expr.A <- as.numeric(expression[idx.A,2:ncol(expression)])
  expr.B <- as.numeric(expression[idx.B,2:ncol(expression)])
  cor <- cor.test(expr.A, expr.B)
  return(cor$estimate)
}

# ==============================================================================
#' Read gene or protein expression from a CSV file
#' 
#' Read gene or protein expression data into \code{flavin} from a CSV file, 
#' where all expression data points for each gene or protein are on the same 
#' row. The first column of the CSV file must be the identifier of the protein.
#' 
#' @param filepath The location of the CSV file containing expression data 
#' 
#' @return The expression data
#' 
#' @export
ReadExpression <- function(filepath) {
  expression <- read.csv(filepath)
  if (ncol(expression) < 4)
    stop("Must have at least three expression points per gene or protein to ", 
         "calculate Pearson correlation coefficient")
  return(expression)
}
