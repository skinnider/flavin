# ==============================================================================
# Calculating enrichment in protein co-expression 
# ==============================================================================
# Read in some protein-protein interactions
interactions <- read.csv("interactions.csv") # from flavin-utils 
# Read formatted expression data 
expression <- ReadExpression("expression.csv") # in flavin-utils e.g. https://github.com/skinnider/pcp-pipeline-analysis/blob/master/data/proteome/HPM_protein_level_expression_matrix_Kim_et_al_052914.csv

# Calculate coexpression for all interacting proteins
coexpression <- CalculateCoexpression(interactions, expression)
# Calculate enrichment over non-interacting proteins 
enrichment <- CalculateCoexpressionEnrichment(coexpression, interactions, 
                                              expression)

# Plot enrichment over non-interacting proteins 
plot <- PlotCoexpressionEnrichment(enrichment)
plot(plot)

# Plot the PCC distribution for interacting proteins and a random sample of
# non-interacting proteins 
plot <- PlotRandomCoexpression(enrichment)
plot(plot)

# Calculate correlation between coexpression and interaction score 
CalculateCoexpressionScoreCorrelation(coexpression)
# Now, plot the relationship 
plot <- PlotCoexpressionScoreCorrelation(coexpression, bins = 10)
plot(plot)

# ==============================================================================
# Calculating enrichment in subcellular localization 
# ==============================================================================

# Read in some protein-protein interactions
interactions <- read.csv("interactions.csv")
# Read formatted subcellular localization data 
localization <- read.csv("subcellular_localization.csv")

# Calculate shared subcellular localization for all interacting proteins 
shared.loc <- CalculateCoexpression(interactions, expression)

