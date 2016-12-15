devtools::use_package("data.table")
devtools::use_package("ggplot2")

# Set default ggplot theme
flavin_theme <- theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        panel.grid=element_blank(),
        panel.border = element_blank(), 
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"), 
        axis.ticks = element_line(colour="grey50"),
        legend.position = "bottom")

# Set colours
flavin_palette <- c("#f58d62", "#8f6ca9", "#5ec5e2", "#1482b2")