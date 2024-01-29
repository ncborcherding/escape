
#' @importFrom stats dist hclust

heatmapEnrichment <- function(enriched, 
                              assay = NULL,
                              group.by =NULL, 
                              gene.set.use = "all", 
                              order.by = NULL,
                              cluster.rows = FALSE,
                              cluster.columns = FALSE,
                              scale = FALSE, 
                              facet.by = NULL, 
                              palette = "inferno") {
  
  
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  enriched <- .prepData(input.data, assay, gene.set.use, group.by, NULL, facet.by) 
  
  if(length(gene.set) == 1 && gene.set == "all") {
    gene.set <- gene.set[gene.set %!in% c(group.by, facet.by)]
  }
  
  enriched.summary <- enriched %>%
                        group_by(enriched[,group.by], enriched[,facet.by]) %>%
                        summarise(across(which(colnames(enriched) %in% gene.set), mean)) %>%
                        as.data.frame()
  colnames(enriched.summary)[1] <- group.by
                          
  
  if(scale) {
    enriched.summary[,gene.set] <- apply(enriched.summary[,gene.set], 2, scale)
  }
  
  reformated.enriched <-  suppressMessages(melt(enriched.summary))
  
  if(cluster.rows) {
    row.order <- gene.set[hclust(dist(t(enriched.summary[,gene.set]), method = "euclidean"), method = "ward.D2")$order]
    reformated.enriched[,2] <- factor(reformated.enriched[,2], levels = row.order)
  }
  
  if(cluster.columns) {
    column.order <- enriched.summary[,1][hclust(dist(enriched.summary[,gene.set], method = "euclidean"), method = "ward.D2")$order]
    reformated.enriched[,1] <- factor(reformated.enriched[,1], levels = as.vector(column.order))
  }
  
  
  #TODO Add facet.by functionality

    plot <- ggplot(reformated.enriched,
                   mapping = aes(x = reformated.enriched[,group.by],
                                 y = reformated.enriched[,2],
                                 fill = reformated.enriched[,3])) +
                  geom_tile(color = "black", linewidth = 0.5) +
                  scale_y_discrete(expand = c(0, 0)) +
                  scale_x_discrete(expand = c(0, 0)) +
                  labs(fill = "Enrichment Scoure") +
                  guides(fill = guide_colorbar(title.position = "top",
                                                title.hjust = 0.5)) + 
                  coord_equal() + 
                  scale_fill_gradientn(colors = .colorizer(palette, 11)) +
                  theme_classic() + 
                  theme(axis.title = element_blank(),
                        axis.ticks = element_blank(), 
                        legend.direction = "horizontal", 
                        legend.position = "bottom")
    
    if (!is.null(facet.by)) {
      plot <- plot + 
                facet_grid(as.formula(paste('. ~', facet.by))) 
    }
    return(plot)
  
}