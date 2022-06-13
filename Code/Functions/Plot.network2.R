#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#


Plot.network.sim.empire2 <- function(ruler, sim.thresh.plot = 1, net.stats, sim.thresh.stats, weight = TRUE, top.reg, clust = FALSE){
  # Set sim.thresh default
  sim.thresh = sim.thresh.plot
  
  # Selet presence absence matrix for specified empire
  if(ruler == "Great_Britain") {COL = "#4777EFFF"} # former#0066CC
  if(ruler == "Spain") {COL = "#DB3A07FF"} # former #CC3300
  if(ruler == "Portugal") {COL = "#30123BFF"} # former "#660066"
  if(ruler == "Netherlands") {COL = "#FE9B2DFF"} # former #FF9933
  
  if(ruler == "Great_Britain") {PA <- pa.GBR}
  if(ruler == "Spain") {PA <- pa.ESP}
  if(ruler == "Portugal") {PA <- pa.PRT}
  if(ruler == "Netherlands") {PA <- pa.NED}
  
  
  # Calculate beta diversity metrics
  beta.div <- beta.pair(PA, "sorensen") # Nestedness component
  ##### beta.div <- beta.pair(PA, "jaccard") # Turnover component
  
  # Create adjacency matrix of dissimilarities using Simpson dissimilarity
  beta.div.simp <- as.matrix(beta.div$beta.sor)
  ##### beta.div.simp <- as.matrix(beta.div$beta.jtu) # Turnover component
  
  # Create long format list of pairwise distances
  beta.div.simp.long <- t(combn(colnames(beta.div.simp), 2))
  beta.div.simp.long <- data.frame(beta.div.simp.long, dist=beta.div.simp[beta.div.simp.long])
  
  # Change dissimilarity to similarity (1-dissimilarity)
  beta.div.simp.long$dist <- 1-beta.div.simp.long$dist
  
  # Select pairwise dissimilarities below a specific threshold
  beta.div.simp.long2 <- beta.div.simp.long[beta.div.simp.long$dist >= sim.thresh.plot,]
  
  
  
  # Create nodes file with geographic coordinates
  nodes <- pred %>%
    filter(OBJIDsic %in% unique(c(beta.div.simp.long2$X1,beta.div.simp.long2$X2))) %>%
    select(OBJIDsic, LON, LAT, name)
  nodes$OBJIDsic <- as.factor(nodes$OBJIDsic)
  

  if(clust == T){  
  # subset of most central nodes based on clusters
  nodes_top_cent <- net.stats %>%
    group_by(cluster.mod) %>%
    top_n(top.reg, cent.eig) %>%
    arrange(cluster.mod) %>%
  left_join(nodes %>% select(OBJIDsic, LON, LAT, name), by = "OBJIDsic")
  }
  
  if(clust == F){
    top.reg = 10
    # subset most central nodes of entire dataset
    nodes_top_cent <- net.stats %>%
      top_n(top.reg, cent.eig) %>%
      left_join(nodes %>% select(OBJIDsic, LON, LAT, name), by = "OBJIDsic")
  }
  
  
  # Create edges file with links and weights
  edges <- beta.div.simp.long2
  colnames(edges) <- c("from", "to", "weight")
  edges$from <- as.factor(edges$from)
  edges$to <- as.factor(edges$to)
  
  # Build graph element
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  edges_for_plot <- edges %>%
    inner_join(nodes %>% 
                 select(OBJIDsic, LON, LAT), by = c('from' = 'OBJIDsic')) %>%
    rename(x = LON, y = LAT) %>%
    inner_join(nodes %>% select(OBJIDsic, LON, LAT), by = c('to' = 'OBJIDsic')) %>%
    rename(xend = LON, yend = LAT)
  assert_that(nrow(edges_for_plot) == nrow(edges))
  
  
  #============================================================================================================#
  
  
  nodes <- nodes%>%
    left_join(net.stats, by = "OBJIDsic")
  
  #============================================================================================================#
  if(weight == TRUE) {
  nodes$weight = nodes$cent.eig
  } else(nodes$weight = 1)
  
  
  # Build plot theme
  maptheme <- theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(legend.position = "bottom") +
    theme(panel.grid = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))
  
  # Get worldmap shape
  country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
                                 data = map_data('world'),
                                 fill = "#CECECE", color = "#515151",
                                 size = 0.15)
  mapcoords <- coord_fixed(xlim = c(-150, 180), ylim = c(-55, 80))
  
  
  #=========================#
  # HACKY PLOT
  
  if(ruler == "Great_Britain") {title.emp <- paste ("British empire")}
  if(ruler == "Spain") {title.emp <- paste("Spanish empire")}
  if(ruler == "Portugal") {title.emp <- paste("Portuguese empire")}
  if(ruler == "Netherlands") {title.emp <- paste("Dutch empire")}
  
  theme_transp_overlay <- theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )
  
  # the base plot showing only the world map
  p_base <- ggplot() + country_shapes + mapcoords + maptheme+
    #ggtitle(paste0(ruler," - similarity >", sim.thresh.stats, " - links shown >", sim.thresh.plot)) ## Title that givs power, statistics threshold and threshold for plotting
    ggtitle(title.emp)
  
  
  
  # first overlay: edges as arcs
  p_edges <- ggplot(edges_for_plot) +
    geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     # draw edges as arcs
                   size = weight), colour = COL, 
               curvature = 0.33, alpha = 0.1) +
    scale_size_continuous(guide = FALSE, range = c(0.5, 2)) +  # scale for edge widths
    mapcoords + maptheme + theme_transp_overlay +
    theme(legend.position = "none")
  
  
  # second overlay: nodes as points
  colo <- nodes$cluster.mod
  colo <- colo %>% 
    recode("1" = viridis_pal(option = "turbo")(8)[3]) %>% # yellow "#edae29"
    recode("2" = viridis_pal(option = "turbo")(8)[4]) %>% # green "#66a182"
    recode("3" = viridis_pal(option = "turbo")(8)[5]) %>% # navy blue "#2e4057"
    recode("4" = viridis_pal(option = "turbo")(8)[8]) # grey "#8d96a3"
  
  p_nodes <- ggplot(nodes) +
    geom_point(aes(x = LON, y = LAT, size = weight),
               shape = 21, fill = colo, color = "black",    # draw nodes
               stroke = 0.5) +
    scale_size_continuous(guide = FALSE, range = c(1, 4)) +    # scale for node size
    mapcoords + maptheme + theme_transp_overlay
  
  # third overlay: node names of the top central nodes
  
  colo2 <- nodes_top_cent$cluster.mod
  colo2 <- colo2 %>% 
    recode("1" = viridis_pal(option = "turbo")(9)[3]) %>% # yellow "#edae29"
    recode("2" = viridis_pal(option = "turbo")(9)[4]) %>% # green "#66a182"
    recode("3" = viridis_pal(option = "turbo")(9)[8]) %>% # navy blue "#2e4057"
    recode("4" = viridis_pal(option = "turbo")(9)[9]) # grey "#8d96a3"
  
  p_cent_nodes <- ggplot(nodes_top_cent) +
    geom_label_repel(aes(x = LON, y = LAT, label = Region_unique),             # draw text labels
            hjust = 0, nudge_x = 1, nudge_y = 4,
            size = 5, color = colo2, fontface = "bold", max.overlaps = 100, segment.color = "black", segment.size = 1, min.segment.length = 0) +
    mapcoords + maptheme + theme_transp_overlay
  
  # combine the overlays to a full plot
  
  # proper positioning of the grobs can be tedious... I found that
  # using `ymin` works quite well but manual tweeking of the
  # parameter seems necessary
  
  p <- p_base +
    annotation_custom(ggplotGrob(p_edges), ymin = -74) +
    annotation_custom(ggplotGrob(p_nodes), ymin = -74) +
    annotation_custom(ggplotGrob(p_cent_nodes), ymin = -74)
  
}