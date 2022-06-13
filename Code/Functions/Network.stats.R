#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#

# Identify groups of regions in network ----

Network.stats <- function(ruler, sim.thresh.stats = 1){
  
  # Selet presence absence matrix for specified empire
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
  
  # Select pairwise similarities above a specific threshold
  beta.div.simp.long2 <- beta.div.simp.long[beta.div.simp.long$dist >= sim.thresh.stats,]
  
  # Build graph object
  graph.emp <- as_tbl_graph(beta.div.simp.long2, directed = F, edges = beta.div.simp.long2$dist)
  
  # Calculate centrality and cluster measures
  
  ### Definitions Centrality
  
  #### Authority & hub centrality: Hub and authority centarlities are generalization of eigenvector centrality. A high hub node points to many good authorities and a high authority node receives from many good hubs.
  
  #### Betweenness centrality: The betweenness centrality for each nodes is the number of the shortest paths that pass through the nodes.
  
  #### Closeness centrality: Closeness centrality measures how many steps are required to access every other nodes from a given nodes. It describes the distance of a node to all other nodes. The more central a node is, the closer it is to all other nodes.
  
  #### Eigenvector centrality: A node is important if it is linked to by other important nodes. The centrality of each node is proportional to the sum of the centralities of those nodes to which it is connected. In general, nodes with high eigenvector centralities are those which are linked to many other nodes which are, in turn, connected to many others (and so on).
  
  
  
  ### Definition Clusters
  
  #### Infomap community finding. It groups nodes by minimizing the expected description length of a random walker trajectory
  
  #### Community structure detection based on edge betweenness. It groups densely connected nodes
  
  graph.cent.clust <- graph.emp %>%
    activate(nodes) %>%
    mutate(cent.auth = centrality_authority(),
           #cent.betw = centrality_betweenness(weights = dist),
           #cent.close = centrality_closeness(weights = dist),
           #cent.hub = centrality_hub(weights = dist),
           cent.eig = centrality_eigen(weights = dist),
           #cluster.comp = as.factor(group_components()), # Group by connected compenents 
           #cluster.e.betw = as.factor(group_edge_betweenness()), # Group densely connected nodes
           cluster.mod = as.factor(group_fast_greedy(weights = dist)) # Group nodes by optimising modularity - https://journals-aps-org.uaccess.univie.ac.at/pre/pdf/10.1103/PhysRevE.70.066111
           #cluster.info = as.factor(group_infomap()), # Group nodes by minimizing description length
           #cluster.eigen = as.factor(group_leading_eigen()) # Group nodes based on the leading eigenvector of the modularity matrix
           #cluster.mod.opt = as.factor(group_optimal()) # Group nodes by optimising the moldularity score
    )
  
  graph.cent.clust
  
  # Get network measure data
  graph.clust.data <- data.frame(
    OBJIDsic = as.integer(as.character(V(graph.cent.clust)$name)),
    cent.auth = V(graph.cent.clust)$cent.auth,
    #cent.betw = V(graph.cent.clust)$cent.betw,
    #cent.close = V(graph.cent.clust)$ cent.close,
    #cent.hub = V(graph.cent.clust)$cent.hub,
    cent.eig = V(graph.cent.clust)$cent.eig,
    #cluster.comp = V(graph.cent.clust)$cluster.comp, # Group by connected compenents 
    #cluster.e.betw = V(graph.cent.clust)$cluster.e.betw, # Group densely connected nodes
    cluster.mod = V(graph.cent.clust)$cluster.mod # Group nodes by optimising modularity 
    #cluster.info = V(graph.cent.clust)$cluster.info, # Group nodes by minimizing description length
    #cluster.eigen = V(graph.cent.clust)$cluster.eigen # Group nodes based on the leading eigenvector of the modularity matrix
    #cluster.mod.opt = V(graph.cent.clust)$cluster.mod.opt # Group nodes by optimising the moldularity score
  )
  
  graph.clust.data <- graph.clust.data %>%
    left_join(regions[,c("OBJIDsic", "Region_unique", "UN_geospat")], by = "OBJIDsic")
  
  graph.clust.data$OBJIDsic <- as.factor(graph.clust.data$OBJIDsic)
  
  return(graph.clust.data)
}
