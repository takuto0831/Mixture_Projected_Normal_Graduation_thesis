library(igraph)
g1 <- graph.formula(X-+Y,X-+Z)
plot(g1,layout = layout.circle,
     edge.arrow.size = 1.5,
     edge.arrow.width = 1.5,
     edge.width = 4,
     vertex.size = 25)
g2 <- graph.formula(X+-Y,X-+Z)
plot(g2,layout = layout.circle,
     edge.arrow.size = 1.5,
     edge.arrow.width = 1.5,
     edge.width = 4,
     vertex.size = 25)
g3 <- graph.formula(X-+Y,X+-Z)
plot(g3,layout = layout.circle,
     edge.arrow.size = 1.5,
     edge.arrow.width = 1.5,
     edge.width = 4,
     vertex.size = 25)

