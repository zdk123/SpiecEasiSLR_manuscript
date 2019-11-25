nets <- load(list.files(pattern=sprintf('%sNetFits.RData', basename(getwd()))))

nets_li <- lapply(mget(nets), getOpt, sign=TRUE)
edge_dli   <- proxy::dist(nets_li, method=edge_dist2)
degree_dli <- proxy::dist(nets_li, method=degree_cordist)
btw_dli    <- proxy::dist(nets_li, method=btw_cordist)
edge_dense  <- sapply(nets_li, edge_density)
edge_ppos   <- sapply(nets_li, proportion_pos)

save(edge_dli, edge_dense, edge_ppos, degree_dli, btw_dli, file='net_dist.RData')
