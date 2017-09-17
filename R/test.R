#data(sim_dat)
#mm=ClussCluster(sim_dat, nclust=4, ws=c(1.3,2,3,4,6,8,10), verbose = F)
#mmm=ClussCluster(sim_dat, centers =sim_dat[,c(10, 20, 35, 55)], ws=c(1.3,2,3,4,6,8,10), verbose = F)
# ClussCluster_Gap(sim_dat, nclust=4, ws=c(1.3,2,3,4,6,8,10), verbose = F)

#Summary_ClussCluster(mm)
#print_ClussCluster(mm)
#plot_ClussCluster(mm)
#plot_ClussCluster(mm[[4]])

# data(Hou)
# hou.dat <-Hou$x
# run.ft <- filter_gene(hou.dat)
# dat.ft <- run.ft$dat.ft
# run.gap <- ClussCluster_Gap(dat.ft)
# s <- run.gap$best
# hou.test = ClussCluster(dat.ft, nclust=3, ws=s, verbose = F)
# Summary_ClussCluster(hou.test[[1]])
# plot_ClussCluster(hou.test[[1]])
# plot_ClussCluster(hou.test)
