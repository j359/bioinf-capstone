
tfinalcopynumber_filtered = t(finalcopynumber_filtered)
tfinaldnameth_filtered = t(finaldnameth_filtered)

dmeth <- rownames(tfinaldnameth_filtered)
copyn <- rownames(tfinalcopynumber_filtered)
finalgenes <- combined5

dmeth2 <- lapply(dmeth, substr, 16, 50)

fcommongenes <- lapply(commongenes, substr, 6,30)

commongenes <- intersect(finalgenes, rownames(demrna))
tcommongenes <- t(commongenes)

write.csv(tcommongenes , file = ("/Users/tesssenftle/Downloads/commondegenes.csv"))

write.csv(cluster_outcome , file = ("/Users/tesssenftle/Downloads/clusteroutcome.csv"))

commongenes2 <- intersect(rownames(demrna), dmeth2)

commongenes <- read.csv(file = "/Users/tesssenftle/Downloads/degenes.csv")

commongenes <- commongenes[0:2, 0]
tcommongenes <- t(commongenes)





finalmrna_filtered$rn <-rownames(finalmrna_filtered)
finalcopynumber_filtered$rn <-rownames(finalcopynumber_filtered)
finaldnameth_filtered$rn <-rownames(finaldnameth_filtered)





require(plyr)
finalcombined <- join_all(list(finalmrna_filtered,finalcopynumber_filtered,finaldnameth_filtered), by = 'rn', type = 'full')
finalcombined <- t(finalcombined)

write.csv( finalcombined, "/Users/tesssenftle/Downloads/combined.csv")
