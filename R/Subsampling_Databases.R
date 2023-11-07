# Random Subsets ------------------------------------------------
print("Subsampling depthRegul")

tempN <- depthRegul$getActiveSampleNumber()

depthRegul1000 <- depthRegul$clone()
depthRegul1000$addSelectionByRanks(seq(from=1,to=tempN,by=floor(tempN/1000)))
depthRegul1000 = Db_createReduce(depthRegul1000)
depthRegul1000$deleteColumns('NewSel')

depthRegul5000 <- depthRegul$clone()
depthRegul5000$addSelectionByRanks(seq(from=1,to=tempN,by=floor(tempN/5000)))
depthRegul5000 = Db_createReduce(depthRegul5000)
depthRegul5000$deleteColumns('NewSel')

depthRegul10000 <- depthRegul$clone()
depthRegul10000$addSelectionByRanks(seq(from=1,to=tempN,by=floor(tempN/10000)))
depthRegul10000 = Db_createReduce(depthRegul10000)
depthRegul10000$deleteColumns('NewSel')

