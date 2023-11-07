myDf <- read.table("Data/FAROFA123_depth_cleaned_regul",header=T)

# creation DB
myDb <- Db_create()
myDb$addColumns(myDf$longitude,'longitude',ELoc_X(),0)
myDb$addColumns(myDf$latitude,'latitude',ELoc_X(),1)
myDb$addColumns(myDf$depth,'depth',ELoc_Z(),0)
myDb$addColumns(myDf$farofa,'farofa')
myDb$addColumns(myDf$distAcross,'distAcross')
myDb$addColumns(myDf$distAlong,'distAlong')
myDb$addColumns(myDf$dirAcross,'dirAcross')
myDb$addColumns(myDf$trapezeNb,'trapezeNb')

#ggplot()+plot(myDb,name_color="depth")
