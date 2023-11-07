# Reading 1 ping data ------------------------------------------------------------
print("Reading data")
df <- read.csv(file = 'Data/FAROFA123_logSa_1ping.csv')
df = na.omit(df)

# Cleaning dataset : removing 5 isolated points
sel <- df$Longitude > -32.41892 & df$Longitude < -32.41783 & df$Latitude > -3.799620 & df$Latitude < -3.798449
df <- df[!sel,]
names(df) <- c("longitude","latitude","sa70","sa200","depth","time","esdu","farofa")

# Cleaning truncated data due to echosounder's threshold 
tempZ <- unique(df$depth)
tempFreq <- sapply(1:length(tempZ), function(i) sum(df$depth==tempZ[i]))
summary(tempFreq) # More than 50% of the observed values are observed only once
plot(tempZ,tempFreq)
tempZSel <- (tempFreq > 2000) | 
  (tempZ > 90 & tempZ < 110 & tempFreq > 100) | 
  (tempZ > 140 & tempZ < 160 & tempFreq > 300) | 
  (tempZ > 290 & tempZ < 310 & tempFreq > 500) |
  tempZ > 390

tempZKept <- tempZ[!tempZSel]
tempFreqKept <- tempFreq[!tempZSel]
plot(tempZKept,tempFreqKept)

tempZToBeRemoved <- tempZ[tempZSel]

sel <- df$depth %in% tempZToBeRemoved
df2 <- df[!sel,]

write.table(df2,"Data/FAROFA123_logSa_1ping_cleaned",row.names=F)