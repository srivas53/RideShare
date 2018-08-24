# Anthony Niemiec 2018, uSMART research group (Purdue University, Dr. Cai)
# NOTE: this script allows a user to upload a spatial data set (.csv) and a
#	street network (in the form of an ESRI shapefile). It can be used to
#	harvest attributes for all streets in a trip over a list of pickup/
#	dropoff coordinates; for example, this script can take in two street
#	intersections, calculate a trip between them, and return a list of
#	all speed limits or street lengths encountered during that trip. This
#	can be useful in conjunction with other data such as measured trip
#	times in order to create estimations for trip travel times based
#	off of predictor variables over multiple road trips in a network.
# This script deals with an ESRI shapefile with its coordinates stored in
#	UTM instead of lat-long pairs, so there is a portion of the script
#	that deals with that; for most uses, THIS WILL BE COMMENTED OUT if
#	your ESRI shapefile already stores its coordinate pickup and dropoff
#	points in lat-long format.
# You will likely need to modify something where I have placed a comment
#	with two stars ** around it, i.e. to show which column of your
#	shapefile's attributes table holds certain values, or if you wish to
#	parse more values from your shapefile's attribute table than just
#	speed limits and road lengths in a given trip between two nodes.

# to download the necessary packages (once per computer)
install.packages("rgdal") # reading in the shapefile
install.packages("shp2graph") # converts shapefile to igraph-compatible object
install.packages("igraph") # comprehensible graphing and spatial functions, i.e. get.shortest.paths()
install.packages("RANN") # the k-nearest-neighbors function, nn2()
# to link the necessary packages once downloaded (every time you run this)
library(rgdal)
library(shp2graph)
library(igraph)
library(RANN)
library(dplyr)

# reading in the .csv
DataSet <- read.csv("yellow Sample.csv", header=TRUE, sep=",") # **your 'sep' may be different depending on your .csv**
# my dataset: yellow_2014_noAxisHug.csv, a cleaned data set from
#	the NYC taxi bureau's website; covers yellow taxi trips from the month
#	of June, 2014

# reading in the ESRI shapefile and creating useful objects (specifically a shapefile, other network types will most likely not work due to the shapefile-specific functions called later on)
#setwd("C:/Users/Niemiec/Desktop/NYC_Manhattan") # change to be the directory that holds all your shapefile's components (not just the .shp)
MAPshp <- readOGR(dsn=".", layer="TSquare") # "layer" will be the name of your shapefile, i.e. whatever precedes the file extension '.shp'
MAPraw <- readshpnw(MAPshp, ELComputed=FALSE) # ELComputed and other parameters may change; see documentation for shp2graph package, basically used to extract edgelist and nodelist
# elements in MAPraw: (will be used to extract data once routes are computed based off of its igraph rendition)
# index 1: boolean for "Detailed" (if every intersection of two lines/edges should be turned into a new node if it isn't one already)
# index 2: "nodelist" object
# index 3: "edgelist" object
# index 4: lengths of edges if "ELComputed" argument is set to true (defaulted to false)
# index 5: original attribute table from the shapefile (in my NYC shapefile, MAPraw[[5]]$SPEED denotes posted speed limits and MAPraw[[5]]$Shape_Leng denotes street lengths)
# index 6: "nodeXlist," contains all node X-coords
# index 7: "nodeYlist," contains all node Y-coords

costMatrix <- MAPraw[[5]]$FromToCost # ** if your shapefile has a list of road "costs," i.e. a proper edge weight of traveling down a street, you may put the column name here **
n = which(costMatrix==-1) # ** in my data set, one-way streets were given an edge weight of -1 when going down the wrong way; your data set may handle this differently, or not at all **
costMatrix[n] <- MAPraw[[5]]$ToFromCost[n] # ** if your data set does not have such a column, you may improvise by dividing your column of street lengths by your column of speed limits; this way, you get a general cost of average street times **
MAPgraph <- nel2igraph(MAPraw[[2]], MAPraw[[3]], weight=costMatrix) # create an igraph object; this format is used to compute shortest paths between nodes
#plot.igraph(MAPgraph) # **visualize the loaded-in nodes and edges to make sure nothing is wrong with the input data**

# parsing DataSet coordinates, to create a startMatrix and endMatrix of trips
DataXo <- DataSet$pickup_longitude # change the $ column to the column which stores origin XY and destination XY data
DataYo <- DataSet$pickup_latitude # pickup latitude
DataXd <- DataSet$dropoff_longitude # dropoff longitude
DataYd <- DataSet$dropoff_latitude # dropoff latitude
startMatrix = matrix(c(DataXo,DataYo), nrow=length(DataXo), ncol=2) # matrix of all pickup coordinates
endMatrix = matrix(c(DataXd, DataYd), nrow=length(DataXo), ncol=2) # matrix of all dropoff coordinates

MAPxCoords <- matrix(nrow=length(MAPraw[[2]][,1]), ncol=1) # empty matrix x coord for street ("node") intersections
MAPyCoords <- matrix(nrow=length(MAPraw[[2]][,1]), ncol=1) # empty matrix y coord for street ("node") intersections
CoordRaw <- MAPraw[[2]][,2]
for (i in 1:(length(MAPxCoords))) { # gathering coordinates from MAPraw, the loaded-in shapefile
  MAPxCoords[i] <- CoordRaw[[i]][1] # x column of shapefile coordinates
  MAPyCoords[i] <- CoordRaw[[i]][2] # y column of shapefile coordinates
}
streetMatrix = matrix(c(MAPxCoords, MAPyCoords), nrow=length(MAPxCoords), ncol=2) # matrix of all street intersection coordinates

# converting from UTM to lat-long coordinates; **COMMENT THIS OUT IF YOUR SHAPEFILE HAS LAT-LONG COORDINATES**
tempSP <- SpatialPoints(streetMatrix, proj4string=CRS("+proj=utm +zone=18 +datum=WGS84"))
tempGEO <- spTransform(tempSP, CRS("+proj=longlat +datum=WGS84"))
streetMatrix = matrix(c(tempGEO$coords.x1, tempGEO$coords.x2), nrow=length(tempGEO$coords.x1), ncol=2)

# setting up objects to be loaded by the street-parsing algorithm; **you may include whatever you wish to parse from your shapefile attribute table here**
MAPspeeds <- as.numeric(data.matrix(unlist(MAPraw[[5]]$SPEED))) # list of all speeds; **your shapefile's attribute table may have a different column name for this **
MAPlengths <- MAPraw[[5]]$Shape_Leng # list of all street lengths; ** again, you may need to change the name of this to mirror your own shapefile's attribute table **
spdlims <- vector(mode='list', length=1) # column to hold all speed limits over all calculated trips
strlens <- vector(mode='list', length=1) # column to hold all street lengths over all calculated trips
edgeIDs <- vector(mode='list', length=1) # holds edge ID collections, used solely during processing

# estimating the nearest street intersection for each DataSet lat-long coord
tripStartNodeIDs <- nn2(data=streetMatrix, query=startMatrix, k=1, treetype="bd", searchtype="radius", radius=.002, eps=0) # **you may wish to change the radius to specify how far a coordinate pair should search for an intersection before turning up null; this helps filter out coordinate errors/outliers in the process
tripEndNodeIDs <- nn2(data=streetMatrix, query=endMatrix, k=1, treetype="bd", searchtype="radius", radius=.002, eps=0) # remember, while the radius may seem comically small, it is comparing distances between lat-long pairs; a radius of .002 here is approximately 222 meters
DataSet$startID <- c(tripStartNodeIDs$nn.idx) # appends a column of start IDs to data set; to filter out all 0/null values in start and end nodes
DataSet$endID <- c(tripEndNodeIDs$nn.idx)
DataSet <- dplyr::filter(DataSet, startID!=0 & endID!=0) # filtering out any values where the start/end ID didn't have a neighbor within ~200 meters (or whatever you changed the radius to)
tripStartNodeIDs <- DataSet$startID # create individual copies of this column again, so that the algorithm is more streamlined/quicker at retrieval from these columns
tripEndNodeIDs <- DataSet$endID

# the attribute-catching algorithm
attCatch <- function(tS, tE) { # tripStartNodeIDs and tripEndNodeIDs are passed in
  # getting node IDs
  path <- get.shortest.paths(MAPgraph, from=tS, to=tE) # get a list of node IDs on the path
  if (length(path$vpath[[1]]) > 1) { # trips where the origin intersection == destination intersection are discarded as "loop" trips/trips that are too short
    # getting edge IDs (one edge processed at a time each pass of this nested loop)
    x1 <- path$vpath[[1]][1:(length(path$vpath[[1]])-1)]
    x2 <- path$vpath[[1]][2:length(path$vpath[[1]])]
    edgeIDs <- get.edge.ids(MAPgraph, c(rbind(x1,x2)))
    spdlims[[1]] <- MAPspeeds[edgeIDs]
    # harvesting the street lengths
    strlens[[1]] <- MAPlengths[edgeIDs]
    # **this would be a good place to include more attributes
  }
  c(spdlims, strlens,edgeIDs) # returns all speed limit groupings and street length groupings
} ### index of inc isn't changing

# calling the algorithm to iterate over all start and end ID pairs,
#	returning a two-column list of attributes ** you will have an n-column
#	list based on how many attributes you wish to gather (inside the if
#	statement) and return (via c() at the end)
tripAttributes <- mapply(tS = tripStartNodeIDs, tE = tripEndNodeIDs, FUN=attCatch)
# in this example, [1,x] refers to speed limit groupings and [2,x] refers
#	to street length groupings



# write data to files ... ** you may want to use setwd to change where this output will appear, because it will currently output to the folder containing your shapefile (from around the start of the program in order to read the shapefile in) **

# writing the filtered data set out to a file, due to nn2()'s nearest-neighbor search
write.table(DataSet, file="cleanDataSet.csv", sep=" ") # ** you may change the file name and the separator if you wish **

# writing out the captured trip attributes to a file in .csv table format
sink("programOutput.csv") # opens a sink; ** rename to what you want your output to be **

cat("TripNumber", "StreetInTrip", "StreetLength", "SpeedLimit", sep=","); cat('\n') # the header; ** remember to add more names if you extract more attributes than just speed and street lengths **
for (i in 1:(length(tripAttributes[1,])/2)) { # for all groupings in the tripAttributes we just gathered
  if (length(tripAttributes[1,i][[1]]) != 0) { # this occurs with loop trips; to "filter" out, we will simply fail to write them
    j <- 1:length(tripAttributes[1,i][[1]]) # what street number of the trip that we are on (resets for each row)
    cat(rbind( i, ',', j, ',', unlist(tripAttributes[1,i]), ',', unlist(tripAttributes[2,i]), '\n' ), sep="") # print the trip row number and necessary attributes
    # ** you may add more attributes to this list by interting a <',', unlist(tripAttributes[x,i]), > before the closing parentheses rbind() [directly after unlist(tripAttributes[2,i])], where x is the 3rd, 4th, etc. attribute gathered from attCatch **
  }
}

sink() # closes sink

# below is if you wish to generate output that calculates the relative
#	trip length error according to what is recorded in your DataSet;
#	if your data set does not include a recorded trip length, then
#	the latter parts of ErrorReport.txt and the filtering process may
#	not apply to you ** remember to change column names of DataSet to
#	model what the proper column names would be in your data **

# generating relative trip length errors, given the algorithm's approximation and your DataSet's recorded values
relErr <- list() # will hold the relative errors
len <- vector(mode="double", length=length(DataSet$trip_distance)) # hold len's
for (i in 1:(length(tripAttributes)/2)) { # same size as sample of trips
  len[[i]] = sum(tripAttributes[2,i][[1]]) # length recalculated per loop iteration
}
real = DataSet$trip_distance * 1609.34
relErr <- abs((len-real)/real)
relErr <- relErr[relErr!=1] # remove trips completely filtered out (100% off)
m = mean(unlist(relErr))
s = sd(unlist(relErr))
me = median(unlist(relErr))

# removing trips whose relative trip length errors are one standard deviation above the mean relative trip length error
relErrFiltered <- unlist(relErr[relErr < (m+s)])

mFilt <- mean(relErrFiltered)
sFilt <- sd(relErrFiltered)
meFilt <- median(relErrFiltered)

loss <- length(relErr) - length(relErrFiltered)
lossPercent <- loss/length(relErr) * 100

sink("programOutputFiltered.csv") # opens a sink for filtered values; ** you may change this name, as well **

cat("TripNumber", "StreetInTrip", "StreetLength", "SpeedLimit", sep=","); cat('\n') # the header; ** again, remember to add more names if you extract more attributes than just speed and street lengths **
for (i in 1:(length(tripAttributes[1,])/2)) { # for all groupings in the tripAttributes we just gathered
  if (length(tripAttributes[1,i][[1]]) != 0 && relErr[i]<(m+s)) { # will only write out trips with appropriate relErr values ** you may seek to change this if you perform custom filtering **
    j <- 1:length(tripAttributes[1,i][[1]])
    cat(rbind( i, ',', j, ',', unlist(tripAttributes[1,i]), ',', unlist(tripAttributes[2,i]), '\n' ), sep="")
    # ** again, you may add more attributes to this list by interting a <',', unlist(tripAttributes[x,i]), > before the last element of the rbind() statement, where x is the 3rd, 4th, etc. attribute gathered from attCatch **
  }
}

sink()

sink("ErrorReport.txt") # a generally helpful report of how closely the algorithm has approximated what your DataSet's recorded trip lengths were; also presents the potential benefit and cost of removing relative errors that are one standard deviation above the mean relatiev error regarding trip length approximation

cat("Error report for data stored in programOutput.txt:", '\n') # ** change this to the filename of your program output specified above **
cat("Mean relative trip length error, compared to DataSet: ", m*100, "%", '\n', sep="")
cat("Standard deviation of relative errors: ", s*100, "%", '\n', sep="")
cat("Median relative error: ", me*100, "%", '\n', sep="")
cat("---", '\n')
cat("Error report for data stored in programOutputFiltered.txt:", '\n', sep="") # ** change this to be the filename of your program output specified above, plus the word "Filtered" at the end **
cat("Mean relative trip length error, with values above one sd removed: ", mFilt*100, "%", '\n', sep="")
cat("Standard deviation: ", sFilt*100, "%", '\n', sep="")
cat("Median: ", mFilt*100, "%", '\n', sep="")
cat("Number of trips filtered out:", loss, '\n')
cat("% of trips filtered out: ", lossPercent, "%", '\n', sep="")

sink()










# ------------------------------------------------------------------------

# delete all below

len <- vector(mode="double", length=10000)
real = DataSet$trip_distance[1:10000] * 1609.34






ts <- tripStartNodeIDs[1:10000]
te <- tripEndNodeIDs[1:10000]
tA <- mapply(tS = ts, tE = te, FUN=attCatch)


relErr <- list() # will hold the relative errors
len <- vector(mode="double", length=10000) # hold len's
for (i in 1:(length(tA)/2)) { # same size as sample of trips
  len[[i]] = sum(tA[2,i][[1]]) # length recalculated per loop iteration
}
real = DataSet$trip_distance[1:10000] * 1609.34
relErr <- abs((len-real)/real)
m = mean(unlist(relErr))
s = sd(unlist(relErr))
me = median(unlist(relErr))

# removing trips whose relative trip length errors are one standard deviation above the mean relative trip length error
relErrFiltered <- unlist(relErr[relErr < (m+s)])

mFilt <- mean(relErrFiltered)
sFilt <- sd(relErrFiltered)
meFilt <- median(relErrFiltered)

loss <- length(relErr) - length(relErrFiltered)
lossPercent <- loss/length(relErr) * 100

m
mFilt

sum(tA[2,1][[1]])

ts[1:10]
te[1:10]
head(DataSet)

head(streetMatrix)

streetMatrix[ts[1:10],]
streetMatrix[te[1:10],]
DataSet$trip_distance[1:10]*1609.34

sum(tA[2,ts[1]][[1]])
DataSet$trip_distance[1]*1609.34
real[1:10]
len[1:10]

tA[2,1]
DataSet[1]
MAPraw[[2]][1037,2]

ts[1]

sum(tA[2,1][[1]])/1609.34

DataSet$trip_distance[1]*1609.34
sum(tA[2,][[1]])

head(tA[2,])
tA[2,1][[1]]

ts[1:10]
sum(tA[2,ts[1]][[1]])
tA[2,ts[1]][[1]]
tA[2,1037][[1]]

sum(tA[2,1037][[1]])
DataSet$trip_distance[1037]*1609.34

sum(tA[2,6][[1]])
DataSet$trip_distance[6]*2489









#--


sum(t[2,4][[1]])
DataSet$trip_distance[1:10]*1609.34



# the attribute-catching algorithm

attCatch <- function(tS, tE) { # tripStartNodeIDs and tripEndNodeIDs are passed in
  path <- get.shortest.paths(MAPgraph, from=tS, to=tE)
  if (length(path$vpath[[1]]) > 1) {
    x1 <- path$vpath[[1]][1:(length(path$vpath[[1]])-1)]
    x2 <- path$vpath[[1]][2:length(path$vpath[[1]])]
    edgeIDs <- get.edge.ids(MAPgraph, c(rbind(x1,x2)))
    spdlims[[1]] <- MAPspeeds[edgeIDs]
    strlens[[1]] <- MAPlengths[edgeIDs]
  }
  c(spdlims, strlens)
}

max = 10
smallStart <- (tripStartNodeIDs[1:max])
smallEnd <- (tripEndNodeIDs[1:max])

spdlims <- vector(mode='list', length=1) # column to hold all speed limits over all calculated trips
strlens <- vector(mode='list', length=1) # column to hold all street lengths over all calculated trips
edgeIDs <- vector(mode='list', length=1) # holds edge ID collections, used solely during processing, but you may wish to refer to this later for personal use
testC <- vector(mode='list', length=2)
#interval = c(1:max)

time1 <- Sys.time()
t <- mapply(tS = smallStart, tE = smallEnd, FUN=attCatch)
time2 <- Sys.time()
time2-time1

t

head(MAPraw[[5]]$FromToCost, n=20)
head(MAPraw[[5]]$ToFromCost, n=20)
MAPraw[[5]]$FromToCost[MAPraw[[5]]$FromToCost != MAPraw[[5]]$ToFromCost] == -1

MAPraw[[5]]$FromToCost[9892]
MAPraw[[5]]$ToFromCost[9892]


# # # # # # # # # # # # # # # # # # # # # # # # # # # #

########## rel Err equation

relErr <- list() # will hold the relative errors
len <- list() # hold len's
for (i in 1:(length(t)/2)) { # same size as sample of trips
  len[[i]] = sum(t[2,i][[1]]) # length recalculated per loop iteration
  real = DataSet$trip_distance[i] * 1609.34 # converting between miles and meters
  #	print(abs(len[[i]]-real)) # total error
  #	print(abs(len[[i]] - real) / real) # relative error
  #	print("-") # blank space between readings
  relErr[[i]] = abs(len[[i]]-real)/real # store the relative error
}
m = mean(unlist(relErr))
s = sd(unlist(relErr))
m; s
median(unlist(relErr))






t[2,39]



lenFiltered <- list() # holds filtered trips
relErrFiltered <- list() # holds relErr for filtered trips
j = 1
for (i in 1:length(len)) {	# filter out items that are one sd above
  if (relErr[[i]] < (m+s)) {
    #		lenFiltered[[i]] = len[[i]]
    real = DataSet$trip_distance[i] * 1609.34
    relErrFiltered[[j]] = (abs(len[[i]]-real))/real
    j = j + 1
  }
}
mFilt <- mean(unlist(relErrFiltered))
sFilt <- sd(unlist(relErrFiltered))
mFilt; sFilt
median(unlist(relErrFiltered))
length(relErr) - length(relErrFiltered)

length(relErrFiltered)


# # # # # # # # # # # # # # # # # # # # # # # # # # # #



# generating a relative error report

length(costMatrix)
length(MAPraw[[5]]$Shape_Leng)
costMatrix[1:10]
MAPraw[[5]]$Shape_Leng[1:10]

relErr <- list() # will hold the relative errors
len <- list() # hold len's
for (i in 1:(length(t)/2)) { # same size as sample of trips
  len[[i]] = sum(t[2,i][[1]]) # length recalculated per loop iteration
  real = DataSet$trip_distance[i] * 1609.34 # converting between miles and meters
  relErr[[i]] = abs(len[[i]]-real)/real # store the relative error
}
m = mean(unlist(relErr))
s = sd(unlist(relErr))
me = median(unlist(relErr))

lenFiltered <- list() # holds filtered trips
relErrFiltered <- list() # holds relErr for filtered trips
j = 1
for (i in 1:length(len)) {	# filter out items that are one sd above
  if (relErr[[i]] < (m+s)) {
    real = DataSet$trip_distance[i] * 1609.34
    relErrFiltered[[j]] = (abs(len[[i]]-real))/real
    j = j + 1
  }
}
mFilt <- mean(unlist(relErrFiltered))
sFilt <- sd(unlist(relErrFiltered))
meFilt <- median(unlist(relErrFiltered))

loss <- length(relErr) - length(relErrFiltered)
lossPercent <- loss/length(relErr) * 100

sink("ErrorReport.txt")

cat("Error report for data stored in Trips.txt:", '\n')
cat("Mean relative trip length error, compared to DataSet: ", m*100, "%", '\n', sep="")
cat("Standard deviation of relative errors: ", s*100, "%", '\n', sep="")
cat("Median relative error: ", me*100, "%", '\n', sep="")
cat("---", '\n')
cat("Error report for data stored in TripsFiltered.txt:", '\n', sep="")
cat("Mean relative trip length error, with values above one sd removed: ", mFilt*100, "%", '\n', sep="")
cat("Standard deviation: ", sFilt*100, "%", '\n', sep="")
cat("Median: ", mFilt*100, "%", '\n', sep="")
cat("Number of trips filtered out:", loss, '\n')
cat("% of trips filtered out: ", lossPercent, "%", '\n', sep="")

sink()


names(DataSet)
DataSet$tod[816]
head(DataSet,n=10)

#--

testIn <- read.table("programOutputFiltered.csv", header=T, sep=",") # reading it all back in from the file
head(testIn, n=20)





sink("programOutput.csv") # opens a sink; ** rename to what you want your output to be **

cat("TripNumber", "SteetinTrip", "StreetLength", "SpeedLimit", sep=","); cat('\n') # the header; ** remember to add more names if you extract more attributes than just speed and street lengths **
for (i in 1:(length(t[1,])/2)) { # for all groupings in the tripAttributes we just gathered
  if (length(t[1,i][[1]]) != 0) { # this occurs with loop trips; to "filter" out, we will simply fail to write them
    j <- 1:length(t[1,i][[1]]) # what street number of the trip that we are on (resets each row)
    cat(rbind( i, ',', j, ',', unlist(t[1,i]), ',', unlist(t[2,i]), '\n' ), sep="")#; cat('\n') # print the trip row number and necessary attributes
    # ** you may add more attributes to this list by interting a <, unlist(t[x,i]), > before the closing parenthesis of rbind() [directly after unlist(t[2,i])], where x is the 3rd, 4th, etc. attribute gathered from attCatch **
  }
}

sink() # closes sink

n <- which(relErr < (m+s))
real <- DataSet$trip_distance[n] * 1609.34
len <- len[n]
relErrFiltered <- abs(len-real) / real

relErrFiltered <- unlist(relErr[relErr < (m+s)])

m+s


length(len) == length(real)
typeof(real)
length(relErrFiltered)
head(relErrFiltered, n=10)
head(relErrFiltered, n=10)

length(n)
head(len)

for (i in 1:length(len)) {	# filter out items that are one sd above
  if (relErr[[i]] < (m+s)) {
    #		real = DataSet$trip_distance[i] * 1609.34
    relErrFiltered[[j]] = (abs(len[[i]]-real))/real
    j = j + 1
  }
}







relErr <- list() # will hold the relative errors
len <- vector(mode="double", length=length(DataSet$trip_distance)) # hold len's
for (i in 1:(length(t)/2)) { # same size as sample of trips
  len[[i]] = sum(t[2,i][[1]]) # length recalculated per loop iteration
}
real = DataSet$trip_distance * 1609.34
relErr <- abs((len-real)/real)
m = mean(unlist(relErr))
s = sd(unlist(relErr))
me = median(unlist(relErr))

length(tripAttributes[1,])
q <- 1:(length(tripAttributes[1,])/2)
q

n <- which(tripAttributes[1,q][[1]] != 0)
n



length(tripAttributes)
tripAttributes

for (i in 1:(length(tripAttributes[1,])/2)) { # for all groupings in the tripAttributes we just gathered
  if (length(tripAttributes[1,i][[1]]) != 0) { # this occurs with loop trips; to "filter" out, we will simply fail to write them
    cat(i, rbind( unlist(tripAttributes[1,i]), unlist(tripAttributes[2,i]) ), sep=","); cat('\n') # print the trip row number and necessary attributes
    # ** you may add more attributes to this list by interting a <, unlist(tripAttributes[x,i]), > before the closing parenthesis of rbind(), where x is the 3rd, 4th, etc. attribute gathered from attCatch **
  }
}





n <- which(relErr < (m+s))

cat("TripNumber", "StreetLength", "SpeedLimit", sep=","); cat('\n') # the header; ** again, remember to add more names if you extract more attributes than just speed and street lengths **
for (i in 1:(length(tripAttributes[1,])/2)) { # for all groupings in the tripAttributes we just gathered
  if (length(tripAttributes[1,i][[1]]) != 0 && relErr[i]<(m+s)) { # will only write out trips with appropriate relErr values ** you may seek to change this if you perform custom filtering **
    cat(i, rbind( unlist(tripAttributes[1,i]), unlist(tripAttributes[2,i]) ), sep=","); cat('\n')
    # ** again, you may add more attributes to this list by interting a <, unlist(tripAttributes[x,i]), > before the closing parenthesis of rbind(), where x is the 3rd, 4th, etc. attribute gathered from attCatch **
  }
}




rm(relErr); rm(relErrFiltered); gc()






pOut <- read.csv("programOutput.csv", header=T, sep=",")
pOutF <- read.csv("programOutputFiltered.csv", header=T, sep=",")
head(pOutF, n=1000)




xCoords <- matrix(nrow=length(MAPraw[[2]][,1]), ncol=1) # x coord for street ("node") intersections
yCoords <- matrix(nrow=length(MAPraw[[2]][,1]), ncol=1) # y coord for street ("node") intersections
CoordRaw <- MAPraw[[2]][,2]
for (i in (1:10)) { # gathering coordinates from MAPraw, the loaded-in shapefile
  MAPxCoords[i] <- CoordRaw[[i]][1] # x column of shapefile coordinates
  MAPyCoords[i] <- CoordRaw[[i]][2] # y column of shapefile coordinates
}
xCoords[1:10,1] <- unlist(CoordRaw)[(1:10)*2 - 1] # odd entries
yCoords[1:10,1] <- unlist(CoordRaw)[(1:10)*2] # even entries

MAPxCoords == xCoords

MAPxCoords
xCoords
yCoords
zCoords[1:20,1]=c(xCoords,yCoords)
zCoords

c(xCoords,yCoords)

yCoords <- matrix(nrow=10, ncol=1) # x coord for street ("node") intersections

as.numeric(as.character(MAPraw[[5]]$SPEED[1:10]))
typeof(MAPraw[[5]]$SPEED)
cat(MAPraw[[5]]$SPEED[6])


xCoords <- CoordRaw[[]][1]
yCoords <- CoordRaw[[1:10]][2]
xCoords
xCoords <- CoordRaw[1:
                      CoordRaw[1:10][[2:3]]
                    c <- CoordRaw[1:3]
                    c[[1]][1]
                    d <- 1:10
                    unlist(CoordRaw[1:3])
                    
                    
                    testFunc <- function(alpha) {
                      bet <<- alpha
                    }
                    
                    as.numeric(data.matrix(unlist(MAPraw[[5]]$Shape_Leng[1:10])))
                    as.numeric(data.matrix(unlist(MAPraw[[5]]$SPEED[1:10])))
                    
                    
                    testF <- function(beep, bop) {
                      return (beep$bop[1])
                    }
                    typeof(s)
                    s <- "SPEED"
                    MAPraw[[5]]$s
                    testF(MAPraw[[5]], "SPEED")
                    MAPraw[[5]]$SPEED[1]
                    
                    MAPraw[[5]]$SPEED[3081]
                    
                    get.shortest.paths(MAPgraph, from=tripStartNodeIDs[1:10], to=tripEndNodeIDs[1:10])$vpath[[1]]
                    
                    testFunc(alpha=2)
                    alpha
                    bet
                    
                    tripStartNodeIDs[1:10]
                    tripEndNodeIDs[1:10]
                    
                    test <- attCatch(tripStartNodeIDs[1:10], tripEndNodeIDs[1:10])
                    
                    test
                    attCatch(1037,2139)
                    
                    length(CoordRaw)
                    unlist(CoordRaw)
                    head(CoordRaw)
                    sMatrix = matrix(c(MAPxCoords, MAPyCoords), nrow=length(MAPxCoords), ncol=2) # matrix of all street intersection coordinates
                    
                    
                    j <- droplevels(head(unlist(MAPraw[[5]]$SPEED))[2])
                    j <- MAPraw[[5]]$SPEED[[2]]
                    droplevels(j)+12
                    j
                    as.integer(j[[1]])
                    as.numeric(as.character(j))+12
                    as.character(j)
                    j[1]
                    levels(j)
                    
                    ?droplevels
                    
                    # # # # # ALGO EFFICIENCY # # # # #
                    
                    # fix the following:
                    # 1. (also adapt this to all other parts of code dealing with these for-loops)
                    MAPxCoords <- matrix(nrow=length(MAPraw[[2]][,1]), ncol=1) # x coord for street ("node") intersections
                    MAPyCoords <- matrix(nrow=length(MAPraw[[2]][,1]), ncol=1) # y coord for street ("node") intersections
                    #CoordRaw <- MAPraw[[2]][,2]
                    #for (i in 1:(length(MAPxCoords))) { # gathering coordinates from MAPraw, the loaded-in shapefile
                    #	MAPxCoords[i] <- CoordRaw[[i]][1] # x column of shapefile coordinates
                    #	MAPyCoords[i] <- CoordRaw[[i]][2] # y column of shapefile coordinates
                    #}
                    CoordRaw <- unlist(MAPraw[[2]][,2])
                    MAPxCoords <- CoordRaw[(1:length(CoordRaw))*2 - 1] # odd entries
                    MAPyCoords <- CoordRaw[(1:length(CoordRaw))*2] # even entries
                    streetMatrix = matrix(c(MAPxCoords, MAPyCoords), nrow=length(MAPxCoords), ncol=2) # matrix of all street intersection coordinates
                    
                    # 2.
                    # why does nn2() remove 2 million rows of data
                    
                    
                    
                    
                    
                    
                    
sink("programOutput4.csv")
cat("TripNumber", "StreetInTrip", "StreetLength", "SpeedLimit","EdgeID", sep=","); cat('\n')
for (i in 1:(length(tripAttributes[1,]))) {
if (length(tripAttributes[1,i][[1]]) != 0) {
          j <- 1:length(tripAttributes[1,i][[1]])
            cat(rbind( i, ',', j, ',', unlist(tripAttributes[1,i]), ',', unlist(tripAttributes[2,i]), ',',unlist(tripAttributes[3,i], '\n' ), sep=""))
                      }
                    }
                    
                    sink()
                    
                    
                    
                    