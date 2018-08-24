#Importing and cleaning data--------------
progop.df<-read.csv('test_sample.csv')
progop.df<-progop.df[-c(1,2,3)]
colnames(progop.df)[20]="Trip Start Node ID"
colnames(progop.df)[21]="Trip End Node ID"
colnames(progop.df)[1]="Speed Limit"
colnames(progop.df)[2]="Street Length"


countID<-plyr::count(progop.df$Trip.ID)
progop.df<-progop.df[order(progop.df$Trip.ID),]
progop.df["Time on Egde"]<-0


#Feature extaction and addition-----------------
#extract the avg. time on each edge/street(trip time/number of streets)-----
for (i in 1:nrow(progop.df)){
  if (progop.df[i,4]==progop.df[i+1,4]){
    progop.df[i,19]=progop.df[i,14]/countID[(countID$x==progop.df[i,4]),2]
    
  }
  else if (progop.df[i,4]==progop.df[i-1,4]){
    progop.df[i,19]=progop.df[i,14]/countID[(countID$x==progop.df[i,4]),2]
    
  }
}

#Previous and next Edge ID----
progop.df["Prev. Edge ID"]<-0
progop.df["Next Edge ID"]<-0

for (i in 1:nrow(progop.df)){
  if (progop.df[i,4]==progop.df[i+1,4]){
    progop.df[i,20]=progop.df[i+1,3]
    
    
    
  }
  
}


for (i in 2:nrow(progop.df)){
  if (progop.df[i,4]==progop.df[i-1,4]){
    progop.df[i,21]=progop.df[i-1,3]
    
    
    
  }
  
}


#Weekday/Weekend
progop.df$pickup_datetime<-as.character(progop.df$pickup_datetime)
progop.df$dropoff_datetime<-as.character(progop.df$dropoff_datetime)
progop.df$pickup_datetime<-substring(progop.df$pickup_datetime,1,8)
progop.df$dropoff_datetime<-substring(progop.df$dropoff_datetime,1,8)

progop.df<-progop.df[-6]


progop.df$pickup_datetime<-as.Date(progop.df$pickup_datetime,"%d-%m-%y")

progop.df$DayofWeek<-weekdays(progop.df$pickup_datetime)

progop.df$DayofWeek<-as.factor(progop.df$DayofWeek)
progop.df$month<-month(progop.df$pickup_datetime) #No point as its only for June
colnames(progop.df)[5]<-"DATE"
#Speed on every Edge---------

progop.df$Edgespeed<-(progop.df$`Street Length`/progop.df$`Time on Egde`)


#Weather Data--------
weather.df<-read.csv('Weather June.csv')

#keeping only relevant weather data---
weather.df<-weather.df[,c(3,4,6,10,11,12,13)]

weather.df$DATE<-as.Date(weather.df$DATE,"%d-%m-%y")

#Merge with progop-----

speed.df<-merge(x = progop.df, y = weather.df, by = "DATE", all.x = TRUE)

speed.df<-speed.df[order(speed.df$EdgeID),]
which(table(speed.df$EdgeID) == max(table(speed.df$EdgeID))) #2671


#Convert edge speed from m/s to mph
speed.df$Edgespeed<-(speed.df$Edgespeed)*2.23

#extracting a sample for EDA------
sample.df<-speed.df[speed.df$EdgeID==1664,]




#Will have to account for error in trip distance and summation of edge length, which I was not able to do using anthony script


#EDA-----

#effect of weather
library(ggplot2)
sample.df$TMAX<-as.numeric(sample.df$TMAX)     



plot(sample.df$TMAX,sample.df$Edgespeed)

plot(sample.df$PRCP,sample.df$Edgespeed)

#prev and next edge id in this case is constant

#effect of day of the week
plot(sample.df$DayofWeek,sample.df$Edgespeed)

#effect of TOD
sample.df$tod<-as.factor(sample.df$tod)
plot(sample.df$tod,sample.df$Edgespeed)

#effect of startid,endid

plot(sample.df$startID,sample.df$Edgespeed)

plot(sample.df$endID,sample.df$Edgespeed)

#fare amount


plot(sample.df$fare_amount,sample.df$Edgespeed)

#relation bw avg speed of the trip vs actual/estimated edge speed


plot(sample.df$speed,sample.df$Edgespeed) #this kind of validates our feature extraction of edge speed















