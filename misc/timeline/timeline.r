require(googleVis) ## googleVis 0.5.0-3

dat <- data.frame(Room=c("Room 1","Room 2","Room 3"),
                  Language=c("English", "German", "French"),
                  start=as.POSIXct(c("2014-03-14 14:00", 
                                     "2014-03-14 15:00",
                                     "2014-03-14 14:30")),
                  end=as.POSIXct(c("2014-03-14 15:00", 
                                   "2014-03-14 16:00",
                                   "2014-03-14 15:30")))
plot(
  gvisTimeline(data=dat, 
               rowlabel="Room", barlabel="Language", 
               start="start", end="end")
)
