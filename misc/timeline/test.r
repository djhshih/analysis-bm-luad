library(ggplot2)

timeset<-data.frame(year=c(1986,1995,2011,1990,1998,2010),text=c('I was born','Had a nice icecream','Spotted a dodo','First swim','Crashed my bicycle','Bought a helmet'),y=c(runif(3,.5,1),runif(3,-1,-.5)))

g <- ggplot(timeset, aes(x=year, y=0))

g <- g + geom_segment(aes(y=0, yend=y, xend=year))

g <- g + geom_text(aes(y=y,label=text),size=2.5,vjust=-1)

g <- g + geom_point(aes(y=y))

g <- g + scale_y_continuous(limits=c(-2,2))

g <- g + geom_hline(y=0,size=1,color='purple')   #draw a vertical line

g <- g + geom_segment(x=2011.4,xend=2012.2,y=.2,yend=0,color='purple',size=1) + geom_segment(x=2011.4,xend=2012.2,y=-.2,yend=0,color='purple',size=1) #drawing the actual arrow

g <- g + opts(axis.text.y =theme_blank(),title='My timeline')+ylab('')+xlab('')

print(g)
