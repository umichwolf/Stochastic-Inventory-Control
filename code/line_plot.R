#TODO: Written in hurry. Clean up this mess and reuse functions

library(ggplot2)

db_illustration <- read.csv("output/db_illustration.csv")
db_ill_NetInventory <- data.frame(TimePeriod = 0:20, ParameterValue=db_illustration$NetInventory, Group = rep("NetInventory", 21))
db_ill_Cost <- data.frame(TimePeriod = 0:20, ParameterValue=db_illustration$Cost, Group = rep("Cost", 21))
db_ill_Demand <- data.frame(TimePeriod = 0:20, ParameterValue=db_illustration$Demand, Group = rep("Demand", 21))
db_ill_Ordering <- data.frame(TimePeriod = 0:20, ParameterValue=db_illustration$Ordering, Group = rep("Ordering", 21))
db_ill_dat <- rbind(db_ill_NetInventory,db_ill_Cost,db_ill_Demand,db_ill_Ordering)
DualBalancingParameters <- ggplot(db_ill_dat, aes(x=TimePeriod, y=ParameterValue, group=Group, colour=Group)) + geom_line(size=1.5) + ggtitle("Dual Balancing Parameters") + theme(axis.text=element_text(size=6),axis.title=element_text(size=6), title=element_text(size=8),legend.position='bottom',legend.title = element_text(color="#ffffff",size=10),legend.text = element_text(size=6))

db_output01 <- read.csv("output/db_output01.csv")
db_output02 <- read.csv("output/db_output02.csv")
db_output03 <- read.csv("output/db_output03.csv")
db_output04 <- read.csv("output/db_output04.csv")

my_output01 <- read.csv("output/my_output01.csv")
my_output02 <- read.csv("output/my_output02.csv")
my_output03 <- read.csv("output/my_output03.csv")
my_output04 <- read.csv("output/my_output04.csv")

myopic_bad_case <- read.csv("output/myopic_bad_case.csv")

tc_to <- read.csv("output/tc_to.csv")
tc_to_y <- (tc_to$TotalCostN/tc_to$TotalOrderingN)
tc_to_y_2 <- (tc_to$TotalCostIN/tc_to$TotalOrderingIN)
tc_to_x <- (0:19)*20
line_1 <- data.frame(TimePeriod = tc_to_x, TotalCostToOrderRatio = tc_to_y, Group = rep("Independent Normal", length(tc_to_x)))
line_2 <- data.frame(TimePeriod = tc_to_x, TotalCostToOrderRatio = tc_to_y_2, Group = rep("Normal Increment", length(tc_to_x)))
tc_to_dat <- rbind(line_1,line_2)
TotalCostToOrderRatio <- ggplot(tc_to_dat, aes(x=TimePeriod, y=TotalCostToOrderRatio, group=Group, colour=Group)) + geom_line(size=1.5)+ylim(2.75,4) + ggtitle("Cost Demand Ratio") + theme(axis.text=element_text(size=6),axis.title=element_text(size=6), title=element_text(size=8),legend.position='bottom',legend.title = element_text(color="#ffffff",size=10),legend.text = element_text(size=6))

net_inv_1 <- read.csv("output/net_inv_1.csv")
net_inv_2 <- read.csv("output/net_inv_2.csv")

ni_1_line_1 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_1$NetInventory0, Group = rep("x=0", 31))
ni_1_line_2 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_1$NetInventory20, Group = rep("x=20", 31))
ni_1_line_3 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_1$NetInventory40, Group = rep("x=40", 31))
ni_1_line_4 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_1$NetInventory60, Group = rep("x=60", 31))
net_inv_1_dat <- rbind(ni_1_line_1,ni_1_line_2,ni_1_line_3,ni_1_line_4)
NetInventory <- ggplot(net_inv_1_dat, aes(x=TimePeriod, y=NetInventory, group=Group, colour=Group)) + geom_line(size=1) + ggtitle("Net Inventory") + theme(axis.text=element_text(size=6),axis.title=element_text(size=6), title=element_text(size=8),legend.position='bottom',legend.title = element_text(color="#ffffff",size=10),legend.text = element_text(size=6))

ni_2_line_1 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_2$LeadTime2, Group = rep("L=2", 31))
ni_2_line_2 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_2$LeadTime4, Group = rep("L=4", 31))
ni_2_line_3 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_2$LeadTime6, Group = rep("L=6", 31))
ni_2_line_4 <- data.frame(TimePeriod = 0:30, NetInventory = net_inv_2$LeadTime8, Group = rep("L=8", 31))
net_inv_2_dat <- rbind(ni_2_line_1,ni_2_line_2,ni_2_line_3,ni_2_line_4)
LeadTime <- ggplot(net_inv_2_dat, aes(x=TimePeriod, y=NetInventory, group=Group, colour=Group)) + geom_line(size=1) + ggtitle("Lead Time") + theme(axis.text=element_text(size=6),axis.title=element_text(size=6), title=element_text(size=8),legend.position='bottom',legend.title = element_text(color="#ffffff",size=10),legend.text = element_text(size=6))

#dp_output01 <- read.csv("output/dp_output01.csv")
#dp_output02 <- read.csv("output/dp_output02.csv")
#dp_output03 <- read.csv("output/dp_output03.csv")
#dp_output04 <- read.csv("output/dp_output04.csv")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plot_2_graph <- function(y,y_title,y2,y2_title,y_axs_title,nrows,title) {
  line_1 <- data.frame(TimePeriod = 0:nrows, ParamValue = y, Group = rep(y_title, nrows+1))
  line_2 <- data.frame(TimePeriod = 0:nrows, ParamValue = y2, Group = rep(y2_title, nrows+1))
  dat <- rbind(line_1,line_2)
  ggplot(dat, aes(x=TimePeriod, y=ParamValue, group=Group, colour=Group)) + ggtitle(title) + theme(axis.text=element_text(size=6),axis.title=element_text(size=6), title=element_text(size=8),legend.position='bottom',legend.title = element_text(color="#ffffff",size=10),legend.text = element_text(size=6))
}

plot_graph <- function(y,y_title,y2,y2_title,y3,y3_title,y4,y4_title,y_axs_title,nrows,title) {
  line_1 <- data.frame(TimePeriod = 0:nrows, ParamValue = y, Group = rep(y_title, nrows+1))
  line_2 <- data.frame(TimePeriod = 0:nrows, ParamValue = y2, Group = rep(y2_title, nrows+1))
  line_3 <- data.frame(TimePeriod = 0:nrows, ParamValue = y3, Group = rep(y3_title, nrows+1))
  line_4 <- data.frame(TimePeriod = 0:nrows, ParamValue = y4, Group = rep(y4_title, nrows+1))
  dat <- rbind(line_1,line_2,line_3,line_4)
  ggplot(dat, aes(x=TimePeriod, y=ParamValue, group=Group, colour=Group)) + geom_point(size=1.5,alpha=.7) + ggtitle(title) + theme(legend.position = "none",axis.text=element_text(size=6),axis.title=element_text(size=10), title=element_text(size=9))
}

#Independent Normal
AccumulativeDemandAndCost_Normal <- plot_graph(db_output01$AccumulativeDemand, "DualBalancingAccumulativeDemand", db_output01$AccumulativeCost,"DualBalancingAccumulativeCost",my_output01$AccumulativeDemand,"MyopicAccumulativeDemand",my_output01$AccumulativeCost,"MyopicAccumulativeCost","AccumulativeDemandAndCost",100,"Accumulative Demand vs. Accumulative Cost:\n Independent Normal")
#Binomial Increment
AccumulativeDemandAndCost_Binomial <- plot_graph(db_output02$AccumulativeDemand, "DualBalancingAccumulativeDemand", db_output02$AccumulativeCost,"DualBalancingAccumulativeCost",my_output02$AccumulativeDemand,"MyopicAccumulativeDemand",my_output02$AccumulativeCost,"MyopicAccumulativeCost","AccumulativeDemandAndCost",100,"Accumulative Demand vs. Accumulative Cost:\n Binomial Increment")
#Brownian Motion with Drift
AccumulativeDemandAndCost_Brownian <- plot_graph(db_output03$AccumulativeDemand, "DualBalancingAccumulativeDemand", db_output03$AccumulativeCost,"DualBalancingAccumulativeCost",my_output03$AccumulativeDemand,"MyopicAccumulativeDemand",my_output03$AccumulativeCost,"MyopicAccumulativeCost","AccumulativeDemandAndCost",100,"Accumulative Demand vs. Accumulative Cost:\n Brownian Motion with Drift")
#3-period Markov Chain
AccumulativeDemandAndCost_Markov <- plot_graph(db_output04$AccumulativeDemand, "DualBalancingAccumulativeDemand", db_output04$AccumulativeCost,"DualBalancingAccumulativeCost",my_output04$AccumulativeDemand,"MyopicAccumulativeDemand",my_output04$AccumulativeCost,"MyopicAccumulativeCost","AccumulativeDemandAndCost",100,"Accumulative Demand vs. Accumulative Cost:\n 3-period Markov Chain")

#Independent Normal
DemandAndNetInventory_Normal <- plot_graph(db_output01$Demand, "DualBalancingDemand", db_output01$NetInventory,"DualBalancingNetInventory",my_output01$Demand,"MyopicDemand",my_output01$NetInventory,"MyopicNetInventory","DemandAndNetInventory",100,"Demand vs. Net Inventory:\n Independent Normal")
#Binomial Increment
DemandAndNetInventory_Binomial <- plot_graph(db_output02$Demand, "DualBalancingDemand", db_output02$NetInventory,"DualBalancingNetInventory",my_output02$Demand,"MyopicDemand",my_output02$NetInventory,"MyopicNetInventory","DemandAndNetInventory",100,"Demand vs. Net Inventory:\n Binomial Increment")
#Brownian Motion with Drift
DemandAndNetInventory_Brownian <- plot_graph(db_output03$Demand, "DualBalancingDemand", db_output03$NetInventory,"DualBalancingNetInventory",my_output03$Demand,"MyopicDemand",my_output03$NetInventory,"MyopicNetInventory","DemandAndNetInventory",100,"Demand vs. Net Inventory:\n Brownian Motion with Drift")
#3-period Markov Chain
DemandAndNetInventory_Markov <- plot_graph(db_output04$Demand, "DualBalancingDemand", db_output04$NetInventory,"DualBalancingNetInventory",my_output04$Demand,"MyopicDemand",my_output04$NetInventory,"MyopicNetInventory","DemandAndNetInventory",100,"Demand vs. Net Inventory:\n 3-period Markov Chain")

#Myopic Bad: Demand
MyopicBadDemand <- plot_2_graph(myopic_bad_case$DemandMY, "MyopicDemand", myopic_bad_case$DemandDB,"DualBalancingDemand","blah",200,"Myopic Bad Case: Comparing Demand") + geom_point() +ylab("Demand")
#Myopic Bad: Net Inventory
MyopicBadNetInventory <- plot_2_graph(myopic_bad_case$NetInventoryMY, "MyopicNetInventory", myopic_bad_case$NetInventoryDB,"DualBalancingNetInventory","blah",200,"Myopic Bad Case: Comparing Net Inventory") + geom_line() +ylab("Net Inventory")
#Myopic Bad: Cumulative Cost
MyopicBadAccumulativeCost <- plot_2_graph(myopic_bad_case$AccumulativeCostMY, "MyopicAccumulativeCost", myopic_bad_case$AccumulativeCostDB,"DualBalancingAccumulativeCost","blah",200,"Myopic Bad Case:\nComparing Accumulative Cost") + geom_line() +ylab("Accumulative Cost")
#Myopic Bad: Ordering
MyopicBadOrdering <- plot_2_graph(myopic_bad_case$OrderingMY, "MyopicOrdering", myopic_bad_case$OrderingDB,"DualBalancingOrdering","blah",200,"Myopic Bad Case: Comparing Ordering")+xlim(0,30) + geom_point() +ylab("Ordering")

ggsave(MyopicBadDemand,file="figures/MyopicBadDemand.png", scale=.5)
ggsave(MyopicBadNetInventory,file="figures/MyopicBadNetInventory.png", scale=.5)
ggsave(MyopicBadAccumulativeCost,file="figures/MyopicBadAccumulativeCost.png", scale=.5)
ggsave(MyopicBadOrdering,file="figures/MyopicBadOrdering.png", scale=.5)
ggsave(AccumulativeDemandAndCost_Normal,file="figures/AccumulativeDemandAndCost_Normal.png", scale=.5)
ggsave(AccumulativeDemandAndCost_Binomial,file="figures/AccumulativeDemandAndCost_Binomial.png", scale=.5)
ggsave(AccumulativeDemandAndCost_Brownian,file="figures/AccumulativeDemandAndCost_Brownian.png", scale=.5)
ggsave(AccumulativeDemandAndCost_Markov,file="figures/AccumulativeDemandAndCost_Markov.png", scale=.5)
ggsave(DemandAndNetInventory_Normal,file="figures/DemandAndNetInventory_Normal.png", scale=.5)
ggsave(DemandAndNetInventory_Binomial,file="figures/DemandAndNetInventory_Binomial.png", scale=.5)
ggsave(DemandAndNetInventory_Brownian,file="figures/DemandAndNetInventory_Brownian.png", scale=.5)
ggsave(DemandAndNetInventory_Markov,file="figures/DemandAndNetInventory_Markov.png", scale=.5)
ggsave(DualBalancingParameters,file="figures/DualBalancingParameters.png", scale=.5)
ggsave(TotalCostToOrderRatio,file="figures/TotalCostToOrderRatio.png", scale=.5)

ggsave(NetInventory,file="figures/NetInventory.png", scale=.5)
ggsave(LeadTime,file="figures/LeadTime.png", scale=.5)

multiplot(p1,p2,p3,p4,cols=2)
#multiplot(q1,q2,q3,q4,cols=2)