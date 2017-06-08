#Shelli Kesler 12/2/16
#Determines fit of degree distribution with power law
#Compares fit between two groups
#Assumes degree data are listed per subject per region and are in adjacent columns
#NOTE: Attach dataframe before running

#See https://cran.r-project.org/web/packages/poweRlaw/vignettes/b_powerlaw_examples.pdf

#degcols = columns containing degree data
powerL <- function(data,groupvar,degcols){ #e.g. powerL(mydata,Group,19:50)
  library(poweRlaw)
  degcols<<-degcols
  #get group labels and subset groups
  Gs<-as.data.frame(table(groupvar)) 
  Gtemp1<-data[which(groupvar==as.character(Gs[1,1])),];Gtemp1=droplevels(Gtemp1)
  Gtemp2<-data[which(groupvar==as.character(Gs[2,1])),];Gtemp2=droplevels(Gtemp2)
  
  #get mean degree for each region
  G1<-as.integer(colMeans(Gtemp1[,degcols]))
  G2<-as.integer(colMeans(Gtemp2[,degcols]))
  
  G1pl <<- displ$new(G1)#fit discrete power law
  est1 <<- estimate_xmin(G1pl)#estimate lower threshold
  G1pl$setXmin(est1)#update the power law object
  bs_p1 <<- bootstrap_p(G1pl)#goodness of fit test, large p = good fit
  
  G2pl <<- displ$new(G2)
  est2 <<- estimate_xmin(G2pl)
  G2pl$setXmin(est1) #each distribution must have the same lower threshold for comparison
  
  #estimate parameters for Group2 using Group1 threshold
  G2pars <- estimate_pars(G2pl)
  G2pl$setPars(G2pars)
  bs_p2 <<- bootstrap_p(G2pl)
  
  #compare distributions
  comp = compare_distributions(G1pl, G2pl)

par(mfrow=c(2,2))
plot(G1pl,ylab="Log Cumulative Distribution",xlab="Log Degree", 
     main=as.character(Gs[1,1]));lines(G1pl, col=2)
plot(G2pl,ylab="Log Cumulative Distribution",xlab="Log Degree", 
     main=as.character(Gs[2,1]));lines(G2pl,col=2)
results=list("G1fit"=bs_p1$p,"G2fit"=bs_p2$p,"R"=comp$test_statistic,"p"=comp$p_one_sided)
return(results)
}