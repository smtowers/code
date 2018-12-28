

if (0){
   tooltip <- function(x) {
     if (is.null(x)) return(NULL)
     if (is.null(x$ID)) return(NULL)
 
     # Pick out the movie with this ID
     if (0){
     all_movies <- isolate(movies())
     movie <- all_movies[all_movies$ID == x$ID, ]
 
     paste0("<b>", movie$Title, "</b><br>",
       movie$Year, "<br>",
       "$", format(movie$BoxOffice, big.mark = ",", scientific = FALSE)
     )
     }
   }
}
##################################################################################
##################################################################################
##################################################################################
return_p_value_string = function(p,ladd_pequals=T){
  pb = round(p,3)
  pb[p>=0.10] = as.character(round(p[p>=0.10],2))
  pb[p>=0.30] = as.character(round(p[p>=0.30],1))
  l = which(p>=0.99&nchar(pb)==1)
  pb[l] = paste("1.0",sep="")
  l = which(p>=0.10&nchar(pb)==3)
  pb[l] = paste(pb[l],"0",sep="")
  l = which(p<0.10&p>0.001&nchar(pb)==4)
  pb[l] = paste(pb[l],"0",sep="")
  pb[p<0.001] = "<0.001"
  if (ladd_pequals==T){
    pb[p>=0.001] = paste("p=",pb[p>=0.001],"",sep="")
    pb[p<0.001] = paste("p<0.001",sep="")
  }else{
    pb = paste("",pb,"",sep="")
  }
  return(pb)
}


##################################################################################
##################################################################################
# https://stackoverflow.com/questions/6177629/how-to-silence-the-output-from-this-r-package
##################################################################################
shut_up = function(expr) {
  #temp file
  f = file()

  #write output to that file
  sink(file = f)

  #evaluate expr in original environment
  y = eval(expr, envir = parent.frame())

  #close sink
  sink()

  #get rid of file
  close(f)

  y
}

##################################################################################
##################################################################################
# The PMF of the negative binomial distribution
#
# mean of the distribution is m
# variance of the distribution is  m+m^2/alpha
##################################################################################
my_dnegbinom=function(yobs,ypred,alpha){
  #https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0000180
  x = yobs
  m = ypred
  k = 1/alpha
  probb = lgamma(k+x) -lgamma(k) -lfactorial(x) +x*log(m/(m+k)) -k*log(1+m/k)
  probb = exp(probb)
  return(probb)
}

##################################################################################
##################################################################################
# The PMF of the truncated negative binomial distribution
##################################################################################
dtrunc_negbinom=function(yobs,mu,alpha,truncation){
  vdenom = rep(1,length(yobs))
  if (truncation>0){
    for (i in 0:(truncation-1)){
      vdenom = vdenom - my_dnegbinom(i,mu,alpha)
    }
  }
  vprob = my_dnegbinom(yobs,mu,alpha)/vdenom
  vprob[yobs<truncation] = 0
  return(vprob)
}

##################################################################################
##################################################################################
# calcuate the expected value of mu for the negative binomial distribution
##################################################################################
negbinom_fun  = function(par,mylist){
  data = (mylist[[1]])
  trunc = mylist[[2]]
  lfit = mylist[[3]]

  r = par[2]
  if (lfit==2) r = par[2]+par[3]*data$banned
  if (lfit==3) r = par[2]+par[3]*data$banned+par[4]*data$none
  if (lfit==4) r = par[2]+par[3]*data$dateb
  if (lfit==5) r = par[2]+par[3]*data$dateb+par[4]*data$banned
  if (lfit==6) r = par[2]+par[3]*data$dateb+par[4]*data$hcsap+par[5]*data$assault
  if (lfit==7) r = par[2]+par[3]*data$dateb+par[4]*data$banned+par[5]*data$none
  if (lfit==8) r = par[2]+par[3]*data$only_pistol+par[4]*data$only_rifle+par[5]*data$pistol_and_rifle
  #cat(lfit,"  ",data[1,],"  ",r,"\n")
  #cat(lfit,r,"\n")
  p = exp(r)
  return(p)
}

##################################################################################
##################################################################################
# calculate the negative log likelihood of the truncated negative binomial
##################################################################################
negll_trunc_negbinom  = function(par,mylist){

  thedata = mylist[[1]]
  trunc = mylist[[2]]

  mu = negbinom_fun(par,mylist)
  vprob = dtrunc_negbinom(thedata$yobs,mu,exp(par[1]),trunc)

  neg_log_like = sum(-log(vprob))
  return(neg_log_like)
}





##################################################################################
##################################################################################
##################################################################################
# log-series distribution
# https://web.archive.org/web/20110726144520/http://www.math.mcgill.ca/~dstephens/556/Papers/Fisher1943.pdf
#@article{fisher1943relation,
#  title={The relation between the number of species and the number of individuals in a random sample of an animal population},
#  author={Fisher, Ronald A and Corbet, A Steven and Williams, Carrington B},
#  journal={The Journal of Animal Ecology},
#  pages={42--58},
#  year={1943},
#  publisher={JSTOR}
#}
#
#https://en.wikipedia.org/wiki/Logarithmic_distribution
#http://influentialpoints.com/Training/logarithmic_series_distribution.htm
#ftp://ftp.ems.psu.edu/data/pub/geosc/pub/brachio/0.pdf
#
# mean of the distribution is -(1/log(1-p))*(p/(1-p))
#
# if mu and k=1/alpha are the parameters of the NB distribution (see above)
# then as k->0 we have p -> m/(m+k)
# and p/(1-p) -> (m/k)
##################################################################################
my_dlogseries=function(yobs,p){
  probb = -(p)^yobs/(log(1-p)*yobs)
  return(probb)
}

##################################################################################
##################################################################################
##################################################################################
dtrunc_logseries = function(yobs,p,truncation){
  vdenom = rep(1,length(yobs))
  if (truncation>1){
    for (i in 1:(truncation-1)){
      vdenom = vdenom - my_dlogseries(i,p)
    }
  }
  vprob = my_dlogseries(yobs,p)/vdenom
  vprob[yobs<truncation] = 0
  return(vprob)
}

##################################################################################
##################################################################################
##################################################################################
logseries_fun  = function(par,mylist){
  data = mylist[[1]]
  trunc = mylist[[2]]
  lfit = mylist[[3]]

  # r is log(p/(1-p))
  r = par[1]
  if (lfit==2) r = par[1] + par[2]*data$banned
  if (lfit==3) r = par[1] + par[2]*data$banned+par[3]*data$none
  if (lfit==4) r = par[1]+par[2]*data$dateb
  if (lfit==5) r = par[1]+par[2]*data$dateb+par[3]*data$banned
  if (lfit==6) r = par[1]+par[2]*data$dateb+par[3]*data$hcsap+par[4]*data$assault
  if (lfit==7) r = par[1]+par[2]*data$dateb+par[3]*data$banned+par[4]*data$none
  p = exp(r)/(1+exp(r))
  return(p)
}

##################################################################################
##################################################################################
##################################################################################
negll_trunc_logseries  = function(par,mylist){

  thedata = mylist[[1]]
  trunc = mylist[[2]]

  p=logseries_fun(par,mylist)
  vprob = dtrunc_logseries(thedata$yobs,p,trunc)

  neg_log_like = sum(-log(vprob))
  return(neg_log_like)
}


##################################################################################
##################################################################################
##################################################################################


##################################################################################
##################################################################################
# read in data from github repository
##################################################################################
read_in_data_github = function(input){

  fname = "https://raw.githubusercontent.com/smtowers/data/master/Towers_et_al_public_mass_shootings_Sep_1994_to_Dec_2018.csv"
  thetable = read.table(fname,header=T,as.is=T,sep=",")
  return(thetable)
}

select_data = function(thetable,input){

  if (input$excluded_events!=""){
    a = input$excluded_events
    vdate = unlist(strsplit(a,";"))
    for (i in 1:length(vdate)){
      v = as.numeric(unlist(strsplit(vdate[i],"/")))
      if (length(v)==3&!is.na(sum(v))){
        thetable = subset(thetable,!(month==v[1]&day==v[2]&year==v[3]))
      }
    }
  }

  icount = 0
  if (input$included_events!=""){
    a = input$included_events
    vdate = unlist(strsplit(a,";"))
    for (i in 1:length(vdate)){
      v = as.numeric(unlist(strsplit(vdate[i],"/")))
      if (length(v)==3&!is.na(sum(v))){
        subtable = subset(thetable,(month==v[1]&day==v[2]&year==v[3]))
        if (icount==0){
          specially_included_table = subtable
        }else{
          specially_included_table = rbind(specially_included_table,subtable)
        }
        icount = 1
      }
    }
  }


  thetable = subset(thetable,year>=input$year_range[1]&year<=input$year_range[2])
  if (input$num_casualties=="no_sel"){
    thetable = subset(thetable,num_gun_fatalities_public>=input$min_num_casualties)
  }else{
    thetable = subset(thetable,num_shot_public>=input$min_num_casualties)
  }

  #if (input$include_only_meet_definition){
    thetable = subset(thetable,meet_inclusion_criteria=="Yes")
  #}
  if (input$venue=="work"){
    thetable = subset(thetable,grepl("ork",venue))
  }
  if (input$venue=="not_work"){
    thetable = subset(thetable,!grepl("ork",venue))
  }
  if (input$venue=="school"){
    thetable = subset(thetable,grepl("chool",venue))
  }
  if (input$venue=="public"){
    thetable = subset(thetable,!grepl("chool",venue)&!grepl("ork",venue))
  }
  if (input$perp_knew_victims=="all_or_some"){
    thetable = subset(thetable,public_victims_knew_perp=="Yes"|public_victims_knew_perp=="Partly")
  }
  if (input$perp_knew_victims=="none"){
    thetable = subset(thetable,public_victims_knew_perp=="No")
  }

  #######################################################################
  # include the specifically included events
  #######################################################################
  if (icount>0) thetable = rbind(thetable,specially_included_table)

   return(thetable)
}


##################################################################################
##################################################################################
# aggregate data by day
##################################################################################
aggregate_data_by_day = function(thetable,year_min,year_max){
  wjulian_day = seq(julian(1,01,year_min),julian(12,31,year_max))
  wnum = rep(0,length(wjulian_day))
  wnum_known_if_FAWB_weapon_involved = rep(0,length(wjulian_day))
  wnum_FAWB_weapon_involved = rep(0,length(wjulian_day))
  wpopulation = 200180208 + 7312.273*wjulian_day

  g = aggregate(rep(1,nrow(thetable)),by=list(thetable$julian_day),FUN="sum")
  i = which(wjulian_day%in%g[[1]])
  wnum[i] = g[[2]][match(wjulian_day[i],g[[1]])]

  b = subset(thetable,known_if_high_capacity_involved_or_not==1)
  g = aggregate(rep(1,nrow(b)),by=list(b$julian_day),FUN="sum")
  i = which(wjulian_day%in%g[[1]])
  wnum_known_if_FAWB_weapon_involved[i] = g[[2]][match(wjulian_day[i],g[[1]])]

  b = subset(thetable,involved_FAWB_banned_firearm==1)
  g = aggregate(rep(1,nrow(b)),by=list(b$julian_day),FUN="sum")
  i = which(wjulian_day%in%g[[1]])
  wnum_FAWB_weapon_involved[i] = g[[2]][match(wjulian_day[i],g[[1]])]

  date=as.Date(wjulian_day,origin="1970-01-01")
  wyear=as.numeric(format(date,"%Y"))
  wjour=as.numeric(format(date,"%j"))
  wdate = wyear+(wjour-0.5)/365
  lgood = wyear%%4==0
  wdate[lgood] = wyear[lgood]+(wjour[lgood]-0.5)/366

  wdat = data.frame(date=wdate,julian_day=wjulian_day,year=wyear,jour=wjour,num=wnum,num_known_if_FAWB_weapon_involved=wnum_known_if_FAWB_weapon_involved,num_FAWB_weapon_involved=wnum_FAWB_weapon_involved,population=wpopulation)
  wdat$x = wdat$julian_day-julian(9,13,1994)
  wdat$xb = wdat$julian_day-julian(9,13,2004)
  wdat$xb[wdat$xb<0] = 0

  return(wdat)
}

##################################################################################
##################################################################################
# aggregate data by year
##################################################################################
aggregate_data_by_year=function(thetable,year_min,year_max,wdat){
  zyear = seq(year_min,year_max)
  znum = rep(0,length(zyear))
  znum_known_if_FAWB_weapon_involved = rep(0,length(zyear))
  znum_FAWB_weapon_involved = rep(0,length(zyear))
  zpopulation = aggregate(wdat$population,by=list(wdat$year),FUN="mean")[[2]]

  g = aggregate(rep(1,nrow(thetable)),by=list(thetable$year),FUN="sum")
  i = which(zyear%in%g[[1]])
  znum[i] = g[[2]][match(zyear[i],g[[1]])]

  b = subset(thetable,known_if_high_capacity_involved_or_not==1)
  g = aggregate(rep(1,nrow(b)),by=list(b$year),FUN="sum")
  i = which(zyear%in%g[[1]])
  znum_known_if_FAWB_weapon_involved[i] = g[[2]][match(zyear[i],g[[1]])]

  b = subset(thetable,involved_FAWB_banned_firearm==1)
  g = aggregate(rep(1,nrow(b)),by=list(b$year),FUN="sum")
  i = which(zyear%in%g[[1]])
  znum_FAWB_weapon_involved[i] = g[[2]][match(zyear[i],g[[1]])]

  zdat = data.frame(year=(zyear+0.5),num=znum,num_known_if_FAWB_weapon_involved=znum_known_if_FAWB_weapon_involved,num_FAWB_weapon_involved=znum_FAWB_weapon_involved,population=zpopulation)
  zdat$x = zdat$year-0.5-1994
  zdat$xb = zdat$year-0.5-2004
  zdat$xb[zdat$xb<0] = 0

  return(zdat)
}

##################################################################################
##################################################################################
##################################################################################
mycolors = function(){
  a = data.frame(data_color="black"
                ,fit_color="red3"
                ,background_color="cornsilk"
                ,notification_color="purple4"
                ,alternate_color1="deepskyblue4"
                ,alternate_color2="darkgreen"
                #,alternate_color3="magenta3"
                ,alternate_color3="darkorange2"
                ,notification_cex=1.2
                ,apch=20
                ,acex=4
                ,alwd=9
                ,cex_lab=1.5
                ,cex_main=1.5
                ,stringsAsFactors=F)
  return(a)
}

##################################################################################
##################################################################################
# fit to number of incidents per day
##################################################################################
fit_to_number_incidents_per_day = function(wdat,lover_dispersion=T){

  if (lover_dispersion){
    ##############################################################################
    # for non-overdispersed data, Negative Binomial fit has difficulty fitting
    # for over-dispersion coefficient, and covariance matrix from fit
    # is not trustworthy
    # If the AIC of the Negative Binomial fit is greater than that of a Poisson
    # likelihood fit, over-dispersion is not an issue, and the Poisson likelihood
    # is preferred.
    ##############################################################################
    options(warn=(-1))
    mylog=capture.output({fit_pop_time = suppressWarnings(gamlss(num~offset(log(population))+x,family=NBI,data=wdat,control=gamlss.control(trace=F),trace=F))})
    mylog=capture.output({fit_pop_time_pois = glm(num~offset(log(population))+x,family="poisson",data=wdat)})
    options(warn=0)
    if (AIC(fit_pop_time)<AIC(fit_pop_time_pois)){
      sfit=summary(fit_pop_time, save=TRUE)
      inc_per_year = sfit$mu.coef.table[2,1]*365.25
      p_pop_time=sfit$mu.coef.table[2,4]
    }else{
      inc_per_year = fit_pop_time_pois$coef[2]*365.25
      p_pop_time=summary(fit_pop_time_pois)$coef[2,4]
    }
  }else{
    mylog=capture.output({fit_pop_time = glm(num~offset(log(population))+x,family="poisson",data=wdat)})
    inc_per_year = fit_pop_time$coef[2]*365.25
    p_pop_time=summary(fit_pop_time)$coef[2,4]
  }
  percent_rate_inc_per_year = ((exp(inc_per_year)-1))*100
  ypred_per_year = 365.25*predict(fit_pop_time,type="response")

  return(list(percent_rate_inc_per_year=percent_rate_inc_per_year
             ,pvalue_percent_rate_inc_per_year=p_pop_time
             ,date=wdat$date
             ,ypred_per_year=ypred_per_year
             ,fit_pop_time=fit_pop_time
             ))

}

##################################################################################
##################################################################################
# plot number of incidents over time
##################################################################################
plot_incidents_over_time=function(zdat
                                 ,myfit
                                 ,thecolors
                                 ,lprint=F){

  plot(zdat$year,zdat$num,col=thecolors$data_color,cex=thecolors$acex,xlab="Date",ylab="\043 incidents per year",main="\043 incidents per year",pch=thecolors$apch,ylim=c(0,max(zdat$num+3)),cex.lab=thecolors$cex_lab,cex.main=thecolors$cex_main)
  u = par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)
  points(zdat$year,zdat$num,col=thecolors$data_color,cex=thecolors$acex,pch=thecolors$apch)
  lines(myfit$date,myfit$ypred_per_year,col=thecolors$fit_color,lwd=thecolors$alwd)
  
  if (lprint){
    xrange = u[2]-u[1]
    yrange = u[4]-u[3]
    xmin = u[1]
    ymin = u[3]
    p = myfit$pvalue_percent_rate_inc_per_year
    add_string = ""
    #if (p<0.05) add_string = "*"
    #if (p<0.01) add_string = "**"
    #if (p<0.001) add_string = "***"
    x = xmin + 0.05*xrange
    y = ymin + 0.80*yrange 
    text(x,y,paste("Percentage increase per year: ",round(myfit$percent_rate_inc_per_year,0),"%",add_string,sep=""),adj=c(0,1),cex=thecolors$notification_cex,col=thecolors$notification_color,font=2)

    p = return_p_value_string(myfit$pvalue_percent_rate_inc_per_year)
    y = ymin + 0.75*yrange 
    text(x,y,paste("Significance: ",p,sep=""),adj=c(0,1),cex=thecolors$notification_cex,col=thecolors$notification_color,font=2)

  }

  return()

}

##################################################################################
##################################################################################
# fit to fraction involving FAWB banned weaponry before and after ban
##################################################################################
fit_to_fraction_involving_banned_weaponry = function(wdat){

  wdat$ban = rep(1,nrow(wdat))
  wdat$ban[wdat$julian_day>=julian(09,13,2004)] = 0
  wdat$banb = 1-wdat$ban

  ldone = 0
  increase_prob = 0
  pvalue_increase_prob = 1

  if (sum(wdat$ban)!=nrow(wdat)&sum(wdat$ban)!=0&sum(wdat$ban*wdat$num_known_if_FAWB_weapon_involved)>1){
    ldone = 1
    myfit = glm(cbind(wdat$num_FAWB_weapon_involved,wdat$num_known_if_FAWB_weapon_involved-wdat$num_FAWB_weapon_involved)~factor(wdat$banb),family="binomial")

    ypred = predict(myfit,type="response")
    increase_prob = (max(ypred)-min(ypred))*100
    pvalue_increase_prob = summary(myfit)$coef[2,4]
  }else{
    myfit = glm(cbind(wdat$num_FAWB_weapon_involved,wdat$num_known_if_FAWB_weapon_involved-wdat$num_FAWB_weapon_involved)~1,family="binomial")

    ypred = predict(myfit,type="response")
  }
  
  return(list(myfit=myfit
             ,date=wdat$date
             ,ypred=ypred
             ,ldone=ldone
             ,increase_prob=increase_prob
             ,pvalue_increase_prob=pvalue_increase_prob
             ))
}



##################################################################################
##################################################################################
# plot fraction involving FAWB banned weaponry before and after ban
##################################################################################
plot_fraction_involving_banned_weaponry = function(zdat
                                                  ,wdat
                                                  ,myfit
                                                  ,thecolors
                                                  ,lprint=T
                                                  ){

  
  wdat$ban = rep(0,nrow(wdat))
  wdat$ban[wdat$xb>0] = 1

  plot(1,1,xlim=c(-1.0,2.0),ylim=c(0,100),xlab="",ylab="Percentage involving weaponry banned during FAWB",col=0,main="Percentage involving weaponry banned during\n Federal Assault Weapons Ban (FAWB)",xaxt="n",cex.main=thecolors$cex_main,cex.lab=thecolors$cex_lab,yaxs="i")
  u = par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)

  i = which(wdat$ban==0)
  f1 =sum(wdat$num_FAWB_weapon_involved[i])/sum(wdat$num_known_if_FAWB_weapon_involved[i])
  ib = which(wdat$ban==1)
  f2 =sum(wdat$num_FAWB_weapon_involved[ib])/sum(wdat$num_known_if_FAWB_weapon_involved[ib])
  f = c(f1,f2)*100
  acol = thecolors$alternate_color1
  #acol = c(acol,thecolors$notification_color)
  acol = c(acol,thecolors$alternate_color3)

  p = return_p_value_string(myfit$pvalue_increase_prob)
  for (j in 1:2){
    x1 = (j-1)-0.5
    x2 = (j-1)+0.5
    y1 = 0
    y2 = f[j]
    rect(x1,y1,x2,y2
        ,col=acol[j]
        ,border=acol[j]
        )
  }
  adiff = round(f[2]-f[1],0)
  if (!is.na(adiff)){
    if (adiff>0) adiff = paste("+",adiff,sep="")
    text(0.5,95,paste("Change after lapse of FAWB=",adiff,"% (",p,")",sep=""),font=2,cex=thecolors$notification_cex)
  }
  wname = paste("During\n FAWB\n","(N=",sum(wdat$num_known_if_FAWB_weapon_involved[i]),")",sep="")
  wname = c(wname,paste("After\n Lapse\n FAWB\n (N=",sum(wdat$num_known_if_FAWB_weapon_involved[ib]),")",sep=""))
  axis(1,las=2,labels=wname,at=seq(0,1),cex.axis=thecolors$cex_lab,cex.lab=thecolors$cex_lab)




  if (0){
  plot(zdat$year,zdat$num_FAWB_weapon_involved/zdat$num_known_if_FAWB_weapon_involved,col=thecolors$data_color,cex=thecolors$acex,xlab="Date",ylab="Fraction involving weaponry banned during FAWB",main="Fraction involving weaponry\n banned during FAWB",pch=thecolors$apch,ylim=c(0,1))

  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)

  points(zdat$year,zdat$num_FAWB_weapon_involved/zdat$num_known_if_FAWB_weapon_involved,col=thecolors$data_color,cex=thecolors$acex,pch=thecolors$apch)

  lines(myfit$date,myfit$ypred,col=thecolors$fit_color,lwd=thecolors$alwd)
  if (myfit$ldone==0){
    year_mid = mean(zdat$year)
    text(year_mid,0.22,"(Note: insufficient data to compare periods before and after FAWB)",cex=1.0,font=2,col=thecolors$notification_color)

  }
  }

}

##################################################################################
##################################################################################
# fit to the temporal trends in the number of casualties
##################################################################################
fit_temporal_trends_casualties = function(temp
                                         ,min_num_casualties
                                         ,luse_logseries=F){

  mydata = data.frame(yobs = temp$num_casualties)
  mydata$dateb = temp$date-min(temp$date)
  mydata$date = temp$date
  
  mydata$after = rep(0,nrow(mydata))
  mydata$after[temp$julian_day>=julian(9,13,2004)] = 1
  mydata$banned=as.numeric(temp$involved_FAWB_banned_firearm)
  mydata$none=as.numeric(temp$public_victims_knew_perp=="No")
  mydata$assault=as.numeric(temp$involved_assault_rifle)
  mydata$hcsap=as.numeric(temp$involved_high_capacity_handgun)
  mydata = mydata[order(mydata$date),]

  ################################################################################
  ################################################################################
  # log series distribution is limit of very over-dispersed negative binomial
  ################################################################################
  mylist = list()
  mylist[[1]] = mydata
  mylist[[2]] = min_num_casualties
  lfit = 1 # just constant fit
  lfit = 2 # constant+factor(banned)
  lfit = 3 # constant+factor(banned)+factor(none)
  lfit = 4 # constant+date
  lfit = 5 # constant+date+factor(banned)
  lfit = 6 # constant+date+factor(hcsap)+factor(assault)
  lfit = 7 # constant+date+factor(banned)+factor(none)
  
  myfit = list()
  myA = list()
  vlike = numeric(0)

  myfitb = list()
  myAb = list()
  vlikeb = numeric(0)
  for (lfit in 4:4){
    #cat(names(mydata),"\n")
    mylist[[3]] = lfit
    npar = 1
    if (lfit==2) npar = 2
    if (lfit==3) npar = 3
    if (lfit==4) npar = 2
    if (lfit==5) npar = 3
    if (lfit==6) npar = 4
    if (lfit==7) npar = 4


    if (luse_logseries){
      mylog=capture.output({rb = optim(par=rep(0,npar),fn=negll_trunc_logseries,hessian=T,mylist=mylist,control=list(maxit=1000))})
      vlikeb = c(vlikeb,rb$value)
      Ab = solve(rb$hessian)
      myfitb[[lfit]] = rb
      myAb[[lfit]] = Ab
      if (lfit==4){
        mydata$p = logseries_fun(rb$par,mylist)
      }
    }else{
      #cat("Got here\n")
      mylog=capture.output({r = optim(par=rep(0.1,(npar+1)),fn=negll_trunc_negbinom,hessian=T,mylist=mylist,control=list(maxit=1000))})
      vlike = c(vlike,r$value)
      A = solve(r$hessian)
      myfit[[lfit]] = r
      myA[[lfit]] = A
      if (lfit==4){
        mydata$mu = negbinom_fun(r$par,mylist)
      }
    }
  }
  if (luse_logseries){
    pvalue = 1-pchisq((myfitb[[4]]$par[2]/sqrt(myAb[[4]][2,2]))^2,1)
  }else{
    pvalue = 1-pchisq((myfit[[4]]$par[3]/sqrt(myA[[4]][3,3]))^2,1)
  }

  mydata$ypred = rep(0,nrow(mydata))
  x = seq(min_num_casualties,1000)
  for (i in 1:nrow(mydata)){
    if (luse_logseries){
      y = dtrunc_logseries(x,mydata$p[i],min_num_casualties)
    }else{
      y = dtrunc_negbinom(x,mydata$mu[i],exp(myfit[[4]]$par[1]),min_num_casualties)
    }
    ymean = weighted.mean(x,y)
    mydata$ypred[i] = ymean
  }

  # work out percentage increase by year
  yend= mydata$ypred[1]
  ybeg= mydata$ypred[nrow(mydata)]
  diff = mydata$date[1]-mydata$date[nrow(mydata)]
  frac_increase_killed = exp(log(yend/ybeg)/diff)-1
  per_increase_casualties = (frac_increase_killed)*100

 
  return(list(myfit=myfit
             ,myA=myA
             ,myfitb=myfitb
             ,myAb=myAb
             ,date=mydata$date
             ,yobs=mydata$yobs
             ,ypred=mydata$ypred
             ,per_increase_casualties=per_increase_casualties
             ,p_value = pvalue
             )
         )
}

##################################################################################
##################################################################################
##################################################################################
plot_number_casualties_over_time = function(thetable
                                           ,myfit
                                           ,input
                                           ,thecolors
                                           ,lprint=F
                                           ){

  xmin = input$year_range[1]
  xmax = input$year_range[2]+1
  plot(thetable$date,thetable$num_casualties,xlab="Date",ylab="\043 casualties per incident",cex=thecolors$acex,pch=thecolors$apch,col=thecolors$data_color,ylim=c(0,max(thetable$num_casualties+5)),xlim=c(xmin,xmax),main="\043 casualties per incident over time",cex.lab=thecolors$cex_lab,cex.main=thecolors$cex_main)

  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)
  points(thetable$date,thetable$num_casualties,cex=thecolors$acex,pch=thecolors$apch,col=thecolors$data_color)
  #lines(myfit$date,myfit$ypred,col=thecolors$fit_color,lwd=thecolors$alwd)
  i = order(myfit$date)
  a = smooth.spline(myfit$date[i],myfit$ypred[i],all.knots=T,df=length(myfit$date))
  x = seq(xmin,xmax,length=100)
  b = predict(a,x)
  lines(b$x,b$y,col=thecolors$fit_color,lwd=thecolors$alwd)
  #lines(b$x,b$y,col="blue",lwd=2)


  if (lprint){
    xrange = u[2]-u[1]
    yrange = u[4]-u[3]
    xmin = u[1]
    ymin = u[3]
    p = myfit$p_value
    add_string = ""
    #if (p<0.05) add_string = "*"
    #if (p<0.01) add_string = "**"
    #if (p<0.001) add_string = "***"
    x = xmin + 0.05*xrange
    y = ymin + 0.80*yrange
    text(x,y,paste("Percentage increase per year: ",round(myfit$per_increase_casualties,0),"%",add_string,sep=""),adj=c(0,1),cex=thecolors$notification_cex,col=thecolors$notification_color,font=2)

    p = return_p_value_string(myfit$p_value)
    y = ymin + 0.75*yrange 
    text(x,y,paste("Significance: ",p,sep=""),adj=c(0,1),cex=thecolors$notification_cex,col=thecolors$notification_color,font=2)

  }

}


##################################################################################
##################################################################################
##################################################################################
violin_plot_casualties = function(thetable
                                  ,thecolors
                                  ,lprint=F
                                  ){

  m = subset(thetable,known_if_high_capacity_involved_or_not==1)
  m$x = m$num_casualties

  m$lperiod = rep(0,nrow(m))
  m$lperiodb = rep(0,nrow(m))
  i = which(m$involved_FAWB_banned_firearm==0)
  m$lperiod[i] = 0
  i = which(m$involved_FAWB_banned_firearm==1)
  m$lperiod[i] = 1
  i = which(m$involved_FAWB_banned_firearm==1&m$involved_high_capacity_handgun==0&m$involved_assault_rifle==1)
  m$lperiodb[i] = 2
  i = which(m$involved_FAWB_banned_firearm==1&m$involved_high_capacity_handgun==1&m$involved_assault_rifle==0)
  m$lperiodb[i] = 3

  plot(1,1,xlim=c(-0.5,3.5),ylim=c(0,max(m$x+5)),xlab="",ylab="Distribution \043 casualties per incident\n (bar width represents \043 incidents w/ that many casualties)",col=0,main="Distribution \043 casualties per incident",xaxt="n",cex.main=thecolors$cex_main,cex.lab=thecolors$cex_lab,yaxs="i")
  u = par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)
  lcustom_h = 1
  if (lcustom_h){
    ah = 0.5
    vioplot(m$x[m$lperiod==0],at=0,add=T,col=thecolors$fit_color,h=ah,wex=0.1)
    vioplot(m$x[m$lperiodb==3],at=1,add=T,col=thecolors$alternate_color2,h=ah)
    vioplot(m$x[m$lperiod==1],at=2,add=T,col=thecolors$notification_color,h=ah)
    vioplot(m$x[m$lperiodb==2],at=3,add=T,col=thecolors$alternate_color1,h=ah)
  }else{
    vioplot(m$x[m$lperiod==0],at=0,add=T,col=thecolors$fit_color)
    vioplot(m$x[m$lperiodb==3],at=1,add=T,col=thecolors$alternate_color2)
    vioplot(m$x[m$lperiod==1],at=2,add=T,col=thecolors$notification_color)
    vioplot(m$x[m$lperiodb==2],at=3,add=T,col=thecolors$alternate_color1)
  }
  vname = "No banned weaponry"
  vname = c(vname,"High-capacity pistol")
  vname = c(vname,"High-capacity rifle")
  vname = c(vname,"Any banned weaponry")
  acol = character(0)
  acol = c(acol,thecolors$fit_color)
  acol = c(acol,thecolors$alternate_color1)
  acol = c(acol,thecolors$notification_color)
  acol = c(acol,thecolors$alternate_color2)
  legend("topleft",legend=vname,col=acol,lwd=7,bty="n",cex=1.0)
  #axis(1,las=2,labels=vname,at=seq(0,3),cex.axis=0.8)
}

##################################################################################
##################################################################################
##################################################################################
battleship_plot_casualties = function(thetable
                                     ,thecolors
                                     ,input
                                     ,lprint=F
                                     ){

  m = subset(thetable,known_if_high_capacity_involved_or_not==1)
  m$y = m$num_casualties

  ################################################################################
  ################################################################################
  ################################################################################
  m$only_pistol = rep(0,nrow(m))
  m$only_rifle = rep(0,nrow(m))
  m$pistol_and_rifle = rep(0,nrow(m))
  m$none_involved = rep(0,nrow(m))
  m$any = rep(0,nrow(m))

  m$none_involved[m$involved_FAWB_banned_firearm==0] = 1
  m$only_pistol[m$involved_FAWB_banned_firearm==1&
                m$involved_high_capacity_handgun==1&
                m$involved_assault_rifle==0] = 1
  m$only_rifle[m$involved_FAWB_banned_firearm==1&
               m$involved_high_capacity_handgun==0&
               m$involved_assault_rifle==1] = 1
  m$pistol_and_rifle[m$involved_FAWB_banned_firearm==1&
                     m$involved_high_capacity_handgun==1&
                     m$involved_assault_rifle==1] = 1
  m$any[m$involved_FAWB_banned_firearm==1] = 1

  vres = numeric(0)
  evres = numeric(0)
  vp = numeric(0)
  viter = numeric(0)
  for (iter in 1:5){
    mydata = data.frame(only_pistol = m$only_pistol
                       ,only_rifle = m$only_rifle
                       ,pistol_and_rifle = m$pistol_and_rifle
                       ,yobs=m$y
                       )
    if (iter==1) mydata = subset(mydata,only_pistol==0&only_rifle==0&pistol_and_rifle==0)
    if (iter==1) mydata_null = mydata
    if (iter==2) mydata = subset(mydata,only_pistol==1&only_rifle==0&pistol_and_rifle==0)
    if (iter==3) mydata = subset(mydata,only_pistol==0&only_rifle==1&pistol_and_rifle==0)
    if (iter==4) mydata = subset(mydata,only_pistol==0&only_rifle==0&pistol_and_rifle==1)
    if (iter==5) mydata = subset(mydata,only_pistol==1|only_rifle==1|pistol_and_rifle==1)

    mylist = list()
    mylist[[1]] = mydata
    mylist[[2]] = input$min_num_casualties
    mylist[[3]] = 1 # type of fit
    npar = 2
  
    if (nrow(mydata)>1){
      #if (lfit==8) r = par[2]+par[3]*data$only_pistol+par[4]*data$only_rifle+par[5]*data$pistol_and_rifle
      #mylog=capture.output({r = optim(par=rep(0.1,npar),fn=negll_trunc_negbinom,hessian=T,mylist=mylist,control=list(maxit=1000))})
      #mydata$mu = negbinom_fun(r$par,mylist)
    
      #negll_trunc_logseries(rep(0.1,(npar-1)),mylist)
      mylog=capture.output({r = optim(par=rep(0.1,(npar-1)),fn=negll_trunc_logseries,hessian=T,mylist=mylist,control=list(maxit=1000),method="Brent",lower=(0),upper=(+20))})
      #mylog=capture.output({r = optimize(f=negll_trunc_logseries,mylist=mylist,lower=(-0),upper=(+10))})
      A = solve(r$hessian)
      mydata$mu = logseries_fun(r$par,mylist)
      vres = c(vres,r$par[1])
      evres = c(evres,sqrt(A[1,1]))
      viter = c(viter,iter)
      #a = ks.test(mydata$yobs,mydata_null$yobs,alternative="less")
      #vp = c(vp,a$p.value)
    }
  }
  z = (vres-vres[viter==1])/sqrt(evres^2+evres[viter==1]^2)
  p = 1-pchisq(z^2,1)
  p = return_p_value_string(p)
  if (length(p)>0) vp = return_p_value_string(vp)

  ################################################################################
  ################################################################################
  ################################################################################
  breaks = seq((input$min_num_casualties-0.5),max(m$y+0.5),1)
  a = hist(m$y[m$none_involved==1],plot=F,breaks=breaks)
  b = hist(m$y[m$only_pistol==1],plot=F,breaks=breaks)
  c = hist(m$y[m$only_rifle==1],plot=F,breaks=breaks)
  d = hist(m$y[m$pistol_and_rifle==1],plot=F,breaks=breaks)
  e = hist(m$y[m$any==1],plot=F,breaks=breaks)
  vN = c(sum(a$counts),sum(b$counts),sum(c$counts),sum(d$counts),sum(e$counts))
  vN = paste("(N=",vN,")",sep="")
  lmax = rep(0,4)
  if (sum(a$counts)>0) lmax[1] = max(breaks[which(a$counts>0)])
  if (sum(b$counts)>0) lmax[2] = max(breaks[which(b$counts>0)])
  if (sum(c$counts)>0) lmax[3] = max(breaks[which(c$counts>0)])
  if (sum(d$counts)>0) lmax[4] = max(breaks[which(d$counts>0)])
  if (sum(e$counts)>0) lmax[5] = max(breaks[which(e$counts>0)])
  ymax = max(c(a$counts,b$counts,c$counts,d$counts,e$counts))
  a$counts = 0.90*a$counts/ymax
  b$counts = 0.90*b$counts/ymax
  c$counts = 0.90*c$counts/ymax
  d$counts = 0.90*d$counts/ymax
  e$counts = 0.90*e$counts/ymax

  acol = character(0)
  acol = c(acol,thecolors$fit_color)
  acol = c(acol,thecolors$alternate_color1)
  acol = c(acol,thecolors$notification_color)
  acol = c(acol,thecolors$alternate_color2)
  acol = c(acol,thecolors$alternate_color3)

  shift = 0.35*max(m$y)
  bshift = 0.05*max(m$y)
  if (shift<5) shift=5
  if (bshift<2) bshift=2
  #cat(shift,bshift,"\n")
  plot(1,1,xlim=c(-0.5,4.5),ylim=c((input$min_num_casualties-0.0),max(m$y+shift)),xlab="",ylab="Distribution \043 casualties per incident\n (width of bars represents \043 of incidents w/ that many casualties)",col=0,main="Distribution \043 casualties per incident\n and dependence on use of FAWB banned weaponry",xaxt="n",cex.main=thecolors$cex_main,cex.lab=thecolors$cex_lab,yaxs="i")
  u = par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)

  for (j in 1:5){
    if (j==1) temp=a
    if (j==2) temp=b
    if (j==3) temp=c
    if (j==4) temp=d
    if (j==5) temp=e
    amax = lmax[j]+0.5
    mycol = acol[j]
    lapply(seq(1,length(a$mids))
     ,function(i,amax,temp,mycol){
        if (temp$counts[i]>=0&temp$mids[i]<=amax){
          rect((j-1)-0.5*temp$counts[i]
              ,temp$mids[i]
              ,(j-1)+0.5*temp$counts[i]
              ,temp$mids[i]+1.0
              ,col=mycol
              ,border=mycol
              )
        }
      }
     ,amax=amax,temp=temp,mycol=mycol)
    if (j>1&sum(viter==j)==1){
      text((j-1),lmax[j]+bshift,p[viter==j])
    }
  }

  vname = "No banned weaponry"
  vname = c(vname,"High-capacity pistol")
  vname = c(vname,"High-capacity rifle")
  vname = c(vname,"High-capacity pistol and rifle")
  vname = c(vname,"Any banned high-capacity (including shotguns)")
  
  legend("topleft",legend=vname,col=acol,lwd=7,bty="n",cex=1.0)
  wname = paste("No\n banned\n weaponry\n",vN[1])
  wname = c(wname,paste("High-cap\n Pistol\n",vN[2]))
  wname = c(wname,paste("High-cap\n Rifle\n",vN[3]))
  wname = c(wname,paste("High-cap\n Pistol and\n Rifle\n",vN[4]))
  wname = c(wname,paste("Any\n banned\n",vN[5]))
  axis(1,las=2,labels=wname,at=seq(0,4),cex.axis=1.2)

  return()

}

##################################################################################
##################################################################################
##################################################################################
histogram_plot_casualties = function(thetable
                                    ,thecolors
                                    ,input
                                    ,lprint=F
                                    ){

  m = subset(thetable,known_if_high_capacity_involved_or_not==1)
  m$y = m$num_casualties

  ################################################################################
  ################################################################################
  ################################################################################
  m$only_pistol = rep(0,nrow(m))
  m$only_rifle = rep(0,nrow(m))
  m$pistol_and_rifle = rep(0,nrow(m))
  m$none_involved = rep(0,nrow(m))
  m$any = rep(0,nrow(m))

  m$none_involved[m$involved_FAWB_banned_firearm==0] = 1
  m$only_pistol[m$involved_FAWB_banned_firearm==1&
                m$involved_high_capacity_handgun==1&
                m$involved_assault_rifle==0] = 1
  m$only_rifle[m$involved_FAWB_banned_firearm==1&
               m$involved_high_capacity_handgun==0&
               m$involved_assault_rifle==1] = 1
  m$pistol_and_rifle[m$involved_FAWB_banned_firearm==1&
                     m$involved_high_capacity_handgun==1&
                     m$involved_assault_rifle==1] = 1
  m$any[m$involved_FAWB_banned_firearm==1] = 1

  vres = numeric(0)
  evres = numeric(0)
  vp = numeric(0)
  viter = numeric(0)
  for (iter in 1:5){
    mydata = data.frame(only_pistol = m$only_pistol
                       ,only_rifle = m$only_rifle
                       ,pistol_and_rifle = m$pistol_and_rifle
                       ,yobs=m$y
                       )
    if (iter==1) mydata = subset(mydata,only_pistol==0&only_rifle==0&pistol_and_rifle==0)
    if (iter==1) mydata_null = mydata
    if (iter==2) mydata = subset(mydata,only_pistol==1&only_rifle==0&pistol_and_rifle==0)
    if (iter==3) mydata = subset(mydata,only_pistol==0&only_rifle==1&pistol_and_rifle==0)
    if (iter==4) mydata = subset(mydata,only_pistol==0&only_rifle==0&pistol_and_rifle==1)
    if (iter==5) mydata = subset(mydata,only_pistol==1|only_rifle==1|pistol_and_rifle==1)

    mylist = list()
    mylist[[1]] = mydata
    mylist[[2]] = input$min_num_casualties
    mylist[[3]] = 1 # type of fit
    npar = 2
  
    if (nrow(mydata)>1){
      #if (lfit==8) r = par[2]+par[3]*data$only_pistol+par[4]*data$only_rifle+par[5]*data$pistol_and_rifle
      #mylog=capture.output({r = optim(par=rep(0.1,npar),fn=negll_trunc_negbinom,hessian=T,mylist=mylist,control=list(maxit=1000))})
      #mydata$mu = negbinom_fun(r$par,mylist)
    
      #negll_trunc_logseries(rep(0.1,(npar-1)),mylist)
      mylog=capture.output({r = optim(par=rep(0.1,(npar-1)),fn=negll_trunc_logseries,hessian=T,mylist=mylist,control=list(maxit=1000),method="Brent",lower=(0),upper=(+20))})
      #mylog=capture.output({r = optimize(f=negll_trunc_logseries,mylist=mylist,lower=(-0),upper=(+10))})
      A = solve(r$hessian)
      mydata$mu = logseries_fun(r$par,mylist)
      vres = c(vres,r$par[1])
      evres = c(evres,sqrt(A[1,1]))
      viter = c(viter,iter)
      a = ks.test(mydata$yobs,mydata_null$yobs,alternative="less")
      vp = c(vp,a$p.value)
    }
  }
  z = (vres-vres[viter==1])/sqrt(evres^2+evres[viter==1]^2)
  p = 1-pchisq(z^2,1)
  p = return_p_value_string(p)
  vp = return_p_value_string(vp)

  ################################################################################
  ################################################################################
  ################################################################################
  breaks = seq((input$min_num_casualties+0.5),max(m$y+0.5),1)
  a = hist(m$y[m$none_involved==1],plot=F,breaks=breaks,freq=T)
  b = hist(m$y[m$only_pistol==1],plot=F,breaks=breaks,freq=T)
  c = hist(m$y[m$only_rifle==1],plot=F,breaks=breaks,freq=T)
  d = hist(m$y[m$pistol_and_rifle==1],plot=F,breaks=breaks,freq=T)
  e = hist(m$y[m$any==1],plot=F,breaks=breaks,freq=T)
  vN = c(sum(a$counts),sum(b$counts),sum(c$counts),sum(d$counts),sum(e$counts))
  vN = paste("(N=",vN,")",sep="")

  a$counts = 1.00*a$counts/sum(a$counts)
  b$counts = 1.00*b$counts/sum(b$counts)
  c$counts = 1.00*c$counts/sum(c$counts)
  d$counts = 1.00*d$counts/sum(d$counts)
  e$counts = 1.00*e$counts/sum(e$counts)

  acol = character(0)
  acol = c(acol,thecolors$fit_color)
  acol = c(acol,thecolors$alternate_color1)
  acol = c(acol,thecolors$notification_color)
  acol = c(acol,thecolors$alternate_color2)
  acol = c(acol,thecolors$alternate_color3)

  bcol = character(0)
  for (i in 1:length(acol)){
    mycol=rgb(t(col2rgb(acol[i])),max=255,alpha=125)
    bcol = c(bcol,mycol)
  }
  bcol = acol
  plot(a$mids,a$counts,col=bcol[1],type="l",lwd=5)
  u = par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)
  lines(a$mids,a$counts,col=bcol[1],lwd=5)
  lines(b$mids,b$counts,col=bcol[2],lwd=5)
  lines(c$mids,c$counts,col=bcol[3],lwd=5)
  lines(d$mids,d$counts,col=bcol[4],lwd=5)
  lines(e$mids,e$counts,col=bcol[5],lwd=5)

  vname = "No banned weaponry"
  vname = c(vname,"High-capacity pistol")
  vname = c(vname,"High-capacity rifle")
  vname = c(vname,"High-capacity pistol and rifle")
  vname = c(vname,"Any banned high-capacity (including shotguns)")
  
  legend("topleft",legend=vname,col=acol,lwd=7,bty="n",cex=0.9)
  wname = paste("No\n banned\n weaponry\n",vN[1])
  wname = c(wname,paste("High-cap\n Pistol\n",vN[2]))
  wname = c(wname,paste("High-cap\n Rifle\n",vN[3]))
  wname = c(wname,paste("High-cap\n Pistol and\n Rifle\n",vN[4]))
  wname = c(wname,paste("Any\n banned\n",vN[5]))
  axis(1,las=2,labels=wname,at=seq(0,4),cex.axis=1.2)
  return()


}





