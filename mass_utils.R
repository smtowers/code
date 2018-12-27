
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
  data = mylist[[1]]
  trunc = mylist[[2]]
  lfit = mylist[[3]]

  r = par[2]
  if (lfit==2) r = par[2]+par[3]*data$banned
  if (lfit==3) r = par[2]+par[3]*data$banned+par[4]*data$none
  if (lfit==4) r = par[2]+par[3]*data$dateb
  if (lfit==5) r = par[2]+par[3]*data$dateb+par[4]*data$banned
  if (lfit==6) r = par[2]+par[3]*data$dateb+par[4]*data$hcsap+par[5]*data$assault
  if (lfit==7) r = par[2]+par[3]*data$dateb+par[4]*data$banned+par[5]*data$none
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
    thetable = subset(thetable,num_gun_fatalities_public>=input$min_num_killed)
  }else{
    thetable = subset(thetable,num_shot_public>=input$min_num_killed)
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
                ,apch=20
                ,acex=4
                ,alwd=9
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
                                 ,thecolors){

  plot(zdat$year,zdat$num,col=thecolors$data_color,cex=thecolors$acex,xlab="Date",ylab="\043 incidents per year",main="\043 incidents per year",pch=thecolors$apch,ylim=c(0,max(zdat$num+2)))
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col = thecolors$background_color, border = thecolors$background_color)
  points(zdat$year,zdat$num,col=thecolors$data_color,cex=thecolors$acex,pch=thecolors$apch)
  lines(myfit$date,myfit$ypred_per_year,col=thecolors$fit_color,lwd=thecolors$alwd)

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
                                                  ,myfit
                                                  ,thecolors
                                                  ){

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

##################################################################################
##################################################################################
# fit to the temporal trends in the number of casualties
##################################################################################
fit_temporal_trends_casualties = function(temp
                                         ,min_num_killed){

  mydata = data.frame(yobs = temp$num_gun_fatalities_public)
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
  mylist[[2]] = min_num_killed
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
    cat(names(mydata),"\n")
    mylist[[3]] = lfit
    npar = 1
    if (lfit==2) npar = 2
    if (lfit==3) npar = 3
    if (lfit==4) npar = 2
    if (lfit==5) npar = 3
    if (lfit==6) npar = 4
    if (lfit==7) npar = 4

    mylog=capture.output({r = optim(par=rep(0.1,(npar+1)),fn=negll_trunc_negbinom,hessian=T,mylist=mylist,control=list(maxit=1000))})
    vlike = c(vlike,r$value)
    A = solve(r$hessian)
    myfit[[lfit]] = r
    myA[[lfit]] = A
    if (lfit==4){
      mydata$mu = negbinom_fun(r$par,mylist)
    }

    if (0){
      mylog=capture.output({rb = optim(par=rep(0,npar),fn=negll_trunc_logseries,hessian=T,mylist=mylist,control=list(maxit=1000))})
      vlikeb = c(vlikeb,rb$value)
      Ab = solve(rb$hessian)
      myfitb[[lfit]] = rb
      myAb[[lfit]] = Ab
      if (lfit==4){
        mydata$p = logseries_fun(rb$par,mylist)
      }
    }
  }
  pvalue = 1-pchisq((myfit[[4]]$par[3]/sqrt(myA[[4]][3,3]))^2,1)

  mydata$ypred = rep(0,nrow(mydata))
  x = seq(min_num_killed,1000)
  for (i in 1:nrow(mydata)){
    if (0){
      y = dtrunc_logseries(x,mydata$p[i],min_num_killed)
    }else{
      y = dtrunc_negbinom(x,mydata$mu[i],exp(myfit[[4]]$par[1]),min_num_killed)
    }
    ymean = weighted.mean(x,y)
    mydata$ypred[i] = ymean
  }

  # work out percentage increase by year
  yend= mydata$ypred[1]
  ybeg= mydata$ypred[nrow(mydata)]
  diff = mydata$date[1]-mydata$date[nrow(mydata)]
  frac_increase_killed = exp(log(yend/ybeg)/diff)-1
  per_increase_killed = (frac_increase_killed)*100

 
  return(list(myfit=myfit
             ,myA=myA
             ,myfitb=myfitb
             ,myAb=myAb
             ,date=mydata$date
             ,yobs=mydata$yobs
             ,ypred=mydata$ypred
             ,per_increase_killed=per_increase_killed
             ,p_value = pvalue
             )
         )


}




