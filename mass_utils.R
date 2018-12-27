

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
# fit to number of incidents per day
##################################################################################
fit_to_number_incidents_per_day = function(wdat,lover_dispersion){
  mylog=capture.output({pois_pop <- glm(num~offset(log(population)),family="poisson",data=wdat)})
  mylog=capture.output({pois_pop_time <- glm(num~offset(log(population))+x,family="poisson",data=wdat)})
  mylog=capture.output({pois_pop_time_ban <- glm(num~offset(log(population))+x+xb,family="poisson",data=wdat)})
  inc_per_year = pois_pop_time$coef[2]*365.25
  p_pop_time=summary(pois_pop_time)$coef[2,4]
  p_pop_time_ban1=summary(pois_pop_time_ban)$coef[2,4]
  p_pop_time_ban2=summary(pois_pop_time_ban)$coef[3,4]
  y_pop_time = predict(pois_pop_time,type="response")

  if (lover_dispersion){
    mylog=capture.output({negbinom_pop <- gamlss(num~offset(log(population)),family=NBI,data=wdat)})
    mylog=capture.output({negbinom_pop_time <- gamlss(num~offset(log(population))+x,family=NBI,data=wdat)})
    mylog=capture.output({negbinom_pop_time_ban <- gamlss(num~offset(log(population))+x+xb,family=NBI,data=wdat)})
    sfit=summary(negbinom_pop_time, save=TRUE)
    p_pop_time=sfit$mu.coef.table[2,4]
    inc_per_year = sfit$mu.coef.table[2,1]*365.25
    sfit=summary(negbinom_pop_time_ban, save=TRUE)
    p_pop_time_ban1=sfit$mu.coef.table[2,4]
    p_pop_time_ban2=sfit$mu.coef.table[3,4]
    y_pop_time = predict(negbinom_pop_time,type="response")
  }
  #diff_time = (max(wdat$x)-min(wdat$x))/365.25
  #ybeg = y_pop_time[which.min(wdat$x)]
  #yend = y_pop_time[which.max(wdat$x)]
  #rate_inc = exp(log(yend/ybeg)/diff_time)-1
  rate_inc = (exp(inc_per_year)-1)

  return(list(rate_inc=rate_inc
             ,p_pop_time=p_pop_time
             ,pois_pop=pois_pop
             ,pois_pop_time=pois_pop_time
             ,pois_pop_time_ban=pois_pop_time_ban
             ,negbinom_pop=negbinom_pop
             ,negbinom_pop_time=negbinom_pop_time
             ,negbinom_pop_time_ban=negbinom_pop_time_ban))

}


##################################################################################
##################################################################################
# fit to the temporal trends in the number of casualties
##################################################################################
fit_temporal_trends_casualties = function(temp,min_num_killed){

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
  myfitb = list()
  myAb = list()
  vlike = numeric(0)
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
    mylog=capture.output({rb = optim(par=rep(0,npar),fn=negll_trunc_logseries,hessian=T,mylist=mylist,control=list(maxit=1000))})
    vlike = c(vlike,r$value)
    vlikeb = c(vlikeb,rb$value)
    A = solve(r$hessian)
    Ab = solve(rb$hessian)
    myfit[[lfit]] = r
    myA[[lfit]] = A
    myfitb[[lfit]] = rb
    myAb[[lfit]] = Ab
    if (lfit==4){
      mydata$p = logseries_fun(rb$par,mylist)
      mydata$mu = negbinom_fun(r$par,mylist)
    }
  }
  pvalue = 1-pchisq((myfit[[4]]$par[3]/sqrt(myA[[4]][3,3]))^2,1)

  mydata$ypred = rep(0,nrow(mydata))
  x = seq(min_num_killed,1000)
  for (i in 1:nrow(mydata)){
    if (1){
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




