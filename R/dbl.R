# Lift Chart Function
# This function produces single and double lift charts
# @param ds dataframe used to develop lift chart
# @param champ incumbent prediction (set variable to uniformly 1.00 for single lift chart)
# @param chall prediction from revised models
# @param target response variable
# @param wate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param bkt number of quantile to display on lift chart
# @param index_tf 0 - no rescaling, 1 = y-axis variables scaled to average to 1
# @param graph_tf 0 - return table, 1 = return plot
# @return contains either a table or a plot containing the lift chart content
dbl <- function(ds , champ, chall, target, wate,bkt,indx_tf,graph_tf){
  a<-subset(ds,select=c(champ,chall,target,wate))%>%rename( prd0=1, prd1=2, tar=3, wgt=4) %>% mutate( chg = as.numeric(round(prd1/prd0,3)), grp=1)
  a<-a%>%arrange(grp,chg)%>%group_by(grp)%>%mutate(qnt = ceiling( bkt * (cumsum(wgt)/sum(wgt)) ), prd0=prd0*wgt, prd1=prd1*wgt,tar=tar*wgt) %>% group_by(qnt) %>% summarise_at(c('wgt','prd0','prd1','tar'),sum)  %>% mutate( prd0=round(prd0/wgt,3) , prd1=round(prd1/wgt,3), tar=round(tar/wgt,3) )
  if (indx_tf == 1) { a<-a%>%mutate( prd0=sum(wgt)*prd0/crossprod(wgt,prd0) , prd1=sum(wgt)*prd1/ crossprod(wgt,prd1) , tar=sum(wgt)*tar/crossprod(wgt,tar) ) }
  a<-a%>%mutate(chg=prd1/prd0)%>%select(qnt,wgt,chg,prd0,prd1,tar)
  a1<-a%>%select(qnt,prd0)%>%rename(prd=2)%>%mutate(metric='champion',prd=round(prd,2))
  a2<-a%>%select(qnt,prd1)%>%rename(prd=2)%>%mutate(metric='challenger',prd=round(prd,2))
  a3<-a%>%select(qnt,tar)%>%rename(prd=2)%>%mutate(metric='target',prd=round(prd,2))
  a0<-bind_rows(a1,a2,a3)%>%mutate(qnt=sprintf("%03d",qnt))
  aa<-ggplot(data=a0, aes(x=qnt, y=prd, fill=metric,color=metric)) + scale_colour_manual(values = c('black','purple','blue')) + scale_fill_manual(values = c('black','purple','blue')) + geom_bar(data=filter(a0,metric %in% c("target")), aes(x = qnt, y = prd) , stat ="identity", position="dodge") + geom_text(aes(label=prd), position=position_dodge(width=0.9), hjust=1, vjust=-0.25, size=3)  + geom_line(data=filter(a0,metric %in% c("champion", "challenger")),aes(x = qnt, y = prd, group=metric, color=metric), size=2) +   theme(axis.text.x = element_text(angle = 45),legend.position="bottom") + ggtitle("double lift chart") +  theme(plot.title = element_text(hjust = 0.5)) + xlab("quantile") + ylab("target")
  if (graph_tf == 1) { return(aa)} else {return(a)}
  gc() }

# Aggregated Component Plus Residual Plot
# This function component plus residual plot type construct
# @param ds dataframe used to develop lift chart
# @param var component variable value
# @param chall prediction from GLM
# @param target response variable
# @param wate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param bkt number of quantile to display on CPR type plot
# @param graph_tf 0 - return table, 1 = return plot
# @return contains either a table or a plot containing the CPR type plot
myCpr<-function(ds,var,target,chall,wate,beta,bkt,graph_tf){
  d<-subset(ds,select=c(var,target,chall,wate))%>%select(v=var,t=target,p=chall,w=wate)%>%arrange(v)
  d<-d%>%mutate(v=exp(v)*w, t=t*w, p=p*w) %>%mutate(c=ceiling( (cumsum(w)/sum(w)) / (1/bkt)))%>%group_by(c)%>%summarise(v=sum(v)/sum(w),t=sum(t)/sum(w),p=sum(p)/sum(w),w=sum(w))
  d <-d%>%mutate(y=log(t)-log(p)+beta*log(v), x=log(v), f=beta*log(v))
  if (graph_tf ==1) { d<- ggplot(data=d, aes(x=x, y=y, group=1)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + geom_line(aes(x=x,y=beta*x)) + theme_bw() + xlab(var) + ylab(target)  +   theme(axis.text.x = element_text(angle = 45),legend.position="bottom")   }
  return(d) }

# Bucketing function for EDA
# This function converts a continuous variable into N equally weighted buckets
# @param ds dataframe containing variable
# @param var variable to put into buckets
# @param newvar name for bucketed variable
# @param wate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param bkt number of quantiles to place variable into
# @return dataframe with column appended
contbkt<-function(ds,var,newvar,wate,bkt) {
  x<- subset( ds, select=c(var,wate))  %>% rename(v=1,w=2)
  x <- x  %>% filter(!is.na(v)) %>% arrange(v,w) %>% mutate(w=pmax(.00001, w)) %>% mutate( b = ceiling( (cumsum(w)/sum(w)) / (1/bkt) ) , vw=w*v)
  x2 <- x %>% group_by(v,b) %>% summarise(w=sum(w)) %>% arrange(v,desc(w)) %>% group_by(v) %>% mutate(ct=1,ct=cumsum(ct)) %>% filter(ct==1) %>% select(v,b)
  x <- left_join(x%>%select(-b),x2,by='v') %>% group_by(b) %>% mutate(vw=sum(vw)/sum(w)) %>% group_by( v, vw) %>% summarise
  x$vw<- as.character(x$vw)
  x<-setnames(x,old=c('v'),new=c(var))
  y<-left_join(ds,x,by=var) %>% mutate (vw=ifelse(is.na(vw), 'z', vw) )
  y<-setnames(y,old=c('vw'),new=c( newvar ) )
  return(y)
  gc() }



# Univariate function for continuous
# This function converts a continuous variable into N equally weighted buckets and produces univariate statistics
# @param ds dataframe containing variable
# @param v variable to put into buckets
# @param w weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param t target variable
# @param q number of quantiles to place variable into
# @return dataframe with n buckets and average target per bucket

uni<-function(ds,v,w,t,q) {
  x<-ds%>%select(v,w,t)%>%rename(V=1,W=2,T=3)
  x1<-x%>%filter(is.na(V))%>%mutate(Q=9999)%>%group_by(Q)%>%summarise(P=sum(W),T=crossprod(T,W)/sum(W))
  x2<-x%>%filter(!is.na(V))%>%arrange(V)%>%mutate(Q=ceiling(cumsum(W)/sum(W)/(1/q)))%>%group_by(Q)%>%summarise(P=sum(W),L=min(V),H=max(V),V=crossprod(V,W)/sum(W),T=crossprod(T,W)/sum(W))
  x<-bind_rows(x2,x1)
  x$T <- round( x$T / as.vector((crossprod(x$P,x$T)/sum(x$P) )), 3)
  return(x) }


# Univariate function for categorical
# This function produces univariate statistics for categorical variable
# @param ds dataframe containing variable
# @param v variable to put into buckets
# @param w weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param t target variable
# @return dataframe with n buckets and average target per bucket

unic<-function(ds,v,w,t) {
  x<-ds%>%select(v,w,t)%>%rename(V=1,W=2,T=3)
  x1<-x%>%filter(is.na(V))%>%mutate(V='Z_NA')%>%group_by(V)%>%summarise(P=sum(W),T=crossprod(T,W)/sum(W))
  x2<-x%>%filter(!is.na(V))%>%mutate(V=as.character(V))%>%group_by(V)%>%summarise(P=sum(W),T=crossprod(T,W)/sum(W))
  x<-bind_rows(x2,x1)
  x$T <- round( x$T / as.vector((crossprod(x$P,x$T)/sum(x$P) )), 3)
  return(x) }


# Multiple variable univariate function for continuous
# This function intakes a list of continuous variables and produces univariate statistics in equally weighted buckets
# @param myds dataframe containing variable
# @param mylist list of continuous variables to put into buckets
# @param wywate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param mytar target variable
# @param mybkt number of quantiles to place variable into
# @return dataframe with n buckets and average target per bucket

uniloop <- function(myds, mylist, mywate, mytar, mybkt) {
  a<-mylist
  d<-myds
  for (ctr in 1:length(a)) {
    x<-a[ctr]
    v<-uni(d,x,mywate,mytar,mybkt)
    v$var<-x
    if (ctr==1) {vv<-v} else {vv<-bind_rows(vv,v)}
    rm(x,v) }
  vv<-vv%>%mutate(L=round(L,3),H=round(H,3),V=round(V,3),T=round(T,3))
  return(vv) }


# Multiple variable univariate function for categorical
# This function intakes a list of categorical variables and produces univariate statistics for any meetign minimum weight threshold
# @param myds dataframe containing variable
# @param mylist list of continuous variables to put into buckets
# @param wywate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param mytar target variable
# @param mythresh cutoff weight for what to display
# @return dataframe with average target per bucket

cateloop <- function(myds, mylist, mywate, mytar, mythresh) {
  a<-mylist
  d<-myds
  for (ctr in 1:length(a)) {
    x<-a[ctr]
    v<-unic(d,x,mywate,mytar)
    v$V <- ifelse(v$P/sum(v$P)<mythresh, 'z_oth', v$V)
    v<-v%>%group_by(V)%>%summarise(T=round(crossprod(P,T)/sum(P),3),P=sum(P))
    v$var<-x
    if (ctr==1) {vv<-v} else {vv<-bind_rows(vv,v)}
    rm(x,v) }
  return(vv) }


# Actual vs Expected function for categorical
# This function produces univariate actual and expected (predicted) statistics for categorical variable
# More generally, it produces weighted averages for two numeric variables (t,p) by a third, categorical variable (v) 
# @param ds dataframe containing variable
# @param v variable to put into buckets
# @param w weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param t target variable
# @param p expected/predicted variable
# @return dataframe with n buckets and average target per bucket

avec<-function(ds,v,w,t,p) {
  x<-ds%>%select(v,w,t,p)%>%rename(V=1,W=2,T=3,P=4)
  x1<-x%>%filter(is.na(V))%>%mutate(V='Z_NA')%>%group_by(V)%>%summarise(Wt=sum(W),T=crossprod(T,W)/sum(W),P=crossprod(P,W)/sum(W))
  x2<-x%>%filter(!is.na(V))%>%mutate(V=as.character(V))%>%group_by(V)%>%summarise(Wt=sum(W),T=crossprod(T,W)/sum(W),P=crossprod(P,W)/sum(W))
  x<-bind_rows(x2,x1)
  x$T <- round( x$T / as.vector((crossprod(x$Wt,x$T)/sum(x$Wt) )), 3)
  x$P <- round( x$P / as.vector((crossprod(x$Wt,x$P)/sum(x$Wt) )), 3)
  return(x) }


# Multiple Variable Actual vs Expected function for categorical
# This function produces univariate actual (target) and expected (predicted) statistics for multiple categorical variables
# More generally, it produces weighted averages for two numeric variables (mytar,mypred) by each categorical variable in list (mylist) 
# @param myds dataframe containing variables
# @param mylist list of categorical variables to put into buckets
# @param wywate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param mytar target (actual) variable
# @param mypred prediction (expected) variable
# @param mythresh cutoff weight for what to display
# @return dataframe with average target and prediction per bucket for each categorical in mylist

avecloop <- function(myds, mylist, mywate, mytar, mypred, mythresh) {
  a<-mylist
  d<-myds
  for (ctr in 1:length(a)) {
    x<-a[ctr]
    v<-avec(d,x,mywate,mytar,mypred)
    v$V <- ifelse(v$Wt/sum(v$Wt)<mythresh, 'z_oth', v$V)
    v<-v%>%group_by(V)%>%summarise(T=round(crossprod(Wt,T)/sum(Wt),3),P=round(crossprod(Wt,P)/sum(Wt),3),Wt=sum(Wt))
    v$var<-x
    if (ctr==1) {vv<-v} else {vv<-bind_rows(vv,v)}
    rm(x,v) }
  return(vv) }

# Actual vs Expected function for continuous
# This function produces univariate actual and expected (predicted) statistics for continuous variable
# More generally, it produces weighted averages for two numeric variables (t,p) by a third, numeric variable (v) 
# @param ds dataframe containing variable
# @param v variable to put into buckets
# @param w weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param t target variable
# @param p expected/predicted variable
# @param q number of quantiles to place variable into
# @return dataframe with n buckets and average target and prediction per bucket

ave<-function(ds,v,w,t,p,q) {
  x<-ds%>%select(v,w,t,p)%>%rename(V=1,W=2,T=3,P=4)
  x1<-x%>%filter(is.na(V))%>%mutate(Q=9999)%>%group_by(Q)%>%summarise(Wt=sum(W),T=crossprod(T,W)/sum(W),P=crossprod(P,W)/sum(W))
  x2<-x%>%filter(!is.na(V))%>%arrange(V)%>%mutate(Q=ceiling(cumsum(W)/sum(W)/(1/q)))%>%group_by(Q)%>%summarise(Wt=sum(W),L=min(V),H=max(V),V=crossprod(V,W)/sum(W),T=crossprod(T,W)/sum(W),P=crossprod(P,W)/sum(W))
  x<-bind_rows(x2,x1)
  x$T <- round( x$T / as.vector((crossprod(x$Wt,x$T)/sum(x$Wt) )), 3)
  x$P <- round( x$P / as.vector((crossprod(x$Wt,x$P)/sum(x$Wt) )), 3)
  return(x) }


# Multiple Variable Actual vs Expected function for continuous
# This function produces univariate actual (target) and expected (predicted) statistics for multiple continuous variables
# More generally, it produces weighted averages for two numeric variables (mytar,mypred) by each numeric variable in list (mylist) 
# @param myds dataframe containing variables
# @param mylist list of continuous variables to put into buckets
# @param wywate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param mytar target (actual) variable
# @param mypred prediction (expected) variable
# @param mybkt number of quantiles to place variable into
# @return dataframe with average target and prediction per bucket for each continuous variable in mylist

aveloop <- function(myds, mylist, mywate, mytar, mypred, mybkt) {
  a<-mylist
  d<-myds
  for (ctr in 1:length(a)) {
    x<-a[ctr]
    v<-ave(d,x,mywate,mytar,mypred,mybkt)
    v$var<-x
    if (ctr==1) {vv<-v} else {vv<-bind_rows(vv,v)}
    rm(x,v) }
  vv<-vv%>%mutate(L=round(L,3),H=round(H,3),V=round(V,3),T=round(T,3))
  return(vv) }


# Relative Net Lift Function
# This function produces "relative net lift" calculations for a specified dataset.
# That is, it calculates the maximum possible Gini gain from a naive to an ideal model, and what % of that gain a candidate model consumes.
# @param myds dataframe containing candidate model(s) and target
# @param myfloor99 naive model
# @param mypred99 candidate model
# @param mytarget99 ideal model
# @param myriem99 number of rectangles to use in riemann sum for area under curve
# @param mymemo99 description to append to output, typically of candidate model
# @return dataframe with naive, candidate, and ideal ginis and resulting opportunity consumption

netlift<-function(ds,myfloor99,mypred99,mytarget99,myriem99,mymemo99){
  f=subset(ds,select=c(myfloor99))%>%rename(f=1)     
  m=subset(ds,select=c(mypred99))%>%rename(m=1)
  c=subset(ds,select=c(mytarget99))%>%rename(c=1)
  a<-cbind(f,m,c)%>%mutate(x=1,f=pmax(0.00001,f) )
  af<- a%>% group_by(x) %>% arrange( f ) %>% mutate( q= ceiling(cumsum( f)/sum(f) * myriem99 ) ) %>% group_by( x,q ) %>% summarise( wf=sum(f), tf=sum(c ) ) %>% arrange(x,q) %>%group_by(x) %>% mutate(wf=cumsum(wf)/sum(wf), tf=cumsum(tf)/sum(tf) ,df=wf-tf)
  am<- a%>% group_by(x) %>% arrange( m ) %>% mutate( q= ceiling(cumsum( f)/sum(f) * myriem99 ) ) %>% group_by( x,q ) %>% summarise( wm=sum(f), tm=sum(c ) ) %>% arrange(x,q) %>%group_by(x) %>% mutate(wm=cumsum(wm)/sum(wm), tm=cumsum(tm)/sum( tm ), dm=wm-tm )
  ac<- a%>% group_by( x) %>% arrange( c ) %>% mutate( q= ceiling(cumsum( f)/sum(f) * myriem99 ) ) %>% group_by( x,q ) %>% summarise( wc=sum(f), tc=sum(c ) ) %>% arrange(x,q) %>%group_by(x) %>% mutate(wc=cumsum(wc)/sum(wc), tc=cumsum(tc)/sum(tc) , dc=wc - tc )
  a<-af%>%inner_join(am)%>%inner_join(ac) %>% mutate( rectangles=1) %>% group_by( x) %>% summarise( dummy=sum( df )/sum(rectangles), model=sum( dm)/sum(rectangles), genius= sum ( dc )/sum(rectangles) )
  a<-a%>%mutate(opportunity=genius-dummy, achievement=model-dummy, relnetlift= achievement/opportunity) %>% ungroup() %>% select(-x) %>% mutate( nametag=mymemo99)
  rm(m,c,f,af,am,ac)
return(a) }



# Tweedie forward selection based on AIC
# This function steps through a list of variables and selects them based on AIC
# At each iteration, it selects the variables that improves AIC the most, holding constant previously selected variable list
# @param myds99 dataframe containing variables, targets, and weights
# @param myvar99 list of variables to step through in standard R list format
# @param mytar99 target variable
# @param wywate99 weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param mypower99 P parameter to use for Tweedie GLM fit
# @return dataframe with order variables are selected and resulting AIC improvement

fwdTweedie<-function( myds99, myvar99, mytar99, mywate99, mypower99 ) {
# weights vector
mywt99<-subset(myds99%>%ungroup(),select=c(mywate99))%>%rename(w=1)
mylen99<-length(myvar99)
for(ctr0 in 1:length(myvar99)) {

# create lists of selected vs. unselected variables to consider each iteration
if (ctr0==1) {myvar90<-myvar99} else
{ myvar90<-as.data.frame(myvar99)%>%rename(myvari=1)%>%  inner_join(myout90%>%filter(iter==(ctr0-1)&rnk!=1)%>%select(myvari))    
     names(myvar90)<−NULL
     myvar90<-as.list(t(myvar90)) 
 myvar9x<-as.data.frame(myvar99)%>%rename(myvari=1)%>%  inner_join(myout90%>%filter(rnk==1)%>%select(myvari))    
     names(myvar9x)<−NULL
     myvar9x<-as.list(t(myvar9x)) }
 for (ctrx in 1:length(myvar9x)) {
        if (ctrx==1) { mystr9x<-as.character(myvar9x[ctrx]) } else
       { mystr9x<-paste(mystr9x,' + ', as.character(myvar9x[ctrx]) ) }
   }

# step through variables one at a time to find minimum AIC
      for (ctr1 in 1:length(myvar90)) {
            myvar91<-as.character(myvar90[ctr1])
            if (ctr0==1) { mystr91<-paste(  mytar99 , ' ~ ', myvar91 ) }
            if (ctr0!=1) { mystr91<-paste(  mytar99 , ' ~ ', mystr9x, ' + ', myvar91 )  }
            myglm91 <-glm(mystr91, myds99 , mywt99$w ,family=tweedie(var.power=mypower99, link.power=0))
            myaic91<-AICtweedie(myglm91)
            myaic91<-as.data.frame(myaic91)%>%rename(aic=1)%>%mutate(myvari=myvar91)
             if (ctr1==1) {myout91<-myaic91} else {myout91<-bind_rows(myout91,myaic91)}
          myout91<-myout91%>%arrange(aic)%>%mutate(rnk=1,rnk=cumsum(rnk),iter=ctr0)
           }
     if (ctr0==1) {myout90<-myout91} else {myout90<-bind_rows(myout90,myout91)}
     if (ctr0==1) {mysel90<-myout90%>%filter(rnk==1)%>%mutate(improve=1.0000)} else {mysel90<-myout90%>%filter(rnk==1)%>%left_join(mysel90%>%mutate(iter=iter+1)%>%select(iter,aic)%>%rename(aic0=aic))%>%mutate(improve=round(aic0/aic,4))%>%select(-aic0) }
gc() }
rm(myaic91,myglm91,myout91,myout90,myvar90,mysel90)
return(mysel90) }


