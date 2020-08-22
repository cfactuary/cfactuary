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
  a<-a%>%arrange(grp,chg)%>%group_by(grp)%>%mutate(qnt = ceiling( bkt * (cumsum(wgt)/sum(wgt)) ), prd0=prd0*wgt, prd1=prd1*wgt,tar=tar*wgt) %>% group_by(qnt) %>% summarize_at(c('wgt','prd0','prd1','tar'),sum)  %>% mutate( prd0=round(prd0/wgt,3) , prd1=round(prd1/wgt,3), tar=round(tar/wgt,3) )
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
  d<-d%>%mutate(v=exp(v)*w, t=t*w, p=p*w) %>%mutate(c=ceiling( (cumsum(w)/sum(w)) / (1/bkt)))%>%group_by(c)%>%summarize(v=sum(v)/sum(w),t=sum(t)/sum(w),p=sum(p)/sum(w),w=sum(w))
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
  x2 <- x %>% group_by(v,b) %>% summarize(w=sum(w)) %>% arrange(v,desc(w)) %>% group_by(v) %>% mutate(ct=1,ct=cumsum(ct)) %>% filter(ct==1) %>% select(v,b)
  x <- left_join(x%>%select(-b),x2,by='v') %>% group_by(b) %>% mutate(vw=sum(vw)/sum(w)) %>% group_by( v, vw) %>% summarize
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
  x1<-x%>%filter(is.na(V))%>%mutate(Q=9999)%>%group_by(Q)%>%summarize(P=sum(W),T=crossprod(T,W)/sum(W))
  x2<-x%>%filter(!is.na(V))%>%arrange(V)%>%mutate(Q=ceiling(cumsum(W)/sum(W)/(1/q)))%>%group_by(Q)%>%summarize(P=sum(W),L=min(V),H=max(V),V=crossprod(V,W)/sum(W),T=crossprod(T,W)/sum(W))
  x<-bind_rows(x2,x1)
  x$T <- round( x$T / (crossprod(x$P,x$T)/sum(x$P) ), 3)
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
  x1<-x%>%filter(is.na(V))%>%mutate(V='Z_NA')%>%group_by(V)%>%summarize(P=sum(W),T=crossprod(T,W)/sum(W))
  x2<-x%>%filter(!is.na(V))%>%mutate(V=as.character(V))%>%group_by(V)%>%summarize(P=sum(W),T=crossprod(T,W)/sum(W))
  x<-bind_rows(x2,x1)
  x$T <- round( x$T / (crossprod(x$P,x$T)/sum(x$P) ), 3)
  return(x) }




