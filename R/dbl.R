# Lift Chart Macro
#
# This function produces single and double lift charts
# @param champion incumbent prediction (set variable to uniformly 1.00 for single lift chart)
# @param challenger prediction from revised models
# @param target response variable
# @param wate weighting variable (set variable to uniformly 1 for equally weighted observations)
# @param bkt number of quantile to display on lift chart
# @param index_tf 0 - no rescaling, 1 = y-axis variables scaled to average to 1
# @param graph_tf 0 - return table, 1 = return plot
# @return contains either a table or a plot containing the lift chart content
dbl <- function(ds , champ, chall, target, wate,bkt,indx_tf,graph_tf){
##  library(ggplot2)
##  library(dplyr)
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
