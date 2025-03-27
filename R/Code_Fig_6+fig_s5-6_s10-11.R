
######0.0预测增温2度对N2O和N2排放的影响#########
{pacman::p_load("ggplot2","ggpmisc","ggpubr","ggprism","doBy","car","readr","crayon","rstatix",
                "lme4","MCMCglmm","tidyverse","sjPlot","directlabels","agricolae","lubridate",
                "patchwork","reshape2","dplyr","animation")  #加载多个包
  Sys.setlocale("LC_TIME","English")  
  
  windowsFonts(Arial=windowsFont("Arial"),
               Times=windowsFont("Times New Roman"))
  
  mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
                axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = .5,vjust =1),
                axis.ticks = element_line(linewidth = .5),
                axis.ticks.length = unit(2,"mm"),
                prism.ticks.length = unit(1,"mm"),
                axis.text.y = element_text(colour = "black",size = 21,hjust=.5), 
                axis.title.x =element_text(size=21), axis.title.y=element_text(colour = "black",size=21),
                legend.text = element_text(size=21), legend.title =element_text(size=21),
                legend.margin = margin(unit(c(0.05,3,5,3),'mm')),
                # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
                panel.background = element_rect(fill =NA, linewidth=0.6,colour = "black", linetype = "solid",
                                                inherit.blank=T),
                panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
                panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
                plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
                legend.key = element_rect(fill = NA), 
                legend.background = element_rect(fill = NA), 
                plot.margin = unit(c(5,5,2,2), 'mm'),   #调整画图区域的间距，从上右下左调整
                strip.background = element_rect(fill = 'snow',color='black'), 
                strip.text = element_text(colour = "black",size = 14,hjust=.5),
                legend.position = "none")
  
  ###计算数据均值，范围，se
  datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                        range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                        n=length(na.omit(x)),
                        se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                        sd=round(sd(na.omit(x)),3))
  }
}

########1.0 底物和功能基因的季节变化及增温响应_Extended Data###########
list.files('D:/Workspace/book/Qingyuan/soil_data/C_N_2023/底物和功能基因数据')

dat<-read.csv('D:/Workspace/book/Qingyuan/soil_data/C_N_2023/底物和功能基因数据/soil dat_2019-2023.csv',header = T)

str(dat)
# dat$plot<-substr(dat$plot_ab,1,1);unique(dat$plot)
dat$treatment<-ifelse(dat$plot %in% c(1,4,5),'warmed','control')
names(dat)[c(8,10)]<-c("EOC","EON")

dat$AOB_AOA=dat$AOB/dat$AOA

dat$AOA_AOB=dat$AOA/dat$AOB

dat$nosZ_nirSK<-dat$nosZ/dat$nirSK

dat$TDN_MBN<-dat$TDN/dat$MBN

dat$EOC_MBC<-dat$EOC/dat$MBC

dat$EON_MBN<-dat$EON/dat$MBN

dat$TIN_MBN<-(dat$NH4_N+dat$NO3_N)/dat$MBN

dat$AOB=dat$AOB*(1+dat$soil.water_dry)
dat$AOA=dat$AOA*(1+dat$soil.water_dry)
dat$amoA=dat$AOA+dat$AOB

dat$nirS=dat$nirS*(1+dat$soil.water_dry)
dat$nirK=dat$nirK*(1+dat$soil.water_dry)
dat$nirSK=dat$nirSK*(1+dat$soil.water_dry)
dat$nosZ=dat$nosZ*(1+dat$soil.water_dry)

str(dat)
unique(dat$soillayer)

dat$month<-as.numeric(format(as.Date(dat$Date),'%m'))

temp<-read.csv("D:/Workspace/book/Qingyuan/soil_data/warming system/Soiltemp-2 degree_comity.csv",header=T)
str(temp)
temp$Temp30<-as.numeric(temp$Temp30)
temp$date<-as.Date(as.POSIXct(temp$MCGS_Time,format="%Y/%m/%d %H"),
                   tz="Asia/Taipei");str(temp)

temp<-aggregate(temp[,c(3:8)], by=list(temp$date),mean, na.rm = TRUE)
names(temp)[1]<-'Date'
temp<-reshape2::melt(temp, id="Date", measure=c(2:7), value.factor=F, na.rm=TRUE)
temp$plot<-substr(temp$variable,5,5)

dat$temp<-NA
for(i in unique(dat$Date)){
  for(j in 1:6){
  dat$temp[dat$Date==i & dat$plot==j]<-
    temp$value[temp$Date==i & temp$plot==j]
  }
}



#MCMCglmm是基于贝叶斯的包, 也可以做混合线性模型分析:

library(MCMCglmm)
M1<-MCMCglmm(NH4_N~treatment+factor(Date),random = ~plot_ab,
             data=dat[which(dat$soillayer=='O' & dat$month %in% seq(5,11,1)),])
summary(M1)

library(ggplot2)

myfun<-function(x){c(m=mean(na.omit(x)),se=sd(na.omit(x))/sqrt(length(na.omit(x))))}

library(doBy)
dat.treat<-summaryBy(data=dat,FUN=myfun,
                     AOA+AOB+AOB_AOA+AOA_AOB+nirK+nirS+nirSK+nosZ+MBC+MBN+EOC+EON+NH4_N+NO3_N+EON_MBN+EOC_MBC+TDN_MBN+TIN_MBN+soil.water_dry~treatment+Date+soillayer)

str(dat.treat)

###计算显著性
dat.res<-data.frame();k=1
for(i in unique(dat$Date)){
  for(j in unique(dat$soillayer)){
    dat_ij<-dat[dat$Date==i & dat$soillayer==j,]
    dat.res[k,1]=i
    dat.res[k,2]=j
    dat.res[k,3]='warmed'
    if(i=='2021/5/22'){
      lm<-t.test(NH4_N~treatment,data=dat_ij);lm
      dat.res[k,16]=mean(dat_ij$NH4_N[dat_ij$treatment=='warmed'])   #NH4_N均值
      dat.res[k,17]=lm$p.value    #P值NH4_N
      
      lm<-t.test(NO3_N~treatment,data=dat_ij);lm
      dat.res[k,18]=mean(dat_ij$NO3_N[dat_ij$treatment=='warmed'])   #NO3_N均值
      dat.res[k,19]=lm$p.value    #P值NO3_N
      
      lm<-t.test(soil.water_dry~treatment,data=dat_ij);lm
      dat.res[k,20]=mean(dat_ij$soil.water_dry[dat_ij$treatment=='warmed'])   #moisture均值
      dat.res[k,21]=lm$p.value  
      
      dat.res[k,40]=mean(dat_ij$NH4_N[dat_ij$treatment=='warmed'])/
        mean(dat_ij$NH4_N[dat_ij$treatment=='control']) #NH4_N均值比值
      
      dat.res[k,41]=mean(dat_ij$NO3_N[dat_ij$treatment=='warmed'])/
        mean(dat_ij$NO3_N[dat_ij$treatment=='control']) #NO3_N均值比值
      
    }else{
      
    lm<-t.test(log(AOA,10)~treatment,data=dat_ij);lm
    dat.res[k,4]=mean(log(dat_ij$AOA[dat_ij$treatment=='warmed'],10))   #AOA均值
    dat.res[k,5]=lm$p.value    #P值AOA
    
    lm<-t.test(log(AOB,10)~treatment,data=dat_ij);lm
    dat.res[k,6]=mean(log(dat_ij$AOB[dat_ij$treatment=='warmed'],10))   #AOB均值
    dat.res[k,7]=lm$p.value    #P值AOB
    
    lm<-t.test(log(nirS,10)~treatment,data=dat_ij);lm
    dat.res[k,8]=mean(log(dat_ij$nirS[dat_ij$treatment=='warmed'],10))   #nirS均值
    dat.res[k,9]=lm$p.value    #P值AOA
    
    lm<-t.test(log(nirK,10)~treatment,data=dat_ij);lm
    dat.res[k,10]=mean(log(dat_ij$nirK[dat_ij$treatment=='warmed'],10))   #nirK均值
    dat.res[k,11]=lm$p.value    #P值nirK
    
    lm<-t.test(log(nosZ,10)~treatment,data=dat_ij);lm
    dat.res[k,12]=mean(log(dat_ij$nosZ[dat_ij$treatment=='warmed'],10))   #nosZ均值
    dat.res[k,13]=lm$p.value    #P值nosZ
    
    lm<-t.test(log(EOC,10)~treatment,data=dat_ij);lm
    dat.res[k,14]=mean(log(dat_ij$EOC[dat_ij$treatment=='warmed'],10))   #EOC均值
    dat.res[k,15]=lm$p.value    #P值EOC
    
    lm<-t.test(NH4_N~treatment,data=dat_ij);lm
    dat.res[k,16]=mean(dat_ij$NH4_N[dat_ij$treatment=='warmed'])   #NH4_N均值
    dat.res[k,17]=lm$p.value    #P值NH4_N
    
    lm<-t.test(NO3_N~treatment,data=dat_ij);lm
    dat.res[k,18]=mean(dat_ij$NO3_N[dat_ij$treatment=='warmed'])   #NO3_N均值
    dat.res[k,19]=lm$p.value    #P值NO3_N
    
    lm<-t.test(soil.water_dry~treatment,data=dat_ij);lm
    dat.res[k,20]=mean(dat_ij$soil.water_dry[dat_ij$treatment=='warmed'])   #moisture均值
    dat.res[k,21]=lm$p.value    #P值moisture 
    
    
    lm<-t.test(nirSK~treatment,data=dat_ij);lm
    dat.res[k,22]=mean(dat_ij$nirSK[dat_ij$treatment=='warmed'])   #nirSK均值
    dat.res[k,23]=lm$p.value    #P值nirSK
    
    lm<-t.test(log(AOA_AOB,10)~treatment,data=dat_ij);lm
    dat.res[k,24]=mean(log(dat_ij$AOA_AOB[dat_ij$treatment=='warmed'],10))   #AOB/AOA均值
    dat.res[k,25]=lm$p.value    #P值AOB/AOA
    
    lm<-t.test(log(EON,10)~treatment,data=dat_ij);lm
    dat.res[k,26]=mean(log(dat_ij$EON[dat_ij$treatment=='warmed'],10))   #EON均值
    dat.res[k,27]=lm$p.value    #P值EON
    
    lm<-t.test(log(MBC,10)~treatment,data=dat_ij);lm
    dat.res[k,28]=mean(log(dat_ij$MBC[dat_ij$treatment=='warmed'],10))   #MBC均值
    dat.res[k,29]=lm$p.value    #P值MBC
    
    lm<-t.test(log(MBN,10)~treatment,data=dat_ij);lm
    dat.res[k,30]=mean(log(dat_ij$MBN[dat_ij$treatment=='warmed'],10))   #MBN均值
    dat.res[k,31]=lm$p.value    #P值MBN
    
    lm<-t.test(log(EOC_MBC,10)~treatment,data=dat_ij);lm
    dat.res[k,32]=mean(log(dat_ij$EOC_MBC[dat_ij$treatment=='warmed'],10))   #EOC_MBC均值
    dat.res[k,33]=lm$p.value    #P值EOC_MBC
    
    lm<-t.test(log(EON_MBN,10)~treatment,data=dat_ij);lm
    dat.res[k,34]=mean(log(dat_ij$EON_MBN[dat_ij$treatment=='warmed'],10))   #EON_MBN均值
    dat.res[k,35]=lm$p.value    #P值EON_MBN
    
    
    lm<-t.test(log(TDN_MBN,10)~treatment,data=dat_ij);lm
    dat.res[k,36]=mean(log(dat_ij$TDN_MBN[dat_ij$treatment=='warmed'],10))   #TDN_MBN均值
    dat.res[k,37]=lm$p.value    #P值TDN_MBN
    
    
    lm<-t.test(log(TIN_MBN,10)~treatment,data=dat_ij);lm
    dat.res[k,38]=mean(log(dat_ij$TIN_MBN[dat_ij$treatment=='warmed'],10))   #TIN_MBN均值
    dat.res[k,39]=lm$p.value    #P值TIN_MBN
    
    dat.res[k,40]=mean(dat_ij$NH4_N[dat_ij$treatment=='warmed'])/
      mean(dat_ij$NH4_N[dat_ij$treatment=='control']) #NH4_N均值比值
    
    dat.res[k,41]=mean(dat_ij$NO3_N[dat_ij$treatment=='warmed'])/
      mean(dat_ij$NO3_N[dat_ij$treatment=='control']) #NO3_N均值比值
    
    dat.res[k,42]=mean(dat_ij$MBC[dat_ij$treatment=='warmed'])/
      mean(dat_ij$MBC[dat_ij$treatment=='control']) #MBN均值比值
    
    dat.res[k,43]=mean(dat_ij$MBN[dat_ij$treatment=='warmed'])/
      mean(dat_ij$MBN[dat_ij$treatment=='control']) #MBN均值比值
    
    dat.res[k,44]=mean(dat_ij$EOC[dat_ij$treatment=='warmed'])/
      mean(dat_ij$EOC[dat_ij$treatment=='control']) #DON均值比值
    
    dat.res[k,45]=mean(dat_ij$EON[dat_ij$treatment=='warmed'])/
      mean(dat_ij$EON[dat_ij$treatment=='control']) #DON均值比值
    
    dat.res[k,46]=mean(dat_ij$AOA[dat_ij$treatment=='warmed'])/
      mean(dat_ij$AOA[dat_ij$treatment=='control']) #AOA均值比值
    
    dat.res[k,47]=mean(dat_ij$AOB[dat_ij$treatment=='warmed'])/
      mean(dat_ij$AOB[dat_ij$treatment=='control']) #AOB均值比值
    
    dat.res[k,48]=mean(dat_ij$nirS[dat_ij$treatment=='warmed'])/
      mean(dat_ij$nirS[dat_ij$treatment=='control']) #nirS均值比值
    
    dat.res[k,49]=mean(dat_ij$nirK[dat_ij$treatment=='warmed'])/
      mean(dat_ij$nirK[dat_ij$treatment=='control']) #nirK均值比值
    
    }
    k=k+1  
  }
}


names(dat.res)<-c('Date','soillayer','treatment','AOA.m','P_AOA','AOB.m','P_AOB','nirS.m','P_nirS','nirK.m','P_nirK','nosZ.m','P_nosZ',
                  'EOC.m','P_EOC','NH4_N.m','P_NH4_N','NO3_N.m','P_NO3_N','soil.water_dry.m','P_soil.water_dry','nirSK.m','P_nirSK',
                  'AOA_AOB.m','P_AOA_AOB','EON.m','P_EON','MBC.m','P_MBC','MBN.m','P_MBN',
                  'EOC_MBC.m','P_EOC_MBC','EON_MBN.m','P_EON_MBN','TDN_MBN.m','P_TDN_MBN','TIN_MBN.m','P_TIN_MBN',
                  'NH4.R','NO3.R','MBC.R','MBN.R','EOC.R','EON.R','AOA.R','AOB.R','nirS.R','nirK.R')

str(dat.treat)

datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                      range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                      n=length(na.omit(x)),
                      se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                      sd=round(sd(na.omit(x)),3))}


unique(dat.res$Date)

datFUN((dat.res$NH4.R[dat.res$soillayer=='O' ]-1))
datFUN((dat.res$NH4.R[dat.res$soillayer=='O' & dat.res$NH4.R<=1]-1))

datFUN((dat.res$NO3.R[dat.res$soillayer=='O']-1))
datFUN((dat.res$NO3.R[dat.res$soillayer=='O' & dat.res$NO3.R>=1]-1))

datFUN((dat.res$MBN.R[dat.res$soillayer=='O']-1))
datFUN((dat.res$MBN.R[dat.res$soillayer=='0-10cm']-1))

datFUN((dat.res$EON.R[dat.res$soillayer=='O']-1))
datFUN((dat.res$EON.R[dat.res$soillayer=='0-10cm']-1))

datFUN((dat.res$AOA.R[dat.res$soillayer=='O']-1))
datFUN((dat.res$AOB.R[dat.res$soillayer=='O']-1))

datFUN((dat.res$nirK.R[dat.res$soillayer=='O']-1))

datFUN((dat.res$nirS.R[dat.res$soillayer=='O']-1))

datFUN((dat.res$nirS.R[dat.res$soillayer=='0-10cm']-1))

datFUN((dat.res$EON.R[dat.res$soillayer=='0-10cm']-1))

library(ggplot2);library(ggprism)
mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
              axis.text.x = element_text(colour = "black",size=18,angle = 0,hjust = .5,vjust =1),
              axis.ticks = element_line(linewidth = .5),
              axis.ticks.length = unit(2,"mm"),
              prism.ticks.length = unit(1,"mm"),
              axis.text.y = element_text(colour = "black",size = 18,hjust=.5), 
              axis.title.x =element_text(size=18), axis.title.y=element_text(colour = "black",size=18,hjust=.5),
              legend.text = element_text(size=18,hjust=0), legend.title =element_text(size=18),
              panel.background = element_rect(fill =NA, linewidth=0.3,colour = "black", linetype = "solid",
                                              inherit.blank=T),
              panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
              panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
              plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
              legend.key = element_rect(fill = NA), 
              legend.background = element_rect(fill = NA), 
              plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   #调整画图区域的间距，从上右下左调整
              strip.background = element_rect(fill = 'snow2',color='black'), 
              strip.text = element_text(colour = "black",size = 14,hjust=.5),
              legend.position = "none")

str(dat.treat)

dat.treat$group<-paste(substr(dat.treat$Date,1,4),dat.treat$treatment,'-')

dat.treat<-within(dat.treat,soillayer<-factor(soillayer,levels = c('O','0-10cm')))
dat.res<-within(dat.res,soillayer<-factor(soillayer,levels = c('O','0-10cm')))

unique(dat.treat$soillayer)

################## Ext.Dat_Fig.3 #################
labels <- c("O" = "O layer","0-10cm"="0-10 cm","10-20cm"="10-20 cm","20-40cm"="20-40 cm")

unique(dat.treat$Date)

pw<-ggplot(dat.treat,aes(Date,soil.water_dry.m*100,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=(soil.water_dry.m-soil.water_dry.se)*100,
                    ymax=(soil.water_dry.m+soil.water_dry.se)*100),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=soil.water_dry.m*100-15,
                label=ifelse(P_soil.water_dry<=0.001,'***',
                             ifelse(P_soil.water_dry<=0.01,'**',
                                    ifelse(P_soil.water_dry<=0.05,'*',
                                           ifelse(P_soil.water_dry<=0.10,paste('P = ',round(P_soil.water_dry,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(limits=c(0,120),n.breaks=4,guide = "prism_offset_minor",expand = c(0, 0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        # panel.grid.minor = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = c(0.06,0.85))+
  labs(x='Date',y=expression(atop(paste('Soil moisture'),paste('(%)'))));pw

p1<-ggplot(dat.treat,aes(Date,NH4_N.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=(NH4_N.m-NH4_N.se),
                    ymax=(NH4_N.m+NH4_N.se)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=NH4_N.m+5,
                label=ifelse(P_NH4_N<=0.001,'***',
                             ifelse(P_NH4_N<=0.01,'**',
                                    ifelse(P_NH4_N<=0.05,'*',
                                           ifelse(P_NH4_N<=0.059,paste('P = ',round(P_NH4_N,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(limits=c(0,30),n.breaks=4,guide = "prism_offset_minor",expand = c(0, 0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('NH'[4]^'+'),paste('(mg N kg'^'-1',')'))));p1

p2<-ggplot(dat.treat,aes(Date,NO3_N.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=(NO3_N.m-NO3_N.se),
                    ymax=(NO3_N.m+NO3_N.se)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=NO3_N.m+5,
                label=ifelse(P_NO3_N<=0.001,'***',
                             ifelse(P_NO3_N<=0.01,'**',
                                    ifelse(P_NO3_N<=0.05,'*',
                                           ifelse(P_NO3_N<=0.059,paste('P = ',round(P_NO3_N,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(limits=c(0,60),breaks=seq(0,60,20),guide = "prism_offset_minor",expand = c(0, 0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('NO'[3]^'-'),paste('(mg N kg'^'-1',')'))));p2

p3<-ggplot(dat.treat,aes(Date,EOC.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=(EOC.m-EOC.se),
                    ymax=(EOC.m+EOC.se)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=10^EOC.m+100,
                label=ifelse(P_EOC<=0.001,'***',
                             ifelse(P_EOC<=0.01,'**',
                                    ifelse(P_EOC<=0.055,'*',
                                           ifelse(P_EOC<=0.059,paste('P = ',round(P_EOC,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(limits=c(0,600),breaks=seq(0,600,200),guide = "prism_offset_minor",expand = c(0, 0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('EOC'),paste('(mg C kg'^'-1',')'))));p3

p4<-ggplot(dat.treat,aes(Date,EON.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  geom_errorbar(aes(ymin=(EON.m-EON.se),
                    ymax=(EON.m+EON.se)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=10^EON.m+15,
                label=ifelse(P_EON<=0.001,'***',
                             ifelse(P_EON<=0.01,'**',
                                    ifelse(P_EON<=0.055,'*',
                                           ifelse(P_EON<=0.059,paste('P = ',round(P_EON,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(limits=c(0,90),breaks=seq(0,90,30),guide = "prism_offset_minor",expand = c(0, 0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('EON'),paste('(mg N kg'^'-1',')'))));p4

p5<-ggplot(dat.treat,aes(Date,log(AOA.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(AOA.m-AOA.se,10),
                    ymax=log(AOA.m+AOA.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=AOA.m+0.5,
                label=ifelse(P_AOA<=0.001,'***',
                             ifelse(P_AOA<=0.01,'**',
                                    ifelse(P_AOA<=0.05,'*',
                                           ifelse(P_AOA<=0.059,paste('P = ',round(P_AOA,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(6,12,2),limits=c(6,12),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^6)),expression(paste('10'^8)),expression(paste('10'^10)),expression(paste('10'^12))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('AOA'),paste('(copies g'^-1,')'))));p5

p6<-ggplot(dat.treat,aes(Date,log(AOB.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(AOB.m-AOB.se,10),
                    ymax=log(AOB.m+AOB.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=AOB.m+0.28,
                label=ifelse(P_AOB<=0.001,'***',
                             ifelse(P_AOB<=0.01,'**',
                                    ifelse(P_AOB<0.055,'*',
                                           ifelse(P_AOB<=0.059,paste('P = ',round(P_AOB,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(5,8,1),limits=c(5,8),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^5)),expression(paste('10'^6)),expression(paste('10'^7)),expression(paste('10'^8))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('AOB'),paste('(copies g'^-1,')'))));p6


str(dat.treat)
p7<-ggplot(dat.treat,aes(Date,log(AOA_AOB.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(AOA_AOB.m-AOA_AOB.se,10),
                    ymax=log(AOA_AOB.m+AOA_AOB.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=0.7+AOA_AOB.m,
                label=ifelse(P_AOA_AOB<=0.001,'***',
                             ifelse(P_AOA_AOB<=0.01,'**',
                                    ifelse(P_AOA_AOB<0.055,'*',
                                           ifelse(P_AOA_AOB<=0.059,paste('P = ',round(P_AOA_AOB,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(0,6,2),limits=c(0,6),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^-1)),expression(paste('10')),expression(paste('10'^3)),
                              expression(paste('10'^6))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('AOA/AOB'),paste(''))));p7


p8<-ggplot(dat.treat,aes(Date,log(nirK.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(nirK.m-nirK.se,10),
                    ymax=log(nirK.m+nirK.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=nirK.m+0.5,
                label=ifelse(P_nirK<=0.001,'***',
                             ifelse(P_nirK<=0.01,'**',
                                    ifelse(P_nirK<=0.05,'*',
                                           ifelse(P_nirK<=0.059,paste('P = ',round(P_nirK,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(6,9,1),limits=c(6,9),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^6)),expression(paste('10'^7)),expression(paste('10'^8)),expression(paste('10'^9))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('nirK'),paste('(copies g'^-1,')'))));p8



p9<-ggplot(dat.treat,aes(Date,log(nirS.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(nirS.m-nirS.se,10),
                    ymax=log(nirS.m+nirS.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=nirS.m-0.75,
                label=ifelse(P_nirS<=0.001,'***',
                             ifelse(P_nirS<=0.01,'**',
                                    ifelse(P_nirS<=0.05,'*','')))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(5,9,2),limits=c(5,9),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^5)),expression(paste('10'^7)),expression(paste('10'^9))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('nirS'),paste('(copies g'^-1,')'))));p9



 p10<-ggplot(dat.treat,aes(Date,log(nirSK.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(nirSK.m-nirSK.se,10),
                    ymax=log(nirSK.m+nirSK.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=nirSK.m+0.5,
                label=ifelse(P_nirSK<=0.001,'***',
                             ifelse(P_nirSK<=0.01,'**',
                                    ifelse(P_nirSK<=0.05,'*',
                                           ifelse(P_nirSK<=0.059,paste('P = ',round(P_nirSK,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(6,9,1),limits=c(6,9),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^6)),expression(paste('10'^7)),expression(paste('10'^8)),expression(paste('10'^9))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                             "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                    labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                             "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='Date')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 60,hjust = 1,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='Date',y=expression(atop(paste('nirS+nirK'),paste('(copies g'^-1,')'))));p10


blank_year<-data.frame(date=rep(c('2019/7/21','2020/9/21','2021/6/24','2022/8/12','2023/7/7'),2),
                       soillayer=c(rep('O',5),rep('0-10cm',5)),
                       treatment=c('control'),
                       labs=rep(c('2019','2020','2021','2022','2023'),2),
                       group=c('control-2019'))

blank_year<-within(blank_year,soillayer<-factor(soillayer,levels=c('O','0-10cm')))

p11<-ggplot(dat.treat,aes(Date,log(nosZ.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(nosZ.m-nosZ.se,10),
                    ymax=log(nosZ.m+nosZ.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=nosZ.m+0.5,
                label=ifelse(P_nosZ<=0.001,'***',
                             ifelse(P_nosZ<=0.01,'**',
                                    ifelse(P_nosZ<=0.05,'*',
                                           ifelse(P_nosZ<=0.059,paste('P = ',round(P_nosZ,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(7,9,1),minor_breaks=seq(7,9,0.5),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression(paste('10'^7)),expression(paste('10'^8)),expression(paste('10'^9))))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"),name='')+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.06,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  # labs(x='',y=expression(atop(paste('nosZ'),paste('(log'[10],', copies g'^-1,')'))))+
  labs(x='',y=expression(atop(paste('nosZ'),paste('(copies g'^-1,')'))))+
  coord_cartesian(clip="off",ylim=c(7,9))+   #在plot绘图区域外面标注的关键
  geom_text(data = blank_year[blank_year$labs %in% c('2019','2020','2021'),],aes(x=date, y=5.75,label=labs),
            vjust=0,color='black',size=4)+
  geom_text(data = blank_year[blank_year$date %in% c('2022/8/12'),],aes(x=10.5, y=5.75,label='2022'),vjust=0,hjust=0.5,
            color='black',size=4)+
  geom_text(data = blank_year[blank_year$date %in% c('2023/7/7'),],aes(x=12.5, y=5.75,label='2023'),vjust=0,hjust=0.5,
            color='black',size=4)+
  geom_segment(data = blank_year,aes(x="2019/7/3", xend="2019/9/22",y=6.15,yend=6.15),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2020/7/17", xend="2020/11/5",y=6.15,yend=6.15),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2021/5/22", xend="2021/9/19",y=6.15,yend=6.15),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2022/8/12", xend="2022/10/1",y=6.15,yend=6.15),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2023/7/7", xend="2023/9/23",y=6.15,yend=6.15),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F);p11

#######
library(patchwork);library(ggpubr)


blank_year<-data.frame(date=rep(c('2019/7/21','2020/9/21','2021/6/24','2022/8/12','2023/7/7'),2),
                       soillayer=c(rep('O',5),rep('0-10cm',5)),
                       treatment=c('control'),
                       labs=rep(c('2019','2020','2021','2022','2023'),2),
                       group=c('control-2019'))

blank_year<-within(blank_year,soillayer<-factor(soillayer,levels=c('O','0-10cm')))

m1<-ggplot(dat.treat,aes(Date,log(MBC.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(MBC.m-MBC.se,10),
                    ymax=log(MBC.m+MBC.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=MBC.m+0.15,
                label=ifelse(P_MBC<=0.001,'***',
                             ifelse(P_MBC<=0.01,'**',
                                    ifelse(P_MBC<=0.05,'*',
                                           ifelse(P_MBC<=0.059,paste('P = ',round(P_MBC,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(2,4,1),guide = "prism_offset_minor",expand = c(0,0),
                     labels=c(expression('10'^2),expression('10'^3),expression('10'^4)))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  # labs(x='',y=expression(atop(paste('nosZ'),paste('(log'[10],', copies g'^-1,')'))))+
  labs(x='',y=expression(atop(paste('MBC'),paste('(mg C kg'^-1,')'))))+
  coord_cartesian(clip="off",ylim=c(2,4));m1

m2<-ggplot(dat.treat,aes(Date,log(MBN.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(MBN.m-MBN.se,10),
                    ymax=log(MBN.m+MBN.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=MBN.m+0.15,
                label=ifelse(P_MBN<=0.001,'***',
                             ifelse(P_MBN<=0.01,'**',
                                    ifelse(P_MBN<=0.05,'*',
                                           ifelse(P_MBN<=0.059,paste('P = ',round(P_MBN,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(2,3,0.5),labels=c(expression('10'^2.0),expression('10'^2.5),expression('10'^3.0)),guide = "prism_offset_minor",expand = c(0,0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/5/22","2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "5/22","6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  # labs(x='',y=expression(atop(paste('nosZ'),paste('(log'[10],', copies g'^-1,')'))))+
  labs(x='',y=expression(atop(paste('MBN'),paste('(mg N kg'^-1,')'))))+
  coord_cartesian(clip="off",ylim=c(2,3))+   #在plot绘图区域外面标注的关键
  geom_text(data = blank_year[blank_year$labs %in% c('2019','2020'),],aes(x=date, y=1.4,label=labs),
            vjust=0,color='black',size=4)+
  geom_text(data = blank_year[blank_year$date %in% c('2021/6/24'),],aes(x='2021/6/24', y=1.4,label='2021'),vjust=0,hjust=0.5,
            color='black',size=4)+
  geom_text(data = blank_year[blank_year$date %in% c('2022/8/12'),],aes(x=10, y=1.4,label='2022'),vjust=0,hjust=0,
            color='black',size=4)+
  geom_text(data = blank_year[blank_year$date %in% c('2023/7/7'),],aes(x=12, y=1.4,label='2023'),vjust=0,hjust=0,
            color='black',size=4)+
  geom_segment(data = blank_year,aes(x="2019/7/3", xend="2019/9/22",y=1.56,yend=1.56),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2020/7/17", xend="2020/11/5",y=1.56,yend=1.56),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2021/5/22", xend="2021/9/19",y=1.56,yend=1.56),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2022/8/12", xend="2022/10/1",y=1.56,yend=1.56),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2023/7/7", xend="2023/9/23",y=1.56,yend=1.56),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F);m2



(pw+rremove("xlab")+rremove("x.text")+
    theme(plot.margin = unit(c(0.2,0.3,0.25,0.3),'cm'),
          legend.text = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5),
          axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+
    (p1+rremove("xlab")+rremove("x.text"))+
    theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.08,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+
    (p2+rremove("xlab")+rremove("x.text"))+
    theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.08,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+
    (p3+rremove("xlab")+rremove("x.text"))+
    theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+
    (p4+rremove("xlab")+rremove("x.text"))+
    theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.14,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=12,hjust = 0,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+
    (m1+rremove("xlab")+rremove("x.text"))+
    theme(plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          strip.background = element_blank(),strip.text = element_blank(),
          axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+  
    (m2+theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.8,0.3),'cm'),
              axis.text.x = element_text(colour = "black",size=12,hjust = 1,vjust =1,angle=30),
              axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
              axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5)))+
    plot_layout(ncol = 1, byrow = T)+
    plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
    theme(plot.tag.position = c(0, 0.98),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
          plot.tag = element_text(size = 14,vjust = 0,hjust=0,face="bold")))->figure;figure

getwd()
ggsave("Seasonal character of substrate+microbes_2019-2023.pdf", figure, width =10, height =10,
       device=cairo_pdf)


################## Ext.Dat_Fig.4 #################

((p5+rremove("xlab")+rremove("x.text"))+
   theme(plot.margin = unit(c(0.3,0.3,0.25,0.3),'cm'),
         axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
         axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+
   (p6+rremove("xlab")+rremove("x.text"))+
   theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
         axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
         axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+    
   (p7+rremove("xlab")+rremove("x.text"))+
   theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
         axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
         axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+  
   (p8+rremove("xlab")+rremove("x.text"))+
   theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
         axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
         axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+  
   (p9+rremove("xlab")+rremove("x.text"))+
   theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
         axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
         axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+  
   (p10+rremove("xlab")+rremove("x.text"))+
   theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
         axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
         axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5))+  
   (p11+theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.8,0.3),'cm'),
              axis.text.x = element_text(colour = "black",size=12,hjust = 1,vjust =1,angle = 30),
              axis.text.y = element_text(colour = "black",size=12,hjust = 1,vjust =0.5),
              axis.title.y = element_text(colour = "black",size=12,hjust = 0.5,vjust =0.5)))+
   plot_layout(ncol = 1, byrow = T)+
   plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
   theme(plot.tag.position = c(0, 0.98),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
         plot.tag = element_text(size = 14,vjust = 0,hjust=0,face="bold")))->figure;figure

setwd('D:/工作目录/202409/Manuscript_kai/Talk_20241113/Data and Code/Fig_6+Ext Dat fig_3_4_7_8')

ggsave("Seasonal character of functional genes_final_N2019-2023.pdf", figure, width =10, height =10,
       device=cairo_pdf)


###在每个Date采样下都存在样方的误差！！！20240821 by kai
###treatment 组间因子，Date组内因子，plot样方差异
summary(aov(soil.water_dry~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(soil.water_dry~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='0-10cm'),]))

summary(aov(NH4_N~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(NO3_N~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(EOC~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(EON~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(EON~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='0-10cm'),]))
summary(aov(MBC~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(MBN~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))

summary(aov(MBC~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='0-10cm'),]))
summary(aov(MBN~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='0-10cm'),]))

summary(aov(AOA~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(AOB~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(log(AOA_AOB,10)~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))

summary(aov(nirS~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))
summary(aov(nirS~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='0-10cm'),]))
summary(aov(nirK~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='0-10cm'),]))
summary(aov(nosZ~treatment*factor(Date)+Error(plot/Date),data=dat[which(dat$soillayer=='O'),]))



###微生物与氮的比值
str(dat.treat)

ggplot(dat.treat,aes(Date,TDN_MBN.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_errorbar(aes(ymin=TDN_MBN.m-TDN_MBN.se,
                    ymax=TDN_MBN.m+TDN_MBN.se),linewidth=0.27,width=0.1)+
  geom_text(data=dat.res,
            aes(x=Date,y=10^TDN_MBN.m+0.065,
                label=ifelse(P_TDN_MBN<=0.001,'***',
                             ifelse(P_TDN_MBN<=0.01,'**',
                                    ifelse(P_TDN_MBN<=0.05,'*',
                                           ifelse(P_TDN_MBN<=0.059,paste('P = ',round(P_TDN_MBN,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(0,0.3,0.1),minor_breaks = seq(0,0.3,0.05),guide = "prism_offset_minor",expand = c(0,0),
                     labels=seq(0,0.3,0.1))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='',y=expression(atop(paste('TDN/MBN ratio'),paste(''))))+
  coord_cartesian(clip="off",ylim=c(0,0.3))+   #在plot绘图区域外面标注的关键
  geom_text(data = blank_year[blank_year$labs %in% c('2019','2020'),],aes(x=date, y=-0.075,label=labs),vjust=0,color='black',size=6)+
  geom_text(data = blank_year[blank_year$date %in% c('2021/6/24'),],aes(x='2021/6/24', y=-0.075,label='2021'),vjust=0,hjust=-0.03,
            color='black',size=6)+
  geom_text(data = blank_year[blank_year$date %in% c('2022/8/12'),],aes(x=9.5, y=-0.075,label='2022'),vjust=0,hjust=0.5,
            color='black',size=6)+
  geom_text(data = blank_year[blank_year$date %in% c('2023/7/7'),],aes(x=11.5, y=-0.075,label='2023'),vjust=0,hjust=0.5,
            color='black',size=6)+
  geom_segment(data = blank_year,aes(x="2019/7/3", xend="2019/9/22",y=-0.047,yend=-0.047),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2020/7/17", xend="2020/11/5",y=-0.047,yend=-0.047),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2021/6/24", xend="2021/9/19",y=-0.047,yend=-0.047),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2022/8/12", xend="2022/10/1",y=-0.047,yend=-0.047),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2023/7/7", xend="2023/9/23",y=-0.047,yend=-0.047),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)->m3;m3

  

ggplot(dat.treat,aes(Date,EON_MBN.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_errorbar(aes(ymin=EON_MBN.m-EON_MBN.se,
                    ymax=EON_MBN.m+EON_MBN.se),linewidth=0.27,width=0.1)+
  geom_text(data=dat.res,
            aes(x=Date,y=10^EON_MBN.m+0.025,
                label=ifelse(P_EON_MBN<=0.001,'***',
                             ifelse(P_EON_MBN<=0.01,'**',
                                    ifelse(P_EON_MBN<=0.05,'*',
                                           ifelse(P_EON_MBN<=0.059,paste('P = ',round(P_EON_MBN,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(0,0.15,0.05),minor_breaks = seq(0,0.15,0.025),guide = "prism_offset_minor",expand = c(0,0),
                     labels=seq(0,0.15,0.05),limits=c(0,0.15))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='',y=expression(atop(paste('EON/MBN ratio'),paste(''))))->m4;m4


str(dat.treat)
ggplot(dat.treat,aes(Date,TIN_MBN.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_errorbar(aes(ymin=TIN_MBN.m-TIN_MBN.se,
                    ymax=TIN_MBN.m+TIN_MBN.se),linewidth=0.27,width=0.1)+
  geom_text(data=dat.res,
            aes(x=Date,y=10^TIN_MBN.m+0.025,
                label=ifelse(P_TIN_MBN<=0.001,'***',
                             ifelse(P_TIN_MBN<=0.01,'**',
                                    ifelse(P_TIN_MBN<=0.05,'*',
                                           ifelse(P_TIN_MBN<=0.059,paste('P = ',round(P_TIN_MBN,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(0,0.15,0.05),minor_breaks = seq(0,0.15,0.025),guide = "prism_offset_minor",expand = c(0,0),
                     labels=seq(0,0.15,0.05),limits=c(0,0.15))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='',y=expression(atop(paste('TIN/MBN ratio'),paste(''))))->m5;m5

ggplot(dat.treat,aes(Date,EOC_MBC.m,col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_errorbar(aes(ymin=EOC_MBC.m-EOC_MBC.se,
                    ymax=EOC_MBC.m+EOC_MBC.se),linewidth=0.27,width=0.1)+
  geom_text(data=dat.res,
            aes(x=Date,y=10^EOC_MBC.m+0.025,
                label=ifelse(P_EOC_MBC<=0.001,'***',
                             ifelse(P_EOC_MBC<=0.01,'**',
                                    ifelse(P_EOC_MBC<=0.05,'*',
                                           ifelse(P_EOC_MBC<=0.059,paste('P = ',round(P_EOC_MBC,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(0,1.2,0.4),minor_breaks = seq(0,1.2,0.2),guide = "prism_offset_minor",expand = c(0,0),
                     labels=seq(0,1.2,0.4),limits=c(0,1.2))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  labs(x='',y=expression(atop(paste('EOC/MBC ratio'),paste(''))))->m6;m6


((m4+rremove("xlab")+rremove("x.text"))+
    theme(plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5))+  
    (m5+rremove("xlab")+rremove("x.text"))+
    theme(plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          strip.background = element_blank(),strip.text = element_blank(),
          axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5))+  
    (m3+theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.8,0.3),'cm'),
               axis.text.x = element_text(colour = "black",size=18,hjust = 0.5,vjust =1),
               axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
               axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5)))+
    plot_layout(ncol = 1, byrow = T)+
    plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
    theme(plot.tag.position = c(0, 0.95),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
          plot.tag = element_text(size = 19,vjust = 0,hjust=0,face="bold")))->figure;figure

getwd()
ggsave("Seasonal character of microbes_final2.pdf", figure, width =16, height =9.5,
       device=cairo_pdf)






m2<-ggplot(dat.treat,aes(Date,log(MBN.m,10),col=treatment))+
  geom_point(size=0.9,pch=1,stroke=1)+geom_line()+facet_wrap(soillayer~., labeller=labeller(soillayer = labels))+
  geom_point(color='white',size=0.75)+geom_line(aes(group=group))+
  # geom_blank(data = blank_data, aes(x = date, y = soilwater_dry.m, group=lab))+
  geom_errorbar(aes(ymin=log(MBN.m-MBN.se,10),
                    ymax=log(MBN.m+MBN.se,10)),linewidth=0.27,width=0.1)+
  scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                      labels=c("control","warmed"))+
  geom_text(data=dat.res,
            aes(x=Date,y=MBN.m+0.45,
                label=ifelse(P_MBN<=0.001,'***',
                             ifelse(P_MBN<=0.01,'**',
                                    ifelse(P_MBN<=0.05,'*',
                                           ifelse(P_MBN<=0.059,paste('P = ',round(P_MBN,2)),''))))),
            color='black',size=7.45)+
  scale_y_continuous(breaks=seq(2,4,1),labels=c(expression('10'^2),expression('10'^3),expression('10'^4)),
                     guide = "prism_offset_minor",expand = c(0,0))+
  scale_x_discrete(limits=c("2019/7/3","2019/7/21","2019/9/22","2020/7/17","2020/9/21","2020/11/5",
                            "2021/6/24","2021/9/19","2022/8/12","2022/10/1","2023/7/7","2023/9/23"),
                   labels=c("7/3","7/21","9/22","7/17","9/21","11/5",
                            "6/24","9/19","8/12","10/1","7/7","9/23"))+mythem+
  theme(axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = 0.5,vjust =1),
        panel.grid.major = element_line(linetype=2,size=0.15,colour = 'grey'),
        panel.border = element_rect(size=.36,fill=NA,colour = "black"),
        plot.margin = unit(c(0.4,0.3,0.8,0.3),'cm'),  
        panel.spacing.y = unit(5, 'mm'),legend.position = 'none')+
  # labs(x='',y=expression(atop(paste('nosZ'),paste('(log'[10],', copies g'^-1,')'))))+
  labs(x='',y=expression(atop(paste('MBN'),paste('(mg N kg'^-1,')'))))+
  coord_cartesian(clip="off",ylim=c(2,4))+   #在plot绘图区域外面标注的关键
  geom_text(data = blank_year[blank_year$labs %in% c('2019','2020'),],aes(x=date, y=1.05,label=labs),vjust=0,color='black',size=6)+
  geom_text(data = blank_year[blank_year$date %in% c('2021/6/24'),],aes(x='2021/6/24', y=1.05,label='2021'),vjust=0,hjust=-0.03,
            color='black',size=6)+
  geom_text(data = blank_year[blank_year$date %in% c('2022/8/12'),],aes(x=9.5, y=1.05,label='2022'),vjust=0,hjust=0.5,
            color='black',size=6)+
  geom_text(data = blank_year[blank_year$date %in% c('2023/7/7'),],aes(x=11.5, y=1.05,label='2023'),vjust=0,hjust=0.5,
            color='black',size=6)+
  geom_segment(data = blank_year,aes(x="2019/7/3", xend="2019/9/22",y=1.46,yend=1.46),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2020/7/17", xend="2020/11/5",y=1.46,yend=1.46),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2021/6/24", xend="2021/9/19",y=1.46,yend=1.46),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2022/8/12", xend="2022/10/1",y=1.46,yend=1.46),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data = blank_year,aes(x="2023/7/7", xend="2023/9/23",y=1.46,yend=1.46),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F);m2

((p1+rremove("xlab")+rremove("x.text"))+
    theme(plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5))+  
    (p2+rremove("xlab")+rremove("x.text"))+
    theme(strip.background = element_blank(),strip.text = element_blank(),
          plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5))+  
    (m1+rremove("xlab")+rremove("x.text"))+
    theme(strip.background = element_blank(),strip.text = element_blank(),
          plot.margin = unit(c(0.1,0.3,0.25,0.3),'cm'),
          axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
          axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5))+  
    (m2+theme(strip.background = element_blank(),strip.text = element_blank(),plot.margin = unit(c(0.1,0.3,0.8,0.3),'cm'),
              axis.text.x = element_text(colour = "black",size=18,hjust = 0.5,vjust =1),
              axis.text.y = element_text(colour = "black",size=18,hjust = 1,vjust =0.5),
              axis.title.y = element_text(colour = "black",size=18,hjust = 0.5,vjust =0.5)))+
    plot_layout(ncol = 1, byrow = T)+
    plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
    theme(plot.tag.position = c(0, 0.95),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
          plot.tag = element_text(size = 19,vjust = 0,hjust=0,face="bold")))->figure;figure

getwd()
ggsave("Seasonal character of IN and microbes_final.pdf", figure, width =16, height =9,
       device=cairo_pdf)


(mean(na.omit(dat$NH4_N[dat$soillayer=='O' & dat$treatment=='warmed']))-
    mean(na.omit(dat$NH4_N[dat$soillayer=='O' & dat$treatment=='control'])))/mean(na.omit(dat$NH4_N[dat$soillayer=='O' & dat$treatment=='control']))

(mean(na.omit(dat$NO3_N[dat$soillayer=='O' & dat$treatment=='warmed']))-
    mean(na.omit(dat$NO3_N[dat$soillayer=='O' & dat$treatment=='control'])))/mean(na.omit(dat$NO3_N[dat$soillayer=='O' & dat$treatment=='control']))

datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                      range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                      n=length(na.omit(x)),
                      se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                      sd=round(sd(na.omit(x)),3))}

datFUN(na.omit(dat$NH4_N[dat$soillayer=='O' & dat$treatment=='control']))
datFUN(na.omit(dat$NO3_N[dat$soillayer=='O' & dat$treatment=='control']))






(mean(na.omit(dat$EON[dat$soillayer=='O' & dat$treatment=='warmed']))-
    mean(na.omit(dat$EON[dat$soillayer=='O' & dat$treatment=='control'])))/mean(na.omit(dat$EON[dat$soillayer=='O' & dat$treatment=='control']))


(mean(na.omit(dat$EON[dat$soillayer=='0-10cm' & dat$treatment=='warmed']))-
    mean(na.omit(dat$EON[dat$soillayer=='0-10cm' & dat$treatment=='control'])))/mean(na.omit(dat$EON[dat$soillayer=='0-10cm' & dat$treatment=='control']))

(mean(na.omit(dat$MBN[dat$soillayer=='0-10cm' & dat$treatment=='warmed']))-
    mean(na.omit(dat$MBN[dat$soillayer=='0-10cm' & dat$treatment=='control'])))/mean(na.omit(dat$MBN[dat$soillayer=='0-10cm' & dat$treatment=='control']))

(mean(na.omit(dat$MBN[dat$soillayer=='O' & dat$treatment=='warmed']))-
    mean(na.omit(dat$MBN[dat$soillayer=='O' & dat$treatment=='control'])))/mean(na.omit(dat$MBN[dat$soillayer=='O' & dat$treatment=='control']))

(mean(na.omit(dat$AOB[dat$soillayer=='O' & dat$treatment=='warmed']))-
    mean(na.omit(dat$AOB[dat$soillayer=='O' & dat$treatment=='control'])))/mean(na.omit(dat$AOB[dat$soillayer=='O' & dat$treatment=='control']))


(mean(na.omit(dat$nirK[dat$soillayer=='0-10cm' & dat$treatment=='warmed']))-
    mean(na.omit(dat$nirK[dat$soillayer=='0-10cm' & dat$treatment=='control'])))/mean(na.omit(dat$nirK[dat$soillayer=='0-10cm' & dat$treatment=='control']))


dat$N_perM<-(dat$NH4_N+dat$NO3_N+dat$EON)/dat$MBN   #微生物对N的竞争
dat$IN_ON<-(dat$NH4_N+dat$NO3_N)/dat$EON  #微生物对N的竞争
dat$MBC_MBN<-dat$MBC/dat$MBN  #微生物对N的竞争

summary(aov(N_perM~treatment*Date+Error(plot/Date),data=dat[dat$soillayer=='O',]))
summary(aov(N_perM~treatment*Date+Error(plot/Date),data=dat[dat$soillayer=='0-10cm',]))

summary(aov(IN_ON~treatment*Date+Error(plot/Date),data=dat[dat$soillayer=='O',]))
summary(aov(IN_ON~treatment*Date+Error(plot/Date),data=dat[dat$soillayer=='0-10cm',]))

summary(aov(MBC_MBN~treatment*Date+Error(plot_ab/Date),data=dat[dat$soillayer=='O',]))
summary(aov(MBC_MBN~treatment*Date+Error(plot_ab/Date),data=dat[dat$soillayer=='0-10cm',]))

###总结增温对微生物碳氮含量的影响
dat$year=substr(dat$Date,1,4)
  
str(dat)
dat_res<-data.frame()
k=1
for(i in 2019:2022){
  for(j in unique(dat$soillayer)){
    for(t in unique(dat$treatment)){
    dat_ij=dat[dat$year==i & dat$soillayer==j & dat$treatment==t,]
    dat_res[k,1]=i
    dat_res[k,2]=j
    dat_res[k,3]=t
    dat_res[k,4]=mean(na.omit(dat_ij$MBN))
    dat_res[k,5]=mean(na.omit(dat_ij$MBC))
    
    dat_res[k,6]=mean(na.omit(dat_ij$EOC))
    dat_res[k,7]=mean(na.omit(dat_ij$EON))
    dat_res[k,8]=mean(na.omit(dat_ij$C_con))
    dat_res[k,9]=mean(na.omit(dat_ij$N_con))
    
    dat_res[k,10]=mean(na.omit(dat_ij$NH4_N))
    dat_res[k,11]=mean(na.omit(dat_ij$NO3_N))
    
    dat_res[k,12]=mean(na.omit(dat_ij$nosZ))
    dat_res[k,13]=mean(na.omit(dat_ij$nirS))
    dat_res[k,14]=mean(na.omit(dat_ij$nirK))
    dat_res[k,15]=mean(na.omit(dat_ij$AOA))
    dat_res[k,16]=mean(na.omit(dat_ij$AOB))
    dat_res[k,17]=mean(na.omit(dat_ij$soil.water_dry))
    k=k+1
    }
  }
}


names(dat_res)<-c('year','soillayer','treatment','MBN','MBC','EOC','EON','TC','TN','NH4_N','NO3_N','nosZ','nirS','nirK','AOA','AOB','soil.water_dry')

write.csv(dat_res,'土壤理化指标按年划分.csv',row.names = F)

################# Fig.6  #######################
pacman::p_load("ggplot2","ggpmisc","ggpubr","ggprism","doBy","car","readr","crayon","rstatix",
               "lme4","MCMCglmm","tidyverse","sjPlot","directlabels","agricolae","lubridate",
               "patchwork","reshape2","dplyr","animation")  #加载多个包
Sys.setlocale("LC_TIME","English")  

windowsFonts(Arial=windowsFont("Arial"),
             Times=windowsFont("Times New Roman"))

mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
              axis.text.x = element_text(colour = "black",size=18,angle = 0,hjust = .5,vjust =1),
              axis.ticks = element_line(linewidth = .5),
              axis.ticks.length = unit(2,"mm"),
              prism.ticks.length = unit(1,"mm"),
              axis.text.y = element_text(colour = "black",size = 18,hjust=.5), 
              axis.title.x =element_text(size=18), axis.title.y=element_text(colour = "black",size=18,hjust=.5),
              legend.text = element_text(size=18,hjust=0), legend.title =element_text(size=18),
              # legend.margin = margin(unit(c(0.05,3,1,3),'mm')),
              # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
              panel.background = element_rect(fill =NA, linewidth=0.6,colour = "black", linetype = "solid",
                                              inherit.blank=T),
              panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
              panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
              plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
              legend.key = element_rect(fill = NA), 
              legend.background = element_rect(fill = NA), 
              plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   #调整画图区域的间距，从上右下左调整
              strip.background = element_rect(fill = NA,color='black'), 
              strip.text = element_text(colour = "black",size = 18,hjust=.5),
              legend.position = "none")



############## 2) Extended fig-净氮矿化硝化实验 #################
setwd("D:/Workspace/book/Qingyuan/soil_data/C_N_2023/Rplots")

library(readxl)
s <- read.csv("D:/Workspace/book/Qingyuan/soil_data/C_N_2023/20231103 2020-2023 input net mineralization and nitrification lab incubation.csv",
                  header=T)

head(s) 
str(s)
names(s)
levels(as.factor(s$date))
s$date[which(s$date=="2020/11/5")]<-"2020/9/21"
s$date[which(s$date=="2021/9/24")]<-"2021/9/19"
names(s)[c(1,9,10)]<-c("date","ammonium","nitrate")
s$sample<-paste(s$date,s$subplot,s$treatment,s$layer,sep = "_")
#####################
unique(s$incubation.time)

s0<-s[s$incubation.time%in% c("0","1"),]
head(s0)
names(s0)
str(s0)
names(s0)[c(9,10)]<-c("ammonium_ini","nitrate_ini")####initial
s0$layer<-factor(s0$layer,ordered = T,levels = c("Oa+e","0-10 cm","10-20 cm","20-40 cm"))
s0$ammonium_ini<-as.numeric(s0$ammonium_ini)
s0$nitrate_ini<-as.numeric(s0$nitrate_ini)
###
library(ggplot2)
ggplot(s0,aes(treatment,ammonium_ini,fill=treatment))+geom_boxplot()+
  facet_grid(layer~date,scales = "free")

ggplot(s0,aes(treatment,nitrate_ini,fill=treatment))+geom_boxplot()+
  facet_grid(layer~date,scales = "free")

############################
s7<-s[s$incubation.time=="7",]###培养7天
head(s7)
s7$year<-format(as.Date(s7$date),"%Y")

s7$ammonium_ini<-NA
s7$nitrate_ini<-NA

unique(s7$year)
for(i in 1:length(s7$sample)){
  if(s7$year[i]=='2024' & s7$layer[i] %in% c('10-20 cm','20-40 cm')){
  s7$date[i]
  s7$sample[i]
  s7$temperature[i]
  s0$ammonium_ini[s0$sample==s7$sample[i]]
  
  s7$ammonium_ini[i]=s0$ammonium_ini[s0$sample==s7$sample[i] & s0$temperature==s7$temperature[i]]
  s7$nitrate_ini[i]=s0$nitrate_ini[s0$sample==s7$sample[i] & s0$temperature==s7$temperature[i]]
  
  s7$ammonium_prod[i]<-(s7$ammonium[i]-s7$ammonium_ini[i])/6
  s7$nitrate_prod[i]<-(s7$nitrate[i]-s7$nitrate_ini[i])/6
  
  }else{
    s7$ammonium_ini[i]=s0$ammonium_ini[s0$sample==s7$sample[i]]
    s7$nitrate_ini[i]=s0$nitrate_ini[s0$sample==s7$sample[i]]
    
    s7$ammonium_prod[i]<-(s7$ammonium[i]-s7$ammonium_ini[i])/7
    s7$nitrate_prod[i]<-(s7$nitrate[i]-s7$nitrate_ini[i])/7
    
    
}
}



# start
sm<-s7   # N transformation during the experimental incubation

head(sm)
str(sm)
sm$temperature<-as.integer(sm$temperature)
sm$ammonium<-as.numeric(sm$ammonium)
sm$nitrate<-as.numeric(sm$nitrate)
sm$ammonium_ini<-as.numeric(sm$ammonium_ini)
sm$nitrate_ini<-as.numeric(sm$nitrate_ini)

### N-min
sm$mineralization<-sm$ammonium_prod+sm$nitrate_prod
sm$layer<-factor(sm$layer,ordered = T,levels = c("Oa+e","0-10 cm","10-20 cm","20-40 cm"))
head(sm);names(sm)

#####################################
###calculate Q10 for each sample
###################################
str(sm)
sm$sample<-as.factor(sm$sample)
sample_name<-levels(sm$sample)
n<-length(sample_name)
n
######
k=1
q10<-data.frame(sample=character(0),
                R0_min=numeric(0),
                b_min=numeric(0),
                R0_nit=numeric(0),
                b_nit=numeric(0))

for(i in sample_name){
  datai<-sm[which(sm$sample==i&sm$mineralization>0),]###sm$mineralization>0，否则计算停止
  if (nrow(datai)>2){
    nlm<-nls(mineralization~R0*exp(b*temperature),
             data = datai,start=list(R0=0.2,b=0.01),trace = TRUE)
    a<-summary(nlm)
    q10[k,1]=i
    q10[k,2]=coef(nlm)[1]   #R0_min
    q10[k,3]=coef(nlm)[2]   #b_min
    

    nlm_i<-nls(nitrate_prod~R0*exp(b*temperature),data = datai,start=list(R0=0.2,b=0.01),
             control=nls.control(maxiter=200),trace = TRUE)
    b<-summary(nlm_i)
    
    q10[k,4]=coef(nlm_i)[1]   # R0_nitrification 
    q10[k,5]=coef(nlm_i)[2]   #b_nitrification
    
    k=k+1
  }
}

library(stringr)
sample_s<-str_split_fixed(q10$sample,"_",n=4)
str(sample_s)
head(q10)
q10$date<-sample_s[,1]
q10$subplot<-sample_s[,2]
q10$treatment<-sample_s[,3]
q10$layer<-sample_s[,4]

q10;rm(nlm,nlm_i,s0,s7,s,datai,a,q10.res,q10all,Q10n,dat_ij,b,R15,R15N,res,Ron,sample_s)


############## cal the Q10,R10,R15,R17
q10$Q10_min<-exp(q10$b_min*10)
q10$R10_min<-q10$R0_min*exp(q10$b_min*10)
q10$R15_min<-q10$R0_min*exp(q10$b_min*15)
q10$R17_min<-q10$R0_min*exp(q10$b_min*17)


q10$Q10_nit<-exp(q10$b_nit*10)
q10$R10_nit<-q10$R0_nit*exp(q10$b_nit*10)
q10$R15_nit<-q10$R0_nit*exp(q10$b_nit*15)
q10$R17_nit<-q10$R0_nit*exp(q10$b_nit*17)

write.csv(q10,"20231104 output Q10 for net N min and nit lab.csv")


##plotting
q10<-read.csv("20231104 output Q10 for net N min and nit lab.csv",header = T)
str(q10)
leveQ10_nitlevels(q10$layer)

q10$layer<-factor(q10$layer,ordered = T,levels = c("Oa+e","0-10 cm","10-20 cm","20-40 cm"));library(data.table)
levels(q10$layer)
q10$year <- factor(substr(q10$date,1,4),ordered = T,levels = c("2020","2021","2022","2023"))
levels(q10$year)

####sm方差分析####
names(q10)

q10.res<-data.frame()
k=1
for(i in unique(q10$year)){
  for(j in unique(q10$layer)){
    dat_ij=q10[which(q10$year== i & q10$layer== j), ]
    if(length(dat_ij$treatment)>=6){
    q10.res[k,1]=i
    q10.res[k,2]=j
    
    res<-t.test(Q10_min~treatment,data=dat_ij)
    q10.res[k,3]=mean(dat_ij$Q10_min)
    q10.res[k,4]=res$p.value
    
    
    res<-t.test(R0_min~treatment,data=dat_ij)
    q10.res[k,5]=mean(dat_ij$R0_min)
    q10.res[k,6]=res$p.value
    
    res<-t.test(R15_min~treatment,data=dat_ij)
    q10.res[k,7]=mean(dat_ij$R15_min)
    q10.res[k,8]=res$p.value
    
    ###nit
    res<-t.test(Q10_nit~treatment,data=dat_ij)
    q10.res[k,9]=mean(dat_ij$Q10_nit)
    q10.res[k,10]=res$p.value
    
    
    res<-t.test(R0_nit~treatment,data=dat_ij)
    q10.res[k,11]=mean(dat_ij$R0_nit)
    q10.res[k,12]=res$p.value
    
    res<-t.test(R15_nit~treatment,data=dat_ij)
    q10.res[k,13]=mean(dat_ij$R15_nit)
    q10.res[k,14]=res$p.value
    
    k=k+1
    }
  }
}


names(q10.res)<-c('year','layer','Q10_min','P_Q10_min','R0_min','P_R0_min','R15_min','P_R15_min',
                                 'Q10_nit','P_Q10_nit','R0_nit','P_R0_nit','R15_nit','P_R15_nit')

q10.res$treatment<-'control'

q10.res<-within(q10.res,layer<-factor(layer,levels=c('Oa+e','0-10 cm','10-20 cm','20-40 cm')))
q10$layer<-factor(q10$layer,ordered = T,levels = c("Oa+e","0-10 cm","10-20 cm","20-40 cm"))

Q10n <- ggplot(q10,aes(year,Q10_min,fill=treatment))+geom_boxplot(width=0.5,outlier.shape = 1,position=position_dodge(0.6))+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(.~layer)+
  ylim(c(0,7))+
  scale_fill_manual(values = c("#1776EB","#F06B49"),labels=c('Control','Warmed'),name="")+
  labs(x=expression('Year'),
       y=expression(paste('Q'[10], ' of soil net N mineralization')))+
  geom_text(data=q10.res[q10.res$P_Q10_min<=0.05,],
            aes(x=year,y=Q10_min+1.5,
                label=ifelse(P_Q10_min<=0.001,'***',
                             ifelse(P_Q10_min<=0.01,'**',
                                    ifelse(P_Q10_min<=0.05,'*',
                                           ifelse(P_Q10_min<=0.10,paste('P = ',round(P_Q10_min,2)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_Q10_min>0.05,],
            aes(x=year,y=Q10_min+2.2,
                label=ifelse(P_Q10_min<=0.001,'***',
                             ifelse(P_Q10_min<=0.01,'**',
                                    ifelse(P_Q10_min<=0.05,'*',
                                           ifelse(P_Q10_min<=0.10,paste('P=',round(P_Q10_min,2)),''))))),
            color='black',size=3,hjust=0.4)+
  theme_bw()+
  theme(legend.position = c(9/10,9.3/10),legend.background =element_rect(fill =NA, colour = NA, linetype = "solid"));Q10n

Ron <- ggplot(q10,aes(year,R0_min,fill=treatment))+geom_boxplot(width=0.5,outlier.shape = 1,position=position_dodge(0.6))+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(.~layer)+
  scale_fill_manual(values = c("#1776EB","#F06B49"),name="Treatment")+
  labs(x=expression('Year'),
       y=expression(atop(paste('R'[0],' of soil net N mineralization'),paste( ' (mg N',' kg soil'^{-1},' day'^{-1},')'))))+
  theme_bw()+ ylim(c(0,3))+
  geom_text(data=q10.res[q10.res$P_R0_min<=0.05,],
            aes(x=year,y=R0_min+.5,
                label=ifelse(P_R0_min<=0.001,'***',
                             ifelse(P_R0_min<=0.01,'**',
                                    ifelse(P_R0_min<=0.05,'*',
                                           ifelse(P_R0_min<=0.10,paste('P = ',round(P_R0_min,2)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_R0_min>0.05,],
            aes(x=year,y=R0_min+.7,
                label=ifelse(P_R0_min<=0.001,'***',
                             ifelse(P_R0_min<=0.01,'**',
                                    ifelse(P_R0_min<=0.05,'*',
                                           ifelse(P_R0_min<=0.10,paste('P=',round(P_R0_min,2)),''))))),
            color='black',size=3,hjust=0.4)+
  theme(legend.position = 'none',
        legend.key =element_rect(fill =NA),
        legend.background =element_rect(fill =NA, colour = NA, linetype = "solid"));Ron


###模拟增温样地土壤矿化减少了

(mean(q10$R17_min[q10$layer=='Oa+e' & q10$treatment=='warmed'])-
  mean(q10$R15_min[q10$layer=='Oa+e' & q10$treatment=='control']))/
  mean(q10$R15_min[q10$layer=='Oa+e' & q10$treatment=='control'])


########### R15_control vs R15_warmed shown as the following
R15N <- ggplot(q10,aes(year,R15_min,fill=treatment))+geom_boxplot(width=0.5,outlier.shape = 1,position=position_dodge(0.6))+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(.~layer)+
  scale_fill_manual(values = c("#1776EB","#F06B49"),name="Treatment")+
  labs(x=expression('Year'),
       y=expression(atop(paste('R'[15],' of soil net N mineralization'),paste( ' (mg N',' kg soil'^{-1},' day'^{-1},')'))))+
  theme_bw()+ ylim(c(0,6))+
  geom_text(data=q10.res[q10.res$P_R15_min<=0.05,],
            aes(x=year,y=R15_min+1.2,
                label=ifelse(P_R15_min<=0.001,'***',
                             ifelse(P_R15_min<=0.01,'**',
                                    ifelse(P_R15_min<=0.05,'*',
                                           ifelse(P_R15_min<=0.10,paste('P = ',round(P_R15_min,2)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_R15_min>0.05,],
            aes(x=year,y=R15_min+1.5,
                label=ifelse(P_R15_min<=0.001,'***',
                             ifelse(P_R15_min<=0.01,'**',
                                    ifelse(P_R15_min<=0.05,'*',
                                           ifelse(P_R15_min<=0.1,paste('P=',round(P_R15_min,2)),''))))),
            color='black',size=3,hjust=0.4)+
  theme(legend.position = 'none',
        legend.key =element_rect(fill =NA),
        legend.background =element_rect(fill =NA, colour = NA, linetype = "solid"));R15N



str(q10)

#######制作R15 vs R17的数据 dat_N
dat_N=q10
dat_N$R15_min[dat_N$treatment=='warmed']=dat_N$R17_min[dat_N$treatment=='warmed']  #将增温的R15换成R17

###统计分析再跑一一下啊
q10.res<-data.frame()
k=1
for(i in unique(q10$year)){
  for(j in unique(q10$layer)){
    dat_ij=dat_N[which(dat_N$year== i & dat_N$layer== j), ]
    if(length(dat_ij$treatment)>=6){
      q10.res[k,1]=i
      q10.res[k,2]=j
      
      ### min
      res<-t.test(Q10_min~treatment,data=dat_ij)
      q10.res[k,3]=mean(dat_ij$Q10_min)
      q10.res[k,4]=res$p.value
      
      
      res<-t.test(R0_min~treatment,data=dat_ij)
      q10.res[k,5]=mean(dat_ij$R0_min)
      q10.res[k,6]=res$p.value
      
      res<-t.test(R15_min~treatment,data=dat_ij)
      q10.res[k,7]=mean(dat_ij$R15_min)
      q10.res[k,8]=res$p.value
      
      ### nit
      res<-t.test(Q10_nit~treatment,data=dat_ij)
      q10.res[k,9]=mean(dat_ij$Q10_nit)
      q10.res[k,10]=res$p.value
      
      
      res<-t.test(R0_nit~treatment,data=dat_ij)
      q10.res[k,11]=mean(dat_ij$R0_nit)
      q10.res[k,12]=res$p.value
      
      res<-t.test(R15_nit~treatment,data=dat_ij)
      q10.res[k,13]=mean(dat_ij$R15_nit)
      q10.res[k,14]=res$p.value
      
      k=k+1
    }
  }
}




names(q10.res)<-c('year','layer','Q10_min','P_Q10_min','R0_min','P_R0_min','R15_min','P_R15_min',
                                 'Q10_nit','P_Q10_nit','R0_nit','P_R0_nit','R15_nit','P_R15_nit')

q10.res$treatment<-'control'

q10.res<-within(q10.res,layer<-factor(layer,levels=c('Oa+e','0-10 cm','10-20 cm','20-40 cm')))
dat_N$layer<-factor(dat_N$layer,ordered = T,levels = c("Oa+e","0-10 cm","10-20 cm","20-40 cm"))




########### R15_control vs R17_warmed shown as the following

R1517N <- ggplot(dat_N,aes(year,R15_min,fill=treatment))+geom_boxplot(width=0.5,outlier.shape = 1,position=position_dodge(0.6))+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(.~layer)+
  scale_fill_manual(values = c("#1776EB","#F06B49"),name="Treatment")+
  labs(x=expression('Year'),
       y=expression(atop(paste('Soil net N mineralization'),paste( ' (mg N',' kg soil'^{-1},' day'^{-1},')'))))+
  theme_bw()+ ylim(c(0,6))+
  geom_text(data=q10.res[q10.res$P_R15_min<=0.05,],
            aes(x=year,y=R15_min+1.2,
                label=ifelse(P_R15_min<=0.001,'***',
                             ifelse(P_R15_min<=0.01,'**',
                                    ifelse(P_R15_min<=0.05,'*',
                                           ifelse(P_R15_min<=0.10,paste('P = ',round(P_R15_min,2)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_R15_min>0.05,],
            aes(x=year,y=R15_min+1.5,
                label=ifelse(P_R15_min<=0.001,'***',
                             ifelse(P_R15_min<=0.01,'**',
                                    ifelse(P_R15_min<=0.05,'*',
                                           ifelse(P_R15_min<=0.1,paste('P=',round(P_R15_min,2)),''))))),
            color='black',size=3,hjust=0.4)+
  theme(legend.position = 'none',
        legend.key =element_rect(fill =NA),
        legend.background =element_rect(fill =NA, colour = NA, linetype = "solid"));R1517N

R1517_nitri  <- ggplot(dat_N,aes(year,R15_nit,fill=treatment))+
  geom_boxplot(width=0.5,outlier.shape = 1,position=position_dodge(0.6))+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(.~layer)+ylim(c(0,6))+
  scale_fill_manual(values = c("#1776EB","#F06B49"),name="Treatment")+
  labs(x=expression('Year'),
       y=expression(atop(paste('Soil net N nitrification'),paste( ' (mg N',' kg soil'^{-1},' day'^{-1},')'))))+
  theme_bw()+
  geom_text(data=q10.res[q10.res$P_R15_nit<=0.05,],
            aes(x=year,y=R15_nit+1.2,
                label=ifelse(P_R15_nit<=0.001,'***',
                             ifelse(P_R15_nit<=0.01,'**',
                                    ifelse(P_R15_nit<=0.05,'*',
                                           ifelse(P_R15_nit<=0.10,paste('P = ',round(P_R15_nit,2)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_R15_nit>0.05,],
            aes(x=year,y=R15_nit+1.5,
                label=ifelse(P_R15_nit<=0.001,'***',
                             ifelse(P_R15_nit<=0.01,'**',
                                    ifelse(P_R15_nit<=0.05,'*',
                                           ifelse(P_R15_nit<=0.10,paste('P = ',round(P_R15_nit,2)),''))))),
            color='black',size=3,hjust=0.4)+
  theme(legend.position = 'none',
        legend.key =element_rect(fill =NA),
        legend.background =element_rect(fill =NA, colour = NA, linetype = "solid"));R1517_nitri 


########Fig 6.ab_<Rmin in 15 and 17>#############
R1517N <- ggplot(dat_N[dat_N$layer %in% c('Oa+e','0-10 cm'),],aes(year,R15_min,fill=treatment))+
  geom_boxplot(width=0.35,outlier.shape = 1,position=position_dodge(0.6),size=0.4)+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(layer~.,labeller = labeller(layer=c('Oa+e'='O layer','0-10 cm'='0-10 cm')))+
  scale_fill_manual(values = c("#1776EB","#F06B50"),name="Treatment")+
  labs(x=expression('Year'),
       y=expression(atop(paste('Soil net N mineralization'),
                         paste('(mg N',' kg soil'^{-1},' day'^{-1},')'))))+
  scale_y_continuous(limits=c(0,6),breaks=seq(0,6,2),expand=c(0,0))+
  geom_text(data=dat_N[dat_N$layer %in% c('Oa+e') & dat_N$R15_min==max(dat_N$R15_min),],col='black',
            aes(y=5.4,x='2023'),label=expression(paste('Treat: ', italic('P'),'= 0.018')),
            position=position_dodge(0.8),size=4,hjust=.5)+
  geom_text(data=dat_N[dat_N$layer %in% c('Oa+e') & dat_N$R15_min==max(dat_N$R15_min),],col='black',
            aes(y=4.7,x='2023'),label=expression(paste('Year: ', italic('P'),'< 0.001')),
            position=position_dodge(0.8),size=4,hjust=.5)+
  geom_text(data=dat_N[dat_N$layer %in% c('Oa+e') & dat_N$R15_min==max(dat_N$R15_min),],col='black',
            aes(y=4.0,x='2023'),label='Treat x Year:  ns',
            position=position_dodge(0.8),size=4,hjust=0.45)+
  geom_text(data=q10.res[q10.res$P_R15_min<=0.05 & q10.res$layer %in% c('Oa+e','0-10 cm'),],
            aes(x=year,y=R15_min+1.5,
                label=ifelse(P_R15_min<=0.001,'***',
                             ifelse(P_R15_min<=0.01,'**',
                                    ifelse(P_R15_min<=0.05,'*',
                                           ifelse(P_R15_min<=0.10,paste('P = ',round(P_R15_min,3)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_R15_min>0.05 & q10.res$layer %in% c('Oa+e','0-10 cm'),],
            aes(x=year,y=R15_min+1.7,
                label=ifelse(P_R15_min<=0.001,'***',
                             ifelse(P_R15_min<=0.01,'**',
                                    ifelse(P_R15_min<=0.05,'*',
                                           ifelse(P_R15_min<=0.1,paste('P=',round(P_R15_min,3)),''))))),
            color='black',size=3.5,hjust=0.4)+
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(colour = "black",size=12,angle = 0,hjust = .5,vjust =1),
        axis.ticks = element_line(linewidth = .5),
        axis.ticks.length = unit(2,"mm"),
        prism.ticks.length = unit(1,"mm"),
        axis.text.y = element_text(colour = "black",size = 12,hjust=.5), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(colour = "black",size=12,hjust=.5),
        legend.text = element_text(size=12,hjust=0), 
        panel.background = element_rect(fill =NA, linewidth=0.26,colour = "black", linetype = "solid",inherit.blank=T),
        panel.spacing.y = unit(6,'mm'),
        panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
        panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
        plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
        legend.key = element_rect(fill = NA), legend.title=element_blank(),
        legend.background = element_rect(fill = NA), 
        plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   #调整画图区域的间距，从上右下左调整
        strip.background = element_rect(fill = NA,color='black'), 
        strip.text = element_text(colour = "black",size = 10,hjust=.5),
        legend.position = c(0.13,0.925));R1517N

R1517_nitri  <- ggplot(dat_N[dat_N$layer %in% c('Oa+e','0-10 cm'),],aes(year,R15_nit,fill=treatment))+
  geom_boxplot(width=0.35,outlier.shape = 1,position=position_dodge(0.6),size=0.4)+
  geom_point(pch=1,position=position_dodge(0.6))+
  facet_grid(layer~.,labeller = labeller(layer=c('Oa+e'='O layer','0-10 cm'='0-10 cm')))+
  scale_y_continuous(limits=c(0,6),breaks=seq(0,6,2),expand=c(0,0))+
  scale_fill_manual(values = c("#1776EB","#F06B50"),name="Treatment")+
  labs(x=expression('Year'),
       y=expression(atop(paste('Soil net N nitrification'),
                         paste( '(mg N',' kg soil'^{-1},' day'^{-1},')'))))+
  geom_text(data=dat_N[dat_N$layer %in% c('Oa+e') & dat_N$R15_nit==max(dat_N$R15_nit),],col='black',
            aes(y=5.4,x='2023'),label=expression(paste('Treat: ', italic('P'),'= 0.006')),
            position=position_dodge(0.8),size=4,hjust=.5)+
  geom_text(data=dat_N[dat_N$layer %in% c('Oa+e') & dat_N$R15_nit==max(dat_N$R15_nit),],col='black',
            aes(y=4.7,x='2023'),label=expression(paste('Year: ', italic('P'),'< 0.001')),
            position=position_dodge(0.8),size=4,hjust=.5)+
  geom_text(data=dat_N[dat_N$layer %in% c('Oa+e') & dat_N$R15_nit==max(dat_N$R15_nit),],col='black',
            aes(y=4.0,x='2023'),label='Treat x Year:  ns',
            position=position_dodge(0.8),size=4,hjust=0.45)+
  geom_text(data=q10.res[q10.res$P_R15_nit<=0.05 & q10.res$layer %in% c('Oa+e','0-10 cm'),],
            aes(x=year,y=R15_nit+1.5,
                label=ifelse(P_R15_nit<=0.001,'***',
                             ifelse(P_R15_nit<=0.01,'**',
                                    ifelse(P_R15_nit<=0.05,'*',
                                           ifelse(P_R15_nit<=0.10,paste('P = ',round(P_R15_nit,3)),''))))),
            color='black',size=6)+
  geom_text(data=q10.res[q10.res$P_R15_nit>0.05 & q10.res$layer %in% c('Oa+e','0-10 cm'),],
            aes(x=year,y=R15_nit+1.7,
                label=ifelse(P_R15_nit<=0.001,'***',
                             ifelse(P_R15_nit<=0.01,'**',
                                    ifelse(P_R15_nit<=0.05,'*',
                                           ifelse(P_R15_nit<=0.10,paste('P = ',round(P_R15_nit,3)),''))))),
            color='black',size=3.5,hjust=0.4)+
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(colour = "black",size=12,angle = 0,hjust = .5,vjust =1),
        axis.ticks = element_line(linewidth = .5),
        axis.ticks.length = unit(2,"mm"),
        prism.ticks.length = unit(1,"mm"),
        axis.text.y = element_text(colour = "black",size = 12,hjust=.5), 
        axis.title.x =element_text(size=12), axis.title.y=element_text(colour = "black",size=12,hjust=.5),
        legend.text = element_text(size=12,hjust=0), 
        panel.background = element_rect(fill =NA, linewidth=0.26,colour = "black", linetype = "solid",inherit.blank=T),
        panel.spacing.y = unit(6,'mm'),
        panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
        panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
        plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
        legend.key = element_rect(fill = NA), legend.title=element_blank(),
        legend.background = element_rect(fill = NA), 
        plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   #调整画图区域的间距，从上右下左调整
        strip.background = element_rect(fill = NA,color='black'), 
        strip.text = element_text(colour = "black",size = 10,hjust=.5),
        legend.position = 'none');R1517_nitri 


R1517N+R1517_nitri+
  plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
  theme(plot.tag.position = c(0, 0.98),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
        plot.tag = element_text(size = 12,vjust = 0,hjust=0,face="bold"))->Fig_6ab;Fig_6ab

setwd('D:/工作目录/202409/Manuscript_kai/Talk_20241113/Data and Code/Fig_6+Ext Dat fig_3_4_7_8')
ggsave('Fig_6ab_Rmin+Rnit in soils_new.pdf',Fig_6ab, width = 12, height = 6,device=cairo_pdf)


###重复测量方差分析
#R15_min
dat_N$plot<-substr(dat_N$subplot,1,1)
str(dat_N)

summary(aov(R15_min~treatment*year+Error(plot/year), 
            data=dat_N[dat_N$layer=='Oa+e',]))
t.test(R15_min~treatment,data=dat_N[dat_N$layer=='Oa+e',])


summary(aov(R15_min~treatment*year+Error(plot/year), 
            data=dat_N[dat_N$layer=='0-10 cm',]))
t.test(R15_min~treatment,data=dat_N[dat_N$layer=='0-10 cm',])



#R15_nit
summary(aov(R15_nit~treatment*year+Error(plot/year), 
            data=dat_N[dat_N$layer=='Oa+e',]))
t.test(R15_nit~treatment,data=dat_N[dat_N$layer=='Oa+e',])


summary(aov(R15_nit~treatment*year+Error(plot/year), 
            data=dat_N[dat_N$layer=='0-10 cm',]))
t.test(R15_nit~treatment,data=dat_N[dat_N$layer=='0-10 cm',])


#表层O层变化
datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                      range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                      n=length(na.omit(x)),
                      se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                      sd=round(sd(na.omit(x)),3))}
names(dat_N)
dat_N.yr<-aggregate(dat_N[,c('R15_min','R15_nit')],mean, na.rm = TRUE, 
                    by=list(dat_N$treatment,dat_N$year,dat_N$layer)) 
names(dat_N.yr)[1:3]<-c('treatment','year','layer')

unique(dat_N.yr$layer)

###年尺度Rmin的mean和se  from 2020 to 2021
datFUN((dat_N.yr$R15_min[dat_N.yr$layer=='Oa+e' & dat_N.yr$treatment=='warmed']-
          dat_N.yr$R15_min[dat_N.yr$layer=='Oa+e'  & dat_N.yr$treatment=='control'])/
         dat_N.yr$R15_min[dat_N.yr$layer=='Oa+e'  & dat_N.yr$treatment=='control'])

datFUN((dat_N.yr$R15_min[dat_N.yr$layer=='0-10 cm' & dat_N.yr$treatment=='warmed']-
          dat_N.yr$R15_min[dat_N.yr$layer=='0-10 cm' & dat_N.yr$treatment=='control'])/
         dat_N.yr$R15_min[dat_N.yr$layer=='0-10 cm' & dat_N.yr$treatment=='control'])

###年尺度Rnit的mean和se
datFUN((dat_N.yr$R15_nit[dat_N.yr$layer=='Oa+e' & dat_N.yr$treatment=='warmed']-
          dat_N.yr$R15_nit[dat_N.yr$layer=='Oa+e'  & dat_N.yr$treatment=='control'])/
         dat_N.yr$R15_nit[dat_N.yr$layer=='Oa+e'  & dat_N.yr$treatment=='control'])

datFUN((dat_N.yr$R15_nit[dat_N.yr$layer=='0-10 cm' & dat_N.yr$treatment=='warmed']-
          dat_N.yr$R15_nit[dat_N.yr$layer=='0-10 cm' & dat_N.yr$treatment=='control'])/
         dat_N.yr$R15_nit[dat_N.yr$layer=='0-10 cm' & dat_N.yr$treatment=='control'])


rm(dat_ij,Q10n,R1517N,R1517_nitri,res,Ron,i,j,k,n,sample_name)

#########利用样地的土壤温度计算年氮的矿化速率
##1)读取土壤BD数据
dat_soil<-read_excel('D:/Workspace/book/Qingyuan/soil_data/C_N_2023/土壤水分数据/Copy of soil weight moisture_2019-2023.xlsx',
                     sheet=1,col_names=T)

names(dat_soil)[6]<-'soil.water_dry'
dat_soil<-dat_soil[,c(1:6,8)]
str(dat_soil)

dat_soil$slayer<-ifelse(dat_soil$layer=='O','Oa+e',
                        ifelse(dat_soil$layer=='0-10cm','0-10 cm',
                               ifelse(dat_soil$layer=='10-20cm','10-20 cm','20-40 cm')))

#2)读取土壤温度数据
ghgs<-read.csv('D:/Workspace/book/Qingyuan/QY_Daily GHGs+NO flux_2019-2023.csv',header=T)
str(ghgs)
ghgs$month<-as.numeric(format(as.Date(ghgs$date),'%m'))
ghgs$year<-as.numeric(format(as.Date(ghgs$date),'%Y'))

ghgs<-ghgs[ghgs$month>= 5 & ghgs$month<=10,]

###q10all是min矿化,,,q10nitrifiction是nit硝化
dat_min.nit<-data.frame()

for(i in unique(q10$subplot)){
  for(j in unique(q10$year)){
    for(t in c('Oa+e','0-10 cm')){
      dat_ijt=ghgs[ghgs$year==j & ghgs$plot==substr(i,1,1),] 
      da_ijt=q10[q10$subplot==i & q10$year==j & q10$layer==t,]
      ds_ijt=dat_soil[dat_soil$subplot==i & dat_soil$year==j & dat_soil$slayer==t,]
      
      if(length(da_ijt$R0_min)==1){
        dat_ijt$subplot=i
        dat_ijt$layer=t
        
        dat_ijt$N.min=da_ijt$R0_min*exp((da_ijt$b_min)*dat_ijt$soil_temp)  # N_min per day per sqm2 
        dat_ijt$N.nit=da_ijt$R0_nit*exp((da_ijt$b_nit)*dat_ijt$soil_temp) # N_nit per day per sqm2 
        
        dat_ijt$N.min_m2=dat_ijt$N.min*ds_ijt$BD   #每平方米每日的产量mg N
        dat_ijt$N.nit_m2=dat_ijt$N.nit*ds_ijt$BD   #每平方米每日的产量mg N
        
        
        dat_min.nit<-rbind(dat_min.nit,dat_ijt)
      }
    }
  }
}

str(dat_min.nit)

dat_min.nit_yr<-aggregate(dat_min.nit[,c('N.min_m2','N.nit_m2')],FUN = sum, na.rm=T, 
                       by=list(dat_min.nit$subplot,dat_min.nit$layer,dat_min.nit$year,dat_min.nit$treatment))

names(dat_min.nit_yr)[1:4]<-c('subplot','layer','year','treatment')

myfun<-function(x){c(m=mean(na.omit(x)),se=sd(na.omit(x))/sqrt(length(na.omit(x))))}

### Extended Data Table 1
dat_min.nit_yr_treat<-summaryBy(data=dat_min.nit_yr,FUN = myfun,
                                dat_min.nit_yr[,c(5:6)]~treatment+year+layer)


###计算O，0-10之和
k=length(dat_min.nit_yr_treat$N.min_m2.m)

for(i in 2020:2023){
  for(j in c('control','warmed')){
    k=k+1
    dat_ij=dat_min.nit_yr_treat[dat_min.nit_yr_treat$treatment==j & dat_min.nit_yr_treat$year==i,]
    dat_min.nit_yr_treat[k,]<-NA   # add rows
    dat_min.nit_yr_treat$treatment[k]=j
    dat_min.nit_yr_treat$year[k]=i
    dat_min.nit_yr_treat$layer[k]='all'
    dat_min.nit_yr_treat$N.min_m2.m[k]= sum(dat_ij$N.min_m2.m)
    dat_min.nit_yr_treat$N.min_m2.se[k]= NA
    dat_min.nit_yr_treat$N.nit_m2.m[k]= sum(dat_ij$N.nit_m2.m)
    dat_min.nit_yr_treat$N.nit_m2.se[k]= NA
    
  }
}


# min from soils
datFUN((dat_min.nit_yr_treat$N.min_m2.m[dat_min.nit_yr_treat$treatment=='warmed' & dat_min.nit_yr_treat$layer=='all']-
  dat_min.nit_yr_treat$N.min_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='all'])/
  dat_min.nit_yr_treat$N.min_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='all'])


# nit from soils
datFUN((dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='warmed' & dat_min.nit_yr_treat$layer=='all']-
          dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='all'])/
         dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='all'])


# Res in the paper
# min from Oa+e
datFUN((dat_min.nit_yr_treat$N.min_m2.m[dat_min.nit_yr_treat$treatment=='warmed' & dat_min.nit_yr_treat$layer=='Oa+e']-
          dat_min.nit_yr_treat$N.min_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='Oa+e'])/
         dat_min.nit_yr_treat$N.min_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='Oa+e'])

# nit from Oa+e
datFUN((dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='warmed' & dat_min.nit_yr_treat$layer=='Oa+e']-
          dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='Oa+e'])/
         dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='Oa+e'])




# nit from soils
datFUN((dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='warmed' & dat_min.nit_yr_treat$layer=='all']-
          dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='all'])/
         dat_min.nit_yr_treat$N.nit_m2.m[dat_min.nit_yr_treat$treatment=='control' & dat_min.nit_yr_treat$layer=='all'])


remove(da_ijt,dn_ijt,ds_ijt,rep.aov,res)
###数值除以100为 kg N ha-1,,*10^-6*10^4



#########Fig_6.d_<N fluxes with N transformation functional genes>#########
###20230805huangkai
setwd("D:/Workspace/book/Qingyuan/soil_data")
list.files("D:/Workspace/book/Qingyuan/soil_data")
mydata<-read.csv("soil dat_2019-2023.csv",header=T);str(mydata)
mydata$Date<-gsub("-","/",mydata$Date)


.libPaths("D:/Workspace/Rtrial/R library 2024")
pacman::p_load("ggplot2","ggpmisc","ggpubr","ggprism","doBy","car",
               "patchwork","reshape2","dplyr","animation")  #加载多个包

mydata$Date<-as.character(as.Date(mydata$Date))
data<-mydata[,c(1:26)]
str(data)


#读取容重数据
names(dat_soil)
names(dat_soil)[c(1,3,4)]<-c("Date","soillayer","plot_ab")
which(names(dat_soil) %in% c('Date','plot_ab','soillayer','BD'))
dat_soil$Date<-as.character(dat_soil$Date)

str(dat_soil)
str(data)
data<-merge(data,dat_soil[,c(1,3,4,7)],by=c('Date','plot_ab','soillayer'),all.x = T)

str(data) #MBC MBN unit: mg C or N /kg

flux<-read.csv('D:/Workspace/book/Qingyuan/soil_data/Daily flux_plotAB in QY 2018-2023.csv',header = T)
flux$Date<-as.character(as.Date(flux$date))
str(flux)
flux$NON2O<-flux$no_N+flux$n2o_N
flux$WFPS<-flux$vwc/(1-0.7/2.65)

flux$n2_N<-flux$n2o_N*(flux$WFPS*0.13-1.10)

unique(flux$plot_ab)

data$treatment<-ifelse(substr(data$plot_ab,1,1) %in% c(1,4,5),'warmed','control')
unique(data$soillayer)

str(data)

######将1，2，3日土壤温湿度和排放数据都列举出来
data$soil_temp1=NA
data$WFPS1=NA
data$NO_N1=NA
data$N2O_N1=NA
data$N2_N1=NA
data$CO2_C1=NA

data$soil_temp2=NA
data$WFPS2=NA
data$NO_N2=NA
data$N2O_N2=NA
data$N2_N2=NA
data$CO2_C2=NA

data$soil_temp3=NA
data$WFPS3=NA
data$NO_N3=NA
data$N2O_N3=NA
data$N2_N3=NA
data$CO2_C3=NA

###3日平均
data$soil_temp=NA
data$WFPS=NA
data$NO_N=NA
data$N2O_N=NA
data$N2_N=NA
data$CO2_C=NA

str(flux)


unique(data$Date)
###3日均值
for(k in 1:length(data$Date)){
  i=as.Date(data$Date[k])
  j=data$plot_ab[k]
  
  ###第1日
  data$soil_temp1[k]=mean(na.omit(flux$soil.temp[flux$date %in% as.character(c(i)) & flux$plot_ab==j]))
  data$WFPS1[k]=mean(na.omit(flux$WFPS[flux$date %in% as.character(c(i)) & flux$plot_ab==j]))
  data$NO_N1[k]=mean(na.omit(flux$no_N[flux$date %in% as.character(c(i)) & flux$plot_ab==j]))
  data$N2O_N1[k]=mean(na.omit(flux$n2o_N[flux$date %in% as.character(c(i)) & flux$plot_ab==j]))
  data$N2_N1[k]=mean(na.omit(flux$n2_N[flux$date %in% as.character(c(i)) & flux$plot_ab==j]))
  data$CO2_C1[k]=mean(na.omit(flux$co2_C[flux$date %in% as.character(c(i)) & flux$plot_ab==j]))
  
  ###第2日
  data$soil_temp2[k]=mean(na.omit(flux$soil.temp[flux$date %in% as.character(c(i-1)) & flux$plot_ab==j]))
  data$WFPS2[k]=mean(na.omit(flux$WFPS[flux$date %in% as.character(c(i-1)) & flux$plot_ab==j]))
  data$NO_N2[k]=mean(na.omit(flux$no_N[flux$date %in% as.character(c(i-1)) & flux$plot_ab==j]))
  data$N2O_N2[k]=mean(na.omit(flux$n2o_N[flux$date %in% as.character(c(i-1)) & flux$plot_ab==j]))
  data$N2_N2[k]=mean(na.omit(flux$n2_N[flux$date %in% as.character(c(i-1)) & flux$plot_ab==j]))
  data$CO2_C2[k]=mean(na.omit(flux$co2_C[flux$date %in% as.character(c(i-1)) & flux$plot_ab==j]))
  
  ###第3日
  data$soil_temp3[k]=mean(na.omit(flux$soil.temp[flux$date %in% as.character(c(i-2)) & flux$plot_ab==j]))
  data$WFPS3[k]=mean(na.omit(flux$WFPS[flux$date %in% as.character(c(i-2)) & flux$plot_ab==j]))
  data$NO_N3[k]=mean(na.omit(flux$no_N[flux$date %in% as.character(c(i-2)) & flux$plot_ab==j]))
  data$N2O_N3[k]=mean(na.omit(flux$n2o_N[flux$date %in% as.character(c(i-2)) & flux$plot_ab==j]))
  data$N2_N3[k]=mean(na.omit(flux$n2_N[flux$date %in% as.character(c(i-2)) & flux$plot_ab==j]))
  data$CO2_C3[k]=mean(na.omit(flux$co2_C[flux$date %in% as.character(c(i-2)) & flux$plot_ab==j]))
  
  ###3日均值
  data$soil_temp[k]=mean(na.omit(flux$soil.temp[flux$date %in% as.character(c(i,i-1,i-2)) & flux$plot_ab==j]))
  data$WFPS[k]=mean(na.omit(flux$WFPS[flux$date %in% as.character(c(i,i-1,i-2)) & flux$plot_ab==j]))
  data$NO_N[k]=mean(na.omit(flux$no_N[flux$date %in% as.character(c(i,i-1,i-2)) & flux$plot_ab==j]))
  data$N2O_N[k]=mean(na.omit(flux$n2o_N[flux$date %in% as.character(c(i,i-1,i-2)) & flux$plot_ab==j]))
  data$N2_N[k]=mean(na.omit(flux$n2_N[flux$date %in% as.character(c(i,i-1,i-2)) & flux$plot_ab==j]))
  data$CO2_C[k]=mean(na.omit(flux$co2_C[flux$date %in% as.character(c(i,i-1,i-2)) & flux$plot_ab==j]))
  
}


###计算容重加权后的
res_dat<-data.frame()
k=1
for(i in unique(data$plot_ab)){
  for(j in unique(data$Date)){
    dat_ij<-data[data$plot_ab==i & data$Date==j,]
    res_dat[k,1]=j
    res_dat[k,2]=i
    res_dat[k,3]=sum(dat_ij$MBC*dat_ij$BD)  #mg C m-2
    res_dat[k,4]=sum(dat_ij$MBN*dat_ij$BD)  #mg C m-2
    res_dat[k,5]=sum(dat_ij$AOA*dat_ij$BD)  #mg C m-2
    res_dat[k,6]=sum(dat_ij$AOB*dat_ij$BD)  #mg C m-2
    res_dat[k,7]=sum(dat_ij$nirK*dat_ij$BD)  #mg C m-2
    res_dat[k,8]=sum(dat_ij$nirS*dat_ij$BD)  #mg C m-2
    res_dat[k,9]=sum(dat_ij$nosZ*dat_ij$BD)  #mg C m-2
    res_dat[k,10]=sum(dat_ij$DOC*dat_ij$BD)  #mg C m-2
    res_dat[k,11]=ifelse(substr(i,1,1) %in% c(1,4,5),'warmed','control')  #mg C m-2
    k=k+1
  }
}
names(res_dat)<-c('Date','plot_ab','MBC','MBN','AOA','AOB','nirK','nirS','nosZ','DOC','treatment')

###换算成干土丰度
data$AOB<-data$AOB*(1+data$soil.water_dry)
data$AOA<-data$AOA*(1+data$soil.water_dry)
data$nirS<-data$nirS*(1+data$soil.water_dry)
data$nirK<-data$nirK*(1+data$soil.water_dry)
data$nirSK<-data$nirSK*(1+data$soil.water_dry)
data$nosZ<-data$nosZ*(1+data$soil.water_dry)

###描述性文字信息
datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                      range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                      n=length(na.omit(x)),
                      se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                      sd=round(sd(na.omit(x)),3))}
datFUN(data$AOA)
datFUN(data$nirK)
datFUN(data$nirS)
datFUN(data$nosZ)

datFUN(data$MBC[data$soillayer=='O'])
datFUN(data$MBN[data$soillayer=='O'])

###表层土壤MBC和MBN和DOC的变化
datFUN(res_dat$MBC[res_dat$treatment=='control'])
(mean(na.omit(res_dat$MBC[res_dat$treatment=='warmed']))-
    mean(na.omit(res_dat$MBC[res_dat$treatment=='control'])))/
mean(na.omit(res_dat$MBC[res_dat$treatment=='control']))

datFUN(res_dat$MBN[res_dat$treatment=='control'])
(mean(na.omit(res_dat$MBN[res_dat$treatment=='warmed']))-
    mean(na.omit(res_dat$MBN[res_dat$treatment=='control'])))/
  mean(na.omit(res_dat$MBN[res_dat$treatment=='control']))

datFUN(res_dat$DOC[res_dat$treatment=='control'])
(mean(na.omit(res_dat$DOC[res_dat$treatment=='warmed']))-
    mean(na.omit(res_dat$DOC[res_dat$treatment=='control'])))/
  mean(na.omit(res_dat$DOC[res_dat$treatment=='control']))

(mean(na.omit(res_dat$AOB[res_dat$treatment=='warmed']))-
    mean(na.omit(res_dat$AOB[res_dat$treatment=='control'])))/
  mean(na.omit(res_dat$AOB[res_dat$treatment=='control']))


(mean(na.omit(res_dat$AOA[res_dat$treatment=='warmed']))-
    mean(na.omit(res_dat$AOA[res_dat$treatment=='control'])))/
  mean(na.omit(res_dat$AOA[res_dat$treatment=='control']))



###不同土层含量不同（mg C or N kg-1）_______:>20231101
datFUN(data$MBC[data$soillayer=='O' & data$treatment=='control'])
datFUN(data$MBN[data$soillayer=='O'& data$treatment=='control'])
datFUN(data$DOC[data$soillayer=='O'& data$treatment=='control'])

datFUN(data$MBC[data$soillayer=='0-10cm' & data$treatment=='control'])
datFUN(data$MBN[data$soillayer=='0-10cm'& data$treatment=='control'])
datFUN(data$DOC[data$soillayer=='0-10cm'& data$treatment=='control'])

(mean(na.omit(data$DOC[data$soillayer=='O'& data$treatment=='warmed']))-mean(na.omit(data$DOC[data$soillayer=='O'& data$treatment=='control'])))/
  mean(na.omit(data$DOC[data$soillayer=='O'& data$treatment=='control']))

(mean(na.omit(data$DOC[data$soillayer=='0-10cm'& data$treatment=='warmed']))-mean(na.omit(data$DOC[data$soillayer=='0-10cm'& data$treatment=='control'])))/
  mean(na.omit(data$DOC[data$soillayer=='0-10cm'& data$treatment=='control']))


(mean(na.omit(data$DON[data$soillayer=='O'& data$treatment=='warmed']))-mean(na.omit(data$DON[data$soillayer=='O'& data$treatment=='control'])))/
  mean(na.omit(data$DON[data$soillayer=='O'& data$treatment=='control']))

(mean(na.omit(data$DON[data$soillayer=='0-10cm'& data$treatment=='warmed']))-mean(na.omit(data$DON[data$soillayer=='0-10cm'& data$treatment=='control'])))/
  mean(na.omit(data$DON[data$soillayer=='0-10cm'& data$treatment=='control']))


(mean(na.omit(data$AOA[data$soillayer=='O' & data$treatment=='warmed']))-
    mean(na.omit(data$AOA[data$soillayer=='O' & data$treatment=='control'])))/
  mean(na.omit(data$AOA[data$soillayer=='O' & data$treatment=='control']))

(mean(na.omit(data$AOB[data$soillayer=='O' & data$treatment=='warmed']))-
    mean(na.omit(data$AOB[data$soillayer=='O' & data$treatment=='control'])))/
  mean(na.omit(data$AOB[data$soillayer=='O' & data$treatment=='control']))


(mean(na.omit(data$AOB[data$soillayer=='0-10cm' & data$treatment=='warmed']))-
    mean(na.omit(data$AOB[data$soillayer=='0-10cm' & data$treatment=='control'])))/
  mean(na.omit(data$AOB[data$soillayer=='0-10cm' & data$treatment=='control']))


(mean(na.omit(data$nirS[data$soillayer=='O' & data$treatment=='warmed']))-
    mean(na.omit(data$nirS[data$soillayer=='O' & data$treatment=='control'])))/
  mean(na.omit(data$nirS[data$soillayer=='O' & data$treatment=='control']))

(mean(na.omit(data$MBC[data$soillayer=='O' & data$treatment=='warmed']))-
    mean(na.omit(data$MBC[data$soillayer=='O' & data$treatment=='control'])))/
  mean(na.omit(data$MBC[data$soillayer=='O' & data$treatment=='control']))

(mean(na.omit(data$MBN[data$soillayer=='O' & data$treatment=='warmed']))-
    mean(na.omit(data$MBN[data$soillayer=='O' & data$treatment=='control'])))/
  mean(na.omit(data$MBN[data$soillayer=='O' & data$treatment=='control']))


(mean(na.omit(data$MBC[data$soillayer=='0-10cm' & data$treatment=='warmed']))-
    mean(na.omit(data$MBC[data$soillayer=='0-10cm' & data$treatment=='control'])))/
  mean(na.omit(data$MBC[data$soillayer=='0-10cm' & data$treatment=='control']))

(mean(na.omit(data$MBN[data$soillayer=='0-10cm' & data$treatment=='warmed']))-
    mean(na.omit(data$MBN[data$soillayer=='0-10cm' & data$treatment=='control'])))/
  mean(na.omit(data$MBN[data$soillayer=='0-10cm' & data$treatment=='control']))

######有机层和0-10cm EOC含量的变化（+土壤烘干含水率）
names(data)[c(8,10)]<-c('EOC','EON')

ggplot(data,aes(Date,EOC,col=treatment))+geom_point(position = position_dodge(0.8),pch=1,size=1.6)+
  facet_wrap(.~soillayer,nrow = 2,scales = 'free_y')+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
  geom_boxplot(position = position_dodge(0.8),width=0.6,lwd=0.3,outlier.shape = 1,alpha=0.65,outlier.size=0.5)+ #geom_violin(fill=NA,width=0.25)+
  stat_mean(position = position_dodge(0.8),pch=16,size=1)+
  stat_compare_means(label = "p.signif", 
                     method = "t.test",hide.ns = F,label.y.npc =0.9,label.x =c(1,2),size=4)

names(data)
ggplot(data,aes(Date,log(AOB/AOA,10),col=treatment))+geom_point(position = position_dodge(0.8),pch=1,size=1.6)+
  facet_wrap(.~soillayer,nrow = 2,scales = 'free_y')+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
  geom_boxplot(position = position_dodge(0.8),width=0.6,lwd=0.3,outlier.shape = 1,alpha=0.65,outlier.size=0.5)+
  stat_mean(position = position_dodge(0.8),pch=16,size=1)+
  stat_compare_means(label = "p.signif", 
                     method = "t.test",hide.ns = F,label.y.npc =0.9,label.x =c(1,2),size=4)

str(data)
ggplot(data,aes(Date,log(nirS,10),col=treatment))+geom_point(position = position_dodge(0.8),pch=1,size=1.6)+
  facet_wrap(.~soillayer,nrow = 2,scales = 'free_y')+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
  geom_boxplot(position = position_dodge(0.8),width=0.6,lwd=0.3,outlier.shape = 1,alpha=0.65,outlier.size=0.5)+
  stat_mean(position = position_dodge(0.8),pch=16,size=1)+
  stat_compare_means(label = "p.signif", 
                     method = "t.test",hide.ns = F,label.y.npc =0.9,label.x =c(1,2),size=4)

ggplot(data,aes(log(AOB/AOA),log(NO_N),col=treatment))+geom_point(position = position_dodge(0.8),pch=1,size=1.6)+
  facet_wrap(.~soillayer,nrow = 2,scales = 'free_y')+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  stat_poly_eq(formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
               eq.with.lhs = "italic(y)~`=`~",   #????y"????ʽ         
               aes(label = paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
               parse = TRUE,label.x.npc = "left", label.y.npc = c(0.47,0.87),size = 5)

ggplot(data,aes(log(AOB/AOA),log(NO_N),col=treatment))+geom_point(position = position_dodge(0.8),pch=1,size=1.6)+
  facet_wrap(.~soillayer,nrow = 2,scales = 'free_y')+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  stat_poly_eq(formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
               eq.with.lhs = "italic(y)~`=`~",   #????y"????ʽ         
               aes(label = paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
               parse = TRUE,label.x.npc = "left", label.y.npc = c(0.47,0.87),size = 5)









ggplot(data,aes(Date,log(MBC,10),col=treatment))+geom_point(position = position_dodge(0.8),pch=1,size=1.6)+
  facet_wrap(.~soillayer,nrow = 2,scales = 'free_y')+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
  geom_boxplot(position = position_dodge(0.8),width=0.6,lwd=0.3,outlier.shape = 1,alpha=0.65,outlier.size=0.5)+
  stat_mean(position = position_dodge(0.8),pch=16,size=1)+
  stat_compare_means(label = "p.signif", 
                     method = "t.test",hide.ns = F,label.y.npc =0.9,label.x =c(1,2),size=4)

###将1，2，3日的通量数据都摘出
names(data)

dat1<-data[,c("Date","plot_ab","soillayer",
              "AOB","AOA","nirK","nirS","nirSK","nosZ","nosZ_nirSK",
              "MBC","MBN","EOC","TDN","EON",
              "NH4_N","NO3_N","soil.water_dry","pH",
              "delta_13C","C_con","delta_15N","N_con","C_N","BD","treatment",
              "soil_temp1","WFPS1","NO_N1","N2O_N1","N2_N1","CO2_C1")]

names(dat1)[27:32]<-c("soil_temp","WFPS","NO_N","N2O_N","N2_N","CO2_C")
dat1$labs<-'Day1'

dat2<-data[,c("Date","plot_ab","soillayer","AOB","AOA","nirK","nirS","nirSK","nosZ","nosZ_nirSK","MBC","MBN","EOC","TDN","EON",
              "NH4_N","NO3_N","soil.water_dry","pH","delta_13C","C_con","delta_15N","N_con","C_N","BD","treatment",
              "soil_temp2","WFPS2","NO_N2","N2O_N2","N2_N2","CO2_C2")]
names(dat2)[27:32]<-c("soil_temp","WFPS","NO_N","N2O_N","N2_N","CO2_C")
dat2$labs<-'Day2'

dat3<-data[,c("Date","plot_ab","soillayer","AOB","AOA","nirK","nirS","nirSK","nosZ","nosZ_nirSK","MBC","MBN","EOC","TDN","EON",
              "NH4_N","NO3_N","soil.water_dry","pH","delta_13C","C_con","delta_15N","N_con","C_N","BD","treatment",
              "soil_temp3","WFPS3","NO_N3","N2O_N3","N2_N3","CO2_C3")]
names(dat3)[27:32]<-c("soil_temp","WFPS","NO_N","N2O_N","N2_N","CO2_C")
dat3$labs<-'Day3'


data<-data[,c("Date","plot_ab","soillayer","AOB","AOA","nirK","nirS","nirSK","nosZ","nosZ_nirSK","MBC","MBN","EOC","TDN","EON",
              "NH4_N","NO3_N","soil.water_dry","pH","delta_13C","C_con","delta_15N","N_con","C_N","BD","treatment",
              "soil_temp","WFPS","NO_N","N2O_N","N2_N","CO2_C")]

names(dat1)
names(dat2)
names(dat3)

library (tidyverse)

#put all data frames into list
df_list <- list(dat1, dat2, dat3)

#merge all data frames into list
Reduce(function(x, y) merge(x, y, all= TRUE ), df_list)->dat

###linear regression

data$nirSK_nosZ<-data$nirSK/data$nosZ
data$NH4_NO3<-data$NH4_N/data$NO3_N
data$MBC_MBN<-data$MBC/data$MBN
data$NH4NO3<-data$NH4_N+data$NO3_N

dat$nirSK_nosZ<-dat$nirSK/dat$nosZ
dat$NH4_NO3<-dat$NH4_N/dat$NO3_N
dat$MBC_MBN<-dat$MBC/dat$MBN
dat$NH4NO3<-dat$NH4_N+dat$NO3_N

mythem_sci<-theme(plot.subtitle = element_text(vjust = 0.5,hjust = 0.5,size=8.5), plot.caption = element_text(vjust = 1), 
                  plot.title = element_text(vjust = 0.5,hjust = 0.5,size=8.5), 
                  axis.text.x = element_text(colour = "black",size=8.5,angle = 0,hjust = .5,vjust =1),
                  axis.ticks = element_line(linewidth = .15),
                  axis.ticks.length = unit(1,"mm"),
                  prism.ticks.length = unit(0.5,"mm"),
                  axis.text.y = element_text(colour = "black",size = 8.5), 
                  axis.title.x =element_text(size=8.5), axis.title.y=element_text(colour = "black",size=8.5,vjust =-1.5),
                  legend.text = element_text(size=8.5,vjust = 0.5), legend.title =element_text(size=8.5,vjust = 0.5),
                  # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
                  panel.background = element_rect(fill =NA, linewidth=0.05,colour = "black", linetype = "solid",
                                                  inherit.blank=T),
                  panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
                  panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
                  plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
                  legend.key = element_rect(fill = NA), 
                  legend.background = element_rect(fill = NA), 
                  plot.margin = unit(c(2,2,1,1),'mm'),   #调整画图区域的间距，从上右下左调整
                  strip.background = element_rect(fill = NA,color='black',linewidth=0.05), 
                  strip.text = element_text(colour = "black",size = 7,hjust=.5),
                  legend.position = "none")

dat0<-dat   #dat0前3天的数据


###########绘制person correlation_3日均值排放
#corplot
data$amoA_nirSK=(data$AOB+data$AOA)/data$nirSK
data$AOB_AOA=data$AOB/data$AOA
getwd()   # check the save place

#define the save olace
setwd('D:/工作目录/202409/Manuscript_kai/Talk_20241113/Data and Code/Fig_6+Ext Dat fig_3_4_7_8')

#### corplot
{dat=data
  
  ###将AOA，AOB, nirS, nirK, nosZ对数化！！！
  dat$AOB=log(dat$AOB,10)
  dat$AOA=log(dat$AOA,10)
  dat$nirS=log(dat$nirS,10)
  dat$nirK=log(dat$nirK,10)
  dat$nirSK=log(dat$nirSK,10)
  dat$nosZ=log(dat$nosZ,10)
  
  str(dat)
  dat$treatment<-ifelse(dat$treatment=='control',0,1)
  dat$treatment<-as.numeric(dat$treatment)
  str(dat)
  
  dat$soillayer<-as.numeric(ifelse(dat$soillayer=='O',0,1))
  
  dat$NON2O<-dat$NO_N+dat$N2O_N
  
  
  names(dat)
  
  #####O层
  dati<-dat[dat$soillayer==0,c("soil_temp",'soil.water_dry',"EOC","EON","NH4_N","NO3_N",
                               "MBC","MBN","AOB","AOA","nirS","nirK" ,"nosZ","nirSK",'AOB_AOA',"amoA_nirSK","nirSK_nosZ","CO2_C",        
                               "NO_N","N2O_N","N2_N","treatment")]
  
  names(dati)[c(1:2,5:6,14:22)]<-c("Temp","Moisture","=paste(NH[4]^'+')","=paste(NO[3]^' -')",
                                   "nirK+nirS","AOB/AOA","amoA/(nirK+nirS)","(nirK+nirS)/nosZ","=paste('Flux',' '[CO[2]])",
                                   "=paste('Flux',' '[NO])","=paste('Flux',' '[N[2]][O])","=paste('Flux',' '[N[2]])","Warming")
  mcor <- cor(dati, use="complete.obs")   #选择需要查看的因子
  
  ######计算p值
  library(vegan)
  r=mcor#计算相关系数矩阵
  r1=r[r>=0] #找出大于等于0的相关系数
  r2=r[r<0]  #找出小于0的相关系数
  n=nrow(dati) #算样本量
  p1=(1-pt(r1*sqrt((n-2)/(1-r1^2)),n-2))*2 #计算大于等于0的相关系数的p值
  p2=pt(r2*sqrt((n-2)/(1-r2^2)),n-2)*2 #计算小0的相关系数的p值
  pvalue=r  #构建pvalue矩阵，让它首先等于相关系数矩阵
  pvalue[r>=0]=p1 #pvalue矩阵中r大于等于0地方等于p1
  pvalue[r<0]=p2  #pvalue矩阵中r小于0地方等于p2
  pvalue
  round(pvalue,2)   #p值矩阵
  
  
  ###0-10cm
  dati1<-dat[dat$soillayer==1,c("soil_temp",'soil.water_dry',"EOC","EON","NH4_N","NO3_N",
                                "MBC","MBN","AOB","AOA","nirS","nirK" ,"nosZ","nirSK",'AOB_AOA',"amoA_nirSK","nirSK_nosZ","CO2_C",        
                                "NO_N","N2O_N","N2_N","treatment")]
  
  names(dati1)[c(1:2,5:6,14:22)]<-c("Temp","Moisture","=paste(NH[4]^'+')","=paste(NO[3]^' -')",
                                   "nirK+nirS","AOB/AOA","amoA/(nirK+nirS)","(nirK+nirS)/nosZ","=paste('Flux',' '[CO[2]])",
                                   "=paste('Flux',' '[NO])","=paste('Flux',' '[N[2]][O])","=paste('Flux',' '[N[2]])","Warming")
  mcor1 <- cor(dati1, use="complete.obs")   #选择需要查看的因子
  
  ######计算p值
  library(vegan)
  r=mcor1#计算相关系数矩阵
  r1=r[r>=0] #找出大于等于0的相关系数
  r2=r[r<0]  #找出小于0的相关系数
  n=nrow(dati1) #算样本量
  p1=(1-pt(r1*sqrt((n-2)/(1-r1^2)),n-2))*2 #计算大于等于0的相关系数的p值
  p2=pt(r2*sqrt((n-2)/(1-r2^2)),n-2)*2 #计算小0的相关系数的p值
  pvalue1=r  #构建pvalue矩阵，让它首先等于相关系数矩阵
  pvalue1[r>=0]=p1 #pvalue矩阵中r大于等于0地方等于p1
  pvalue1[r<0]=p2  #pvalue矩阵中r小于0地方等于p2
  pvalue1
  round(pvalue1,2)   #p值矩阵
  
  mcor
  library(corrplot)
  
  #创建颜色梯度"blue","white","red"
  col3<-colorRampPalette(c("dodgerblue4","white","brown3"))
  
  
  pdf("corplot_2019-2023_all.pdf",width=10,height=12)
  # par(mfrow=c(2,2))#散点图边距，字体为myFont 只能设置成2行2列
  layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(2,2),height=c(2,2))
  
  cor.plot<-corrplot(mar=c(0.5,0.1,0.1,0.1),corr=mcor,p.mat = pvalue,insig = "label_sig",
                     col = col3(10),sig.level = c(0.001,0.01,0.05),pch.cex =1.1, #mcor是相关系数矩阵
                     method = "circle",cl.cex=1.0,  #cl.pos="n"无图例,r右图例
                     tl.col="black",tl.cex = 0.8,diag=FALSE,   #diag去除自相关
                     number.cex = 16,number.font = 8,
                     type = "lower")
  mtext("(a) Organic layer",side=3,line = -1.5,cex =1.0,adj=0.57,padj = 0.4,col = 1 )#输出文本，*在下边
  mtext("Person's R correlations",side=1,line = 3.5,cex =1.2,adj=0.6,padj = 0.4,col = 1 )#输出文本，*在下边
  
  
  cor.plot1<-corrplot(mar=c(0.5,0.1,0.1,0.1),corr=mcor1,p.mat = pvalue1,insig = "label_sig",
                      sig.level = c(0.001,0.01,0.05),pch.cex =1.1, #mcor是相关系数矩阵
                      col = col3(10),method = "circle",cl.cex=1.0,   #cl.pos="n"无图例,r右图例
                      tl.col="black",tl.cex = 0.8,diag=FALSE,   #diag去除自相关
                      number.cex = 16,number.font = 8,
                      type = "lower")
  mtext("(b) 0-10 cm",side=3,line = -1.5,cex =1.0,adj=0.57,padj = 0.4,col = 1 )#输出文本，*在下边
  mtext("Person's R correlations",side=1,line = 3.5,cex =1.0,adj=0.6,padj = 0.4,col = 1 )#输出文本，*在下边
  
  
  # mtext("*",side=1,line = 3.6,cex =2,adj=0.42,padj = 0.4,col = 1 )#输出文本，*在下边
  # mtext("P.value<0.05",side=1,line = 3.6,cex =1,adj=0.52,padj = 0.4,col = 1 )#输出文本，"表示p<0.05"在下边
  dev.off()
}

#### Day1, Day2（当天）

{ dat=dat0[dat0$labs %in%c('Day1','Day2'),]
  # names(dat)
  # dat<-aggregate(dat[,c(4:37,40:41)], by=list(dat$Date,dat$soillayer,dat$plot,dat$year,dat$treatment),mean, na.rm = TRUE)
  # 
  # names(dat)[1:5]<-c('Date','soillayer','plot','year','treatment')
  str(dat)
  dat$AOB_AOA=log(dat$AOB/dat$AOA,10)
  dat$nirSK_nosZ=log((dat$nirS+dat$nirK)/dat$nosZ,10)
  dat$NH4_NO3=log(dat$NH4_N/dat$NO3_N,10)
  
  dat$AOB=log(dat$AOB,10)
  dat$AOA=log(dat$AOA,10)
  dat$nirS=log(dat$nirS,10)
  dat$nirK=log(dat$nirK,10)
  dat$nirSK=log(dat$nirSK,10)
  dat$nosZ=log(dat$nosZ,10)
  dat$EOC=log(dat$EOC,10)
  # dat$DON=log(dat$DON,10)
  dat$TDN=log(dat$TDN,10)
  
  dat$MBC=log(dat$MBC,10)
  dat$MBN=log(dat$MBN,10)
  dat$MBC_MBN=log(dat$MBC_MBN,10)
  dat$EON_MBN=log(dat$EON/dat$MBN,10)
  dat$TIN_MBN=log(dat$NH4NO3/dat$MBN,10)
  
  
  dat$NH4_N=log(dat$NH4_N,10)
  dat$NO3_N=log(dat$NO3_N,10)
  dat$NH4NO3=log(dat$NH4NO3,10)
  
  
  dat$soil.water_dry=log(dat$soil.water_dry,10)
  dat$soil_temp=log(dat$soil_temp,10)
  
  # dat$min=log(dat$min,10)
  # dat$nit=log(dat$nit,10)
  dat$CO2_C=log(dat$CO2_C,10)
  # dat$NO_N=log(dat$NO_N,10)
  # dat$N2O_N=log(dat$N2O_N+20,10)
  # dat$N2_N=log(dat$N2_N,10)

  unique(dat$soillayer)
  ###df_c0 
  unique(dat$Date)
  names(dat)
  
  names(dat)
  dat_c0<-as.matrix(dat[dat$soillayer=="O" & dat$treatment=='control',
                        c("AOB","AOA","nirK","nirS","nirSK","nosZ","MBC", "MBN","MBC_MBN","EOC","EON","EON_MBN",
                          "NH4_N","NO3_N","soil.water_dry","soil_temp","NO_N","N2O_N","N2_N","CO2_C","NH4_NO3",
                          "NH4NO3","AOB_AOA","nirSK_nosZ")]) #利用as.matrix()将所需数据集转换为matrix格式，才可在corrplot中跑
  
  dat_c0=data.frame(scale(dat_c0))#数据标准化
  
  # dat_c0<-cor(dat_c0,method="spearman") 
  
  names(dat_c0)
  
  library(psych);library(magrittr);library(grid)
  
  pp<-corr.test(dat_c0[,c('soil_temp','soil.water_dry','EOC','EON','NH4_N','NO3_N',
                          'AOB','AOA','nirK','nirS','nosZ',"AOB_AOA","nirSK_nosZ",
                          'MBN','MBC_MBN','EON_MBN','CO2_C')],
                dat_c0[,c("NO_N","N2O_N")],method="pearson",adjust = "fdr")  
  
  
  
  cor <- pp$r
  pvalue <- pp$p
  
  df_c0 <- melt(cor) %>% mutate(pvalue=melt(pvalue)[,3],
                                p_signif=symnum(pvalue, corr = FALSE, na = FALSE,  
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                symbols = c("***", "**", "*", "+", " "))) %>% 
    set_colnames(c("env","obs","r","p","p_signif"))
  
  str(df_c0)
  
  names(dat_c0)[c(16:18)]
  df_c0<-within(df_c0,obs<-factor(obs,levels=c('gasN','N2_N','N2O_N','NO_N')))
  
  names(dat_c0)[c(1:6,8:15,19:20)]
  # df_c0<-within(df_c0,env<-factor(env,levels=c('soil_temp','soil.water_dry','DOC','NH4_N','NO3_N','NH4NO3',
  #                                              'AOA','nirK','nirS','nirSK','nosZ',"nirSK_nosZ",
  #                                              'MBC','MBN','CO2_C',"min","nit")))
  
  ###df_w0
  dat_w0<-as.matrix(dat[dat$soillayer=='O' & dat$treatment=='warmed',
                        c("AOB","AOA","nirK","nirS","nosZ", "MBN","EOC","EON","TDN",
                          "NH4_N","NO3_N","soil.water_dry","soil_temp","NO_N","N2O_N","N2_N","CO2_C","NH4NO3","MBC_MBN","EON_MBN",
                          "AOB_AOA","nirSK_nosZ")])  #利用as.matrix()将所需数据集转换为matrix格式，才可在corrplot中跑
  
  dat_w0=data.frame(scale(dat_w0))#数据标准化
  
  # dat_c0<-cor(dat_c0,method="spearman") 
  
  names(dat_w0)
  
  pp<-corr.test(dat_w0[,c('soil_temp','soil.water_dry','EOC','EON','NH4_N','NO3_N',
                          'AOB','AOA','nirK','nirS','nosZ',"AOB_AOA","nirSK_nosZ",
                          'MBN','MBC_MBN','EON_MBN','CO2_C')],
                dat_w0[,c("NO_N","N2O_N")],method="pearson",adjust = "fdr")  
  cor <- pp$r
  pvalue <- pp$p
  
  df_w0 <- melt(cor) %>% mutate(pvalue=melt(pvalue)[,3],
                                p_signif=symnum(pvalue, corr = FALSE, na = FALSE,  
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                symbols = c("***", "**", "*", "+", " "))) %>% 
    set_colnames(c("env","obs","r","p","p_signif"))
  
  str(df_w0)
  
  names(dat_w0)[c(16:18)]
  df_w0<-within(df_w0,obs<-factor(obs,levels=c('N2O_N','NO_N')))
  
  names(dat_c0)[c(1:15,19:20)]
  # df_w0<-within(df_w0,env<-factor(env,levels=c('soil_temp','soil.water_dry','DOC','NH4_N','NO3_N','NH4NO3',
  #                                              'AOA','nirK','nirS','nirSK','nosZ',"nirSK_nosZ",
  #                                              'MBC','MBN','CO2_C',"min","nit")))
  
  
  df_w0$treatment='warmed'
  df_c0$treatment='control'
  
  df<-rbind(df_c0,df_w0)
  
  unique(df$env)
  
  df$obs_tr<-paste(df$obs,df$treatment,sep='-')
  
  unique(df$obs_tr)
  df<-within(df,obs_tr<-factor(obs_tr,levels=c("N2O_N-warmed","N2O_N-control",
                                               "NO_N-warmed","NO_N-control")))
  
  unique(df$env)
  
  df<-within(df,env<-factor(env,levels=c('soil_temp','soil.water_dry','NH4_N','NO3_N','NH4NO3','EOC','EON',
                                               'AOB','AOA','nirK','nirS','nirSK','nosZ','AOB_AOA',"nirSK_nosZ",
                                               'MBC','MBN','MBC_MBN','EON_MBN','CO2_C')))
  
  
  unique(df$env)
  unique(df$obs)
  
  # df<-df[df$obs_tr %in% c("N2_N-warmed","N2_N-control","N2O_N-warmed","N2O_N-control","NO_N-warmed","NO_N-control"),]
  p<-ggplot(df,aes(env,obs_tr,col=r,fill=r,group=treatment))+
    geom_tile(color="grey80",fill="white",size=0.3)+
    geom_point(aes(size =abs(r)),shape=21)+facet_grid('O layer'~.)+
    geom_text(data=df[df$p_signif!='+',],aes(label=p_signif),size=6,color="white",hjust=0.5,vjust=0.7)+
    # geom_text(data=df[df$p_signif=='+',],aes(label=p_signif),size=3,color="white",hjust=0.5,vjust=0.5)+
    labs(x = NULL,y = NULL,color=NULL,fill=NULL)+
    scale_color_gradientn(limits=c(-1,1),colours = rev(RColorBrewer::brewer.pal(11,"RdBu")),name=expression(atop("Person's R","correlations")))+
    scale_fill_gradientn(limits=c(-1,1),colours = rev(RColorBrewer::brewer.pal(11,"RdBu")),name=expression(atop("Person's R","correlations")))+
    scale_x_discrete(expand=c(0,0),position = 'top',
                     labels=c('AOA'='AOA','nirS'='nirS','nirK'='nirK','nirSK_nosZ'='nirSK/nosZ',
                              'MBC'='MBC','MBN'='MBN','MBC_MBN'='MBC/MBN','EON_MBN'='EON/MBN','EOC'='EOC','EON'='EON','AOB_AOA'='AOB/AOA',
                              "NH4_N"=expression(paste(NH[4]^+'','-N')),"NO3_N"=expression(paste(NO[3]^-'','-N')),'CO2_C'=expression(paste(italic(F),''[CO2])),
                              'soil_temp'='temp.','soil.water_dry'=expression(paste('moisture'))))+
    scale_y_discrete(expand=c(0,0),labels=c('NO_N-control'='C',#expression(paste(italic(F),''[NO])),
                                            'NO_N-warmed'='W',#expression(paste(italic(F),''[NO])),
                                            'N2O_N-control'='C',#expression(paste(italic(F),''[N2O])),
                                            'N2O_N-warmed'='W',#expression(paste(italic(F),''[N2O])),
                                            position = 'left')) +
    theme(axis.text.x=element_text(angle =0,hjust =0.5,vjust =0.5,
                                   color="black",face = "bold",size = 10),
          axis.text.y=element_text(color="black",face = "bold",size =10,hjust=1),
          strip.background = element_rect(colour='black',fill=NA,linewidth = .3),
          strip.text = element_text(color="black",face = "bold",size =10,hjust=0.5),
          axis.ticks= element_blank(),
          legend.position = 'right',
          panel.border = element_rect(colour='black',fill=NA,linewidth = .3),
          plot.margin = unit(c(0.6,0.05,0.1,0.65),units="inches"),
          panel.spacing.y = unit(0,"cm"))+
    scale_size(range=c(1,10),guide=NULL)+
    guides(color=guide_colorbar(direction="vertical",reverse=F,barwidth=unit(.5,"cm"),
                                barheight=unit(4,"cm")))+
    coord_cartesian(clip='off')+
    #对y进行修饰
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=-0.35,xmax=0.5,
                      ymin=0.5,ymax=2.5)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=-0.35,xmax=0.5,
                      ymin=2.5,ymax=4.5)+
    #对x进行修饰
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=0.5,xmax=6.5,
                      ymin=4.5,ymax=5.2)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=6.5,xmax=13.5,
                      ymin=4.5,ymax=5.2)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=13.5,xmax=17.5,
                      ymin=4.5,ymax=5.2)+
    annotation_custom(grob = grid::textGrob(label = "Environmental factors & Substrates", hjust=0.5,gp=gpar(col="black",cex=1)),
                      xmin=3,xmax=3,ymin=5.0,ymax=5.0)+
    annotation_custom(grob = grid::textGrob(label = "Functional genes", hjust=0.5,gp=gpar(col="black",cex=1)),
                      xmin=10,xmax=10,ymin=5.0,ymax=5.0)+
    annotation_custom(grob = grid::textGrob(label = "Microbial biomass", hjust=0.5,gp=gpar(col="black",cex=1,angle=90)),
                      xmin=15.5,xmax=15.5,ymin=5.0,ymax=5.0)+
    
    annotation_custom(grob = grid::textGrob(label = expression(paste(italic(F),''[NO])), hjust=0.5, rot = 90, gp=gpar(col="black",cex=1)),
                      xmin=-0.1,xmax=-0.1,ymin=3.5,ymax=3.5)+
    annotation_custom(grob = grid::textGrob(label = expression(paste(italic(F),''[N2O])), hjust=0.5,rot = 90, gp=gpar(col="black",cex=1)),
                      xmin=-0.1,xmax=-0.1,ymin=1.5,ymax=1.5)+
    geom_vline(xintercept = 0.09,col='black',size=0.3)+
    geom_hline(yintercept = 4.9,col='black',size=0.3);p
  
  ggsave(file="heatmap_Ol.pdf",width=11.3,height=4.5,units="in",dpi=300)
  
  
  
  ###df_c1
  dat_c0<-as.matrix(dat[dat$soillayer=="0-10cm" & dat$treatment=='control',
                        c("AOB","AOA","nirK","nirS","nirSK","nosZ","MBC", "MBN","EON_MBN","EOC","TDN","EON",
                          "NH4_N","NO3_N","soil.water_dry","soil_temp","NO_N","N2O_N","N2_N","CO2_C","NH4_NO3","MBC_MBN",'EON_MBN',
                          "NH4NO3","AOB_AOA","nirSK_nosZ")]) #利用as.matrix()将所需数据集转换为matrix格式，才可在corrplot中跑
  names(dat_c0)
  
  dat_c0=data.frame(scale(dat_c0))#数据标准化
  
  # dat_c0<-cor(dat_c0,method="spearman") 
  
  names(dat_c0)
  
  pp<-corr.test(dat_c0[,c('soil_temp','soil.water_dry','EOC','EON','NH4_N','NO3_N',
                          'AOB','AOA','nirK','nirS','nosZ',"AOB_AOA",'nirSK_nosZ',
                          'MBN','MBC_MBN','EON_MBN','CO2_C')],
                dat_c0[,c("NO_N","N2O_N")],method="pearson",adjust = "fdr")
  
  
  
  cor <- pp$r
  pvalue <- pp$p
  
  df_c0 <- melt(cor) %>% mutate(pvalue=melt(pvalue)[,3],
                                p_signif=symnum(pvalue, corr = FALSE, na = FALSE,  
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                symbols = c("***", "**", "*", "+", " "))) %>% 
    set_colnames(c("env","obs","r","p","p_signif"))
  
  str(df_c0)
  
  names(dat_c0)[c(16:18)]
  df_c0<-within(df_c0,obs<-factor(obs,levels=c('gasN','N2_N','N2O_N','NO_N')))
  
  names(dat_c0)[c(1:6,7:15,19:20)]
  # df_c0<-within(df_c0,env<-factor(env,levels=c('soil_temp','soil.water_dry','DOC','NH4_N','NO3_N','NH4NO3',
  #                                              'AOA','nirK','nirS','nirSK','nosZ',"nirSK_nosZ",
  #                                              'MBC','MBN','CO2_C',"min","nit")))
  
  ###df_w0
  dat_w0<-as.matrix(dat[dat$soillayer=='0-10cm' & dat$treatment=='warmed',
                        c("AOB","AOA","nirK","nirS","nirSK","nosZ","MBC", "MBN",'MBC_MBN',"EOC","EON","EON_MBN",
                          "NH4_N","NO3_N","soil.water_dry","soil_temp","NO_N","N2O_N","N2_N","CO2_C","NH4_NO3","MBC_MBN",
                          "NH4NO3","AOB_AOA","nirSK_nosZ")])  #利用as.matrix()将所需数据集转换为matrix格式，才可在corrplot中跑
  
  dat_w0=data.frame(scale(dat_w0))#数据标准化
  
  # dat_c0<-cor(dat_c0,method="spearman") 
  
  library(psych);library(magrittr);library(grid)
  names(dat_w0)
  
  pp<-corr.test(dat_w0[,c('soil_temp','soil.water_dry','EOC','EON','NH4_N','NO3_N',
                          'AOB','AOA','nirK','nirS','nosZ',"AOB_AOA",'nirSK_nosZ',
                          'MBN','MBC_MBN','EON_MBN','CO2_C')],
                dat_w0[,c("NO_N","N2O_N")],method="pearson",adjust = "fdr")  
  cor <- pp$r
  pvalue <- pp$p
  
  df_w0 <- melt(cor) %>% mutate(pvalue=melt(pvalue)[,3],
                                p_signif=symnum(pvalue, corr = FALSE, na = FALSE,  
                                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                symbols = c("***", "**", "*", "+", " "))) %>% 
    set_colnames(c("env","obs","r","p","p_signif"))
  
  str(df_w0)
  
  names(dat_w0)[c(16:18)]
  df_w0<-within(df_w0,obs<-factor(obs,levels=c('N2_N','N2O_N','NO_N')))
  
  names(dat_c0)[c(1:15,19:20)]
  # df_w0<-within(df_w0,env<-factor(env,levels=c('soil_temp','soil.water_dry','DOC','NH4_N','NO3_N','NH4NO3',
  #                                              'AOA','nirK','nirS','nirSK','nosZ',"nirSK_nosZ",
  #                                              'MBC','MBN','CO2_C',"min","nit")))
  
  
  df_w0$treatment='warmed'
  df_c0$treatment='control'
  
  df1<-rbind(df_c0,df_w0)
  
  df1$obs_tr<-paste(df1$obs,df1$treatment,sep='-')
  
  unique(df1$obs_tr)
  df1<-within(df1,obs_tr<-factor(obs_tr,levels=c("N2O_N-warmed","N2O_N-control",
                                                 "NO_N-warmed","NO_N-control")))
  
  
  
  # df1<-df1[df1$obs_tr %in% c("N2_N-warmed","N2_N-control","N2O_N-warmed","N2O_N-control","NO_N-warmed","NO_N-control"),]
  
  p1<-ggplot(df1,aes(env,obs_tr,col=r,fill=r,group=treatment))+
    geom_tile(color="grey80",fill="white",size=0.3)+
    geom_point(aes(size =abs(r)),shape=21)+facet_grid('0-10 cm'~.)+
    geom_text(data=df1[df1$p_signif!='+',],aes(label=p_signif),size=6,color="white",hjust=0.5,vjust=0.7)+
    # geom_text(data=df[df$p_signif=='+',],aes(label=p_signif),size=3,color="white",hjust=0.5,vjust=0.5)+
    labs(x = NULL,y = NULL,color=NULL,fill=NULL)+
    scale_color_gradientn(limits=c(-1,1),colours = rev(RColorBrewer::brewer.pal(11,"RdBu")),name=expression(atop("Person's R","correlations")))+
    scale_fill_gradientn(limits=c(-1,1),colours = rev(RColorBrewer::brewer.pal(11,"RdBu")),name=expression(atop("Person's R","correlations")))+
    scale_x_discrete(expand=c(0,0),position = 'top',
                     labels=c('AOA'='AOA','nirS'='nirS','nirK'='nirK','nirSK'='nirS+nirK','nirSK_nosZ'='nirSK/nosZ',
                              'MBC'='MBC','MBN'='MBN','MBC_MBN'='MBC/MBN','EOC'='EOC','AOB_AOA'='AOB/AOA','EON_MBN'='EON/MBN',
                              "NH4_N"=expression(paste(NH[4]^+'','-N')),'CO2_C'=expression(paste(italic(F),''[CO2])),
                              "NO3_N"=expression(paste(NO[3]^-'','-N')),"NH4NO3"=expression('TIN'),
                              'soil_temp'='Temp','soil.water_dry'=expression(paste('Moisture'))))+
    scale_y_discrete(expand=c(0,0),labels=c('NO_N-control'='C',#expression(paste(italic(F),''[NO])),
                                            'NO_N-warmed'='W',#expression(paste(italic(F),''[NO])),
                                            'N2O_N-control'='C',#expression(paste(italic(F),''[N2O])),
                                            'N2O_N-warmed'='W',#expression(paste(italic(F),''[N2O])),
                                            'N2_N-control'='C',#expression(paste(italic(F),''[NO+N2O])),
                                            'N2_N-warmed'='W',#expression(paste(italic(F),''[NO+N2O]))),
                                            position = 'left')) +
    theme(axis.text.x=element_text(angle =0,hjust =0.5,vjust =0.5,
                                   color="black",face = "bold",size = 10),
          axis.text.y=element_text(color="black",face = "bold",size =10,hjust=1),
          strip.background = element_rect(colour='black',fill=NA,linewidth = .3),
          strip.text = element_text(color="black",face = "bold",size =10,hjust=0.5),
          axis.ticks= element_blank(),
          legend.position = 'right',
          panel.border = element_rect(colour='black',fill=NA,linewidth = .3),
          plot.margin = unit(c(0.6,0.05,0.1,0.65),units="inches"),
          panel.spacing.y = unit(0,"cm"))+
    scale_size(range=c(1,10),guide=NULL)+
    guides(color=guide_colorbar(direction="vertical",reverse=F,barwidth=unit(.5,"cm"),
                                barheight=unit(4,"cm")))+
    coord_cartesian(clip='off')+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=-0.35,xmax=0.5,
                      ymin=0.5,ymax=2.5)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=-0.35,xmax=0.5,
                      ymin=2.5,ymax=4.5)+
    #对x进行修饰
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=0.5,xmax=6.5,
                      ymin=4.5,ymax=5.2)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=6.5,xmax=13.5,
                      ymin=4.5,ymax=5.2)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=13.5,xmax=17.5,
                      ymin=4.5,ymax=5.2)+
    annotation_custom(grob = grid::textGrob(label = "Environmental factors & Substrates", hjust=0.5,gp=gpar(col="black",cex=1)),
                      xmin=3,xmax=3,ymin=5.0,ymax=5.0)+
    annotation_custom(grob = grid::textGrob(label = "Functional genes", hjust=0.5,gp=gpar(col="black",cex=1)),
                      xmin=10,xmax=10,ymin=5.0,ymax=5.0)+
    annotation_custom(grob = grid::textGrob(label = "Microbial biomass", hjust=0.5,gp=gpar(col="black",cex=1,angle=90)),
                      xmin=15.5,xmax=15.5,ymin=5.0,ymax=5.0)+
    
    annotation_custom(grob = grid::textGrob(label = expression(paste(italic(F),''[NO])), hjust=0.5, rot = 90, gp=gpar(col="black",cex=1)),
                      xmin=-0.1,xmax=-0.1,ymin=3.5,ymax=3.5)+
    annotation_custom(grob = grid::textGrob(label = expression(paste(italic(F),''[N2O])), hjust=0.5,rot = 90, gp=gpar(col="black",cex=1)),
                      xmin=-0.1,xmax=-0.1,ymin=1.5,ymax=1.5)+
    geom_vline(xintercept = 0.09,col='black',size=0.3)+
    geom_hline(yintercept = 4.9,col='black',size=0.3);p1
  
  ggsave(file="heatmap_1l.pdf",width=11.3,height=4.5,units="in",dpi=300)
  
  
  df$layer='O layer'
  df1$layer='0-10 cm'
  
  df<-rbind(df,df1)
  
  annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
  {
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
  }
  
  unique(df$obs_tr)
  
  df<-within(df,layer<-factor(layer,levels = c('O layer','0-10 cm')))
  p1<-ggplot(df,aes(env,obs_tr,col=r,fill=r,group=treatment))+
    geom_tile(color="grey80",fill="white",size=0.3)+
    geom_point(aes(size =abs(r)),shape=21)+facet_grid(layer~.)+
    geom_text(data=df[df$p_signif!='+',],aes(label=p_signif),size=6,color="black",hjust=0.5,vjust=0.7)+
    labs(x = NULL,y = NULL,color=NULL,fill=NULL)+
    scale_color_gradientn(limits=c(-1,1),colours = rev(RColorBrewer::brewer.pal(11,"RdBu")),name=expression("Person's R correlations"))+
    scale_fill_gradientn(limits=c(-1,1),colours = rev(RColorBrewer::brewer.pal(11,"RdBu")),name=expression("Person's R correlations"))+
    scale_x_discrete(expand=c(0,0),position = 'top',
                     labels=c('AOA'=expression(paste('AOA')),'nirS'=expression(paste('nirS')),'nirK'=expression(paste('nirK')),
                              'nirSK'=expression(paste('nirSK')),'nirSK_nosZ'=expression('nirSK/nosZ'),'AOB_AOA'='AOB/AOA',
                              'MBC'=expression(paste('MBC')),'MBN'=expression(paste('MBN')),'EOC'=expression(paste('EOC')),'MBC_MBN'='MBC/MBN','EON_MBN'='EON/MBN','TIN_MBN'='TIN/MBN',
                              "NH4_N"=expression(paste(NH[4]^+'')),'CO2_C'=expression(paste('CO'[2],' flux')),
                              "NO3_N"=expression(paste(NO[3]^-'')),"NH4NO3"=expression('TIN'),'min'=expression('N-min'),'nit'=expression('N-nit'),
                              'soil_temp'=expression(paste('temp.')),'soil.water_dry'=expression(paste('moisture'))))+
    scale_y_discrete(expand=c(0,0),labels=c('NO_N-control'='C',#expression(paste(italic(F),''[NO])),
                                            'NO_N-warmed'='W',#expression(paste(italic(F),''[NO])),
                                            'N2O_N-control'='C',#expression(paste(italic(F),''[N2O])),
                                            'N2O_N-warmed'='W',#expression(paste(italic(F),''[N2O])),
                                            'N2_N-control'='C',#expression(paste(italic(F),''[NO+N2O])),
                                            'N2_N-warmed'='W',#expression(paste(italic(F),''[NO+N2O]))),
                                            position = 'left')) +
    theme(axis.text.x=element_text(angle =0,hjust =0.5,vjust =-2,
                                   color="black",face = "bold",size = 8),
          axis.text.y=element_text(color="black",face = "bold",size =8,hjust=1),
          strip.background = element_rect(colour='black',fill=NA,linewidth = .59),
          strip.text = element_text(color="black",face = "bold",size =10,hjust=0.5,
                                    margin=margin(c(1, 5, 1, 5), 'cm')),
          axis.ticks= element_blank(),legend.title = element_text(face = "bold",hjust =0.5,vjust =0.9,color="black",size = 10),
          legend.position = 'bottom',legend.text = element_text(hjust =0.5,vjust =0,color="black",size = 10),
          panel.border = element_rect(colour='black',fill=NA,linewidth = .3),
          plot.margin = unit(c(0.35,0.05,0.1,0.4),units="inches"),
          panel.spacing.y = unit(0.15,"cm"))+
    scale_size(range=c(1,10),guide=NULL)+
    guides(color=guide_colorbar(direction="horizontal",reverse=F,barwidth=unit(12,"cm"),
                                barheight=unit(.4,"cm")))+
    coord_cartesian(clip='off')+
    #对y进行修饰
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=-0.35,xmax=0.5,
                      ymin=0.5,ymax=2.5)+
    annotation_custom(grob=rectGrob(gp=gpar(col="black",fill="transparent")),xmin=-0.35,xmax=0.5,
                      ymin=2.5,ymax=4.5)+
    #对x进行修饰
    annotation_custom2(data = df %>% filter(layer == 'O layer'),grob=rectGrob(gp=gpar(col="black",fill="blue",alpha=0.15)),
                       xmin=0.5,xmax=6.5,ymin=4.5,ymax=5.4)+
    annotation_custom2(data = df %>% filter(layer == 'O layer'),grob=rectGrob(gp=gpar(col="black",fill="red",alpha=0.15)),
                       xmin=6.5,xmax=13.5,ymin=4.5,ymax=5.4)+
    annotation_custom2(data = df %>% filter(layer == 'O layer'),grob=rectGrob(gp=gpar(col="black",fill="grey",alpha=0.15)),
                       xmin=13.5,xmax=17.5,ymin=4.5,ymax=5.4)+
    annotation_custom2(data = df %>% filter(layer == 'O layer'),grob=rectGrob(gp=gpar(col="black",fill="transparent",cex=0.05)),
                       xmin=0.5,xmax=17.5,ymin=4.5,ymax=5.4)+
    
    
    
    annotation_custom2(data = df %>% filter(layer == 'O layer'),xmin=3.5,xmax=3.5,ymin=5.15,ymax=5.15,
                       grob = grid::textGrob(label = "Environmental factors & Substrates", hjust=0.5,gp=gpar(col="black",cex=0.8,face = "bold")))+
    annotation_custom2(data = df %>% filter(layer == 'O layer'),xmin=10,xmax=10,ymin=5.15,ymax=5.15,
                       grob = grid::textGrob(label = "Functional genes", hjust=0.5,gp=gpar(col="black",cex=0.8,face = "bold")))+
    annotation_custom2(data = df %>% filter(layer == 'O layer'),xmin=15.5,xmax=15.5,ymin=5.15,ymax=5.15,
                       grob = grid::textGrob(label = "Microbial biomass", hjust=0.5,gp=gpar(col="black",cex=0.8,face = "bold")))+
    
    annotation_custom(grob = grid::textGrob(label = expression(paste('NO flux')), hjust=0.5, rot = 90, gp=gpar(col="black",cex=0.8)),xmin=-0.1,xmax=-0.1,ymin=3.5,ymax=3.5)+
    annotation_custom(grob = grid::textGrob(label = expression(paste('N'[2],'O flux')), hjust=0.5,rot = 90, gp=gpar(col="black",cex=0.8)),xmin=-0.1,xmax=-0.1,ymin=1.5,ymax=1.5)+
    # geom_hline(data = df[df$layer == 'O layer',],aes(yintercept = 4.9),col='black',size=0.3)+
    geom_vline(xintercept = 0.09,col='black',size=0.3);p1
  
  ggsave(file="heatmap_all_t48hr.pdf",p1,width=10.5,height=7,units="in",dpi=300)
  

  (Fig_6ab)/p1+
    plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
    theme(plot.tag.position = c(0, 0.98),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
          plot.tag = element_text(size = 12,vjust = 0,hjust=0,face="bold"))->Fig_6;Fig_6   
  ####先运行出R1517N+R1517_nitri
  

  ggsave(file="Fig_6abd_heatmap_all_t48hr_Rmin+Rnit.pdf",Fig_6,width=11,height=11,units="in",dpi=300)
}


########## Fig_6.c_<Mean effect across the sample dates> #############
dat=dat0[dat0$labs %in%c('Day1','Day2'),]

dat$doy<-as.numeric(format(as.Date(dat$Date),"%j"))

dat$pair<-ifelse(substr(dat$plot_ab,1,1) %in% c(1,2),'p1',
                 ifelse(substr(dat$plot_ab,1,1) %in% c(3,4),'p2','p3'))   
                 
#计算treat的比值warmed/control
dat_t<-data.frame();k=1
for(i in unique(dat$Date)){
  for(j in unique(dat$soillayer)){
    # for(t in unique(dat$pair)){
  dat_i=dat[dat$Date==i & dat$soillayer==j, ]
  str(dat_i)
  dat_i$TIN_MBN=dat_i$NH4NO3/dat_i$MBN
  dat_i$EON_MBN=dat_i$EON/dat_i$MBN
  
  dat_t[k,1]=i
  dat_t[k,2]=j
  
  dat_t[k,3]=mean(dat_i$soil.water_dry[dat_i$treatment=='warmed'])/
    mean(dat_i$soil.water_dry[dat_i$treatment=='control'])-1

  dat_t[k,4]=mean(dat_i$NH4_N[dat_i$treatment=='warmed'])/
    mean(dat_i$NH4_N[dat_i$treatment=='control'])-1
  dat_t[k,5]=mean(dat_i$NO3_N[dat_i$treatment=='warmed'])/
    mean(dat_i$NO3_N[dat_i$treatment=='control'])-1
  dat_t[k,6]=mean(dat_i$NH4NO3[dat_i$treatment=='warmed'])/
    mean(dat_i$NH4NO3[dat_i$treatment=='control'])-1

  dat_t[k,7]=mean(dat_i$TDN[dat_i$treatment=='warmed'])/
    mean(dat_i$TDN[dat_i$treatment=='control'])-1
  dat_t[k,8]=mean(dat_i$EOC[dat_i$treatment=='warmed'])/
    mean(dat_i$EOC[dat_i$treatment=='control'])-1
  dat_t[k,9]=mean(dat_i$EON[dat_i$treatment=='warmed'])/
    mean(dat_i$EON[dat_i$treatment=='control'])-1
  
  dat_t[k,10]=mean(dat_i$TIN_MBN[dat_i$treatment=='warmed'])/
    mean(dat_i$TIN_MBN[dat_i$treatment=='control'])-1
  
  dat_t[k,11]=mean(dat_i$EON_MBN[dat_i$treatment=='warmed'])/
    mean(dat_i$EON_MBN[dat_i$treatment=='control'])-1
  
  dat_t[k,12]=mean(dat_i$MBC[dat_i$treatment=='warmed'])/
    mean(dat_i$MBC[dat_i$treatment=='control'])-1
  dat_t[k,13]=mean(dat_i$MBN[dat_i$treatment=='warmed'])/
    mean(dat_i$MBN[dat_i$treatment=='control'])-1
  dat_t[k,14]=mean(dat_i$MBC_MBN[dat_i$treatment=='warmed'])/
    mean(dat_i$MBC_MBN[dat_i$treatment=='control'])-1
  
  dat_t[k,15]=mean(dat_i$AOA[dat_i$treatment=='warmed'])/
    mean(dat_i$AOA[dat_i$treatment=='control'])-1
  dat_t[k,16]=mean(dat_i$AOB[dat_i$treatment=='warmed'])/
    mean(dat_i$AOB[dat_i$treatment=='control'])-1
  dat_t[k,17]=mean(dat_i$nirK[dat_i$treatment=='warmed'])/
    mean(dat_i$nirK[dat_i$treatment=='control'])-1
  dat_t[k,18]=mean(dat_i$nirS[dat_i$treatment=='warmed'])/
    mean(dat_i$nirS[dat_i$treatment=='control'])-1
  dat_t[k,19]=mean(dat_i$nirSK[dat_i$treatment=='warmed'])/
    mean(dat_i$nirSK[dat_i$treatment=='control'])-1
  dat_t[k,20]=mean(dat_i$nosZ[dat_i$treatment=='warmed'])/
    mean(dat_i$nosZ[dat_i$treatment=='control'])-1
  dat_t[k,21]=mean(dat_i$nirSK_nosZ[dat_i$treatment=='warmed'])/
    mean(dat_i$nirSK_nosZ[dat_i$treatment=='control'])-1
  
  dat_t[k,22]=mean(dat_i$NO_N[dat_i$treatment=='warmed'])/
    mean(dat_i$NO_N[dat_i$treatment=='control'])-1
  dat_t[k,23]=mean(dat_i$N2O_N[dat_i$treatment=='warmed'])/
    mean(dat_i$N2O_N[dat_i$treatment=='control'])-1
  k=k+1
  # dat_t[k,24]=t}
  }
}

names(dat_t)[1:23]<-c('Date','layer','moisture',
                      'NH4','NO3','NH4NO3','TDN','EOC','EON',
                      'TIN_MBN','EON_MBN','MBC','MBN','MBC_MBN',
                      'AOA','AOB','nirK','nirS','nirSK','nosZ','nirSK_nosZ',
                      'NO','N2O')

unique(dat$Date)

names(dat_t)
library(reshape2)
dat_M<-reshape2::melt(dat_t,id.vars=c("Date","layer"),variable.name="variable",value.name="effect_size",
            measure.vars=c("moisture","NH4","NO3","NH4NO3","TDN","EOC","EON","TIN_MBN","EON_MBN","MBC","MBN",
                           "MBC_MBN","AOA","AOB","nirK","nirS","nirSK","nosZ","nirSK_nosZ","NO","N2O"))

dat_M<-within(dat_M,layer<-factor(layer,levels=c('O','0-10cm')))

dat_M<-within(dat_M,variable<-
                factor(variable,levels=rev(c("moisture","NH4","NO3","NH4NO3","TDN","EOC","EON","TIN_MBN","EON_MBN","MBC","MBN",
                                "MBC_MBN","AOA","AOB","nirK","nirS","nirSK","nosZ","nirSK_nosZ","NO","N2O"))))

mean(dat_M$effect_size[dat_M$layer=='O' & dat_M$variable=='NH4' & dat_M$Date!='2021-05-22'])  #NH4的效应值
sd(dat_M$effect_size[dat_M$layer=='O' & dat_M$variable=='NH4' & dat_M$Date!='2021-05-22'])/sqrt(12)

mean(dat_M$effect_size[dat_M$layer=='O' & dat_M$variable=='EON' & dat_M$Date!='2021-05-22'])  #NH4的效应值
sd(dat_M$effect_size[dat_M$layer=='O' & dat_M$variable=='EON' & dat_M$Date!='2021-05-22'])/sqrt(12)

mean(dat_M$effect_size[dat_M$layer=='0-10cm' & dat_M$variable=='EON' & dat_M$Date!='2021-05-22'])  #NH4的效应值
sd(dat_M$effect_size[dat_M$layer=='0-10cm' & dat_M$variable=='EON' & dat_M$Date!='2021-05-22'])/sqrt(12)

mean(dat_M$effect_size[dat_M$layer=='0-10cm' & dat_M$variable=='MBN' & dat_M$Date!='2021-05-22'])  #NH4的效应值
sd(dat_M$effect_size[dat_M$layer=='0-10cm' & dat_M$variable=='MBN' & dat_M$Date!='2021-05-22'])/sqrt(12)

unique(dat_M$variable)

'%notin%' <- Negate('%in%')

library(gghalves);library(ggprism);library(grid)
dat_M<-na.omit(dat_M)
unique(dat_M$variable)

ggplot(dat_M[dat_M$variable %notin% c('NO','N2O','TDN'),],
       aes(variable,effect_size))+
  facet_grid(layer~.,labeller = labeller(layer = c('O' = "O layer", '0-10cm' = "0-10 cm")))+
  # geom_half_violin(fill='grey',col='transparent',alpha=0.7, position = position_nudge(x = 0.2, y = 0),
  #                  side = 'R', adjust = 0.8, trim = F,width=2)+
  geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
  geom_vline(xintercept = 6.5,lty=1,linewidth=0.2)+geom_vline(xintercept = 13.5,lty=1,linewidth=0.2)+
  geom_point(position = position_jitter(0.08),pch=1,size=0.6,alpha=0.79,stroke=0.3)+
  geom_boxplot(outlier.shape = NA,width = 0.32,alpha = 1,fill='white',linewidth = 0.35)+
  scale_y_continuous(limits=c(-Inf,1.5),breaks = seq(-1,2,1),expand=c(0,0))+
  scale_x_discrete(position='top',
                   limits=c("moisture","NH4","NO3","NH4NO3","EOC","EON",
                            "AOB","AOA","nirK","nirS","nirSK","nosZ","nirSK_nosZ","MBC","MBN","MBC_MBN","EON_MBN"),
                   labels=c('moisture',expression(paste('NH'[4]^'+')),expression('NO'[3]^'-'),
                            expression(paste('NH'[4]^'+'~'+NO'[3]^'-')),'EOC','EON',
                            'AOB','AOA','nirK','nirS','nirSK','nosZ','nirSK/nosZ','MBC','MBN','MBC/MBN',"EON/MBN"),
                   expand=c(0,0))+
  coord_cartesian(ylim=c(-1,1.5),xlim = c(0.5,17.5),clip = 'off')+
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(colour = "black",size=8,angle = 0,hjust = .5,vjust =0.5),
        axis.ticks = element_line(linewidth = .15),
        axis.ticks.length = unit(1,"mm"),
        prism.ticks.length = unit(0.5,"mm"),
        axis.text.y = element_text(colour = "black",size = 12), 
        axis.title.x =element_text(size=10), 
        axis.title.y=element_text(colour = "black",size=12,vjust =0.5,face='bold'),
        legend.text = element_text(size=10), legend.title =element_text(size=7),
        panel.background = element_rect(fill =NA, linewidth=0.05,colour = "black", linetype = "solid",inherit.blank=T),
        panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
        panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
        plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        plot.margin = unit(c(0.6,0.2,0.1,0.2),'cm'),   #调整画图区域的间距，从上右下左调整
        strip.background = element_rect(fill = NA,color='black'), 
        strip.text = element_text(colour = "black",size = 10,hjust=.5,face='bold'),
        legend.position = "none")+
  labs(x=NULL,y=expression(paste('Effect size')))+
  annotation_custom2(data = dat_M %>% filter(layer == 'O'),grob=rectGrob(gp=gpar(col="black",fill="blue",alpha=0.15)),
                     xmin=0.5,xmax=6.5,ymin=1.5,ymax=2.4)+
  annotation_custom2(data = dat_M %>% filter(layer == 'O'),grob=rectGrob(gp=gpar(col="black",fill="red",alpha=0.15)),
                     xmin=6.5,xmax=13.5,ymin=1.5,ymax=2.4)+
  annotation_custom2(data = dat_M %>% filter(layer == 'O'),grob=rectGrob(gp=gpar(col="black",fill="grey",alpha=0.15)),
                     xmin=13.5,xmax=17.5,ymin=1.5,ymax=2.4)+
  annotation_custom2(data = dat_M %>% filter(layer == 'O'),grob=rectGrob(gp=gpar(col="black",fill="transparent",cex=0.05)),
                     xmin=0.5,xmax=17.5,ymin=1.5,ymax=2.4)+
  geom_text(data = dat_M %>% filter(layer == 'O' & Date == '2019-07-03' & variable=='AOA'),size = 3.5,
            aes(x=3.5,y=1, label = "Environmental factors & Substrates"), vjust = -5.6, hjust = 0.5, angle = 0,fontface='plain')+
  geom_text(data = dat_M %>% filter(layer == 'O'& Date == '2019-07-03' & variable=='AOA'),size = 3.5,
            aes(x=10.5,y=1, label = "Functional genes"), vjust = -5.6, hjust = 0.5, angle = 0,fontface='plain')+
  geom_text(data = dat_M %>% filter(layer == 'O' & Date == '2019-07-03' & variable=='AOA'),size = 3.5,
            aes(x=15.5,y=1, label = "Microbial biomass"), vjust = -5.6, hjust = 0.5, angle = 0,fontface='plain')->p2;p2
  
ggsave("Response ratio_2019-2023.pdf", p2, width =12, height =6,
       device=cairo_pdf)

########## Fig_6_<Combined the panels> #############
(Fig_6ab/p2/(p1+theme(plot.margin = unit(c(1,0.2,0,0.2),'cm'))))+
  plot_layout(heights = c(0.95, 0.95, 1.35))+
  plot_annotation(tag_levels = 'a', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
  theme(plot.tag.position = c(0, 0.98),   #c(0,1)左上角，c(0,0)左下角，其他位置依次类推
        plot.tag = element_text(size = 12,vjust = 0,hjust=0,face="bold"))->Fig_6;Fig_6   

ggsave(file="Fig_6abcd_heatmap_all_t48hr_Rmin+Rnit+effect_size.pdf",Fig_6,width=12.5,height=12,units="in",dpi=300)

##########@@@##########
str(dat);names(dat)
dat$Date<-as.factor(dat$Date)

summary(aov(soil.water_dry~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O',]))

summary(aov(NH4_N~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O',]))

summary(aov(NO3_N~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O',]))

summary(aov(EOC~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(EON~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(AOA~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(AOB~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(nirK~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(nirS~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(nosZ~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(MBC~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(MBN~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='O' & is.na(dat$AOB)!='TRUE',]))

summary(aov(EON~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='0-10cm' & is.na(dat$AOB)!='TRUE',]))

summary(aov(MBC~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='0-10cm' & is.na(dat$AOB)!='TRUE',]))

summary(aov(MBN~treatment*Date+Error(1/plot),data=dat[dat$soillayer=='0-10cm' & is.na(dat$AOB)!='TRUE',]))


###### Ext.Data.Fig.7_SEM对Min+Nit R15 R17的影响 ######
######## O layer
str(dat_N)   #提取R15分析the main factors
dat_N$plot_ab<-dat_N$subplot
dat_N$Date<-as.Date(dat_N$date)

unique(dat$Date)

str(dat)
dat_i=dat[dat$Date %in% c("2020-11-05","2021-09-19","2022-10-01","2023-09-23") & dat$soillayer=='O' & dat$labs=='Day1',
            c('Date','treatment','plot','plot_ab','EOC','TDN','EON','MBN',
              'AOA','AOB','nirS','nirK','nosZ','NH4_N','NO3_N',
              'soil.water_dry','soil_temp')]
names(dat_N)
dat_i$Date<-as.Date(dat_i$Date)
dat_i$year<-format(dat_i$Date,"%Y")


dat_i<-merge(dat_i,dat_N[dat_N$layer=='Oa+e',c('year','plot_ab','R15_min','R15_nit')],by=c('year','plot_ab'),all=T)

dat_i$treatment<-ifelse(dat_i$treatment=='control',0,1)
dat_i$plot<-as.numeric(dat_i$plot)
str(dat_i)

dat_i[,c('treatment','plot','EOC','TDN','EON','MBN','AOA','AOB','nirS','nirK','nosZ','NH4_N','NO3_N',
         'R15_min','R15_nit','soil.water_dry',"soil_temp")]<-
  scale(dat_i[,c('treatment','plot','EOC','TDN','EON','MBN','AOA','AOB','nirS','nirK','nosZ','NH4_N','NO3_N',
                 'R15_min','R15_nit','soil.water_dry',"soil_temp")])

ggplot(dat_i) + 
  geom_point(data=dat_i,aes(MBN,R15_min,col=year))+
  geom_smooth(data=dat_i,aes(MBN,R15_min,col=year),method='lm',se=F)

ggplot(dat_i,aes(EON,R15_min)) + 
  geom_point(data=dat_i,aes(EON,R15_min,col=year))+
  geom_smooth(data=dat_i,aes(EON,R15_min,col=year),method='lm',se=F)+
  geom_smooth(data=dat_i,aes(EON,R15_min),method='lm',se=F,col='black')+
  stat_poly_eq(formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
               eq.with.lhs = "italic(y)~`=`~",   #????y"????ʽ         
               aes(label = paste(stat(adj.rr.label),stat(p.value.label),sep ="*\",\"~~")),
               parse = TRUE,label.x.npc = "left", label.y.npc = 0.87,size = 5)




library(piecewiseSEM);library(lme4);library(nlme)

psem_cat<-psem(
  # lme(CO2_C~soil_temp+WFPS+treatment,random=~1|plot,data=dat,method="ML"),
  lme(R15_min~MBN+EON,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(R15_nit~soil.water_dry+MBN+R15_min,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(EON~soil.water_dry+treatment,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(soil.water_dry~treatment,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(MBN~treatment+soil.water_dry+EON,random=~1|plot,data=na.omit(dat_i),method="ML")
  
  # min_growing%~~%nit_growing  ###表示相关，即不会考虑这条路径，会缺少这条路径，R2和AIC会减少
)

summary(psem_cat)  ###O layer

res<-summary(psem_cat);res 
#fisherC
fisherC(psem_cat)

AIC(psem_cat)



########0-10cm
str(dat_N)   #提取R15分析the main factors
dat_N$plot_ab<-dat_N$subplot
dat_N$Date<-as.Date(dat_N$date)

unique(dat$Date)

str(dat)
dat_i=dat[dat$Date %in% c("2020-11-05","2021-09-19","2022-10-01","2023-09-23") & dat$soillayer=='0-10cm' & dat$labs=='Day1',
          c('Date','treatment','plot','plot_ab','EOC','TDN','EON','MBN',
            'AOA','AOB','nirS','nirK','nosZ','NH4_N','NO3_N',
            'soil.water_dry','soil_temp')]
names(dat_N)
dat_i$Date<-as.Date(dat_i$Date)
dat_i$year<-format(dat_i$Date,"%Y")


dat_i<-merge(dat_i,dat_N[dat_N$layer=='0-10 cm',c('year','plot_ab','R15_min','R15_nit')],by=c('year','plot_ab'),all=T)

dat_i$treatment<-ifelse(dat_i$treatment=='control',0,1)
dat_i$plot<-as.numeric(dat_i$plot)
str(dat_i)

dat_i[,c('treatment','plot','EOC','TDN','EON','MBN','AOA','AOB','nirS','nirK','nosZ','NH4_N','NO3_N',
         'R15_min','R15_nit','soil.water_dry',"soil_temp")]<-
  scale(dat_i[,c('treatment','plot','EOC','TDN','EON','MBN','AOA','AOB','nirS','nirK','nosZ','NH4_N','NO3_N',
                 'R15_min','R15_nit','soil.water_dry',"soil_temp")])

ggplot(dat_i) + 
  geom_point(data=dat_i,aes(MBN,R15_min,col=year))+
  geom_smooth(data=dat_i,aes(MBN,R15_min,col=year),method='lm',se=F)

ggplot(dat_i,aes(EON,R15_min)) + 
  geom_point(data=dat_i,aes(EON,R15_min,col=year))+
  geom_smooth(data=dat_i,aes(EON,R15_min,col=year),method='lm',se=F)+
  geom_smooth(data=dat_i,aes(EON,R15_min),method='lm',se=F,col='black')+
  stat_poly_eq(formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
               eq.with.lhs = "italic(y)~`=`~",   #????y"????ʽ         
               aes(label = paste(stat(adj.rr.label),stat(p.value.label),sep ="*\",\"~~")),
               parse = TRUE,label.x.npc = "left", label.y.npc = 0.87,size = 5)




library(piecewiseSEM);library(lme4);library(nlme)

psem_cat<-psem(
  # lme(CO2_C~soil_temp+WFPS+treatment,random=~1|plot,data=dat,method="ML"),
  lme(R15_min~MBN+EON,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(R15_nit~soil.water_dry+MBN+R15_min,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(EON~soil.water_dry+treatment,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(soil.water_dry~treatment,random=~1|plot,data=na.omit(dat_i),method="ML"),
  lme(MBN~treatment+soil.water_dry+EON,random=~1|plot,data=na.omit(dat_i),method="ML")
  
  # min_growing%~~%nit_growing  ###表示相关，即不会考虑这条路径，会缺少这条路径，R2和AIC会减少
)

summary(psem_cat)  ###O layer

res<-summary(psem_cat);res 
#fisherC
fisherC(psem_cat)

AIC(psem_cat)

#--------------------obtain the coefficients and R2
coefs(psem_cat)
#coefs(keeley.sem2,intercepts=T)#截距也给出
rsquared(psem_cat)




##########计算响应比之间的关系20240203############
#####将6个对照求均值，分别求出3个样方增温的响应比，10*3=30

dat$gasN<-dat$NO_N+dat$N2O_N+dat$N2_N

unique(dat$labs)

dat_res<-dat[,c('Date','treatment','plot','soillayer',"AOB",'AOA','nirK','nirS','nirSK','nosZ','nirSK_nosZ','EOC','EON',
                'NH4_N','NO3_N','NH4NO3','soil_temp','soil.water_dry','MBC','MBN','CO2_C','NO_N','N2O_N','N2_N','gasN')]

dat_res$pair<-ifelse(dat_res$plot %in% c('1','2'),'pair1',
                 ifelse(dat_res$plot %in% c('3','4'),'pair2','pair3'))
str(dat_res)

myres<-data.frame()
str(dat_res)
k=1
for(i in unique(dat_res$Date)){
  for(j in unique(dat_res$soillayer)){
    for(p in unique(dat_res$pair)){
      
      dat_i<-dat_res[dat_res$Date==i & dat_res$soillayer==j & dat_res$pair==p,]
      dat_ij<-dat_res[dat_res$Date==i & dat_res$soillayer==j,]
      
      myres[k,1]=i
      myres[k,2]=j
      myres[k,3]=p
      myres[k,4]=ifelse(as.numeric(format(as.Date(i),"%m"))<=10,'growing','nongrowing')
      
      ###对照的数据
      myres[k,5]=mean(dat_i$AOB[dat_i$treatment=='control'])
      myres[k,6]=mean(dat_i$AOA[dat_i$treatment=='control'])
      myres[k,7]=mean(dat_i$nirK[dat_i$treatment=='control'])
      myres[k,8]=mean(dat_i$nirS[dat_i$treatment=='control'])
      myres[k,9]=mean(dat_i$nirSK[dat_i$treatment=='control'])
      myres[k,10]=mean(dat_i$nosZ[dat_i$treatment=='control'])
      myres[k,11]=mean(dat_i$nirSK_nosZ[dat_i$treatment=='control'])
      myres[k,12]=mean(dat_i$EOC[dat_i$treatment=='control'])
      myres[k,13]=mean(dat_i$NH4_N[dat_i$treatment=='control'])
      myres[k,14]=mean(dat_i$NO3_N[dat_i$treatment=='control'])
      myres[k,15]=mean(dat_i$NH4NO3[dat_i$treatment=='control'])
      myres[k,16]=mean(dat_i$soil_temp[dat_i$treatment=='control'])
      myres[k,17]=mean(dat_i$soil.water_dry[dat_i$treatment=='control'])
      myres[k,18]=mean(dat_i$MBC[dat_i$treatment=='control'])
      myres[k,19]=mean(dat_i$MBN[dat_i$treatment=='control'])
      myres[k,20]=mean(dat_i$NO_N[dat_i$treatment=='control'])
      myres[k,21]=mean(dat_i$N2O_N[dat_i$treatment=='control'])
      myres[k,22]=mean(dat_i$gasN[dat_i$treatment=='control'])
      myres[k,23]=mean(dat_i$CO2_C[dat_i$treatment=='control'])
      
      ###比值数据W/C,对照用3个样方均值
      myres[k,24]=mean(dat_i$AOB[dat_i$treatment=='warmed'])/mean(dat_i$AOB[dat_i$treatment=='control'])
      myres[k,25]=mean(dat_i$AOA[dat_i$treatment=='warmed'])/mean(dat_i$AOA[dat_i$treatment=='control'])
      myres[k,26]=mean(dat_i$nirK[dat_i$treatment=='warmed'])/mean(dat_i$nirK[dat_i$treatment=='control'])
      myres[k,27]=mean(dat_i$nirS[dat_i$treatment=='warmed'])/mean(dat_i$nirS[dat_i$treatment=='control'])
      myres[k,28]=mean(dat_i$nirSK[dat_i$treatment=='warmed'])/mean(dat_i$nirSK[dat_i$treatment=='control'])
      myres[k,29]=(mean(dat_i$AOA[dat_i$treatment=='warmed'])/mean(dat_i$nirSK[dat_i$treatment=='warmed']))/
        (mean(dat_i$AOA[dat_i$treatment=='control'])/mean(dat_i$nirSK[dat_i$treatment=='control']))
      myres[k,30]=mean(dat_i$nosZ[dat_i$treatment=='warmed'])/mean(dat_i$nosZ[dat_i$treatment=='control'])
      myres[k,31]=mean(dat_i$nirSK_nosZ[dat_i$treatment=='warmed'])/mean(dat_i$nirSK_nosZ[dat_i$treatment=='control'])
      myres[k,32]=mean(dat_i$EOC[dat_i$treatment=='warmed'])/mean(dat_i$EOC[dat_i$treatment=='control'])
      myres[k,33]=mean(dat_i$EON[dat_i$treatment=='warmed'])/mean(dat_i$EON[dat_i$treatment=='control'])
      
      myres[k,34]=mean(dat_i$NH4_N[dat_i$treatment=='warmed'])/mean(dat_i$NH4_N[dat_i$treatment=='control'])
      myres[k,35]=mean(dat_i$NO3_N[dat_i$treatment=='warmed'])/mean(dat_i$NO3_N[dat_i$treatment=='control'])
      myres[k,36]=mean(dat_i$NH4NO3[dat_i$treatment=='warmed'])/mean(dat_i$NH4NO3[dat_i$treatment=='control'])
      myres[k,37]=mean(dat_i$soil_temp[dat_i$treatment=='warmed'])/mean(dat_i$soil_temp[dat_i$treatment=='control'])
      myres[k,38]=mean(dat_i$soil.water_dry[dat_i$treatment=='warmed'])/mean(dat_i$soil.water_dry[dat_i$treatment=='control'])
      myres[k,39]=mean(dat_i$MBC[dat_i$treatment=='warmed'])/mean(dat_i$MBC[dat_i$treatment=='control'])
      myres[k,40]=mean(dat_i$MBN[dat_i$treatment=='warmed'])/mean(dat_i$MBN[dat_i$treatment=='control'])
      myres[k,41]=(mean(dat_i$MBC[dat_i$treatment=='warmed'])/mean(dat_i$MBN[dat_i$treatment=='warmed']))/
        (mean(dat_i$MBC[dat_i$treatment=='control'])/mean(dat_i$MBN[dat_i$treatment=='control']))
      
      myres[k,42]=mean(dat_i$NO_N[dat_i$treatment=='warmed'])/mean(dat_i$NO_N[dat_i$treatment=='control'])
      
      ###比值小于0的进行比值x+sqrt(1+x^2)
      ###全部＋10，影响太大
      # myres[k,43]=(mean(dat_i$N2O_N[dat_i$treatment=='warmed'])+10)/(mean(dat_i$N2O_N[dat_i$treatment=='control'])+10)
      # myres[k,43]=ifelse(mean(dat_i$N2O_N[dat_i$treatment=='warmed'])*mean(dat_i$N2O_N[dat_i$treatment=='control'])>0,
      #                    mean(dat_i$N2O_N[dat_i$treatment=='warmed'])/mean(dat_i$N2O_N[dat_i$treatment=='control']),
      #                    mean(dat_i$N2O_N[dat_i$treatment=='warmed'])/mean(dat_i$N2O_N[dat_i$treatment=='control'])+
      #                      sqrt(1+(mean(dat_i$N2O_N[dat_i$treatment=='warmed'])/mean(dat_i$N2O_N[dat_i$treatment=='control']))^2))
      # 
      myres[k,43]=ifelse(mean(dat_i$N2O_N[dat_i$treatment=='warmed'])*mean(dat_i$N2O_N[dat_i$treatment=='control'])>0 & mean(dat_i$N2O_N[dat_i$treatment=='control'])>0,
                         mean(dat_i$N2O_N[dat_i$treatment=='warmed'])/mean(dat_i$N2O_N[dat_i$treatment=='control']),NA)
      
      
      myres[k,44]=mean(dat_i$gasN[dat_i$treatment=='warmed'])/mean(dat_i$gasN[dat_i$treatment=='control'])
      
      myres[k,45]=mean(dat_i$CO2_C[dat_i$treatment=='warmed'])/mean(dat_i$CO2_C[dat_i$treatment=='control'])
      
      myres[k,46]=mean(dat_i$soil.water_dry[dat_i$treatment=='warmed'])-mean(dat_i$soil.water_dry[dat_i$treatment=='control'])
      
      #min
      myres[k,47]=mean(dat_i$min[dat_i$treatment=='warmed'])/mean(dat_i$min[dat_i$treatment=='control'])
      # #nit
      myres[k,48]=mean(dat_i$nit[dat_i$treatment=='warmed'])/mean(dat_i$nit[dat_i$treatment=='control'])
      
      #N2
      myres[k,49]=ifelse(mean(dat_i$N2_N[dat_i$treatment=='warmed'])>0 & mean(dat_i$N2_N[dat_i$treatment=='control'])>0,
                         mean(dat_i$N2_N[dat_i$treatment=='warmed'])/mean(dat_i$N2_N[dat_i$treatment=='control']),
                         mean(dat_i$N2_N[dat_i$treatment=='warmed'])/mean(dat_i$N2_N[dat_i$treatment=='control'])+
                           sqrt(1+(mean(dat_i$N2_N[dat_i$treatment=='warmed'])/mean(dat_i$N2_N[dat_i$treatment=='control']))^2))
      
      myres[k,50]=(mean(dat_i$AOA[dat_i$treatment=='warmed'])+mean(dat_i$AOB[dat_i$treatment=='warmed']))/
        (mean(dat_i$AOA[dat_i$treatment=='control'])+mean(dat_i$AOB[dat_i$treatment=='control']))
      
      k=k+1
    }
  }
}

# myres<-na.omit(myres)

names(myres)<-c('Date','soillayer','plot_ab','season',
                'AOB.c','AOA.c','nirK.c','nirS.c','nirSK.c','nosZ.c','nirSK_nosZ.c','EOC.c','NH4_N.c','NO3_N.c',
                'NH4NO3.c','temp.c','moisture.c','MBC.c','MBN.c','NO_N.c','N2O_N.c','gas_N.c','CO2_C.c',
                'AOB.r','AOA.r','nirK.r','nirS.r','nirSK.r','AOA_nirSK.r','nosZ.r','nirSK_nosZ.r','EOC.r','EON.r','NH4_N.r','NO3_N.r',
                'NH4NO3.r','temp.r','moisture.r','MBC.r','MBN.r','MBC_MBN.r','NO_N.r','N2O_N.r','gas_N.r','CO2_C.r','moisture.d',
                'min.r','nit.r','N2_N.r','amoA.r')


myres$moisture.per<-myres$moisture.d/myres$moisture.c


myres$MBC.lnr<-log(myres$MBC.r)
myres$MBN.lnr<-log(myres$MBN.r)

myres$MBC.r[is.na(myres$MBC.r)]=1
myres$MBN.r[is.na(myres$MBN.r)]=1

myres$MBC.lnr[is.na(myres$MBC.lnr)]=0
myres$MBN.lnr[is.na(myres$MBN.lnr)]=0


myres$EOC.lnr<-log(myres$EOC.r)
myres$EON.lnr<-log(myres$EON.r)

myres$AOA.lnr<-log(myres$AOA.r)
myres$AOB.lnr<-log(myres$AOB.r)
myres$amoA.lnr<-log(myres$amoA.r)


myres$nirS.lnr<-log(myres$nirS.r)
myres$nirK.lnr<-log(myres$nirK.r)
myres$nirSK.lnr<-log(myres$nirSK.r)
myres$nosZ.lnr<-log(myres$nosZ.r)
myres$nirSK_nosZ.lnr<-log(myres$nirSK_nosZ.r)
myres$min.lnr<-log(myres$min.r)
myres$nit.lnr<-log(myres$nit.r)
myres$NH4.lnr<-log(myres$NH4_N.r)
myres$NO3.lnr<-log(myres$NO3_N.r)
myres$NH4NO3.lnr<-log(myres$NH4NO3.r)
myres$moisture.lnr<-log(myres$moisture.per+1)

myres$NO.lnr<-log(myres$NO_N.r)
myres$N2O.lnr<-log(myres$N2O_N.r)
myres$N2.lnr<-log(myres$N2_N.r)

myres$gasN.lnr<-log(myres$gas_N.r)


str(myres)

###NO 0-10cm
model.step<-step(lm(NO_N.r~temp.r+NO3_N.r+moisture.c+moisture.per+moisture.d+NH4_N.r+NO3_N.r+AOA.r+AOB.r+MBC.r+MBN.r,
                    data=myres[myres$soillayer=='0-10cm' & myres$Date!='2021-05-22',]),direction='both')


model.step<-step(lm(NO.lnr~temp.r+moisture.c+moisture.lnr+NH4.lnr+NO3.lnr+EON.lnr+AOA.lnr+AOB.lnr+nirS.lnr+nirK.lnr+nosZ.lnr+MBN.lnr,
                    data=myres[myres$soillayer=='0-10cm' & myres$Date!='2021-05-22',]),direction='both')

summary(lm(NO.lnr ~ moisture.lnr + NH4.lnr + EON.lnr + AOB.lnr,
           data=myres[myres$soillayer=='0-10cm' & myres$Date!='2021-05-22',]))   #NO的响应与

summary(lm(NO.lnr ~ temp.r + moisture.lnr  + EON.lnr + EOC.lnr  + nirS.lnr + nosZ.lnr +NH4.lnr,
           data=myres[myres$soillayer=='0-10cm' & myres$Date!='2021-05-22',]))   #NO的响应与


myres$month<-format(as.Date(myres$Date),"%m")
ggplot(myres[myres$soillayer=='0-10cm' &
               log(myres$NO_N.r)>-0.66,],
       aes(moisture.d,log(NO_N.r)))+geom_point()+ #facet_wrap(month~.)+
  geom_smooth(method='lm',formula = 'y~x')+xlim(-0.4,0.4)+
  stat_poly_eq(formula = y~x, 
               aes(label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 ..p.value.label..,
                                 sep = "~~~")), 
               parse = TRUE)

str(myres)
###0-10 cm
{
  #######NO
  unique(myres$Date)
  myres$month<-substr(myres$Date,6,7)
  # kdat<-myres[myres$NO_N.c>0 & myres$soillayer=='0-10cm',]
  
  myres$date_pair<-paste(myres$Date,myres$plot_ab,sep='_')
  
  # myres$NO_N.r[myres$date_pair %in% c('2022-10-01_pair1')]=1.2
  # myres$NO_N.r[myres$date_pair %in% c('2019-07-21_pair1')]=1
  
  p1<-ggplot(myres[myres$soillayer=='0-10cm',],
             aes(moisture.r,log(NO_N.r),col=soillayer))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(0,1.5),breaks = seq(0,1.5,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-1,1),labels=c('-1','-.5','0','.5','1'),breaks = seq(-1,1,.5),expand=c(0,0))+
    geom_line(aes(y=-3))+
    geom_smooth(linewidth=0.3,alpha=0.1,
                method='lm',formula = y~x,se=T)+
    stat_poly_eq(size=2.6,parse = TRUE,col='black',
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")),label.x.npc = 0.25, label.y.npc = c(0.97))+
    stat_poly_eq(formula = y~x, size=2.6,col='black',
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")),
                 parse = TRUE,label.x.npc = 0.25,label.y.npc = c(0.87))+
    labs(x=expression(paste('Ratio of moisture (W/C)')),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'none');p1
  
  
  #######N2O
  datk<-myres[myres$soillayer=='0-10cm',]
  p2<-ggplot(myres[myres$soillayer=='0-10cm',],
             aes(moisture.per+1,log(N2O_N.r)))+
    geom_point(data=myres[myres$soillayer=='0-10cm',],shape=1,col='darkred')+
    scale_x_continuous(limits = c(0,1.5),breaks = seq(0,1.5,0.5),expand=c(0,0))+    
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm',],
                linewidth=0.3,alpha=0.1,col='darkred',method='lm',formula = 'y~x',se=T)+ #xlim(-0.2,0.2)+ylim(-2,2)+
    stat_poly_eq(data=myres[myres$N2O_N.c>0 & myres$soillayer=='0-10cm',],
                 size=2.6,parse = TRUE,formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")),label.x.npc = 0.05, label.y.npc = c(0.17))+
    stat_poly_eq(data=myres[myres$N2O_N.c>0 & myres$soillayer=='0-10cm',],
                 size=2.6,formula = y~x, aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.07))+
    labs(x=expression(paste('Ratio  of moisture (W/C)')),y=expression(paste('ln RR of N'[2],'O')))+
    mythem_sci;p2
  
  
  
  ######N2
  p3<-ggplot(myres[myres$month %in% c('05','06','07','08','09','10','11') & myres$soillayer=='0-10cm',],
             aes(moisture.per+1,log(N2_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_smooth(data=myres[myres$N2O_N.c>0 & myres$soillayer=='0-10cm'& log(myres$N2O_N.r)> -4,],linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+ 
    stat_poly_eq(data=myres[myres$N2O_N.c>0 & myres$soillayer=='0-10cm' & log(myres$N2O_N.r)> -4,],size=2.6,parse = TRUE,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")),label.x.npc = 0.05, label.y.npc = c(0.17))+
    stat_poly_eq(data=myres[myres$N2O_N.c>0 & myres$soillayer=='0-10cm' & log(myres$N2O_N.r)> -4,],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.07))+
    geom_point(data=myres[myres$soillayer=='0-10cm' & myres$N2O_N.c>0 & log(myres$N2O_N.r)> -4,],shape=1,col='darkred')+
    geom_point(data=myres[myres$soillayer=='0-10cm' & myres$N2O_N.c>0 & log(myres$N2O_N.r)<= -4,],shape=19,col='darkred')+
    scale_x_continuous(limits = c(0,1.5),breaks = seq(0,1.5,0.5),expand=c(0,0))+    
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    labs(x=expression(paste('Ratio of moisture (W/C)')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'right');p3
  
  myres$DOC.r
  p4<-ggplot(myres[myres$month %in% c('05','06','07','08','09','10','11') & myres$soillayer=='0-10cm',],
             aes(log(EOC.r),log(NO_N.r)))+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],col='darkred',
                linewidth=0.3,alpha=0.1,method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                 size=2.6,formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")),label.x.npc = 0.25, label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                 size=2.6,aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,
                 label.x.npc = 0.25,label.y.npc = c(0.87))+
    scale_y_continuous(limits = c(-1,1),labels=c('-1','-.5','0','.5','1'),breaks = seq(-1,1,.5),expand=c(0,0))+
    labs(x=expression(paste('ln RR of EOC')),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'right');p4
  
  
  p5<-ggplot(myres[myres$month %in% c('05','06','07','08','09','10','11') & myres$soillayer=='0-10cm',],
             aes(log(EOC.r),log(N2O_N.r)))+
    geom_point(shape=1,col='darkred')+
    geom_smooth(data=myres[myres$soillayer=='0-10cm',],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+ #xlim(-0.2,0.2)+ylim(-2,2)+
    scale_x_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),expand=c(0,0))+
    scale_y_continuous(breaks = seq(-4,2,2),expand=c(0,0))+
    coord_cartesian(clip="on",ylim=c(-4,2))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm',],
                 size=2.6,formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")),label.x.npc = 0.05, label.y.npc = c(0.17))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm',],
                 size=2.6,aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.07))+
    labs(x=expression(paste('ln RR of EOC')),y=expression(paste('ln RR of N'[2],'O')))+
    mythem_sci+theme(legend.position = 'right');p5
  
  p6<-ggplot(myres[myres$month %in% c('05','06','07','08','09','10','11') & myres$soillayer=='0-10cm',],
             aes(log(EOC.r),log(N2_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+ #xlim(-0.2,0.2)+ylim(-2,2)+
    scale_x_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                 size=2.6, formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.17))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                 size=2.6,aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.07))+
    labs(x=expression(paste('ln RR of EOC')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'right');p6
  
  myres$NH4_N.r
  ####NH4_N
  p7<-ggplot(myres[myres$soillayer=='0-10cm',],aes(log(NH4_N.r),log(NO_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+ #xlim(-0.2,0.2)+ylim(-2,2)+
    scale_x_continuous(limits = c(-2,1),breaks = seq(-2,1,1),expand=c(0,0))+
    scale_y_continuous(limits = c(-1,1),labels=c('-1','-.5','0','.5','1'),breaks = seq(-1,1,.5),expand=c(0,0))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm',],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.25,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm',],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.25,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of NH'[4]^'+','-N')),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'right');p7
  
  myres$N2O_N.c
  p8<-ggplot(myres[myres$soillayer=='0-10cm' ,],aes(log(NH4_N.r),log(N2O_N.r),col=soillayer,fill=soillayer))+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-2,1),breaks = seq(-2,1,1),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    labs(x=expression(paste('ln RR of NH'[4]^'+','-N')),y=expression(paste('ln RR of N'[2],'O')))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+ #xlim(-0.2,0.2)+ylim(-2,2)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm',],size=2.6,col='black',
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.17))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm',],size=2.6,col='black',
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.07))+
    mythem_sci+theme(legend.position = 'right');p8
  
  p9<-ggplot(myres[myres$month %in% c('05','06','07','08','09','10','11') & myres$soillayer=='0-10cm',],
             aes(log(NH4_N.r),log(N2_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-2,1),breaks = seq(-2,1,1),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(size=2.6,formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.17))+
    stat_poly_eq(size=2.6,aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.07))+
    labs(x=expression(paste('ln RR of NH'[4]^'+','-N')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'right');p9
  
  
  p13<-ggplot(myres[myres$soillayer=='0-10cm',],aes(amoA.lnr,log(NO_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,2,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-1,1),labels=c('-1','-.5','0','.5','1'),breaks = seq(-1,1,.5),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x+I(x^2), aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.25,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,formula = y~x+I(x^2),
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.25,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of amoA')),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'none');p13
  
  p14<-ggplot(myres[myres$soillayer=='0-10cm',],aes(amoA.lnr,log(N2O_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    # scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of amoA')),y=expression(paste('ln RR of N'[2],'O')))+
    mythem_sci+theme(legend.position = 'none');p14
  
  
  p15<-ggplot(myres[myres$soillayer=='0-10cm',],aes(amoA.lnr,log(N2_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of amoA')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'none');p15
  
  p16<-ggplot(myres[myres$soillayer=='0-10cm',],aes(log(AOB.r),log(NO_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    # scale_x_continuous(limits = c(-1,3),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-1,1),labels=c('-1','-.5','0','.5','1'),breaks = seq(-1,1,.5),expand=c(0,0))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.25,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.25,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of AOB')),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'none');p16
  
  p17<-ggplot(myres[myres$soillayer=='0-10cm',],aes(log(AOB.r),log(N2O_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    # scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of AOB')),y=expression(paste('ln RR of N'[2],'O')))+
    mythem_sci+theme(legend.position = 'none');p17
  
  
  p18<-ggplot(myres[myres$soillayer=='0-10cm',],aes(log(nirS.r),log(N2_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of AOB')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'none');p18
  
  p19<-ggplot(myres[myres$soillayer=='0-10cm',],aes(amoA.lnr,log(NO_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of AOB')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'none');p19
  
  myres$MBC.lnr
  p20<-ggplot(myres[myres$soillayer=='0-10cm',],aes(MBN.lnr,log(NO_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of MBN')),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'none');p20
  
  
  p21<-ggplot(myres[myres$soillayer=='0-10cm',],aes(MBN.lnr,log(N2O_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of MBN')),y=expression(paste('ln RR of N'[2],'O')))+
    mythem_sci+theme(legend.position = 'none');p21
  
  p22<-ggplot(myres[myres$soillayer=='0-10cm',],aes(MBN.lnr,log(N2_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of MBN')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'none');p22
  
  myres$CO2_C.r
  p23<-ggplot(myres[myres$soillayer=='0-10cm',],aes(log(CO2_C.r),log(NO_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of CO'[2])),y=expression(paste('ln RR of NO')))+
    mythem_sci+theme(legend.position = 'none');p23
  
  p24<-ggplot(myres[myres$soillayer=='0-10cm',],aes(log(CO2_C.r),log(N2O_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of CO'[2])),y=expression(paste('ln RR of N'[2],'O')))+
    mythem_sci+theme(legend.position = 'none');p24
  
  myres$NO3.lnr
  p24<-ggplot(myres[myres$soillayer=='0-10cm',],aes(MBN.lnr,log(N2O_N.r)))+
    # geom_hline(yintercept = 0,lty=2,linewidth=0.2)+
    # geom_vline(xintercept = 1,lty=2,linewidth=0.2)+
    geom_point(shape=1,col='darkred')+
    scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5),expand=c(0,0))+
    scale_y_continuous(limits = c(-4,2),breaks = seq(-4,2,2),expand=c(0,0))+
    geom_smooth(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],
                linewidth=0.3,alpha=0.1,col='darkred',
                method='lm',formula = 'y~x',se=T)+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 formula = y~x, aes(label = paste(..eq.label..,sep = "~~~")), label.x.npc = 0.05,label.y.npc = c(0.97))+
    stat_poly_eq(data=myres[myres$soillayer=='0-10cm' & myres$month %in% c('05','06','07','08','09','10','11'),],size=2.6,
                 aes(label = paste(..rr.label.., ..p.value.label..,sep = "~~~")), parse = TRUE,label.x.npc = 0.05,label.y.npc = c(0.87))+
    labs(x=expression(paste('ln RR of AOB')),y=expression(paste('ln RR of N'[2])))+
    mythem_sci+theme(legend.position = 'none');p24
  
  
  figure<-ggarrange(p1+theme(plot.margin = unit(c(5,0,0,7),'mm'))+rremove("xlab"),
                    p4+theme(plot.margin = unit(c(5,0,0,7),'mm'))+rremove("xlab")+rremove("ylab"),
                    p7+theme(plot.margin = unit(c(5,3,0,7),'mm'))+rremove("xlab")+rremove("ylab"),
                    
                    p2+theme(plot.margin = unit(c(3,0,0,7),'mm'))+rremove("xlab"),
                    p5+theme(plot.margin = unit(c(3,0,0,7),'mm'))+rremove("xlab")+rremove("ylab"),
                    p8+theme(plot.margin = unit(c(3,3,0,7),'mm'))+rremove("xlab")+rremove("ylab"),
                    
                    p3+theme(plot.margin = unit(c(3,0,2,7),'mm')),
                    p6+theme(plot.margin = unit(c(3,0,2,7),'mm'))+rremove("ylab"),
                    p9+theme(plot.margin = unit(c(3,3,2,7),'mm'))+rremove("ylab"),
                    
                    labels = c(paste("(",letters[1:9],")",sep="")),
                    label.x = c(0.05,rep(c(0.001),3)),label.y = 0.98,
                    ncol = 3, nrow = 3,align = "none",   ##"v"竖直对齐
                    font.label = list(size = 10, color ="black"),
                    widths = c(6.2,5.8,5.8), heights = c(5,5,5.4),
                    common.legend = FALSE,legend.grob = get_legend(p4, position = 'bottom'),legend='none');figure
  
  ggsave("0-10cm lnRR with moisture in samples_new2.pdf", figure, width =10, height =7,
         device=cairo_pdf)
  
  getwd() 
}


#########随机森林验证
str(myres)
unique(myres$Date)

set.seed(123)

library(randomForest)
library(rfPermute)

str(myres)



####NO  Olayer

myres<-myres[myres$Date!='2021-05-22',]

unique(myres$soillayer)

unique(myres$Date)

dat_1= myres[myres$soillayer=='O',
             c('moisture.c','moisture.per','MBN.r','EOC.r',
               'AOB.r','AOA.r','nirS.r','nirK.r','nosZ.r',
               'NH4_N.r','NO3_N.r','NO_N.r')]


str(dat_1)

dat_1<-na.omit(dat_1)
NO_rf <- randomForest(log(NO_N.r) ~ ., data= dat_1,
                      importance=TRUE,proximity=TRUE)
NO_rf
# 通过建立500棵决策树，发现模型对目标变量（richness）变化的对的解释量为26.8%。

#通过置换检验，检测该模型的显著性代码如下：
set.seed(123)
library(rfUtilities);library(rfPermute)
names(dat_1)[14]

NO_perm <- rf.significance(NO_rf, dat_1[,-which(names(dat_1)=='NO_N.r')], nperm=99, ntree=501)
NO_perm
#结果表明模型的显著性为0.01，模型的解释量为25.45%。

#最后检测模型中每个变量对目标变量的重要性，代码如下：
set.seed(123)

NO_rfP<- rfPermute(log(NO_N.r) ~ ., data =dat_1, 
                   ntree = 500,na.action = na.omit, nrep = 100,num.cores = 1)

NO_dat <- importance(NO_rfP, sort.by = NULL, decreasing = TRUE)

# gasN_dat$variables<-
library(tidyr)

NO_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))->dai_no

substrRight<- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}  #选取右数最后一个字符到第n个字符

dai_no$x_names=NA
i=1
for(i in 1:length(dai_no$names)){
  dai_no$x_names[i]<-ifelse(substrRight(as.character(dai_no$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_no$names[i]),"[.]")[[1]][1],')',sep=''),
                            ifelse(substrRight(as.character(dai_no$names[i]),1)=='c',strsplit(as.character(dai_no$names[i]),"[.]")[[1]][1],
                                   ifelse(substrRight(as.character(dai_no$names[i]),1)=='d','Dmoisture',
                                          'Percent change of moisture')))
  i=i+1}

levels(dai_no$names)
str(dai_no)

ggplot(dai_no, aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
  geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
  scale_y_continuous(limits=c(-8,16),breaks=seq(-8,16,8),expand=c(0,0))+
  # scale_x_discrete(limits=c("moisture.per","AOA.r","NH4NO3.r","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",  
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"),
  #                  labels=c("R(moisture)","R(AOA)","R(NH4NO3)","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",  
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"))+
  geom_bar(aes(fill = group),stat = "identity")+
  scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
  geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 0.3,X.IncMSE-0.4),label = label))+
  labs(x = "", y = "% IncMSE",title = expression(paste('ln RR of NO')))+
  coord_flip()+mythem_sci->a0;a0


###0-10cm  NO
dat_1= myres[myres$soillayer=='0-10cm',
             c('moisture.c','moisture.per','MBN.r','EOC.r',
               'nirS.r','nirK.r','nosZ.r','AOA.r','AOB.r',
               'NH4_N.r','NO3_N.r','NO_N.r')]


dat_1<-na.omit(dat_1)
NO_rf <- randomForest(log(NO_N.r) ~ ., data= dat_1,
                      importance=TRUE,proximity=TRUE)
NO_rf
# 通过建立500棵决策树，发现模型对目标变量（richness）变化的对的解释量为26.8%。

#通过置换检验，检测该模型的显著性代码如下：
set.seed(123)
library(rfUtilities);library(rfPermute)
names(dat_1)[18]
NO_perm <- rf.significance(NO_rf, dat_1[,-which(names(dat_1)=='NO_N.r')], nperm=99, ntree=501)
NO_perm
#结果表明模型的显著性为0.01，模型的解释量为25.45%。

#最后检测模型中每个变量对目标变量的重要性，代码如下：
set.seed(123)

NO_rfP<- rfPermute(log(NO_N.r) ~ ., data =dat_1, 
                   ntree = 500,na.action = na.omit, nrep = 100,num.cores = 1)

NO_dat <- importance(NO_rfP, sort.by = NULL, decreasing = TRUE)

# gasN_dat$variables<-
library(tidyr)

NO_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))->dai_no

substrRight<- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}  #选取右数最后一个字符到第n个字符

dai_no$x_names=NA
i=1
for(i in 1:length(dai_no$names)){
  dai_no$x_names[i]<-ifelse(substrRight(as.character(dai_no$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_no$names[i]),"[.]")[[1]][1],')',sep=''),
                            ifelse(substrRight(as.character(dai_no$names[i]),1)=='c',strsplit(as.character(dai_no$names[i]),"[.]")[[1]][1],
                                   ifelse(substrRight(as.character(dai_no$names[i]),1)=='d','Dmoisture',
                                          'Percent change of moisture')))
  i=i+1}

levels(dai_no$names)
str(dai_no)

ggplot(dai_no, aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
  geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
  scale_y_continuous(limits=c(-8,16),breaks=seq(-8,16,8),expand=c(0,0))+
  # scale_x_discrete(limits=c("moisture.per","AOA.r","NH4NO3.r","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",  
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"),
  #                  labels=c("R(moisture)","R(AOA)","R(NH4NO3)","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",  
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"))+
  geom_bar(aes(fill = group),stat = "identity")+
  scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
  geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 0.3,X.IncMSE-0.4),label = label))+
  labs(x = "", y = "% IncMSE",title = expression(paste('ln RR of NO')))+
  coord_flip()+mythem_sci->a1;a1


####N2O  Olayer
dat_1= myres[myres$soillayer=='O' & myres$N2O_N.c>0,
             c('moisture.c','moisture.per','MBN.r','EOC.r',
               'nirS.r','nirK.r','nosZ.r','AOA.r','AOB.r',
               'NH4_N.r','NO3_N.r','N2O_N.r')]

dat_1<-dat_1[dat_1$N2O_N.r>0 & log(dat_1$N2O_N.r)>-2,]  #为什么？
dat_1<-na.omit(dat_1)
N2O_rf <- randomForest(log(N2O_N.r) ~ ., data= dat_1,
                       importance=TRUE,proximity=TRUE)
N2O_rf
# 通过建立500棵决策树，发现模型对目标变量（richness）变化的对的解释量为26.8%。

#通过置换检验，检测该模型的显著性代码如下：
set.seed(123)
library(rfUtilities);library(rfPermute)
names(dat_1)[18]
N2O_perm <- rf.significance(N2O_rf,  dat_1[,-which(names(dat_1)=='N2O_N.r')], nperm=99, ntree=501)
N2O_perm
#结果表明模型的显著性为0.01，模型的解释量为25.45%。

#最后检测模型中每个变量对目标变量的重要性，代码如下：
set.seed(123)

N2O_rfP<- rfPermute(log(N2O_N.r) ~ ., data =dat_1, 
                    ntree = 500,na.action = na.omit, nrep = 100,num.cores = 1)

N2O_dat <- importance(N2O_rfP, sort.by = NULL, decreasing = TRUE)

# gasN_dat$variables<-
library(tidyr)

# rownames(N2O_dat)<-c("Percent change of moisture","R (DOC)",expression(paste(Delta,'moisture')),"Soil moisture","NO3_N.r",
#                      "nirSK_nosZ.r","nosZ.r","AOA.r", "NH4NO3.r","MBN.r","temp.c","nirS.r",
#                      "MBC.r","nirSK.r","temp.r","NH4_N.r")

substrRight<- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}  #选取右数最后一个字符到第n个字符
substrRight("gcb or sbb",1)

N2O_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))->dai_n2o

dai_n2o$x_names=NA
i=1
for(i in 1:length(dai_n2o$names)){
  dai_n2o$x_names[i]<-ifelse(substrRight(as.character(dai_n2o$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_n2o$names[i]),"[.]")[[1]][1],')',sep=''),
                             ifelse(substrRight(as.character(dai_n2o$names[i]),1)=='c',strsplit(as.character(dai_n2o$names[i]),"[.]")[[1]][1],
                                    ifelse(substrRight(as.character(dai_n2o$names[i]),1)=='d','Dmoisture',
                                           'Percent change of moisture')))
  i=i+1}

levels(dai_n2o$names)
str(dai_n2o)

ggplot(dai_n2o,aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
  geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
  scale_y_continuous(limits=c(-8,16),breaks=seq(-8,16,8),expand=c(0,0))+
  # scale_x_discrete(limits=c("moisture.per","AOA.r","NH4NO3.r","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"),
  #                  labels=c("R(moisture)","R(AOA)","R(NH4NO3)","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"))+
  geom_bar(aes(fill = group),stat = "identity")+
  scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
  geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 0.3,X.IncMSE-0.4),label = label))+
  labs(x = "", y = "% IncMSE",title = expression(paste('ln RR of N'[2],'O')))+
  coord_flip()+mythem_sci->b0;b0

###0-10cm  N2O

# dat_1= myres[myres$soillayer=='0-10cm' & myres$N2O_N.r>0,
#              c('temp.c','temp.r','moisture.c','moisture.r','MBC.lnr','MBN.lnr','DOC.lnr',
#                'AOA.lnr','nirS.lnr','nirSK.lnr','nosZ.lnr','nirSK_nosZ.lnr',
#                'NH4.lnr','NO3.lnr','NH4NO3.lnr','N2O_N.r')]

dat_1= myres[myres$soillayer=='0-10cm',
             c('moisture.c','moisture.per','MBN.r','EOC.r',
               'nirS.r','nirK.r','nosZ.r','AOA.r','AOB.r',
               'NH4_N.r','NO3_N.r','N2O_N.r')]

# dat_1$moisture.per <- dat_1$moisture.r-1

# dat_1<-dat_1[dat_1$N2O_N.r>0 & log(dat_1$N2O_N.r)>-2,]  #为什么？

dat_1<-na.omit(dat_1)


N2O_rf <- randomForest(log(N2O_N.r) ~ ., data= dat_1,
                       importance=TRUE,proximity=TRUE)
N2O_rf
# 通过建立500棵决策树，发现模型对目标变量（richness）变化的对的解释量为26.8%。

#通过置换检验，检测该模型的显著性代码如下：
set.seed(123)
library(rfUtilities);library(rfPermute)
names(dat_1)

which(names(dat_1)=='N2O_N.r')

N2O_perm <- rf.significance(N2O_rf, dat_1[,-which(names(dat_1)=='N2O_N.r')], nperm=99, ntree=501)
N2O_perm
#结果表明模型的显著性为0.01，模型的解释量为25.45%。

#最后检测模型中每个变量对目标变量的重要性，代码如下：
set.seed(123)

N2O_rfP<- rfPermute(log(N2O_N.r) ~ ., data =dat_1, 
                    ntree = 500,na.action = na.omit, nrep = 100,num.cores = 1)

N2O_dat <- importance(N2O_rfP, sort.by = NULL, decreasing = TRUE)

# gasN_dat$variables<-
library(tidyr)

# rownames(N2O_dat)<-c("Percent change of moisture","R (DOC)",expression(paste(Delta,'moisture')),"Soil moisture","NO3_N.r",
#                      "nirSK_nosZ.r","nosZ.r","AOA.r", "NH4NO3.r","MBN.r","temp.c","nirS.r",
#                      "MBC.r","nirSK.r","temp.r","NH4_N.r")

substrRight<- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}  #选取右数最后一个字符到第n个字符
substrRight("gcb or sbb",1)

N2O_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))->dai_n2o

dai_n2o$x_names=NA
i=1
for(i in 1:length(dai_n2o$names)){
  dai_n2o$x_names[i]<-ifelse(substrRight(as.character(dai_n2o$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_n2o$names[i]),"[.]")[[1]][1],')',sep=''),
                             ifelse(substrRight(as.character(dai_n2o$names[i]),1)=='c',strsplit(as.character(dai_n2o$names[i]),"[.]")[[1]][1],
                                    ifelse(substrRight(as.character(dai_n2o$names[i]),1)=='d','Dmoisture',
                                           'Percent change of moisture')))
  i=i+1}

levels(dai_n2o$names)
str(dai_n2o)

ggplot(dai_n2o,aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
  geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
  scale_y_continuous(limits=c(-8,16),breaks=seq(-8,16,8),expand=c(0,0))+
  # scale_x_discrete(limits=c("moisture.per","AOA.r","NH4NO3.r","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"),
  #                  labels=c("R(moisture)","R(AOA)","R(NH4NO3)","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",
  #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"))+
  geom_bar(aes(fill = group),stat = "identity")+
  scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
  geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 0.3,X.IncMSE-0.4),label = label))+
  labs(x = "", y = "% IncMSE",title = expression(paste('ln RR of N'[2],'O')))+
  coord_flip()+mythem_sci->b1;b1


####N2
unique(myres$soillayer)
dat_1= myres[myres$soillayer=='O',
             c('moisture.c','moisture.per','MBN.r','EOC.r',
               'nirS.r','nirK.r','nosZ.r','AOA.r','AOB.r',
               'NH4_N.r','NO3_N.r','N2_N.r')]

# dat_1 <- dat_1[log(dat_1$N2_N.r)> -2,]
dat_1<-na.omit(dat_1)

gasN_rf <- randomForest(log(N2_N.r) ~ ., data= dat_1,
                        importance=TRUE,proximity=TRUE)
gasN_rf
# 通过建立500棵决策树，发现模型对目标变量（richness）变化的对的解释量为26.8%。

#通过置换检验，检测该模型的显著性代码如下：
set.seed(123)
library(rfUtilities);library(rfPermute)
names(dat_1)
gasN_perm <- rf.significance(gasN_rf, dat_1[,-which(names(dat_1)=='N2_N.r')], nperm=100, ntree=501)
gasN_perm
#结果表明模型的显著性为0.01，模型的解释量为25.45%。

#最后检测模型中每个变量对目标变量的重要性，代码如下：
set.seed(123)

gasN_rfP<- rfPermute(N2_N.r ~ ., data =dat_1, 
                     ntree = 500,na.action = na.omit, nrep = 100,num.cores = 1)

gasN_dat <- importance(gasN_rfP, sort.by = NULL, decreasing = TRUE)

# gasN_dat$variables<-rownames(gasN_dat)

# dai_nxo$x_names=NA
# i=1
# for(i in 1:length(dai_nxo$names)){
#   dai_nxo$x_names[i]<-ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],')',sep=''),
#                              ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='c',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],
#                                     ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='d','Dmoisture',
#                                            'Percent change of moisture')))
#   i=i+1}

gasN_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))->dai_nxo



dai_nxo$x_names=NA
i=1
for(i in 1:length(dai_nxo$names)){
  dai_nxo$x_names[i]<-ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],')',sep=''),
                             ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='c',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],
                                    ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='d','Dmoisture',
                                           'Percent change of moisture')))
  i=i+1}

levels(dai_nxo$names)
str(dai_nxo)  

ggplot(dai_nxo, aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
  geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
  scale_y_continuous(limits=c(-5,15),breaks=seq(-10,15,5),expand=c(0,0))+
  geom_bar(aes(fill = group),stat = "identity")+
  scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
  geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 0.3,X.IncMSE-0.4),label = label))+
  labs(x = NULL, y = "% IncMSE",title = expression(paste('ln RR of N'[2])))+
  coord_flip()+mythem_sci->t0;t0


###0-10cm
str(myres)
dat_1= myres[myres$soillayer=='0-10cm',
             c('moisture.c','moisture.per','MBN.r','EOC.r',
               'nirS.r','nirK.r','nosZ.r','AOA.r','AOB.r',
               'NH4_N.r','NO3_N.r','N2_N.r')]



gasN_rf <- randomForest(N2_N.r ~ ., data= dat_1,
                        importance=TRUE,proximity=TRUE)
gasN_rf
# 通过建立500棵决策树，发现模型对目标变量（richness）变化的对的解释量为26.8%。

#通过置换检验，检测该模型的显著性代码如下：
set.seed(123)
library(rfUtilities);library(rfPermute)
names(dat_1)[18]
gasN_perm <- rf.significance(gasN_rf, dat_1[,-which(names(dat_1)=='N2_N.r')], nperm=100, ntree=501)
gasN_perm
#结果表明模型的显著性为0.01，模型的解释量为17%。

#最后检测模型中每个变量对目标变量的重要性，代码如下：
set.seed(123)

gasN_rfP<- rfPermute(log(N2_N.r) ~ ., data =dat_1, 
                     ntree = 500,na.action = na.omit, nrep = 100,num.cores = 1)

gasN_dat <- importance(gasN_rfP, sort.by = NULL, decreasing = TRUE)

# gasN_dat$variables<-rownames(gasN_dat)

dai_nxo$x_names=NA
i=1
for(i in 1:length(dai_nxo$names)){
  dai_nxo$x_names[i]<-ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],')',sep=''),
                             ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='c',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],
                                    ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='d','Dmoisture',
                                           'Percent change of moisture')))
  i=i+1}

gasN_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names))->dai_nxo



dai_nxo$x_names=NA
i=1
for(i in 1:length(dai_nxo$names)){
  dai_nxo$x_names[i]<-ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='r', paste('R(',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],')',sep=''),
                             ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='c',strsplit(as.character(dai_nxo$names[i]),"[.]")[[1]][1],
                                    ifelse(substrRight(as.character(dai_nxo$names[i]),1)=='d','Dmoisture',
                                           'Percent change of moisture')))
  i=i+1}

levels(dai_nxo$names)
str(dai_nxo)  

ggplot(dai_nxo, aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
  geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
  scale_y_continuous(limits=c(-5,15),breaks=seq(-10,15,5),expand=c(0,0))+
  geom_bar(aes(fill = group),stat = "identity")+
  scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
  geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 0.3,X.IncMSE-0.4),label = label))+
  labs(x = NULL, y = "% IncMSE",title = expression(paste('ln RR of N'[2])))+
  coord_flip()+mythem_sci->t1;t1


# %IncMSE指平均增加的均方误差，也就是说去除该变量后，均方误差的增加量。
# 结果表明，pH对物种的丰富度的重要性更高，而其它变量对丰富度的变化均不重要

library(ggpubr)
figure<-ggarrange(a0+theme(plot.margin = unit(c(5,2,2,2),'mm')),
                  b0+theme(plot.margin = unit(c(5,2,2,2),'mm')),
                  t0+theme(plot.margin = unit(c(5,2,2,2),'mm')),
                  a1+theme(plot.margin = unit(c(5,5,2,2),'mm')),
                  b1+theme(plot.margin = unit(c(5,5,2,2),'mm')),
                  t1+theme(plot.margin = unit(c(5,5,2,2),'mm')),
                  labels = c("(a)","(b)","(c)","(d)","(e)","(f)"),
                  label.x = rep(c(0.080),6),label.y = 1.00,
                  ncol = 3, nrow = 2,align = "none",   ##"v"竖直对齐
                  font.label = list(size = 12, color ="black"),
                  widths = c(6,6,6), heights = c(5,5),legend = "none",
                  common.legend=T);figure
# setwd('D:/工作目录/202409/Manuscript_kai/Talk_20241113/Data and Code/Fig_6+Ext Dat fig_3_4_7_8')
ggsave("lnRR in RandomForest_all1.pdf", figure, width =9, height =6,
       device=cairo_pdf) 


################## Ext.Data.Fig.8_Random forest and liner regressions #############################
###终极一图
{ggplot(dai_no, aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
    geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
    scale_y_continuous(limits=c(-6,18),breaks=seq(-6,18,6),expand=c(0,0))+
    # scale_x_discrete(limits=c("moisture.per","AOA.r","NH4NO3.r","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",  
    #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"),
    #                  labels=c("R(moisture)","R(AOA)","R(NH4NO3)","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",  
    #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"))+
    geom_bar(aes(fill = group),stat = "identity",col='black',width=0.85,linewidth=0.13)+
    scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
    geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 1,X.IncMSE-0.4),label = ifelse(label!='ns',label,''),vjust=0.6))+
    labs(x = "", y = "% IncMSE")+
    coord_flip()+mythem_sci->a1;a1
  
  ggplot(dai_n2o,aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
    geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
    scale_y_continuous(limits=c(-6,18),breaks=seq(-6,18,6),expand=c(0,0))+
    # scale_x_discrete(limits=c("moisture.per","AOA.r","NH4NO3.r","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",
    #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"),
    #                  labels=c("R(moisture)","R(AOA)","R(NH4NO3)","moisture.c","nirSK.r","nirSK_nosZ.r","moisture.d",
    #                           "nirS.r","temp.c","DOC.r","MBN.r","temp.r","MBC.r","nosZ.r","NO3_N.r","NH4_N.r"))+
    geom_bar(aes(fill = group),stat = "identity",col='black',width=0.85,linewidth=0.13)+
    scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
    geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 1,X.IncMSE-0.4),label = ifelse(label!='ns',label,''),vjust=0.6))+
    labs(x = "", y = "% IncMSE")+
    coord_flip()+mythem_sci->b1;b1
  
  ggplot(dai_nxo, aes(x = reorder(x_names,X.IncMSE), y = X.IncMSE))+
    geom_hline(yintercept = 0,lty=1,linewidth=0.3)+
    scale_y_continuous(limits=c(-6,18),breaks=seq(-6,18,6),expand=c(0,0))+
    geom_bar(aes(fill = group),stat = "identity",col='black',width=0.85,linewidth=0.13)+
    scale_fill_manual(limits=c("In_sig","Sig"),values = c('skyblue','orangered'))+
    geom_text(aes(y = ifelse(X.IncMSE>0,X.IncMSE + 1,X.IncMSE-0.4),label = ifelse(label!='ns',label,''),vjust=0.6))+
    labs(x = NULL, y = "% IncMSE")+
    coord_flip()+mythem_sci->t1;t1
  
  figure<-ggarrange(a1+theme(plot.margin = unit(c(15,0,0,2),'mm'))+rremove("xlab"),
                    p1+theme(plot.margin = unit(c(15,0,0,2),'mm'))+rremove("xlab"),
                    p4+theme(plot.margin = unit(c(15,0,0,2),'mm'))+rremove("xlab")+rremove("ylab"),
                    p7+theme(plot.margin = unit(c(15,3,0,2),'mm'))+rremove("xlab")+rremove("ylab"),
                    
                    b1+theme(plot.margin = unit(c(3,0,2,2),'mm')),
                    p2+theme(plot.margin = unit(c(3,0,2,2),'mm')),
                    p5+theme(plot.margin = unit(c(3,0,2,2),'mm'))+rremove("ylab"),
                    p8+theme(plot.margin = unit(c(3,3,2,2),'mm'))+rremove("ylab"),
                    
                    labels = c(paste("(",letters[1:12],")",sep="")),
                    label.x = rep(c(0.001,0.14,0.09,0.09),3),label.y = c(rep(0.78,4),rep(0.95,8)),
                    ncol = 4, nrow = 2,align = "none",   ##"v"竖直对齐
                    font.label = list(size = 10, color ="black"),
                    widths = c(7,6.4,6,6.3), heights = c(6.2,5.4),
                    common.legend = T,legend.grob = get_legend(p9, position = 'none'),legend='none');figure
  
  ggsave("Randomforest lnRR with moisture in samples.pdf", figure, width =10, height =5,
         device=cairo_pdf) 
  
}


