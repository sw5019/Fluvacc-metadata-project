library(data.table)
library(ggplot2)
library(grid)
library(gplots)
library(sjPlot)
library(sjmisc)
library(emmeans)
library(PerformanceAnalytics)
theme_set(theme_sjplot())
setwd('G:/My Drive/Flu vaccine response project/Analysis on UGA1-5 cohorts')

## figure 1.
metadata<-as.data.frame(fread('1 - Seroconversion_prediction_data_frame_1461_entries.txt',header=T))
metadata<-metadata[metadata$PreVacc_status!='Mixed_vacc_status',]
pdf('Figure 1.pdf',width=10,height=9)
pushViewport(viewport(layout = grid.layout(3,4)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

metadata_v1<-metadata
metadata_v1$D28_Titer_composite<-log2(metadata_v1$D28_Titer_H1N1)+log2(metadata_v1$D28_Titer_H3N2)+log2(metadata_v1$D28_Titer_IBV_Yam)+log2(metadata_v1$D28_Titer_IBV_Vic)
metadata_v1_children<-metadata_v1[metadata_v1$Age<18,]
metadata_v1_adult1<-metadata_v1[metadata_v1$Age>=18 & metadata_v1$Age<65,]
metadata_v1_adult2<-metadata_v1[metadata_v1$Age>=65,]
metadata_v1_all<-data.frame(Titer=c(metadata_v1_children$Composite_baseline,metadata_v1_adult1$Composite_baseline,metadata_v1_adult2$Composite_baseline,
                                    metadata_v1_children$D28_Titer_composite,metadata_v1_adult1$D28_Titer_composite,metadata_v1_adult2$D28_Titer_composite),
                            Group=rep(c(rep('Children',dim(metadata_v1_children)[1]),rep('Adult1',dim(metadata_v1_adult1)[1]),rep('Adult2',dim(metadata_v1_adult2)[1])),2),
                            Time=rep(c('D0','D28'),each=dim(metadata_v1)[1]))
metadata_v1_all$Group<-factor(metadata_v1_all$Group,levels=c('Children','Adult1','Adult2'))
p1<-ggplot(metadata_v1_all,aes(x=Group,y=Titer))+stat_boxplot(aes(x=Group,y=Titer,fill=Time),geom='errorbar')+geom_boxplot(aes(fill=Time),outlier.shape = 1,size=.75)+
  scale_fill_manual(values=c('white','grey'))+
  labs(y='Composite HAI titer')+
  # scale_y_continuous(breaks=round(c(log2(10^3),log2(10^6),log2(10^9),log2(10^12)),0))+
  geom_hline(yintercept=round(4*log2(40),1),linetype='dashed')+
  theme_classic()+theme(legend.position = 'top',axis.title.x = element_blank())

metadata_v2<-data.frame(SC=c(metadata_v1_children$Composite_seroconversion,metadata_v1_adult1$Composite_seroconversion,metadata_v1_adult2$Composite_seroconversion),
                        Age=c(metadata_v1_children$Age,metadata_v1_adult1$Age,metadata_v1_adult2$Age),
                        BMI=c(metadata_v1_children$BMI,metadata_v1_adult1$BMI,metadata_v1_adult2$BMI),
                        Group=c(rep('Children',dim(metadata_v1_children)[1]),rep('Adult1',dim(metadata_v1_adult1)[1]),rep('Adult2',dim(metadata_v1_adult2)[1])))
metadata_v2$Group<-factor(metadata_v2$Group,levels=c('Children','Adult1','Adult2'))
p2<-ggplot(metadata_v2,aes(x=Group,y=SC))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  labs(y='Composite seroconversion')+
  scale_y_continuous(limits=c(-6,30))+
  theme_classic()+theme(axis.title.x = element_blank())

p3<-ggplot(metadata_v2,aes(x=Group,y=Age))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  labs(y='Age (yeas)')+
  scale_y_continuous(limits=c(10,90))+
  theme_classic()+theme(axis.title.x = element_blank())

p4<-ggplot(metadata_v2,aes(x=Group,y=BMI))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  labs(y='BMI (kg/m2)')+
  scale_y_continuous(limits=c(14,65))+
  theme_classic()+theme(axis.title.x = element_blank())

gender_frac_children<-table(metadata_v1_children$Gender)/dim(metadata_v1_children)[1]
gender_frac_adult1<-table(metadata_v1_adult1$Gender)/dim(metadata_v1_adult1)[1]
gender_frac_adult2<-table(metadata_v1_adult2$Gender)/dim(metadata_v1_adult2)[1]
race_frac_children<-table(metadata_v1_children$Race)/dim(metadata_v1_children)[1]
race_frac_adult1<-table(metadata_v1_adult1$Race)/dim(metadata_v1_adult1)[1]
race_frac_adult2<-table(metadata_v1_adult2$Race)/dim(metadata_v1_adult2)[1]
comor_frac_children<-table(metadata_v1_children$Comorbidities)/dim(metadata_v1_children)[1]
comor_frac_adult1<-table(metadata_v1_adult1$Comorbidities)/dim(metadata_v1_adult1)[1]
comor_frac_adult2<-table(metadata_v1_adult2$Comorbidities)/dim(metadata_v1_adult2)[1]
preVacc_frac_children<-table(metadata_v1_children$PreVacc_status)/dim(metadata_v1_children)[1]
preVacc_frac_adult1<-table(metadata_v1_adult1$PreVacc_status)/dim(metadata_v1_adult1)[1]
preVacc_frac_adult2<-table(metadata_v1_adult2$PreVacc_status)/dim(metadata_v1_adult2)[1]

metadata_v3<-data.frame(Gender_frac=c(as.numeric(gender_frac_children),as.numeric(gender_frac_adult1),as.numeric(gender_frac_adult2)),
                        Race_frac=c(as.numeric(race_frac_children),as.numeric(race_frac_adult1),as.numeric(race_frac_adult2)),
                        Comor_frac=c(as.numeric(comor_frac_children),as.numeric(comor_frac_adult1),as.numeric(comor_frac_adult2)),
                        PreVacc_frac=c(as.numeric(preVacc_frac_children),as.numeric(preVacc_frac_adult1),as.numeric(preVacc_frac_adult2)),
                        Gender=rep(names(gender_frac_children),3),
                        Race=rep(names(race_frac_children),3),
                        Comorbidities=rep(names(comor_frac_children),3),
                        PreVacc.=rep(names(preVacc_frac_children),3),
                        Group=rep(c('Children','Adult1','Adult2'),each=2))
metadata_v3$Group<-factor(metadata_v3$Group,levels=c('Children','Adult1','Adult2'))
metadata_v3$Gender<-factor(metadata_v3$Gender,levels=c('Male','Female'))
metadata_v3$Race<-factor(metadata_v3$Race,levels=c('White','Non-white'))
metadata_v3$Comorbidities<-factor(metadata_v3$Comorbidities,levels=c('Yes','No'))
metadata_v3$PreVacc.<-factor(metadata_v3$PreVacc.,levels=c('Naive','Prevaccinated'))
p5<-ggplot(metadata_v3,aes(x=Group,y=Gender_frac))+geom_bar(stat='identity',aes(fill=Gender),col='black',width=.75,size=.75)+
  scale_fill_manual(values=c('grey','black'))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(y='Fraction of participants')+
  theme_classic()+theme(legend.position = 'top',axis.title.x = element_blank())
p6<-ggplot(metadata_v3,aes(x=Group,y=Race_frac))+geom_bar(stat='identity',aes(fill=Race),col='black',width=.75,size=.75)+
  scale_fill_manual(values=c('grey','black'))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(y='Fraction of participants')+
  theme_classic()+theme(legend.position = 'top',axis.title.x = element_blank())
p7<-ggplot(metadata_v3,aes(x=Group,y=Comor_frac))+geom_bar(stat='identity',aes(fill=Comorbidities),col='black',width=.75,size=.75)+
  scale_fill_manual(values=c('grey','black'))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(y='Fraction of participants')+
  theme_classic()+theme(legend.position = 'top',axis.title.x = element_blank())
p8<-ggplot(metadata_v3,aes(x=Group,y=PreVacc_frac))+geom_bar(stat='identity',aes(fill=PreVacc.),col='black',width=.75,size=.75)+
  scale_fill_manual(values=c('grey','black'))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(y='Fraction of participants')+
  theme_classic()+theme(legend.position = 'top',axis.title.x = element_blank())

dose_frac_adult2<-table(metadata_v1_adult2$Vaccine_dose)/dim(metadata_v1_adult2)[1]
dose_df<-data.frame(Dose_frac=c(0,1,0,1,as.numeric(dose_frac_adult2)),
                    Dose=rep(c('High','Standard'),3),
                    Group=rep(c('Children','Adult1','Adult2'),each=2))
dose_df$Group<-factor(dose_df$Group,levels=c('Children','Adult1','Adult2'))
dose_df$Dose<-factor(dose_df$Dose,levels=c('Standard','High'))
p9<-ggplot(dose_df,aes(x=Group,y=Dose_frac))+geom_bar(stat='identity',aes(fill=Dose),col='black',width=.75,size=.75)+
  scale_fill_manual(values=c('grey','black'))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(y='Fraction of participants')+
  theme_classic()+theme(legend.position = 'top',axis.title.x = element_blank())

metadata_v1_children$Month_vaccinated<-factor(metadata_v1_children$Month_vaccinated,levels=c('Sep.','Oct.','Nov.','Dec.','Jan.','Feb.'))
metadata_v1_adult1$Month_vaccinated<-factor(metadata_v1_adult1$Month_vaccinated,levels=c('Sep.','Oct.','Nov.','Dec.','Jan.','Feb.'))
metadata_v1_adult2$Month_vaccinated<-factor(metadata_v1_adult2$Month_vaccinated,levels=c('Sep.','Oct.','Nov.','Dec.','Jan.','Feb.'))
p10<-ggplot(metadata_v1_children,aes(x=Month_vaccinated))+geom_bar(aes(y=(..count..)/sum(..count..)),col='black',fill='white',width=.75,size=.75)+
  labs(x='Month of vaccination',y='Fraction of participants',title='Children')+
  scale_y_continuous(limits=c(0,0.45))+
  theme_classic()

p11<-ggplot(metadata_v1_adult1,aes(x=Month_vaccinated))+geom_bar(aes(y=(..count..)/sum(..count..)),col='black',fill='white',width=.75,size=.75)+
  labs(x='Month of vaccination',y='Fraction of participants',title='Adult1')+
  scale_y_continuous(limits=c(0,0.45))+
  theme_classic()

p12<-ggplot(metadata_v1_adult2,aes(x=Month_vaccinated))+geom_bar(aes(y=(..count..)/sum(..count..)),col='black',fill='white',width=.75,size=.75)+
  labs(x='Month of vaccination',y='Fraction of participants',title='Adult2')+
  scale_y_continuous(limits=c(0,0.45))+
  theme_classic()
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(1,4))
print(p5,vp=vplayout(2,1))
print(p6,vp=vplayout(2,2))
print(p7,vp=vplayout(2,3))
print(p8,vp=vplayout(2,4))
print(p9,vp=vplayout(3,1))
print(p10,vp=vplayout(3,2))
print(p11,vp=vplayout(3,3))
print(p12,vp=vplayout(3,4))
dev.off()


## figure 2.
metadata<-as.data.frame(fread('1 - Seroconversion_prediction_data_frame_1461_entries.txt',header=T))
metadata<-metadata[metadata$PreVacc_status!='Mixed_vacc_status',]
pdf('Figure 2.pdf',width=9,height=6)
pushViewport(viewport(layout = grid.layout(2,3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
p1<-ggplot(data=metadata[metadata$Age<18,],aes(x=PreVacc_status,y=Composite_baseline))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  scale_y_continuous(limits=range(metadata$Composite_baseline))+
  labs(x='Prevacc. status',y='Composite baseline',title='Children')+
  # geom_text(x=,y=,label='P = 1.6e-9')+
  theme_classic()
p2<-ggplot(data=metadata[metadata$Age>=18 & metadata$Age<65,],aes(x=PreVacc_status,y=Composite_baseline))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  scale_y_continuous(limits=range(metadata$Composite_baseline))+
  labs(x='Prevacc. status',y='Composite baseline',title='Adult1')+
  # geom_text(x=,y=,label='P = 4.8e-6')+
  theme_classic()
p3<-ggplot(data=metadata[metadata$Age>=65,],aes(x=PreVacc_status,y=Composite_baseline))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  scale_y_continuous(limits=range(metadata$Composite_baseline))+
  labs(x='Prevacc. status',y='Composite baseline',title='Adult2')+
  # geom_text(x=,y=,label='P = 0.73')+
  theme_classic()
p4<-ggplot(data=metadata[metadata$Age<18,],aes(x=log2(Age),y=log2(BMI)))+geom_point(alpha=.5)+
  scale_y_continuous(limits=range(log2(metadata$BMI)))+
  labs(x='log2 (Age, years)',y='log2 (BMI, kg/m2)',title='Children')+
  # geom_text(x=,y=,label='R = 0.18')+
  theme_classic()
p5<-ggplot(data=metadata[metadata$Age>=18 & metadata$Age<65,],aes(x=log2(Age),y=log2(BMI)))+geom_point(alpha=.5)+
  scale_y_continuous(limits=range(log2(metadata$BMI)))+
  labs(x='log2 (Age, years)',y='log2 (BMI, kg/m2)',title='Adult1')+
  # geom_text(x=,y=,label='P = 0.37')+
  theme_classic()
p6<-ggplot(data=metadata[metadata$Age>=65,],aes(x=log2(Age),y=log2(BMI)))+geom_point(alpha=.5)+
  scale_y_continuous(limits=range(log2(metadata$BMI)))+
  labs(x='log2 (Age, years)',y='log2 (BMI, kg/m2)',title='Adult2')+
  # geom_text(x=,y=,label='P = -0.01')+
  theme_classic()
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()


## figure 3 - generated in Adobe Illustrator.


## figure 4.
residuals_SC<-as.data.frame(fread('1 - Observed_vs._predicted_seroconversion_for_1461_entries.txt',header=T))
residuals_SC<-residuals_SC[residuals_SC$PreVacc_status!='Mixed_vacc_status',]
residuals_SC_adult1<-residuals_SC[residuals_SC$Age>round(log2(17),2) & residuals_SC$Age<round(log2(65),2) & residuals_SC$Cohort=='UGA4',]
residuals_SC_adult2<-residuals_SC[residuals_SC$Age>=round(log2(65),2) & residuals_SC$Cohort=='UGA4',]
residuals_SC_children<-residuals_SC[residuals_SC$Age<=round(log2(17),2) & residuals_SC$Cohort=='UGA4',]
residuals_BL<-as.data.frame(fread('2 - Observed_vs._predicted_baseline_for_672_entries.txt',header=T))
residuals_BL<-residuals_BL[residuals_BL$Cohort=='BL3',]
residuals_BL_adult1<-residuals_BL[residuals_BL$Group=='Adult1',]
residuals_BL_adult2<-residuals_BL[residuals_BL$Group=='Adult2',]
residuals_BL_children<-residuals_BL[residuals_BL$Group=='Children',]

pdf('Figure 4.pdf',width=9,height=6)
pushViewport(viewport(layout = grid.layout(2,3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
p1<-ggplot(residuals_SC_children,aes(x=Composite_seroconversion,y=Predicted_composite))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  scale_y_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  labs(x='Observed seroconversion',y='Predicted seroconversion',title='Children')+
  coord_fixed()+theme_classic()

p2<-ggplot(residuals_SC_adult1,aes(x=Composite_seroconversion,y=Predicted_composite))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  scale_y_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  labs(x='Observed seroconversion',y='Predicted seroconversion',title='Adult1')+
  coord_fixed()+theme_classic()

p3<-ggplot(residuals_SC_adult2,aes(x=Composite_seroconversion,y=Predicted_composite))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  scale_y_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  labs(x='Observed seroconversion',y='Predicted seroconversion',title='Adult2')+
  coord_fixed()+theme_classic()

p4<-ggplot(residuals_BL_children,aes(x=Observed_BL,y=Predicted_BL))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  scale_y_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  labs(x='Observed D0 titer in the following year',y='Predicted D0 titer in the following year',title='Children')+
  geom_vline(xintercept=round(4*log2(40),1),linetype='dashed')+
  coord_fixed()+theme_classic()

p5<-ggplot(residuals_BL_adult1,aes(x=Observed_BL,y=Predicted_BL))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  scale_y_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  labs(x='Observed D0 titer in the following year',y='Predicted D0 titer in the following year',title='Adult1')+
  geom_vline(xintercept=round(4*log2(40),1),linetype='dashed')+
  coord_fixed()+theme_classic()

p6<-ggplot(residuals_BL_adult2,aes(x=Observed_BL,y=Predicted_BL))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  scale_y_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  labs(x='Observed D0 titer in the following year',y='Predicted D0 titer in the following year',title='Adult2')+
  geom_vline(xintercept=round(4*log2(40),1),linetype='dashed')+
  coord_fixed()+theme_classic()
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()


## figures 5 and 6.
RC_composite_SC<-read.table('1 - Relative_contribution_of_variables_for_SC_prediction_with_cutoffs - Composite.txt',sep='\t',header=T)
RC_strains_SC<-read.table('2 - Relative_contribution_of_variables_for_SC_prediction_with_cutoffs - Individual strains.txt',sep='\t',header=T)
rownames(RC_composite_SC)<-RC_composite_SC[,1]
rownames(RC_strains_SC)<-RC_strains_SC[,1]
RC_composite_BL<-read.table('3 - Relative_contribution_of_variables_for_BL_prediction_with_cutoffs - Composite.txt',sep='\t',header=T)
RC_strains_BL<-read.table('4 - Relative_contribution_of_variables_for_BL_prediction_with_cutoffs - Individual strains.txt',sep='\t',header=T)
rownames(RC_composite_BL)<-RC_composite_BL[,1]
rownames(RC_strains_BL)<-RC_strains_BL[,1]

pdf('Figures 5.pdf',width=4,height=8)
colors1 = c(seq(-1,1,length=300))
my_palette1 <- colorRampPalette(c("blue" ,"white", "red"))(n = 299)

p1<-heatmap.2(as.matrix(RC_composite_SC[,2:dim(RC_composite_SC)[2]]),col=my_palette1,breaks=colors1,trace='none',
              Rowv=F,Colv=F,dendrogram='none')
p3<-heatmap.2(as.matrix(RC_composite_BL[,2:dim(RC_composite_BL)[2]]),col=my_palette1,breaks=colors1,trace='none',
              Rowv=F,Colv=F,dendrogram='none')
dev.off()

pdf('Figures 6.pdf',width=6,height=8)
colors2 = c(seq(-1,1,length=1200))
my_palette2 <- colorRampPalette(c("blue" ,"white", "red"))(n = 1199)

p2<-heatmap.2(as.matrix(RC_strains_SC[,2:dim(RC_strains_SC)[2]]),col=my_palette2,breaks=colors2,trace='none',
              Rowv=F,Colv=F,dendrogram='none')
p4<-heatmap.2(as.matrix(RC_strains_BL[,2:dim(RC_strains_BL)[2]]),col=my_palette2,breaks=colors2,trace='none',
              Rowv=F,Colv=F,dendrogram='none')
dev.off()


## figure 7.
metadata<-as.data.frame(fread('1 - Seroconversion_prediction_data_frame_1461_entries.txt',header=T))
metadata<-metadata[metadata$PreVacc_status!='Mixed_vacc_status',]
metadata_cols<-metadata[,c(1:4,6:9,11:12,13,15,17,19,37,21:25)]
metadata_cols$Age<-log2(metadata_cols$Age)
metadata_cols$BMI<-log2(metadata_cols$BMI)
metadata_cols$Gender[metadata_cols$Gender=='Female']<-1
metadata_cols$Gender[metadata_cols$Gender=='Male']<-0
metadata_cols$Gender<-as.numeric(metadata_cols$Gender)
metadata_cols$Race[metadata_cols$Race=='Non-white']<-1 ## for H3N2 among children, Non_white is 0.
metadata_cols$Race[metadata_cols$Race=='White']<-0 ## for H3N2 among children, White is 1.
metadata_cols$Race<-as.numeric(metadata_cols$Race)
metadata_cols$Comorbidity[metadata_cols$Comorbidity=='No']<-1
metadata_cols$Comorbidity[metadata_cols$Comorbidity=='Yes']<-0
metadata_cols$Comorbidity<-as.numeric(metadata_cols$Comorbidity)
metadata_cols$PreVacc_status[metadata_cols$PreVacc_status=='Naive']<-1
metadata_cols$PreVacc_status[metadata_cols$PreVacc_status=='Prevaccinated']<-0
metadata_cols$PreVacc_status<-as.numeric(metadata_cols$PreVacc_status)
metadata_cols$Month_vaccinated[metadata_cols$Month_vaccinated=='Sep.']<-0
metadata_cols$Month_vaccinated[metadata_cols$Month_vaccinated=='Oct.']<-1
metadata_cols$Month_vaccinated[metadata_cols$Month_vaccinated=='Nov.']<-2
metadata_cols$Month_vaccinated[metadata_cols$Month_vaccinated=='Dec.']<-3
metadata_cols$Month_vaccinated[metadata_cols$Month_vaccinated=='Jan.']<-4
metadata_cols$Month_vaccinated[metadata_cols$Month_vaccinated=='Feb.']<-5
metadata_cols$Month_vaccinated<-as.numeric(metadata_cols$Month_vaccinated)

pdf('Figure 7.pdf',width=8,height=6)
pushViewport(viewport(layout = grid.layout(2,3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
metadata_children_cols<-metadata_cols[metadata_cols$Age<log2(18),]
inter_composite<-26.6
coef_composite_age<--2.75
coef_composite_preVacc<-5.63
coef_composite_month<-0.67
coef_composite_baseline<--0.52
metadata_children_cols$Predicted_composite<-inter_composite+coef_composite_age*metadata_children_cols$Age+coef_composite_preVacc*metadata_children_cols$PreVacc_status+coef_composite_month*metadata_children_cols$Month_vaccinated+coef_composite_baseline*metadata_children_cols$Composite_baseline
metadata_children_cols$Residual_composite<-metadata_children_cols$Composite_seroconversion-metadata_children_cols$Predicted_composite
metadata_children_cols$Bins_for_age<-cut(metadata_children_cols$Age,breaks=6)
p1<-ggplot(metadata_children_cols,aes(x=Bins_for_age,y=coef_composite_age*Age+Residual_composite))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  # geom_text(size=5,x=4,y=-1,label=paste('r = ',format(as.numeric(cor.test(metadata_children_cols$Age,coef_composite_age*metadata_children_cols$Age+metadata_children_cols$Residual_composite)$estimate),digits=2),sep=''))+
  labs(x='log2 (Age, years)',y='Seroconversion (partial residual for log2 (Age))')+
  scale_y_continuous(limits=c(-20,7.5))+
  theme_classic()
print(p1,vp=vplayout(1,1))

p4<-ggplot(metadata_children_cols,aes(x=as.factor(Month_vaccinated),y=coef_composite_month*Month_vaccinated+Residual_composite))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  # geom_text(size=5,x=4,y=11,label=paste('r = ',format(as.numeric(cor.test(metadata_children_cols$Month_vaccinated,coef_composite_month*metadata_children_cols$Month_vaccinated+metadata_children_cols$Residual_composite)$estimate),digits=2),sep=''))+
  labs(x='Month of vaccination',y='Seroconversion (partial residual for month of vaccination)')+
  scale_y_continuous(limits=c(-12,12))+
  theme_classic()
print(p4,vp=vplayout(2,1))

metadata_adult1_cols<-metadata_cols[metadata_cols$Age>=log2(18) & metadata_cols$Age<log2(65),]
inter_composite<-12.2
coef_composite_age<--0.70
coef_composite_BMI<-0.67
coef_composite_comor<-0.18
coef_composite_preVacc<-5.99
coef_composite_month<-0.32
coef_composite_baseline<--0.40
metadata_adult1_cols$Predicted_composite<-inter_composite+coef_composite_age*metadata_adult1_cols$Age+coef_composite_BMI*metadata_adult1_cols$BMI+coef_composite_comor*metadata_adult1_cols$Comorbidity+coef_composite_preVacc*metadata_adult1_cols$PreVacc_status+coef_composite_month*metadata_adult1_cols$Month_vaccinated+coef_composite_baseline*metadata_adult1_cols$Composite_baseline
metadata_adult1_cols$Residual_composite<-metadata_adult1_cols$Composite_seroconversion-metadata_adult1_cols$Predicted_composite
metadata_adult1_cols$Bins_for_age<-cut(metadata_adult1_cols$Age,breaks=6)
p2<-ggplot(metadata_adult1_cols,aes(x=Bins_for_age,y=coef_composite_age*Age+Residual_composite))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  # geom_text(size=5,x=5,y=7.5,label=paste('r = ',format(as.numeric(cor.test(metadata_adult1_cols$Age,coef_composite_age*metadata_adult1_cols$Age+metadata_adult1_cols$Residual_composite)$estimate),digits=2),sep=''))+
  labs(x='log2 (Age, years)',y='Seroconversion (partial residual for log2 (Age))')+
  scale_y_continuous(limits=c(-20,7.5))+
  theme_classic()
print(p2,vp=vplayout(1,2))

metadata_adult1_cols$Bins_for_BMI<-cut(metadata_adult1_cols$BMI,breaks=6)
p3<-ggplot(metadata_adult1_cols,aes(x=Bins_for_BMI,y=coef_composite_BMI*BMI+Residual_composite))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  # geom_text(size=5,x=5.5,y=15,label=paste('r = ',format(as.numeric(cor.test(metadata_adult1_cols$BMI,coef_composite_BMI*metadata_adult1_cols$BMI+metadata_adult1_cols$Residual_composite)$estimate),digits=2),sep=''))+
  labs(x='log2 (BMI, kg/m2)',y='Seroconversion (partial residual for log2 (BMI))')+
  theme_classic()
print(p3,vp=vplayout(1,3))

p5<-ggplot(metadata_adult1_cols,aes(x=as.factor(Month_vaccinated),y=coef_composite_month*Month_vaccinated+Residual_composite))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  # geom_text(size=5,x=1,y=12.5,label=paste('r = ',format(as.numeric(cor.test(metadata_adult1_cols$Month_vaccinated,coef_composite_month*metadata_adult1_cols$Month_vaccinated+metadata_adult1_cols$Residual_composite)$estimate),digits=2),sep=''))+
  labs(x='Month of vaccination',y='Seroconversion (partial residual for month of vaccination)')+
  scale_y_continuous(limits=c(-12,12))+
  theme_classic()
print(p5,vp=vplayout(2,2))

BL_df<-as.data.frame(fread('2 - Baseline_prediction_data_frame_672_entries.txt',header=T))
BL_df_adult1<-BL_df[BL_df$Age>=log2(18) & BL_df$Age<log2(65),]
BL_df_adult1$Month_vaccinated[BL_df_adult1$Month_vaccinated=='Sep.']<-0
BL_df_adult1$Month_vaccinated[BL_df_adult1$Month_vaccinated=='Oct.']<-1
BL_df_adult1$Month_vaccinated[BL_df_adult1$Month_vaccinated=='Nov.']<-2
BL_df_adult1$Month_vaccinated[BL_df_adult1$Month_vaccinated=='Dec.']<-3
BL_df_adult1$Month_vaccinated[BL_df_adult1$Month_vaccinated=='Jan.']<-4
BL_df_adult1$Month_vaccinated[BL_df_adult1$Month_vaccinated=='Feb.']<-5
BL_df_adult1$Month_vaccinated<-as.numeric(BL_df_adult1$Month_vaccinated)
inter_composite<-12.7
coef_composite_age<--1.65
coef_composite_month<-0.42
coef_composite_baseline<-0.62
BL_df_adult1$Predicted_composite<-inter_composite+coef_composite_age*BL_df_adult1$Age+coef_composite_month*BL_df_adult1$Month_vaccinated+coef_composite_baseline*BL_df_adult1$Composite_titer
BL_df_adult1$Residual_composite<-BL_df_adult1$Next_year_BL_composite-BL_df_adult1$Predicted_composite
p6<-ggplot(BL_df_adult1,aes(x=as.factor(Month_vaccinated),y=coef_composite_month*Month_vaccinated+Residual_composite))+stat_boxplot(geom='errorbar',width=.25,size=.75)+geom_boxplot(outlier.shape = 1,width=.6,size=.75)+
  # geom_text(size=5,x=1,y=12.5,label=paste('r = ',format(as.numeric(cor.test(BL_df_adult1$Month_vaccinated,coef_composite_month*BL_df_adult1$Month_vaccinated+BL_df_adult1$Residual_composite)$estimate),digits=2),sep=''))+
  labs(x='Month of vaccination',y='BaselineSY (partial residual for month of vaccination)')+
  theme_classic()
print(p6,vp=vplayout(2,3))
dev.off()


# figure s1.
metadata<-as.data.frame(fread('1 - Seroconversion_prediction_data_frame_1461_entries.txt',header=T))
metadata_rm_mixed_status<-metadata[metadata$PreVacc_status!='Mixed_vacc_status',]
metadata_rm_mixed_status_low_BL<-metadata_rm_mixed_status[metadata_rm_mixed_status$Baseline_category_Num_SeroPos_strains=='Low',]
pdf('Figure s1.pdf',width=8,height=5)
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
p1<-ggplot(metadata_rm_mixed_status,aes(x=PreVacc_status,y=Composite_baseline))+geom_boxplot(fill='blue',alpha=.25,width=.6)+geom_jitter(alpha=.25)+
  # scale_y_continuous(limits=c(-6,25))+
  labs(y='Composite baseline')+
  theme_classic()
p2<-ggplot(metadata_rm_mixed_status_low_BL,aes(x=PreVacc_status,y=Composite_baseline))+geom_boxplot(fill='blue',alpha=.25,width=.6)+geom_jitter(alpha=.25)+
  # scale_y_continuous(limits=c(-6,25))+
  labs(y='Composite baseline')+
  theme_classic()
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
dev.off()


## figure s2.
metadata<-as.data.frame(fread('1 - Seroconversion_prediction_data_frame_1461_entries.txt',sep='\t',header=T))
metadata<-metadata[metadata$PreVacc_status!='Mixed_vacc_status',]
metadata_subset<-metadata[,c(3:4,6:9,11:12,37)]
metadata_subset$Gender[metadata_subset$Gender=='Male']<-0
metadata_subset$Gender[metadata_subset$Gender=='Female']<-1
metadata_subset$Gender<-as.numeric(metadata_subset$Gender)

metadata_subset$Race[metadata_subset$Race=='White']<-0
metadata_subset$Race[metadata_subset$Race=='Non-white']<-1
metadata_subset$Race<-as.numeric(metadata_subset$Race)

metadata_subset$Comorbidities[metadata_subset$Comorbidities=='No']<-0
metadata_subset$Comorbidities[metadata_subset$Comorbidities=='Yes']<-1
metadata_subset$Comorbidities<-as.numeric(metadata_subset$Comorbidities)

metadata_subset$PreVacc_status[metadata_subset$PreVacc_status=='Naive']<-0
metadata_subset$PreVacc_status[metadata_subset$PreVacc_status=='Prevaccinated']<-1
metadata_subset$PreVacc_status<-as.numeric(metadata_subset$PreVacc_status)

metadata_subset$Month_vaccinated[metadata_subset$Month_vaccinated=='Sep.']<-0
metadata_subset$Month_vaccinated[metadata_subset$Month_vaccinated=='Oct.']<-1
metadata_subset$Month_vaccinated[metadata_subset$Month_vaccinated=='Nov.']<-2
metadata_subset$Month_vaccinated[metadata_subset$Month_vaccinated=='Dec.']<-3
metadata_subset$Month_vaccinated[metadata_subset$Month_vaccinated=='Jan.']<-4
metadata_subset$Month_vaccinated[metadata_subset$Month_vaccinated=='Feb.']<-5
metadata_subset$Month_vaccinated<-as.numeric(metadata_subset$Month_vaccinated)

metadata_subset$Vaccine_dose[metadata_subset$Vaccine_dose=='Standard']<-0
metadata_subset$Vaccine_dose[metadata_subset$Vaccine_dose=='High']<-1
metadata_subset$Vaccine_dose<-as.numeric(metadata_subset$Vaccine_dose)

metadata_subset$Age<-log2(metadata_subset$Age)
metadata_subset$BMI<-log2(metadata_subset$BMI)
names(metadata_subset)[1:2]<-c('log2(Age)','log2(BMI)')

pdf('Figure s2.pdf',width=8,height=8)
chart.Correlation(metadata_subset, histogram = T, pch=19)
dev.off()


## figure s3.
residuals_SC<-as.data.frame(fread('1 - Observed_vs._predicted_seroconversion_for_1461_entries.txt',header=T))
residuals_SC<-residuals_SC[residuals_SC$PreVacc_status!='Mixed_vacc_status',]
residuals_SC_adult1<-residuals_SC[residuals_SC$Age>round(log2(17),2) & residuals_SC$Age<round(log2(65),2),]
residuals_SC_adult2<-residuals_SC[residuals_SC$Age>=round(log2(65),2),]
residuals_SC_children<-residuals_SC[residuals_SC$Age<=round(log2(17),2),]
residuals_BL<-as.data.frame(fread('2 - Observed_vs._predicted_baseline_for_672_entries.txt',header=T))
residuals_BL_adult1<-residuals_BL[residuals_BL$Group=='Adult1',]
residuals_BL_adult2<-residuals_BL[residuals_BL$Group=='Adult2',]
residuals_BL_children<-residuals_BL[residuals_BL$Group=='Children',]

pdf('Figure s3.pdf',width=9,height=6)
pushViewport(viewport(layout = grid.layout(2,3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
p1<-ggplot(residuals_SC_children,aes(x=Composite_seroconversion,y=Predicted_composite))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  scale_y_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  labs(x='Observed seroconversion',y='Predicted seroconversion',title='Children')+
  coord_fixed()+theme_classic()

p2<-ggplot(residuals_SC_adult1,aes(x=Composite_seroconversion,y=Predicted_composite))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  scale_y_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  labs(x='Observed seroconversion',y='Predicted seroconversion',title='Adult1')+
  coord_fixed()+theme_classic()

p3<-ggplot(residuals_SC_adult2,aes(x=Composite_seroconversion,y=Predicted_composite))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  scale_y_continuous(limits=c(min(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite)),max(c(residuals_SC$Composite_seroconversion,residuals_SC$Predicted_composite))))+
  labs(x='Observed seroconversion',y='Predicted seroconversion',title='Adult2')+
  coord_fixed()+theme_classic()

p4<-ggplot(residuals_BL_children,aes(x=Observed_BL,y=Predicted_BL))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  scale_y_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  labs(x='Observed D0 titer in the following year',y='Predicted D0 titer in the following year',title='Children')+
  geom_vline(xintercept=round(4*log2(40),1),linetype='dashed')+
  coord_fixed()+theme_classic()

p5<-ggplot(residuals_BL_adult1,aes(x=Observed_BL,y=Predicted_BL))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  scale_y_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  labs(x='Observed D0 titer in the following year',y='Predicted D0 titer in the following year',title='Adult1')+
  geom_vline(xintercept=round(4*log2(40),1),linetype='dashed')+
  coord_fixed()+theme_classic()

p6<-ggplot(residuals_BL_adult2,aes(x=Observed_BL,y=Predicted_BL))+geom_point(alpha=.5)+
  scale_x_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  scale_y_continuous(limits=c(min(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL)),max(c(residuals_BL$Observed_BL,residuals_BL$Predicted_BL))))+
  labs(x='Observed D0 titer in the following year',y='Predicted D0 titer in the following year',title='Adult2')+
  geom_vline(xintercept=round(4*log2(40),1),linetype='dashed')+
  coord_fixed()+theme_classic()
print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))
print(p4,vp=vplayout(2,1))
print(p5,vp=vplayout(2,2))
print(p6,vp=vplayout(2,3))
dev.off()


