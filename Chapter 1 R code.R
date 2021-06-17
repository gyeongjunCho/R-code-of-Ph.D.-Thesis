# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.10")

# +++++++++++++++++++++++++++++++++++++++++
# Microbiome analysis
# +++++++++++++++++++++++++++++++++++++++++
## This code base on DADA2 Pipline tutorial (ver. 1.12)
## https://benjjneb.github.io/dada2/tutorial.html
##=====================================================================
Packages <- c("dada2", "DECIPHER", "doParallel", "ggplot2", "phyloseq", 
              "iNEXT", "dplyr", "RColorBrewer", "picante", "reshape2", 
              "ggnetwork", "network", "sna", "keras", "tidyr", "tidygraph", 
              "entropart", "gtable", "grid", "RColorBrewer", "picante", 
              "ggExtra", "dplyr", "ggrepel", "DESeq2", "DT", "gridExtra", "agricolae", 
              "progress", "ecodist","ggalt", "ggpubr", "userfriendlyscience", 
              "compiler", "ggtext", "nparcomp", "cowplot"
              )

attachRequiredPackages <-  function() {lapply(Packages, 
                                              FUN = function(Packages){
                                                do.call("require", list(Packages))
                                                })
}

attachRequiredPackages()


qRT <- read.csv("phenotype_recode/Fusarium solani - Quantification Cq Results_0.csv")
qRT
qRT$Sample <- paste(qRT$Target, qRT$Sample, sep="")
qRT$Target
Std <- qRT[which(qRT$Target==""),]
Std <- Std[which(Std$Sample !=""),]
Std

Cyc <- subset(qRT, Target != "")
Cyc$Time <- factor(Cyc$Target, levels=unique(Cyc$Target))
Cyc$Time <- gsub("A","",Cyc$Time)
Cyc$Time <- as.numeric(as.character(Cyc$Time))

Sample <- c("A1C1","A2C1","A3C1","A4C1","A5C1",
            "A6C1","A7C1","A8C1","A9C1","A10C1",
            
            "A1C3","A2C3","A3C3","A4C3","A5C3",
            "A6C3","A7C3","A8C3","A9C3","A10C3",
            
            "A1C5","A2C5","A3C5","A4C5","A5C5",
            "A6C5","A7C5","A8C5","A9C5","A10C5",
            
            "A1C7","A2C7","A3C7","A4C7","A5C7",
            "A6C7","A7C7","A8C7","A9C7","A10C7")



inputsoil <- c(0.76, 0.42, 0.21, 0.16, 0.20,
               0.35, 0.25, 0.32, 0.26, 0.41,
               
               0.37, 0.15, 0.28, 0.35, 0.29, 
               0.28, 0.63, 0.41, 0.49, 0.39,
               
               0.31, 0.25, 0.18, 0.31, 0.32, 
               0.28, 0.47, 0.36, 0.47, 0.27,
               
               0.42, 0.18, 0.16, 0.16, 0.37,
               0.32, 0.64, 0.43, 0.28, 0.29)

Cyc$inputsoil <- as.character(Cyc$Sample)

for(i in 1:40){
  Cyc$inputsoil <- gsub(Sample[i], as.character(inputsoil[i]), Cyc$inputsoil)
}

Cyc$inputsoil <- as.numeric(Cyc$inputsoil)
Cyc$inputsoil

Cyc$times_per_std_soil<- 0.5/Cyc$inputsoil
Std$Target<- c(9.8*10^5, 9.8*10^4, 9.8*10^5, 9.8*10^4, 9.8*0.75*10^5, 9.8*0.75*10^5,
               9.8*0.5*10^5, 9.8*0.5*10^5, 9.8*0.25*10^5, 9.8*0.25*10^5)



std_p <- ggplot(Std[-c(1,3),], aes(Target,Cq))+
  stat_summary(geom="point",fun.y='mean')+
  stat_summary(geom="line", fun.y='mean')+
  stat_summary(geom="errorbar", fun.data = "mean_sdl", width=0.01)+
  ylab(expression(Cycle[threshold]))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(.x)))+
  xlab(expression(paste("Standard ", (log[10](CFU/g)))))+
  geom_smooth(method='lm',se = FALSE)+
  theme_classic()
std_p

Std$Sample <- as.numeric(as.character(Std$Sample))

lm_formular <- lm(Std[-c(1,3),]$Cq ~ log10(Std[-c(1,3),]$Target))

summary(lm_formular)
a <-lm_formular$coefficients[2]
b <-lm_formular$coefficients[1]


CFU_g<- function(x){10^((x-b)/a)}

mathexpression <- bquote(
  y == -4.306893%*%log[10]~bgroup("",
                                 x,
                                 "")+49.5849
)
mathexpression2 <- bquote(R^2 == bgroup("",0.97962,""))

std_p2 <- std_p +
  annotate(geom="text", x=600000, y=28, label=deparse(mathexpression,width.cutoff = 150), parse=T,size=5)+
  annotate(geom="text", x=600000, y=27.8, label=deparse(mathexpression2,width.cutoff = 150), parse=T,size=5)
std_p2
#============

Cyc$CFU_per_soil <- CFU_g(Cyc$Cq)*Cyc$times_per_std_soil

F.solani_CFU_per_soil <- aggregate(CFU_per_soil~ Sample, Cyc, mean)
F.solani_CFU_per_soil <- F.solani_CFU_per_soil[c(5:40,1:4),]
rownames(F.solani_CFU_per_soil) <- NULL
write.csv(F.solani_CFU_per_soil, "phenotype_recode/Fusarium density.csv")

Cyc$pot<- gsub("A[0-9]{1,}", "", Cyc$Sample)
Cyc$pot
Cyc$Time <- factor(as.character(Cyc$Time), levels = as.character(1:10))

for(i in 1:10){
  print(shapiro.test((Cyc$CFU_per_soil[which(Cyc$Time == i)])))
} #some sample Fusarium density are not normality.

Cyc_fusarium_kruskal <- kruskal.test(CFU_per_soil ~ Time, Cyc)
Cyc_fusarium_kruskal$p.value

Cyc_fusarium_post_hoc <- nparcomp(CFU_per_soil ~ Time, Cyc)
Cyc_fusarium_post_hoc_p <- Cyc_fusarium_post_hoc$Analysis
Cyc_fusarium_post_hoc_p <- Cyc_fusarium_post_hoc_p[which(Cyc_fusarium_post_hoc_p$p.Value < 0.05),]
Cyc_fusarium_post_hoc_p$compare1_cycling <- c(1,2,2,2,3,7)
Cyc_fusarium_post_hoc_p$compare2_cycling <- c(10,8,9,10,10,10)

Cyc_fusarium_post_hoc_p$P_char <- ifelse(Cyc_fusarium_post_hoc_p$p.Value < 0.001, "***",
                                         ifelse(Cyc_fusarium_post_hoc_p$p.Value < 0.01, "**","*"))

Cyc_fusarium_post_hoc_p$loc <- 2^(1:nrow(Cyc_fusarium_post_hoc_p))*10^7.2
Cyc_fusarium_post_hoc_p
Cyc_fusarium_kruskal$p.value

Fusarium <- ggplot(Cyc, aes(Time, CFU_per_soil))+
  stat_summary(geom="bar", position=position_dodge(0.2), fun.data = "mean_se", alpha=0.5, color="black", fill="white")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.3, position=position_dodge(0.2), alpha=0.5)+
  stat_summary(geom="point", aes(color=pot), shape=1, width=0.3, position=position_jitterdodge(0.2))+
  xlab("Cycling")+
  ylab(expression(paste(italic(F.)," ",italic(solani)," CFU / Rhizosphere 1g")))+
  coord_cartesian(xlim=c(0.5,10))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand=c(0,0), limit=c(1,1e11))+
  theme_classic()+
  labs(color="pot")+
  annotation_logticks(sides="l")+
  annotate("text", 
           x=1.5,
           y=1e10,
           label=expression(italic(P) ==  7.106%*%10^-5))+
  annotate("text", 
           x=Cyc_fusarium_post_hoc_p$compare1_cycling/2 + Cyc_fusarium_post_hoc_p$compare2_cycling/2,
           y=Cyc_fusarium_post_hoc_p$loc,
           label=Cyc_fusarium_post_hoc_p$P_char)+
  annotate("segment", 
           x=Cyc_fusarium_post_hoc_p$compare1_cycling, 
           xend=Cyc_fusarium_post_hoc_p$compare2_cycling,
           y=Cyc_fusarium_post_hoc_p$loc,
           yend=Cyc_fusarium_post_hoc_p$loc)+
  annotate("segment", 
           x=Cyc_fusarium_post_hoc_p$compare1_cycling, 
           xend=Cyc_fusarium_post_hoc_p$compare1_cycling,
           y=Cyc_fusarium_post_hoc_p$loc,
           yend=Cyc_fusarium_post_hoc_p$loc/1.3)+
  annotate("segment", 
           x=Cyc_fusarium_post_hoc_p$compare2_cycling, 
           xend=Cyc_fusarium_post_hoc_p$compare2_cycling,
           y=Cyc_fusarium_post_hoc_p$loc,
           yend=Cyc_fusarium_post_hoc_p$loc/1.3)


Fusarium

#============================================================================================

DI <- read.csv("phenotype_recode/ginseng_DI.csv")
DI
DI$cycle_ID <- paste("A", DI$A.Cycle.number, sep="")
DI$cycle_pot_ID <- paste(DI$cycle_ID, DI$potID, sep="")

Cyc_Fusarium_mean <-aggregate(CFU_per_soil ~ Sample, Cyc, mean)
colnames(Cyc_Fusarium_mean)[1] <- "cycle_pot_ID"

shoot <- read.csv("phenotype_recode/shoot length.csv")

shoot_mean <-aggregate(length ~ ID, shoot, mean)
shoot_mean
shoot_sd <-aggregate(length ~ ID, shoot, sd)

colnames(shoot_mean) <- c("cycle_pot_ID", "length_mean")
colnames(shoot_sd) <- c("cycle_pot_ID", "length_sd")

DI <- merge(DI, shoot_mean, by="cycle_pot_ID",all.x=T)

for(i in 1:10){
  print(shapiro.test((DI$pot_DI[which(DI$ A.Cycle.number== i)])))
} #some sample DI are not normality.

Cyc_DI_kruskal <- kruskal.test(pot_DI ~ A.Cycle.number, DI)
Cyc_DI_kruskal

Cyc_DI_post_hoc <- nparcomp(pot_DI ~ A.Cycle.number, DI, correlation = FALSE)

Cyc_DI_post_hoc_p <- Cyc_DI_post_hoc$Analysis
Cyc_DI_post_hoc_p$p.Value < 0.05
Cyc_DI_post_hoc_p <- Cyc_DI_post_hoc_p[which(Cyc_DI_post_hoc_p$p.Value < 0.05),]
Cyc_DI_post_hoc_p$compare1_cycling <- c(1, 3)
Cyc_DI_post_hoc_p$compare2_cycling <- c(5, 5)

Cyc_DI_post_hoc_p$P_char <- ifelse(Cyc_DI_post_hoc_p$p.Value < 0.001, "***",
                                         ifelse(Cyc_DI_post_hoc_p$p.Value < 0.01, "**","*"))

Cyc_DI_post_hoc_p
Cyc_DI_post_hoc_p$loc <- 1:nrow(Cyc_DI_post_hoc_p)*6+60

DI_plot <- ggplot(DI,aes(A.Cycle.number,pot_DI))+
  stat_summary(geom="point")+
  stat_summary(geom="errorbar", fun.data="mean_se",width=0.2)+
  stat_summary(geom="line", fun.data="mean_se")+
  geom_point(shape=1, position= position_dodge(0.2), aes(color=potID))+
  geom_line(position= position_dodge(0.2),aes(color=potID), alpha=0.2)+
  scale_x_continuous(breaks=1:10)+
  scale_y_continuous(labels = c("significantly\nno diffrence",0:3*25), breaks=c(-7,0:3*25))+
  labs(y="Disease index")+
  labs(x="Cycling")+  annotate("text", 
                               x=8,
                               y=70,
                               label=expression(italic(P) ==  5.764%*%10^-5))+
  
  labs(color="pot")+
  annotate("text", 
           x=Cyc_DI_post_hoc_p$compare1_cycling/2 + Cyc_DI_post_hoc_p$compare2_cycling/2,
           y=Cyc_DI_post_hoc_p$loc,
           label=Cyc_DI_post_hoc_p$P_char)+
  annotate("segment", 
           x=Cyc_DI_post_hoc_p$compare1_cycling, 
           xend=Cyc_DI_post_hoc_p$compare2_cycling,
           y=Cyc_DI_post_hoc_p$loc,
           yend=Cyc_DI_post_hoc_p$loc)+
  annotate("segment", 
           x=Cyc_DI_post_hoc_p$compare1_cycling, 
           xend=Cyc_DI_post_hoc_p$compare1_cycling,
           y=Cyc_DI_post_hoc_p$loc,
           yend=Cyc_DI_post_hoc_p$loc-1.5)+
  annotate("segment", 
           x=Cyc_DI_post_hoc_p$compare2_cycling, 
           xend=Cyc_DI_post_hoc_p$compare2_cycling,
           y=Cyc_DI_post_hoc_p$loc,
           yend=Cyc_DI_post_hoc_p$loc-1.5)+
  theme_classic()

DI_plot


for(i in 1:10){
  print(shapiro.test((DI$germination.ratio[which(DI$ A.Cycle.number== i)])))
} #some seedling are not normality.

Cyc_seedling_kruskal <- kruskal.test(germination.ratio ~ A.Cycle.number, DI)
Cyc_seedling_kruskal

Cyc_seedling_post_hoc <-  nparcomp(germination.ratio ~ A.Cycle.number, DI)
Cyc_seedling_post_hoc_p <- Cyc_seedling_post_hoc$Analysis
Cyc_seedling_post_hoc_p <- Cyc_seedling_post_hoc_p[which(Cyc_seedling_post_hoc_p$p.Value < 0.05),]

Cyc_seedling_post_hoc_p$compare1_cycling <- as.numeric(
  gsub(" , [0-9]{1,} [)]","",
       gsub("p[(] ","",Cyc_seedling_post_hoc_p$Comparison)
  )
)

Cyc_seedling_post_hoc_p$compare2_cycling <- as.numeric(
  gsub(" [)]","",
       gsub("p[(] [0-9]{1,} , ","",Cyc_seedling_post_hoc_p$Comparison)
  )
)

Cyc_seedling_post_hoc_p$P_char <- ifelse(Cyc_seedling_post_hoc_p$p.Value < 0.001, "***",
                                   ifelse(Cyc_seedling_post_hoc_p$p.Value < 0.01, "**","*"))

Cyc_seedling_post_hoc_p$loc <- 1:nrow(Cyc_seedling_post_hoc_p)*0.04+1
Cyc_seedling_post_hoc_p

budding_ratio <- ggplot(DI,aes(A.Cycle.number, germination.ratio))+
  stat_summary(geom="point")+
  stat_summary(geom="errorbar", fun.data="mean_se",width=0.2)+
  stat_summary(geom="line", fun.data="mean_se")+
  geom_point(shape=1, position = position_dodge(0.2), aes(color=potID))+
  geom_line(position= position_dodge(0.2),aes(color=potID), alpha=0.1)+
  scale_x_continuous(breaks=1:10)+
  scale_y_continuous(labels = c("significantly\nno diffrence", scales::percent(0:5/5)), breaks = c(-0.15,(0:5)/5))+
  labs(x="Cycling")+
  labs(y="Sprouting ratio")+
  labs(color="pot")+
  theme_classic()+
  annotate("text", 
           x=Cyc_seedling_post_hoc_p$compare1_cycling/2 + Cyc_seedling_post_hoc_p$compare2_cycling/2,
           y=Cyc_seedling_post_hoc_p$loc,
           label=Cyc_seedling_post_hoc_p$P_char)+
  annotate("segment", 
           x=Cyc_seedling_post_hoc_p$compare1_cycling, 
           xend=Cyc_seedling_post_hoc_p$compare2_cycling,
           y=Cyc_seedling_post_hoc_p$loc,
           yend=Cyc_seedling_post_hoc_p$loc)+
  annotate("segment", 
           x=Cyc_seedling_post_hoc_p$compare1_cycling, 
           xend=Cyc_seedling_post_hoc_p$compare1_cycling,
           y=Cyc_seedling_post_hoc_p$loc,
           yend=Cyc_seedling_post_hoc_p$loc-0.01)+
  annotate("segment", 
           x=Cyc_seedling_post_hoc_p$compare2_cycling, 
           xend=Cyc_seedling_post_hoc_p$compare2_cycling,
           y=Cyc_seedling_post_hoc_p$loc,
           yend=Cyc_seedling_post_hoc_p$loc-0.01)+
  labs(x="Cycling")+  annotate("text", 
                               x=2,
                               y=2,
                               label=expression(italic(P) ==  9.262%*%10^-10))+
  theme_classic()

budding_ratio




for(i in c(1:3,6:10)){
  print(shapiro.test((DI$length_mean[which(DI$ A.Cycle.number== i)])))
} #all shoot lengths have normality without 4~5 cycling.

bartlett.test(length_mean ~ A.Cycle.number, subset(DI, A.Cycle.number%in%c(1:3,6:10))) # equal var

shoot_length_aov <- aov(length_mean ~ as.character(A.Cycle.number), subset(DI, A.Cycle.number%in%c(1:3,6:10)))

shoot_length_aov

summary_shoot_length_aov <- summary(shoot_length_aov)
unlist(summary_shoot_length_aov)
shoot_length_aov_post_hoc <- TukeyHSD(shoot_length_aov)
shoot_length_aov_post_hoc
shoot_length_aov_post_hoc <- shoot_length_aov_post_hoc$`as.character(A.Cycle.number)`
shoot_length_aov_post_hoc <- as.data.frame(shoot_length_aov_post_hoc)
shoot_length_aov_post_hoc <- subset(shoot_length_aov_post_hoc, `p adj` < 0.05)

shoot_length_aov_post_hoc$P_char <- ifelse(shoot_length_aov_post_hoc$`p adj` < 0.001, "***",
                                         ifelse(shoot_length_aov_post_hoc$`p adj` < 0.01, "**","*"))

shoot_length_aov_post_hoc$compare1_cycling <- 
  t(
    as.data.frame(
      strsplit(
        rownames(shoot_length_aov_post_hoc),
        "-"
      )
    )
  )[,1] %>% as.numeric()

shoot_length_aov_post_hoc$compare2_cycling <- 
  t(
    as.data.frame(
      strsplit(
        rownames(shoot_length_aov_post_hoc),
        "-"
      )
    )
  )[,2] %>% as.numeric()


shoot_length_aov_post_hoc <- shoot_length_aov_post_hoc[c(1:3,5:12,4),]
shoot_length_aov_post_hoc$loc <- 1:nrow(shoot_length_aov_post_hoc)/5+6

shoot_length <- ggplot(DI,aes(A.Cycle.number, length_mean))+
  stat_summary(geom="point")+
  stat_summary(geom="errorbar", fun.data="mean_se",width=0.2)+
  stat_summary(data=subset(DI, A.Cycle.number %in% 6:10),geom="line", fun.data="mean_se")+
  stat_summary(data=subset(DI, A.Cycle.number %in% 1:3),geom="line", fun.data="mean_se")+
  geom_point(shape=1, position= position_dodge(0.2), aes(color=potID))+
  geom_line(position= position_dodge(0.2),aes(color=potID), alpha=0.2)+
  scale_x_continuous(breaks=1:10)+
  labs(x="Cycling")+
  labs(y="Shoot length (cm)")+
  labs(color="pot")+
  theme_classic()+
  labs(x="Cycling")+  annotate("text", 
                               x=2,
                               y=8,
                               label=expression(italic(P) ==  4.979%*%10^-12))+
  annotate("text", 
           x=shoot_length_aov_post_hoc$compare1_cycling/2 + shoot_length_aov_post_hoc$compare2_cycling/2,
           y=shoot_length_aov_post_hoc$loc,
           label=shoot_length_aov_post_hoc$P_char)+
  annotate("segment", 
           x=shoot_length_aov_post_hoc$compare1_cycling, 
           xend=shoot_length_aov_post_hoc$compare2_cycling,
           y=shoot_length_aov_post_hoc$loc,
           yend=shoot_length_aov_post_hoc$loc)+
  annotate("segment", 
           x=shoot_length_aov_post_hoc$compare1_cycling, 
           xend=shoot_length_aov_post_hoc$compare1_cycling,
           y=shoot_length_aov_post_hoc$loc,
           yend=shoot_length_aov_post_hoc$loc-0.1)+
  annotate("segment", 
           x=shoot_length_aov_post_hoc$compare2_cycling, 
           xend=shoot_length_aov_post_hoc$compare2_cycling,
           y=shoot_length_aov_post_hoc$loc,
           yend=shoot_length_aov_post_hoc$loc-0.1)

shoot_length

ggarrange(DI_plot,  budding_ratio, shoot_length, Fusarium, labels=c("A","B","C","D"))

Fusarium_correlation <- merge(DI, Cyc_Fusarium_mean)
Fusarium_correlation$log10_CFU_per_soil <- log10(Fusarium_correlation$CFU_per_soil )

correlation(Fusarium_correlation$log10_CFU_per_soil, Fusarium_correlation$pot_DI)

F_vs_DI <- ggplot(Fusarium_correlation, aes(x=log10_CFU_per_soil,  y=pot_DI, color=potID, label = cycle_pot_ID))+
  geom_point()+
  geom_text(vjust=-0.5, show.legend = FALSE)+
  stat_smooth(method=lm, aes(group=NULL, color=NULL), level = 0, show.legend = FALSE)+
  annotate(geom = "text", x=5.5, y = 60, label=expression(paste(italic(rho)==-0.037, ",  ", italic(P)==0.971)), label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Disease index")+
  xlab(expression(log[10](italic('F. solani')~textstyle(CFU/soil(g)))))+
  scale_x_continuous(breaks = 1:10)+
  guides(color = guide_legend(override.aes = list(size=3.5)))+
  theme(panel.background = element_rect(fill="white"),
        axis.line.x.bottom = element_line(color="black"),
        axis.line.y.left = element_line(color="black"),
        legend.key = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))


F_vs_DI


correlation(Fusarium_correlation$log10_CFU_per_soil, Fusarium_correlation$germination.ratio*100)

F_vs_budding <- ggplot(Fusarium_correlation, aes(x=log10_CFU_per_soil,  y=germination.ratio*100, color=potID, label = cycle_pot_ID))+
  geom_point()+
  geom_text(vjust=-0.5, show.legend = FALSE)+
  stat_smooth(method=lm, aes(group=NULL, color=NULL), level = 0, show.legend = FALSE)+
  annotate(geom = "text", x=5.5, y = 60, label=expression(paste(italic(rho)==0.01, "0,  ", italic(P)==0.953)), label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Sprouting ratio (%)")+
  xlab(expression(log[10](italic('F. solani')~textstyle(CFU/soil(g)))))+
  scale_x_continuous(breaks = 1:10)+
  guides(color = guide_legend(override.aes = list(size=3.5)))+
  theme(panel.background = element_rect(fill="white"),
        axis.line.x.bottom = element_line(color="black"),
        axis.line.y.left = element_line(color="black"),
        legend.key = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))

F_vs_budding

correlation(Fusarium_correlation$length_mean, Fusarium_correlation$log10_CFU_per_soil)

F_vs_length <- ggplot(Fusarium_correlation, aes(x=log10_CFU_per_soil,  y=length_mean, color=potID, label = cycle_pot_ID))+
  geom_point()+
  geom_text(vjust=-0.5, show.legend = FALSE)+
  stat_smooth(method=lm, aes(group=NULL, color=NULL), level = 0, show.legend = FALSE)+
  annotate(geom = "text", x=5.5, y = 5, label=expression(paste(italic(rho)==-0.243, ",  ", italic(P)==0.187)), label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Shoot length (cm)")+
  xlab(expression(log[10](italic('F. solani')~textstyle(CFU/soil(g)))))+
  scale_x_continuous(breaks = 1:10)+
  guides(color = guide_legend(override.aes = list(size=3.5)))+
  theme(panel.background = element_rect(fill="white"),
        axis.line.x.bottom = element_line(color="black"),
        axis.line.y.left = element_line(color="black"),
        legend.key = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))

F_vs_length

correlation(Fusarium_correlation$A.Cycle.number, Fusarium_correlation$log10_CFU_per_soil)

F_vs_cycling <- ggplot(Fusarium_correlation, aes(x=A.Cycle.number,  y=log10_CFU_per_soil, color=potID, label = cycle_pot_ID))+
  geom_point()+
  geom_text(vjust=-0.5, show.legend = FALSE)+
  stat_smooth(method=lm, aes(group=NULL, color=NULL), level = 0, show.legend = FALSE)+
  annotate(geom = "text", x=8, y = 7.5, label=expression(paste(italic(rho)==-0.398, ",  ", italic(P)==0.011)), label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  xlab("Cycling")+
  ylab(expression(log[10](italic('F. solani')~textstyle(CFU/soil(g)))))+
  scale_x_continuous(breaks = 1:10)+
  guides(color = guide_legend(override.aes = list(size=3.5)))+
  theme(panel.background = element_rect(fill="white"),
        axis.line.x.bottom = element_line(color="black"),
        axis.line.y.left = element_line(color="black"),
        legend.key = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12))


F_vs_cycling

ggarrange(F_vs_DI, F_vs_budding, F_vs_length, F_vs_cycling, labels = c("A", "B", "C", "D"), ncol = 2, nrow=2)

#===============
# DADA2 analysis
#===============

# Input raw files directory
path <- "raw_fastq_files"

# Forward reads list
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))

# Reverse reads list
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extraction, sample, names, from, file, list
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample_ord <- c("EA1C1", "EA1C3", "EA1C5", "EA1C7", "EA2C1", 
                "EA2C3", "EA2C5", "EA2C7", "EA3C1", "EA3C3", 
                "EA3C5", "EA3C7", "EA4C1", "EA4C3", "EA4C5",
                "EA4C7", "EA5C1", "EA5C3", "EA5C5", "EA5C7",
                "EA6C1", "EA6C3", "EA6C5", "EA6C7", "EA7C1",
                "EA7C3", "EA7C5", "EA7C7", "EA8C1", "EA8C3",
                "EA8C5", "EA8C7", "EA9C1", "EA9C3", "EA9C5",
                "EA9C7", "EA10C1", "EA10C3", "EA10C5", "EA10C7", 
                "RA1C1", "RA1C3", "RA1C5", "RA1C7", "RA2C1", 
                "RA2C3", "RA2C5", "RA2C7", "RA3C1", "RA3C3", 
                "RA3C5", "RA3C7", "RA4C1", "RA4C3", "RA4C5", 
                "RA4C7", "RA5C1", "RA5C3", "RA5C5", "RA5C7",
                "RA6C1", "RA6C3", "RA6C5", "RA6C7", "RA7C1", 
                "RA7C3", "RA7C5", "RA7C7", "RA8C1", "RA8C3",
                "RA8C5", "RA8C7", "RA9C1", "RA9C3", "RA9C5",
                "RA9C7", "RA10C1", "RA10C3", "RA10C5", "RA10C7")

# Visualization raw forward quality
raw_foward_Q <- plotQualityProfile(fnFs)
raw_foward_Q

# Visualization raw reverse quality
raw_reverse_Q <- plotQualityProfile(fnRs)
raw_reverse_Q

# Extract raw files list
list.files(path) # files in filtered/ subdirectory
filt_path <- "filtered_fastq_files"
filtFs <- file.path(filt_path , paste0(sample.names, "_F_filt.fastq.gz")) # set filtered forward path info
filtFs

filtRs <- file.path(filt_path , paste0(sample.names, "_R_filt.fastq.gz")) # set filtered reverse path info
filtFs

# Set sample names for filtered
names(filtFs) <- sample.names 
names(filtRs) <- sample.names  

filterAndTrim(fwd = fnFs, # forward reads
              filt = filtFs, # filtered forward reads save info
              rev = fnRs, # reverse reads 
              filt.rev = filtRs, # filtered reverse reads save info
              truncLen=c(270,200), # trim to qualify score 30 or higher
              trimLeft = c(nchar("GTGYCAGCMGCCGCGGTAA"),
                           nchar("GGACTACNVGGGTWTCTAAT")), # trim the 515F, 806R primer
              maxN = 0, maxEE = c(1,1), rm.phix = TRUE, # Other statistical quality check criteria
              n=1e8, # sampling reads number
              compress = TRUE, # compress
              multithread = FALSE) # On Windows set multithread=FALSE; An error has been reported when the operation has an unacceptable amount of RAM. In this case, "FALSE" is recommended.


# Visualization filtered forward quality
filt_foward_Q <- plotQualityProfile(filtFs[sample_ord])
filt_foward_Q

# Visualization filtered reverse quality
filt_reverse_Q <- plotQualityProfile(filtRs[sample_ord])
filt_reverse_Q

errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+11)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e+11)
errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errF, nominalQ=TRUE)
errF_plot
errR_plot

# Sample Inference using Divisive Amplicon Denoising Algorithm (DADA)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 200:300]

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


# Assign taxonomy with IDTAXA
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet
load("IDTAXA SILVA DB/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
rm(trainingSet) # remove SILVA DB on ram for ram management

# No-bacterial OTUs removing
taxid_OTUtable <- cbind(data.frame(taxid), as.data.frame(t(seqtab.nochim)[,sample_ord]))
taxid_OTUtable <- subset(taxid_OTUtable, domain == "Bacteria")
taxid_OTUtable <- subset(taxid_OTUtable, order != "Chloroplast")
taxid_OTUtable <- subset(taxid_OTUtable, family != "Mitochondria")

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))

inputF <- raw_foward_Q$layers[[6]]$data$rc
names(inputF) <- raw_foward_Q$layers[[6]]$data$label

inputR <- raw_reverse_Q$layers[[6]]$data$rc
names(inputR) <- raw_reverse_Q$layers[[6]]$data$label

filteredF <- filt_foward_Q$layers[[6]]$data$rc
names(filteredF) <- filt_foward_Q$layers[[6]]$data$label

filteredR <- filt_reverse_Q$layers[[6]]$data$rc
names(filteredR) <- filt_reverse_Q$layers[[6]]$data$label

inoutput <- data.frame(inputF,inputR,filteredF,filteredR)

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
track <- track[sample_ord,]
tax_Bacteria <- colSums(taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
track <- cbind(inputF, inputR, filteredF, filteredR, track, tax_Bacteria)
rownames(track) <- sample_ord

write.csv(track, "metagenome_analysis_results/DADA2_results/Reads_track.csv")
#===========================
# DADA2 results output files
#===========================

registerDoParallel(cores=64) # CPU threads to multithreads

# OTUs ID generation
OTUsID <- foreach(i = 1:nrow(taxid_OTUtable), .combine="rbind")%dopar%{
  OTUsN <- paste("OTU", i, sep="_")
  print(OTUsN)
}
OTUsID <- as.vector(OTUsID)  
OTUsID 

# OTU fasta files
fasta <- foreach(i= 1:nrow(taxid_OTUtable), .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], rownames(taxid_OTUtable)[i])
  print(match)
} 
fasta <- as.vector(fasta)
fasta <- as.data.frame(gsub("^OTU_", ">OTU_",fasta))
colnames(fasta) <- NA
write.table(fasta, file = "metagenome_analysis_results/DADA2_results/OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

endosphere_fasta <- foreach(i = which(rowSums(OTU_counts_table[,2:41]) != 0), .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], rownames(taxid_OTUtable)[i])
  print(match)
}
endosphere_fasta <- as.vector(endosphere_fasta)
endosphere_fasta <- as.data.frame(gsub("^OTU_", ">OTU_", endosphere_fasta))
write.table(endosphere_fasta, file = "metagenome_analysis_results/DADA2_results/Cycling/endosphere_OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

rhizosphere_fasta <- foreach(i = which(rowSums(OTU_counts_table[,42:81]) != 0), .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], rownames(taxid_OTUtable)[i])
  print(match)
}
rhizosphere_fasta <- as.vector(rhizosphere_fasta)
rhizosphere_fasta <- as.data.frame(gsub("^OTU_", ">OTU_", rhizosphere_fasta))
write.table(rhizosphere_fasta, file = "metagenome_analysis_results/DADA2_results/Cycling/rhizosphere_OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

# OTU counts table
OTU_counts_table <- cbind(OTUsID, taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
dim(OTU_counts_table)
row.names(OTU_counts_table) <- NULL
colnames(OTU_counts_table)[1] <- "#OTU ID"
write.table(OTU_counts_table, file = "metagenome_analysis_results/DADA2_results/OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(OTU_counts_table[,1:81], file = "metagenome_analysis_results/DADA2_results/Cycling/OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(OTU_counts_table[which(rowSums(OTU_counts_table[,2:41]) != 0),1:41], file = "metagenome_analysis_results/DADA2_results/Cycling/endosphere_OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(OTU_counts_table[which(rowSums(OTU_counts_table[,42:81]) != 0),c(1,42:81)], file = "metagenome_analysis_results/DADA2_results/Cycling/rhizosphere_OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# OTU tax table
OTU_tax_table <- cbind(OTUsID, taxid_OTUtable[,1:7])
row.names(OTU_tax_table) <- NULL
colnames(OTU_tax_table)[1] <- "#OTU ID"
write.table(OTU_tax_table, file = "metagenome_analysis_results/DADA2_results/OTU_tax_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# sequence tax counts table
OTU_seq_tax_counts_table <- cbind(as.vector(row.names(taxid_OTUtable)), OTUsID, taxid_OTUtable)
row.names(OTU_seq_tax_counts_table) <- NULL
colnames(OTU_seq_tax_counts_table)[1] <- "sequence" 
write.table(OTU_seq_tax_counts_table, file = "metagenome_analysis_results/DADA2_results/OTU_seq_tax_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ++++++++++++++++++++++++++++++++++
# Microbiome transform data anlaysis
# ++++++++++++++++++++++++++++++++++

# OTU
OTU <- otu_table(
  t(
    read.csv("metagenome_analysis_results/DADA2_results/OTU_counts_table.tsv", sep = "\t",row.names = 1)
  ),
  taxa_are_rows = FALSE
)

# Taxa
TAX <- tax_table(
  as.matrix(
    read.csv("metagenome_analysis_results/DADA2_results/OTU_tax_table.tsv", sep = "\t",row.names = 1)
  )
)

# Sample inforamtion
SAM <- sample_data(
  read.csv("phenotype_recode/sample_information.tsv", sep="\t", row.names = 1)
)
SAM
# To make phyloseq obj
# reads base obj
phy_obj_reads <- phyloseq(OTU,TAX,SAM)

# to modified relative OTU phyloseq obj
phy_obj_relative <- phy_obj_reads
phy_obj_relative@otu_table <- phy_obj_reads@otu_table/rowSums(phy_obj_reads@otu_table)
write.csv(t(phy_obj_relative@otu_table[,paste0("OTU_",1:ncol(phy_obj_relative@otu_table))]), "metagenome_analysis_results/DADA2_results/relative_abundance_by_OTU.csv")

#---------
# Coverage 
#----------
OTU_long_table <- data.frame(t(OTU))

# Chao's coverage
chao_coverage <-  unlist(
  foreach(i = 1:dim(OTU_long_table)[2])%dopar%{
    Coverage(OTU_long_table[,i], Estimator = "Chao")
  }
)
names(chao_coverage) <- colnames(OTU_long_table)

# Good's coverage
good_coverage <- unlist(
  foreach(i = 1:dim(OTU_long_table)[2])%dopar%{
    1 - (sum(OTU_long_table[,i] == 1) / sum(OTU_long_table[,i]))
  }
)

# Coverage summary
coverage <- data.frame(good_coverage, chao_coverage)*100
colnames(coverage) <- c("Good's coverage (%)", "Chao's coverage (%)")
coverage
write.csv(coverage, "metagenome_analysis_results/coverage.csv")

#------------------
# Alphyadiversity
#----------------
# Calcuation
alpha_diversity <- plot_richness(phy_obj_reads)$data
alpha_diversity

# visualization
cycling_alpha_diversity <- ggplot(data = subset(alpha_diversity, variable=="Observed"|variable=="Shannon"|variable=="Simpson"), mapping = aes(x = Cycling, y=value))+
  stat_summary(geom = "point", fun.data = "mean_se")+
  stat_summary(geom = "errorbar",size=0.2, fun.data = "mean_se")+
  geom_point(alpha=0.5, aes(shape=PotID, color=PotID))+
  theme_light()+
  scale_x_continuous(breaks=1:10)+
  scale_shape_manual(values = c(1,2,0,3))+
  facet_wrap(Location ~ variable, scales = "free_y", nrow=2)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  ylab("Value")+
  theme(legend.position = "top",
        panel.grid = element_blank())

cycling_alpha_diversity

alpha_diversity$samples2 <- gsub("^[ER]", "", alpha_diversity$samples)

alpha_diversity

endosphere_alpha <- subset(alpha_diversity, Location=="Endosphere")
endosphere_alpha

correlation(subset(endosphere_alpha, variable=="Observed")$value, subset(endosphere_alpha, variable=="Observed")$Pot_DI)

EObserved_DI<- ggscatter(subset(endosphere_alpha, variable=="Observed"), x="Pot_DI", 
                         y="value", color="PotID", label = "samples2",
                         add = "reg.line",
                         add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  annotate(geom = "richtext", x=50, y = 100, label="*&rho;* = 0.066,  *P* = 0.685", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Observed")+
  xlab("Disease index")+
  facet_wrap(Location~variable)+
  guides(colour=guide_legend(override.aes=list(size=4)))+
  theme(strip.background = element_rect(color="white"), legend.position = "none")

EObserved_DI

correlation(subset(endosphere_alpha, variable=="Shannon")$value, subset(endosphere_alpha, variable=="Shannon")$Pot_DI)

EShannon_DI <- ggscatter(subset(endosphere_alpha, variable=="Shannon"), x="Pot_DI", 
                         y="value", color="PotID", label = "samples2",
                         add = "reg.line",
                         add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  annotate(geom = "richtext", x=50, y = 3, label="*&rho;* = 0.104, *P* = 0.522", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Simpson")+
  xlab("Disease index")+
  facet_wrap(Location~variable)+
  guides(colour=guide_legend(override.aes=list(size=4)))+
  theme(strip.background = element_rect(color="white"), legend.position = "none")

EShannon_DI

correlation(subset(endosphere_alpha, variable=="Simpson")$value, subset(endosphere_alpha, variable=="Simpson")$Pot_DI)

ESimpson_DI <- ggscatter(subset(endosphere_alpha, variable=="Simpson"), x="Pot_DI", 
                         y="value", color="PotID", label = "samples2",
                         add = "reg.line",
                         add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  annotate(geom = "richtext", x=50, y = 0.95, label="*&rho;* = 0.173, *P* = 0.286", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Simpson")+
  xlab("Root rot disease index")+
  facet_wrap(Location~variable)+
  guides(colour=guide_legend(override.aes=list(size=4)))+
  theme(strip.background = element_rect(color="white"), legend.position = "none")

ESimpson_DI


rhizosphere_alpha <- subset(alpha_diversity, Location=="Rhizosphere")
rhizosphere_alpha

correlation(subset(rhizosphere_alpha, variable=="Observed")$value, subset(rhizosphere_alpha, variable=="Observed")$Pot_DI)

RObserved_DI<- ggscatter(subset(rhizosphere_alpha, variable=="Observed"), x="Pot_DI", 
          y="value", color="PotID", label = "samples2",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  #stat_cor(method = "pearson", label.x = 3, label.y =350)+
  annotate(geom = "richtext", x=50, y = 400, label="*&rho;* = -0.570, *P* = 1.252 &#10005; 10<sup>-4</sup>", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Observed")+
  xlab("Root rot disease index")+
  facet_wrap(Location~variable)+
  guides(colour=guide_legend(override.aes=list(size=4)))+
  theme(strip.background = element_rect(color="white"), legend.position = "none")

RObserved_DI

correlation(subset(rhizosphere_alpha, variable=="Shannon")$value, subset(rhizosphere_alpha, variable=="Shannon")$Pot_DI)

RShannon_DI <- ggscatter(subset(rhizosphere_alpha, variable=="Shannon"), x="Pot_DI", 
          y="value", color="PotID", label = "samples2",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  #stat_cor(method = "pearson", label.x = 3, label.y =4)+
  annotate(geom = "richtext", x=50, y = 4, label="*&rho;* = -0.547, *P* = 2.624 &#10005; 10<sup>-4</sup>", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Simpson")+
  xlab("Root rot disease index")+
  facet_wrap(Location~variable)+
  guides(colour=guide_legend(override.aes=list(size=4)))+
  theme(strip.background = element_rect(color="white"), legend.position = "none")

RShannon_DI

correlation(subset(rhizosphere_alpha, variable=="Simpson")$value, subset(rhizosphere_alpha, variable=="Simpson")$Pot_DI)

RSimpson_DI <- ggscatter(subset(rhizosphere_alpha, variable=="Simpson"), x="Pot_DI", 
          y="value", color="PotID", label = "samples2",
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  #stat_cor(method = "pearson", label.x = 3, label.y =350)+
  annotate(geom = "richtext", x=50, y = 0.95, label="*&rho;* = -0.491, *P* = 1.292 &#10005; 10<sup>-3</sup>", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  labs(color="Pot")+
  ylab("Simpson")+
  xlab("Root rot disease index")+
  facet_wrap(Location~variable)+
  guides(colour=guide_legend(override.aes=list(size=4)))+
  theme(strip.background = element_rect(color="white"), legend.position = "none")

RSimpson_DI

ggarrange(EObserved_DI, EShannon_DI, ESimpson_DI, RObserved_DI, RShannon_DI, RSimpson_DI)

#=================================================
# Relative abandunce bar graph (beta diversity)
#=================================================
# sorting relative abandunce data  
relative_abundance <- plot_bar(phy_obj_relative)$data
relative_abundance <- relative_abundance[as.character(1:dim(relative_abundance)[1]),]
relative_abundance$Cycling <- gsub("^","Cycling ", relative_abundance$Cycling)


# output relative abandance data and NA replacement
write.table(relative_abundance,  file = "metagenome_analysis_results/relative_abundance.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "Not assigned")
relative_abundance <- read.csv("metagenome_analysis_results/relative_abundance.tsv", sep="\t")

# Sums data by taxa levels
phylum_relative_abundance <- aggregate(Abundance ~ phylum + Sample + Location + Cycling + PotID, relative_abundance, sum)
class_relative_abundance <- aggregate(Abundance ~ class + Sample + Location + Cycling + PotID, relative_abundance, sum)
order_relative_abundance <- aggregate(Abundance ~ order + Sample + Location + Cycling + PotID, relative_abundance, sum)
family_relative_abundance <- aggregate(Abundance ~ family + Sample + Location + Cycling + PotID, relative_abundance, sum)
genus_relative_abundance <- aggregate(Abundance ~ genus + Sample + Location + Cycling + PotID, relative_abundance, sum)


# top 10 sorting
# phylum_level
sum_phylum_Cycling <- aggregate(Abundance ~ phylum, relative_abundance, sum)
sum_phylum_Cycling <- sum_phylum_Cycling %>% 
  arrange(desc(Abundance))
top_phylum_Cycling <- as.vector(sum_phylum_Cycling[1:10,1])


# class_level
sum_class_Cycling <- aggregate(Abundance ~ class, relative_abundance, sum)
sum_class_Cycling <- sum_class_Cycling %>% 
  arrange(desc(Abundance))
top_class_Cycling <- as.vector(sum_class_Cycling[1:10,1])


# order_level
sum_order_Cycling <- aggregate(Abundance ~ order, relative_abundance, sum)
sum_order_Cycling <- sum_order_Cycling %>% 
  arrange(desc(Abundance))
top_order_Cycling <- as.vector(sum_order_Cycling[1:10,1])


# family_level
sum_family_Cycling <- aggregate(Abundance ~ family, relative_abundance, sum)
sum_family_Cycling <- sum_family_Cycling %>% 
  arrange(desc(Abundance))
top_family_Cycling <- as.vector(sum_family_Cycling[1:10,1])


# genus_level
sum_genus_Cycling <- aggregate(Abundance ~ genus, relative_abundance, sum)
sum_genus_Cycling <- sum_genus_Cycling %>% 
  arrange(desc(Abundance))
top_genus_Cycling <- as.vector(sum_genus_Cycling[1:10,1])


# taxa top bar graph
top_bar <- function(df, fill_label="Top taxa", top_toxa_names){
  levels(df$Cycling) <- levels(df$Cycling)[c(1,3:10,2)]
  
  df$top_taxa <- df[,1]
  df$top_taxa <- factor(df$top_taxa, levels = c(top_toxa_names,"others"))
  df$top_taxa[is.na(df$top_taxa)] <- "others"
  
  df$Sample <- factor(df$Sample, levels = sample_ord) 
  
  df$Cycling <- factor(df$Cycling, levels = unique(df$Cycling)[c(1,3:10,2)])
  
  agg_df <- aggregate(Abundance ~ top_taxa + Sample + Location + Cycling +PotID, df, sum)
  
  p <- ggplot(data = agg_df, mapping = aes(x = Sample, y=Abundance, fill=top_taxa))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c(brewer.pal(n=10, name="Paired"),"#999999"))+
    facet_wrap(Location ~ Cycling, scales = "free_x", nrow=2)+
    theme_light()+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill="grey90", color = "white", size=1),
          strip.text = element_text(color="black"))+
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), expand = c(0,0))+ # fill option
    ylab("Relative abundance")+
    xlab("")+
    labs(fill=fill_label)
}

Cycling_top_10_phylum_bar <- top_bar(phylum_relative_abundance, fill_label = "Top 10 phylum", top_toxa_names = top_phylum_Cycling )
Cycling_top_10_phylum_bar 

Cycling_top_10_class_bar <- top_bar(class_relative_abundance,  fill_label = "Top 10 class", top_toxa_names = top_class_Cycling )
Cycling_top_10_class_bar

Cycling_top_10_order_bar <- top_bar(order_relative_abundance, fill_label = "Top 10 order", top_toxa_names = top_order_Cycling )
Cycling_top_10_order_bar

Cycling_top_10_family_bar <- top_bar(family_relative_abundance, fill_label = "Top 10 family", top_toxa_names = top_family_Cycling )
Cycling_top_10_family_bar

Cycling_top_10_genus_bar <- top_bar(genus_relative_abundance, fill_label = "Top 10 genus", top_toxa_names = top_genus_Cycling )
Cycling_top_10_genus_bar

#================
# 2. 
#================
rarefaction_curve_data <- ggiNEXT( 
  iNEXT(as.data.frame(
    t(OTU)
  )
  )
)$data

rarefaction_curve_data
rarefaction_curve_data$Location <- rarefaction_curve_data$site
rarefaction_curve_data$Location <- gsub("^E[A-Z0-9]{1,}","Endosphere",rarefaction_curve_data$Location)
rarefaction_curve_data$Location <- gsub("^R[A-Z0-9]{1,}","Rhizosphere",rarefaction_curve_data$Location)

rarefaction_curve_data$Cycling <- rarefaction_curve_data$site
rarefaction_curve_data$Cycling <- gsub("^[A-Z]{2,}","Cycling ",rarefaction_curve_data$Cycling)
rarefaction_curve_data$Cycling <- gsub("C[1357]$","",rarefaction_curve_data$Cycling)

rarefaction_curve_data$PotID <- rarefaction_curve_data$site
rarefaction_curve_data$PotID <- gsub("^[A-Z]{1,}[0-9]{1,}","",rarefaction_curve_data$PotID)

rarefaction_curve_data <- arrange(transform(rarefaction_curve_data, Cycling = factor(Cycling, levels=unique(rarefaction_curve_data$Cycling))), Cycling)

rarefaction_curve_max_data <- aggregate(x ~ site + lty + Location + Cycling + PotID, subset(rarefaction_curve_data, lty=="interpolated"), max)
rarefaction_curve_max_data$y <- aggregate(y ~ site + lty + Location + Cycling + PotID, subset(rarefaction_curve_data, lty=="interpolated"), max)$y
rarefaction_curve_max_data$Cycling2 <- as.integer(gsub("[A-Za-z ]{1,}","", rarefaction_curve_max_data$Cycling))
rarefaction_curve_max_data <- arrange(transform(rarefaction_curve_max_data, Cycling = factor(Cycling, levels=unique(rarefaction_curve_data$Cycling))), Cycling)
rarefaction_curve_max_data

cycling_rarefaction_curve_plot <- ggplot(data = subset(rarefaction_curve_data, PotID!="HD"&PotID!="LD"), mapping = aes(x=x , y=y, color=PotID, shape=PotID, linetype=lty, label=site))+
  geom_line()+
  geom_point(data = subset(rarefaction_curve_max_data, PotID!="HD"&PotID!="LD"))+
  theme_light()+
  scale_x_continuous(limits=c(0,110000))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  facet_grid(Location ~ Cycling, scales = "free_y")+
  xlab("Reads")+
  ylab("Inferenced OTU number")+
  labs(color="Pot", shape="Pot")+
  labs(linetype="Method")+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "bottom",
        panel.grid = element_blank())
cycling_rarefaction_curve_plot 

#=======================================
# For sparCC levels (top abundance 95% reads)
#=======================================
# To discard 5% abundance of small amount read OTUs per sample
OTU_relative <- as.matrix(read.csv("metagenome_analysis_results/DADA2_results/relative_abundance_by_OTU.csv", row.names = 1))

OTU_cut_cumulative_rel <- foreach(i = 1:ncol(OTU_relative))%dopar%{
  OTU_rel <- OTU_relative[,i]
  OTU_rel <- sort(OTU_rel, decreasing = TRUE)
  cumsum_OTU_rel <- cumsum(c(OTU_rel))
  cut_OTU <- names(cumsum_OTU_rel[cumsum_OTU_rel <= 0.95])
  as.numeric(gsub("^OTU_","",cut_OTU))
}

OTU_cut_cumulative_rel

# To make top 95% phloseq object
cumulative_cut_OTU_long_table <- t(OTU)
cumulative_cut_OTU_long_table
for(i in 1:dim(OTU_relative)[2]){
  cumulative_cut_OTU_long_table[-OTU_cut_cumulative_rel[[i]],i] <- 0
}
colSums(cumulative_cut_OTU_long_table)
OTU_top_0.95 <- t(cumulative_cut_OTU_long_table)
top_0.95_TAX <- TAX[-which(as.vector(colSums(OTU_top_0.95))==0),]
OTU_top_0.95 <- OTU_top_0.95[,-which(as.vector(colSums(OTU_top_0.95))==0)]
phy_obj_reads_top_0.95 <- phyloseq(OTU_top_0.95,top_0.95_TAX,SAM)

# To make top 95% taxa table for phloseq
OTU_counts_data_top_0.95 <- plot_bar(phy_obj_reads_top_0.95)$data
row.names(OTU_counts_data_top_0.95) <- NULL

# family
unclassified_family <- levels(as.factor(OTU_counts_data_top_0.95$OTU[which(is.na(OTU_counts_data_top_0.95$family)==T)]))
levels(OTU_counts_data_top_0.95$family) <- c(levels(OTU_counts_data_top_0.95$family),unclassified_family)

for(i in 1:length(unclassified_family)){
  OTU_counts_data_top_0.95$family[which(OTU_counts_data_top_0.95$OTU == unclassified_family[i])] <- unclassified_family[i]
}     

top_0.95_family_abundance <- aggregate(Abundance ~ family + Sample + Location + Cycling, OTU_counts_data_top_0.95, sum)

for_sparCC_family <- dcast(data = top_0.95_family_abundance, family ~ Sample)[, c("family", row.names(OTU))]
c("family", row.names(OTU))
colnames(for_sparCC_family)[1] <- "#OTU ID"
c("family",row.names(OTU))
for_sparCC_family


write.table(for_sparCC_family[which(rowSums(for_sparCC_family[,c(2:41)])!=0),c(1:41)], 
            "metagenome_analysis_results/fastspar/Cycling_raw/Endosphere_at_family_for_fastspar.tsv", sep="\t", row.names = F, quote = F)
write.table(for_sparCC_family[which(rowSums(for_sparCC_family[,c(42:81)])!=0),c(1,42:81)], 
            "metagenome_analysis_results/fastspar/Cycling_raw/Rhizosphere_at_family_for_fastspar.tsv" , sep = "\t", row.names = F, quote = F)

#======================
# 1. Rhizosphere family SparCC_df
#======================
# handling SparCC data
sparCC_long_df <- function(sparcc_df){
  remove_rep <- row.names(sparcc_df)
  for(i in 1:ncol(sparcc_df)){
    remove_rep <- cbind(remove_rep, c(as.numeric(rep(NA,i)),as.numeric(rep(0,ncol(sparcc_df)-i))))
  }
  remove_rep <- as.data.frame(remove_rep)[,-1]
  for(i in 1:ncol(sparcc_df)){
    remove_rep[,i]<- as.numeric(remove_rep[,i])
  }
  remove_rep <- sparcc_df + remove_rep
  remove_rep <- melt(cbind(row.names(sparcc_df),remove_rep))
  remove_rep <- na.omit(remove_rep)
  colnames(remove_rep) <- c("taxa1", "taxa2", "value")
  print(remove_rep)
}

R_family_sparCC_cor <- sparCC_long_df(read.csv("metagenome_analysis_results/fastspar/rhizosphere_cycling_family/R_family_median_correlation_fastspar.tsv", sep="\t", row.names = 1))
R_family_sparCC_cor$pvalue <- sparCC_long_df(read.csv("metagenome_analysis_results/fastspar/rhizosphere_cycling_family/R_family_pvalues.tsv", sep="\t", row.names = 1))[,3]
R_family_sparCC_cor
colnames(R_family_sparCC_cor)[3] <- "cor"

R_family_sparCC_cor$abs_cor <- abs(R_family_sparCC_cor$cor)
R_family_sparCC_cor$PorN <- ifelse(R_family_sparCC_cor$cor>=0, "Positive", "Negative")

R_family_sparCC_cor_family_cut_off <- subset(R_family_sparCC_cor, pvalue < 0.01)
R_family_sparCC_cor_family_cut_off <- subset(R_family_sparCC_cor, abs_cor > 0.5)

# To make network obj
R_family_sparCC_net <- network(as.matrix(R_family_sparCC_cor_family_cut_off[,1:2]), direct = FALSE)

# vertex
R_family_vertex <- data.frame(family = network.vertex.names(R_family_sparCC_net))

Centrality <- R_family_sparCC_cor_family_cut_off %>%
  as_tbl_graph(directed=FALSE) %>%
  activate(nodes) %>%
  mutate(eigen = centrality_eigen(directed = F),
         pagerank = centrality_pagerank(directed = F))

Centrality <- as.data.frame(Centrality)
colnames(Centrality)[1] <- "family" 

R_family_vertex <- merge(R_family_vertex, Centrality, by= "family")

Rhizo_family_net_group_A <- c("Cellulomonadaceae", "Alcaligenaceae",
                              "Weeksellaceae", "Flavobacteriaceae","Xanthobacteraceae",
                              "Burkholderiaceae", "Rhizobiaceae", "Sphingobacteriaceae",
                              "Microbacteriaceae", "Devosiaceae", "Sphingomonadaceae",
                              "Promicromonosporaceae", "Beijerinckiaceae","Chitinophagaceae",
                              "Xanthomonadaceae","Comamonadaceae",
                              "Mycobacteriaceae", "Oxalobacteraceae", "Streptomycetaceae", "Inquilinaceae", "Labraceae")

keystonetaxa <- Rhizo_family_net_group_A [Rhizo_family_net_group_A %in%
                               Centrality$family[which(Centrality$eigen > 0.6)]]

Centrality

Rhizo_family_net_group_B <- c("Enterobacteriaceae", "Pseudomonadaceae", "Micrococcaceae", "Nocardiaceae", "Staphylococcaceae",
                              "Brevibacteriaceae","Dermabacteraceae","Sanguibacteraceae", "Nocardioidaceae", "Intrasporangiaceae")

R_family_vertex$group <- R_family_vertex$family

for(i in 1:length(Rhizo_family_net_group_A)){
  R_family_vertex$group <- gsub(Rhizo_family_net_group_A[i],"Group A",R_family_vertex$group)
}

for(i in 1:length(Rhizo_family_net_group_B)){
  R_family_vertex$group <- gsub(Rhizo_family_net_group_B[i],"Group B",R_family_vertex$group)
}

R_family_vertex$group[which(R_family_vertex$group!="Group A"&R_family_vertex$group!="Group B")] <- "No groups"

R_family_sparCC_net %v% 
  "The family network group"= R_family_vertex$group

R_family_sparCC_net %v%
  "Eigen" <- R_family_vertex$eigen

R_family_sparCC_net %e%
  "PorN" <- R_family_sparCC_cor_family_cut_off$PorN

R_family_sparCC_net %e%
  "size" <- R_family_sparCC_cor_family_cut_off$abs_cor

node.color.pal <-c("#F8766D","#00BFC4")
names(node.color.pal) <- c("Group A", "No groups")

set.seed(452)
R_family_sparCC_net_df <- ggnetwork(R_family_sparCC_net)
head(R_family_sparCC_net_df)
R_family_sparCC_net_df$PorN <- factor(R_family_sparCC_net_df$PorN, levels = c("Positive", "Negative"))
R_family_sparcc_net_E <- ggplot(R_family_sparCC_net_df, aes(x=x, y=y, xend = xend, yend = yend))+
  geom_edges(aes(size=size, color=PorN), curvature = 0.15, alpha=0.5)+
  scale_size_continuous(range=c(0,4), breaks = c(0.5, 0.75, 1,3), limits = c(0.4,1))+
  theme_blank()+
  labs(size = "SparCC correlation",
       color = "Positive or negative\ncorrelation")+
  scale_color_manual(values = c("#6CAE00","Orange"))+
  theme(legend.box.margin = margin(t = -205,b=0,l=0,r=20))
R_family_sparcc_net_E

R_family_sparcc_net_N <- ggplot(R_family_sparCC_net_df, aes(x=x, y=y, xend = xend, yend = yend))+
  geom_nodes(aes(size=Eigen, color=`The family network group`))+
  geom_nodetext_repel(aes(label=vertex.names, vjust=Eigen+1), size=4)+
  theme_blank()+
  labs(size = "Eigencentrality")+
  scale_color_manual(values=c("#F8766D","#00BFC4","#555555"))+
  scale_size(range=c(3,10))+
  theme(plot.background = element_blank(),
        panel.background = element_blank())+
  theme(legend.box.margin = margin(t = 205,b=0,l=0,r=0))
R_family_sparcc_net_N

plot.new()
net_E_gtable <- ggplot_gtable(ggplot_build(R_family_sparcc_net_E))
net_N_gtable <- ggplot_gtable(ggplot_build(R_family_sparcc_net_N))

net_panel <-  c(subset(net_E_gtable$layout, name == "panel", se = t:r))
net_legend <- c(subset(net_E_gtable$layout, name == "guide-box", se = t:r))
net_E_gtable$grobs
net_EN_gtable <- gtable_add_grob(net_E_gtable, net_N_gtable$grobs[[which(net_N_gtable$layout$name == "panel")]], net_panel$t, 
                                 net_panel$l, net_panel$b, net_panel$l)

net_EN_gtable <- gtable_add_grob(net_EN_gtable, net_N_gtable$grobs[[15]], #15는 net_N_gtable에서 guide_box의 좌표. 
                                 net_legend$t, net_legend$l, net_legend$b, net_legend$r)

grid.draw(net_EN_gtable)
#===================
# 1. Phylogenetic tree
#===================
# Load fast file
fasta_path <- "metagenome_analysis_results/DADA2_results/OTUs_sequence.fasta"
seqs <- readDNAStringSet(fasta_path)
seqs <- OrientNucleotides(seqs)

# Alignment
aligned_seqs <- AlignSeqs(seqs)
endosphere_aligned_seqs <- aligned_seqs[OTUsID[which(rowSums(OTU_counts_table[,-1][,1:40])>0)]]
rhizosphere_aligned_seqs <- aligned_seqs[OTUsID[which(rowSums(OTU_counts_table[,-1][,41:80])>0)]]
BrowseSeqs(aligned_seqs, highlight = 0)
writeXStringSet(aligned_seqs, file="metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs_aligned_sequence.fasta")
writeXStringSet(endosphere_aligned_seqs, file="metagenome_analysis_results/metagenomic sequence pylogenetic tree/endosphere_OTUs_aligned_sequence.fasta")
writeXStringSet(rhizosphere_aligned_seqs, file="metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere_OTUs_aligned_sequence.fasta")

# Aligned sequences into a database
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(aligned_seqs, type = "XStringSet", dbFile = dbConn, "")

Add2DB(myData=data.frame(identifier=OTUsID,
                         stringsAsFactors=FALSE),
       dbConn)

cons <- IdConsensus(dbConn,
                    threshold=0.3,
                    minInformation=0.1)
#-----endosphere
endosphere_dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(endosphere_aligned_seqs, type = "XStringSet", dbFile = endosphere_dbConn, "")

Add2DB(myData=data.frame(identifier=OTUsID[which(rowSums(OTU_counts_table[,-1][,1:40])>0)],
                         stringsAsFactors=FALSE),
       endosphere_dbConn)

endosphere_cons <- IdConsensus(endosphere_dbConn,
                               threshold=0.3,
                               minInformation=0.1)
#-----rhizosphere
rhizosphere_dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(rhizosphere_aligned_seqs, type = "XStringSet", dbFile = rhizosphere_dbConn, "")

Add2DB(myData=data.frame(identifier=OTUsID[which(rowSums(OTU_counts_table[,-1][,41:80])>0)],
                         stringsAsFactors=FALSE),
       rhizosphere_dbConn)

rhizosphere_cons <- IdConsensus(rhizosphere_dbConn,
                                threshold=0.3,
                                minInformation=0.1)

# calculate a maximum like-hood tree
d <- DistanceMatrix(cons, correction="Jukes-Cantor")
dend <- IdClusters(d, method="ML", type="dendrogram", myXStringSet=cons)
dend <- dendrapply(dend,
                   FUN=function(n) {
                     if(is.leaf(n)) 
                       attr(n, "label") <- 
                         as.expression(substitute(italic(leaf),
                                                  list(leaf=attr(n, "label"))))
                     n
                   })
#------endosphere
endosphere_d <- DistanceMatrix(endosphere_cons, correction="Jukes-Cantor")
endosphere_dend <- IdClusters(endosphere_d, method="ML", type="dendrogram", myXStringSet=endosphere_cons)
endosphere_dend <- dendrapply(endosphere_dend,
                              FUN=function(n) {
                                if(is.leaf(n)) 
                                  attr(n, "label") <- 
                                    as.expression(substitute(italic(leaf),
                                                             list(leaf=attr(n, "label"))))
                                n
                              })

#------rhizosphere
rhizosphere_d <- DistanceMatrix(rhizosphere_cons, correction="Jukes-Cantor")
rhizosphere_dend <- IdClusters(rhizosphere_d, method="ML", type="dendrogram", myXStringSet=rhizosphere_cons)
rhizosphere_dend <- dendrapply(rhizosphere_dend,
                               FUN=function(n) {
                                 if(is.leaf(n)) 
                                   attr(n, "label") <- 
                                     as.expression(substitute(italic(leaf),
                                                              list(leaf=attr(n, "label"))))
                                 n
                               })

# display the phylogenetic tree
plot(dend, yaxt="n", horiz=TRUE)
arrows(-0.1, 6, -0.2, 6, angle=90, length=0.05, code=3)
text(-0.15, 6, "0.1", adj=c(0.5, -0.5))

WriteDendrogram(dend, file = "metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum lkelihood tree (HKY85+G4).newick",
                quoteLabels = F,
                internalLabels = F)
newickfile <- read.table("metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum likelihood tree (HKY85+G4).newick", header = F)$V1
newickfile <- gsub('italic[(][\"]',"",newickfile)
newickfile <-gsub('\")',"",newickfile)
write.table(newickfile, quote = F, col.names = F,row.names = F,"metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum likelihood tree (HKY85+G4).newick")

#-----------endosphere
WriteDendrogram(endosphere_dend , file = "metagenome_analysis_results/metagenomic sequence pylogenetic tree/endosphere OTUs maximum likelihood tree (TN93+G4).newick", 
                quoteLabels = F,
                internalLabels = F)
endosphere_newickfile <- read.table("metagenome_analysis_results/metagenomic sequence pylogenetic tree/endosphere OTUs maximum likelihood tree (TN93+G4).newick", header = F)$V1
endosphere_newickfile <- gsub('italic[(][\"]',"",endosphere_newickfile)
endosphere_newickfile <-gsub('\")',"",endosphere_newickfile)
write.table(endosphere_newickfile, quote = F, col.names = F,row.names = F,"metagenome_analysis_results/metagenomic sequence pylogenetic tree/endosphere OTUs maximum likelihood tree (TN93+G4).newick")

#-----------rhizosphere
WriteDendrogram(rhizosphere_dend , file = "metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere OTUs maximum likelihood tree (F81+G4).newick", 
                quoteLabels = F,
                internalLabels = F)
rhizosphere_newickfile <- read.table("metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere OTUs maximum likelihood tree (F81+G4).newick", header = F)$V1
rhizosphere_newickfile <- gsub('italic[(][\"]',"",rhizosphere_newickfile)
rhizosphere_newickfile <-gsub('\")',"",rhizosphere_newickfile)
write.table(rhizosphere_newickfile, quote = F, col.names = F, row.names = F,"metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere OTUs maximum likelihood tree (F81+G4).newick")

#===================================================================
# bNTI and beta RC
#==================================================================
rownames(OTU_counts_table) <- OTU_counts_table$`#OTU ID`

#------------ endosphere
endosphere_phylogeny <- read.tree("metagenome_analysis_results/metagenomic sequence pylogenetic tree/endosphere OTUs maximum likelihood tree (TN93+G4).newick")
endosphere_tip.label <- endosphere_phylogeny$tip.label

endosphere_match.phylo.OTU <- match.phylo.data(endosphere_phylogeny, OTU_counts_table[endosphere_tip.label,-1][1:40])
# calculate empirical betaMNTD
endosphere.beta.mntd.weighted <- as.matrix(comdistnt(t(endosphere_match.phylo.OTU$data),cophenetic(endosphere_match.phylo.OTU$phy),abundance.weighted=T))
write.csv(endosphere.beta.mntd.weighted, "metagenome_analysis_results/RCI and bNTI/endosphere.betaMNTD_weighted.csv",quote=F);

identical(colnames(endosphere_match.phylo.OTU$data),colnames(endosphere.beta.mntd.weighted)) # just a check, should be TRUE

# calculate randomized betaMNTD
reps <- 999; # number of randomizations

endosphere.rand.weighted.bMNTD.comp <- array(c(-999),dim=c(ncol(endosphere_match.phylo.OTU$data),ncol(endosphere_match.phylo.OTU$data), reps))

pb <- progress_bar$new(total = reps)
for (i in 1:reps) {
  endosphere.rand.weighted.bMNTD.comp[,,i] <- as.matrix(comdistnt(t(endosphere_match.phylo.OTU$data),taxaShuffle(cophenetic(endosphere_match.phylo.OTU$phy)),abundance.weighted=T,exclude.conspecifics = F))
  pb$tick()
}

endosphere.weighted.bNTI <- matrix(c(NA),nrow=ncol(endosphere_match.phylo.OTU$data),ncol=ncol(endosphere_match.phylo.OTU$data))

for (columns in 1:(ncol(endosphere_match.phylo.OTU$data)-1)) {
  for (rows in (columns+1):ncol(endosphere_match.phylo.OTU$data)) {
    endosphere.rand.vals <- endosphere.rand.weighted.bMNTD.comp[rows,columns,]
    endosphere.weighted.bNTI[rows,columns] <- (endosphere.beta.mntd.weighted[rows,columns] - mean(endosphere.rand.vals)) / sd(endosphere.rand.vals)
    rm(endosphere.rand.vals)
  }
}
endosphere.weighted.bNTI_mat <- as.matrix(as.dist(endosphere.weighted.bNTI))
colnames(endosphere.weighted.bNTI_mat) <- sample_ord[1:40]
rownames(endosphere.weighted.bNTI_mat) <- sample_ord[1:40]

for(i in 1:40){
endosphere.weighted.bNTI_mat[i,i] <- NA
}

endosphere.weighted.bNTI_ldf <- melt(endosphere.weighted.bNTI_mat)
endosphere.weighted.bNTI_ldf
colnames(endosphere.weighted.bNTI_ldf) <- c("Sample1" ,"Sample2", "weighted_bNTI")
endosphere.weighted.bNTI_ldf$Sample1 <- factor(endosphere.weighted.bNTI_ldf$Sample1, levels = sample_ord)
endosphere.weighted.bNTI_ldf$Sample2 <- factor(endosphere.weighted.bNTI_ldf$Sample2, levels = sample_ord)

endosphere.weighted.bNTI_ldf$weighted_bNTI2 <- ifelse(endosphere.weighted.bNTI_ldf$weighted_bNTI > 2, "betaNTI > 2",
                                                      ifelse(endosphere.weighted.bNTI_ldf$weighted_bNTI > -2, "2 > betaNTI > -2", "-2 > betaNTI"))

endosphere.weighted.bNTI_ldf$weighted_bNTI2 <- factor(endosphere.weighted.bNTI_ldf$weighted_bNTI2 , levels=c("betaNTI > 2", "2 > betaNTI > -2", "-2 > betaNTI"))


Endosphere_weighted_beta_NTI1 <-
  ggplot(endosphere.weighted.bNTI_ldf, aes(x=Sample1,y=Sample2))+
  geom_tile(aes(fill = weighted_bNTI))+
  scale_fill_gradient2(low = "steelblue", mid="grey80", high = "tomato", limits=c(-6,6), breaks=-3:3*2, na.value = "black")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste(Endosphere~beta,NTI)))
Endosphere_weighted_beta_NTI1

Endosphere_weighted_beta_NTI2 <-
  ggplot(endosphere.weighted.bNTI_ldf, aes(x=Sample1,y=Sample2))+
  geom_tile(aes(fill = weighted_bNTI2))+
  scale_fill_manual(values=c("#FF866C","grey80","#1752FF"),na.value="black",
                    labels = c(expression(paste(beta,NTI>2)),
                               expression(paste(2>beta,NTI)>-2),
                               expression(paste(-2>beta,NTI)),
                               expression(self~pairwise)
                    ))+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste(Endosphere~beta,NTI)))
Endosphere_weighted_beta_NTI2

write.csv(endosphere.weighted.bNTI_mat,"metagenome_analysis_results/RCI and bNTI/endosphere.weighted_bNTI.csv",quote=F)

#-----------------------------
rhizosphere_phylogeny <- read.tree("metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere OTUs maximum likelihood tree (F81+G4).newick")
rhizosphere_tip.label <- rhizosphere_phylogeny$tip.label

rhizosphere_match.phylo.OTU <- match.phylo.data(rhizosphere_phylogeny, OTU_counts_table[rhizosphere_tip.label,-1][41:80])
# calculate empirical betaMNTD
rhizosphere.beta.mntd.weighted <- as.matrix(comdistnt(t(rhizosphere_match.phylo.OTU$data),cophenetic(rhizosphere_match.phylo.OTU$phy),abundance.weighted=T))
write.csv(rhizosphere.beta.mntd.weighted, "metagenome_analysis_results/RCI and bNTI/rhizosphere.betaMNTD_weighted.csv",quote=F);

identical(colnames(rhizosphere_match.phylo.OTU$data),colnames(rhizosphere.beta.mntd.weighted)) # just a check, should be TRUE

# calculate randomized betaMNTD
reps <- 999; # number of randomizations

rhizosphere.rand.weighted.bMNTD.comp <- array(c(-999),dim=c(ncol(rhizosphere_match.phylo.OTU$data),ncol(rhizosphere_match.phylo.OTU$data), reps))

pb <- progress_bar$new(total = reps)
for (i in 1:reps) {
  rhizosphere.rand.weighted.bMNTD.comp[,,i] <- as.matrix(comdistnt(t(rhizosphere_match.phylo.OTU$data),taxaShuffle(cophenetic(rhizosphere_match.phylo.OTU$phy)),abundance.weighted=T,exclude.conspecifics = F))
  pb$tick();
}

rhizosphere.weighted.bNTI <- matrix(c(NA),nrow=ncol(rhizosphere_match.phylo.OTU$data),ncol=ncol(rhizosphere_match.phylo.OTU$data))

for (columns in 1:(ncol(rhizosphere_match.phylo.OTU$data)-1)) {
  for (rows in (columns+1):ncol(rhizosphere_match.phylo.OTU$data)) {
    rhizosphere.rand.vals <- rhizosphere.rand.weighted.bMNTD.comp[rows,columns,]
    rhizosphere.weighted.bNTI[rows,columns] <- (rhizosphere.beta.mntd.weighted[rows,columns] - mean(rhizosphere.rand.vals)) / sd(rhizosphere.rand.vals)
    rm(rhizosphere.rand.vals)
  }
}

rhizosphere.weighted.bNTI_mat <- as.matrix(as.dist(rhizosphere.weighted.bNTI))
for(i in 1:40){
  rhizosphere.weighted.bNTI_mat[i,i] <- NA
}

colnames(rhizosphere.weighted.bNTI_mat) <-sample_ord[41:80]
rownames(rhizosphere.weighted.bNTI_mat) <-sample_ord[41:80]
rhizosphere.weighted.bNTI_ldf <- melt(rhizosphere.weighted.bNTI_mat)
colnames(rhizosphere.weighted.bNTI_ldf) <- c("Sample1" ,"Sample2", "weighted_bNTI")

rhizosphere.weighted.bNTI_ldf$Sample1 <- factor(rhizosphere.weighted.bNTI_ldf$Sample1, levels = sample_ord)
rhizosphere.weighted.bNTI_ldf$Sample2 <- factor(rhizosphere.weighted.bNTI_ldf$Sample2, levels = sample_ord)

rhizosphere.weighted.bNTI_ldf$weighted_bNTI2 <- ifelse(rhizosphere.weighted.bNTI_ldf$weighted_bNTI > 2, "betaNTI > 2",
                                                       ifelse(rhizosphere.weighted.bNTI_ldf$weighted_bNTI > -2, "2 > betaNTI > -2", "-2 > betaNTI"))

rhizosphere.weighted.bNTI_ldf$weighted_bNTI2 <- factor(rhizosphere.weighted.bNTI_ldf$weighted_bNTI2 , levels=c("betaNTI > 2", "2 > betaNTI > -2", "-2 > betaNTI"))

Rhizosphere_weighted_beta_NTI1 <-
  ggplot(rhizosphere.weighted.bNTI_ldf, aes(x=Sample1,y=Sample2))+
  geom_tile(aes(fill = weighted_bNTI))+
  scale_fill_gradient2(low = "steelblue", mid="grey80", high = "tomato", limits=c(-6,6.5), breaks=-3:3*2, na.value="black")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste(Rhizosphere~beta,NTI)))
Rhizosphere_weighted_beta_NTI1

Rhizosphere_weighted_beta_NTI2 <-
  ggplot(rhizosphere.weighted.bNTI_ldf, aes(x=Sample1,y=Sample2))+
  geom_tile(aes(fill = weighted_bNTI2))+
  scale_fill_manual(values=c("#FF866C","grey80","#1752FF"), na.value="black",
                    labels = c(expression(paste(beta,NTI>2)),
                               expression(paste(2>beta,NTI)>-2),
                               expression(paste(-2>beta,NTI)),
                               expression(self~pairwise)
                    ))+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste(Rhizosphere~beta,NTI)))
Rhizosphere_weighted_beta_NTI2

write.csv(rhizosphere.weighted.bNTI,"metagenome_analysis_results/RCI and bNTI/rhizosphere.weighted_bNTI.csv",quote=F)

#-----------------------------
# betaRC function Chase et al. 2011 
#-----------------------------
beta_raup_crick <- function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model	
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite, MARGIN=2, FUN=sum)
  
  
  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur<-rep(1,gamma)
  }
  
  
  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(apply(spXsite, MARGIN=1, FUN=sum)))
  
  ##make_null:
  
  ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels. 
  
  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1
  
  ##null_array will hold the actual null distribution values.  Each element of the array corresponds to a null distribution for each combination of alpha values.  The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  Later the function will find the row of alpha_table with the right combination of alpha values.  That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels. 
  
  
  null_array<-list()
  
  ##looping over each combination of alpha levels:
  pb <- progress_bar$new(total = length(alpha_levels)) # modified by Cho. display progress
  for(a1 in 1:length(alpha_levels)){
    pb$tick()  # modified by Cho. display progress
    for(a2 in a1:length(alpha_levels)){
      
      ##build a null distribution of the number of shared species for a pair of alpha values:
      null_shared_spp<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
        
        
        ##same for com2:
        com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
        
        ##how many species are shared in common?
        null_shared_spp[i]<-sum((com1+com2)>1)
        
      }
      
      
      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp
      
      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[a1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[a2]
      
      #increment the counter for the columns of the alpha table/ elements of the null array
      col_count<-col_count+1
      
      
      
    }
    
  }
  
  ##create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")
  
  
  #####################
  ##do the test:
  
  
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  
  ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      
      ##how many species are shared between the two sites:
      n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      
      ## what was the observed richness of each site?
      obs_a1<-sum(spXsite[i,])
      obs_a2<-sum(spXsite[j,])
      
      ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))
      
      ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)
      
      ##how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)
      
      
      
      rc<-(num_greater_in_null)/reps
      
      
      
      
      if(split_ties){
        
        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }
      
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc<-(rc-.5)*2
      }
      
      
      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }
      
      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }
      
      
      ##store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)
      
      
    }
  }
  
  
  if(as.distance.matrix){
    results<-as.dist(results)
  }	
  
  
  return(results)
  }
beta_raup_crick <- cmpfun(beta_raup_crick)

#=================================================================================================
endosphere_rel_OTU_for.betaRC <- cbind(data.frame(sample_name=sample_ord[1:40]), 
                                  t(OTU_relative[which(rownames(OTU_relative)%in%rownames(endosphere_match.phylo.OTU$data)),1:40])
                                  )

endosphere.betaRC <- beta_raup_crick(endosphere_rel_OTU_for.betaRC)
endosphere.betaRC_matrix <- as.matrix(endosphere.betaRC)

for(i in 1:40){
endosphere.betaRC_matrix[i,i] <- NA
}

write.csv(endosphere.betaRC_matrix, "metagenome_analysis_results/RCI and bNTI/endosphere.betaRC.csv", quote=F)

endosphere.betaRC_ldf <- melt(cbind(data.frame(sample_ord[1:40]), endosphere.betaRC_matrix))
colnames(endosphere.betaRC_ldf) <- c("Sample1", "Sample2", "betaRC")
endosphere.betaRC_ldf$betaRC2 <- ifelse(endosphere.betaRC_ldf$betaRC > 0.95, "betaRC > 0.95",
                                  ifelse(endosphere.betaRC_ldf$betaRC > -0.95, "0.95 > betaRC > -0.95","-0.95 > betaRC"))

endosphere.betaRC_ldf$betaRC2 <- factor(endosphere.betaRC_ldf$betaRC2, levels=c("betaRC > 0.95","0.95 > betaRC > -0.95","-0.95 > betaRC"))
hist(endosphere.betaRC_ldf$betaRC)
endosphere.betaRC_ldf$Sample1 <- factor(endosphere.betaRC_ldf$Sample1, levels=sample_ord[1:40])

Endosphere_betaRC_p <-
ggplot(endosphere.betaRC_ldf, aes(x=Sample1,y=Sample2, fill=betaRC))+
  geom_tile()+
  scale_fill_gradient2(low = "steelblue", mid="grey80", high = "tomato", limits=c(-1,1), na.value="black")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste('Endosphere ', beta[RC])))
Endosphere_betaRC_p

Endosphere_betaRC_p2 <-
  ggplot(endosphere.betaRC_ldf, aes(x=Sample1, y=Sample2, fill=betaRC2))+
  geom_tile()+
  scale_fill_manual(values = c("tomato", "grey80", "steelblue"), na.value="black",
                    labels = c(expression(beta[RC]>0.95),
                               expression(paste(0.95>beta[RC])>-0.95),
                               expression(-0.95>beta[RC]),
                               expression(self~pairwise)
                    )
                    )+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste('Endosphere ', beta[RC])))
Endosphere_betaRC_p2

#---------------------------
rhizosphere_rel_OTU_for.betaRC <- cbind(data.frame(sample_name=sample_ord[41:80]), 
                                       t(OTU_relative[which(rownames(OTU_relative)%in%rownames(rhizosphere_match.phylo.OTU$data)),41:80])
)
rhizosphere.betaRC <- beta_raup_crick(rhizosphere_rel_OTU_for.betaRC)
rhizosphere.betaRC_matrix <- as.matrix(rhizosphere.betaRC)

for(i in 1:40){
  rhizosphere.betaRC_matrix[i,i] <- -NA
}

write.csv(rhizosphere.betaRC_matrix, "metagenome_analysis_results/RCI and bNTI/rhizosphere.betaRC.csv", quote=F)

rhizosphere.betaRC_ldf <- melt(cbind(data.frame(sample_ord[41:80]), rhizosphere.betaRC_matrix))

colnames(rhizosphere.betaRC_ldf) <- c("Sample1", "Sample2", "betaRC")

rhizosphere.betaRC_ldf$betaRC2 <- ifelse(rhizosphere.betaRC_ldf$betaRC > 0.95, "betaRC > 0.95",
                                         ifelse(rhizosphere.betaRC_ldf$betaRC > -0.95, "0.95 > betaRC > -0.95","-0.95 > betaRC"))

rhizosphere.betaRC_ldf$betaRC2 <- factor(rhizosphere.betaRC_ldf$betaRC2, levels=c("betaRC > 0.95","0.95 > betaRC > -0.95","-0.95 > betaRC"))
hist(rhizosphere.betaRC_ldf$betaRC)
rhizosphere.betaRC_ldf$Sample1 <- factor(rhizosphere.betaRC_ldf$Sample1, levels=sample_ord[41:80])

Rhizosphere_betaRC_p <-
  ggplot(rhizosphere.betaRC_ldf, aes(x=Sample1,y=Sample2, fill=betaRC))+
  geom_tile()+
  scale_fill_gradient2(low = "steelblue", mid="grey80", high = "tomato", limits=c(-1,1), na.value = "black")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste('Rhizosphere ', beta[RC])))
Rhizosphere_betaRC_p

Rhizosphere_betaRC_p2 <-
  ggplot(rhizosphere.betaRC_ldf, aes(x=Sample1, y=Sample2, fill=betaRC2))+
  geom_tile()+
  scale_fill_manual(values = c("tomato", "grey80", "steelblue"), na.value = "black",
                    labels = c(expression(beta[RC]>0.95),
                               expression(paste(0.95>beta[RC])>-0.95),
                               expression(-0.95>beta[RC]),
                               expression(self~pairwise)
                    )
                    )+
  annotate(geom="segment", x=Inf, xend=-Inf, y=1:9*4+0.5, yend=1:9*4+0.5, size=0.2)+
  annotate(geom="segment", x=1:10*4+0.5, xend=1:10*4+0.5, y=Inf, yend=-Inf, size=0.2)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=2)+
  theme_light()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        panel.grid = element_blank())+
  xlab("")+
  ylab("")+
  labs(fill=expression(paste('Rhizosphere ', beta[RC])))
Rhizosphere_betaRC_p2

endosphere.betaRC_ldf$site <- "endosphere"
rhizosphere.betaRC_ldf$site <- "rhizosphere"

#-------------------------------------
Stochasticity_turnover<- 
rbind(
  cbind(endosphere.weighted.bNTI_ldf, endosphere.betaRC_ldf[,-(1:2)]),
  cbind(rhizosphere.weighted.bNTI_ldf,rhizosphere.betaRC_ldf[,-(1:2)])
)

Stochasticity_turnover

Stochasticity_turnover$Sample1_cycle <- gsub("C[0-9]$", "", Stochasticity_turnover$Sample1)
Stochasticity_turnover$Sample1_cycle <- as.numeric(gsub("[ER]A", "", Stochasticity_turnover$Sample1_cycle))
Stochasticity_turnover$Sample2_cycle <- gsub("C[0-9]$", "", Stochasticity_turnover$Sample2)
Stochasticity_turnover$Sample2_cycle <- as.numeric(gsub("[ER]A", "", Stochasticity_turnover$Sample2_cycle))

Stochasticity_turnover$sig <- paste(Stochasticity_turnover$weighted_bNTI2, Stochasticity_turnover$betaRC2, sep =", ")

Stochasticity_turnover$site <- gsub("endo", "Endo", Stochasticity_turnover$site)
Stochasticity_turnover$site <- gsub("rhi", "Rhi", Stochasticity_turnover$site)

Stochasticity_turnover_p <- 
ggplot(Stochasticity_turnover)+
  geom_bar(aes(x=Sample1_cycle, y=..count../16, fill=sig), width = 0.8, size=0.2,color="black")+
  scale_x_continuous(breaks = 1:10)+
  scale_y_continuous(labels = scales::percent, breaks = 0:5/5,
                     sec.axis = sec_axis(~., name = "compareed cycle"))+
  scale_fill_manual(values = c(brewer.pal(n=10, name="Paired")),
                    labels = c(expression(paste("[-",infinity, " ~ -2]",",\t[-1.00 ~ -0.95]")),
                                expression(paste("[-",infinity, " ~ -2]",",\t [-0.95 ~ 0.95]")),
                                   expression(paste("[-",infinity, " ~ -2]",",\t  [0.95 ~ 1.00]")),
                                   expression(paste("[-2 ~ 2]",",\t[-0.95 ~ -1.00]")),
                                   expression(paste("[-2 ~ 2]",",\t [-0.95 ~ 0.95]")),
                                   expression(paste("[-2 ~ 2]",",\t  [0.95 ~ 1.00]")),
                                   expression(paste("[ 2 ~ ",infinity, "]",",\t[-0.95 ~ -1.00]")),
                                   expression(paste("[ 2 ~ ",infinity, "]",",\t [-0.95 ~ 0.95]")),
                                   expression(paste("[ 2 ~ ",infinity, "]",",\t  [0.95 ~ 1.00]")),
                               "NA (self pairwised)"))+
  labs(x="Cycle", y="Relative importance", fill=expression(paste("\t   ",beta, NTI, "\t\t   ",beta,"RC")))+
  theme_light()+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  facet_grid(Sample2_cycle~site)


Stochasticity_turnover_p2 <- 
  ggplot(subset(Stochasticity_turnover,Sample2_cycle==1))+
  geom_bar(aes(x=Sample1_cycle, y=..count../16, fill=sig), width = 0.8, size=0.2,color="black")+
  scale_x_continuous(breaks = 1:10)+
  scale_y_continuous(labels = scales::percent, breaks = 0:5/5, expand=c(0,0), limits = c(0,1.02))+
  scale_fill_manual(values = c(brewer.pal(n=10, name="Paired")),
                    labels = c(expression(paste("[-",infinity, " ~ -2]",",\t[-1.00 ~ -0.95]")),
                               expression(paste("[-",infinity, " ~ -2]",",\t [-0.95 ~ 0.95]")),
                               expression(paste("[-",infinity, " ~ -2]",",\t  [0.95 ~ 1.00]")),
                               expression(paste("[-2 ~ 2]",",\t[-0.95 ~ -1.00]")),
                               expression(paste("[-2 ~ 2]",",\t [-0.95 ~ 0.95]")),
                               expression(paste("[-2 ~ 2]",",\t  [0.95 ~ 1.00]")),
                               expression(paste("[ 2 ~ ",infinity, "]",",\t[-0.95 ~ -1.00]")),
                               expression(paste("[ 2 ~ ",infinity, "]",",\t [-0.95 ~ 0.95]")),
                               expression(paste("[ 2 ~ ",infinity, "]",",\t  [0.95 ~ 1.00]")),
                               "self pairwised"))+
  labs(x="Cycle", y="Relative importance", fill=expression(paste("\t   ",beta, NTI, "\t\t   ",beta,"RC")))+
  theme_light()+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  facet_grid(. ~ site)+
  annotate(geom="segment", x=c(-Inf,-Inf), xend=c(Inf,-Inf), y=c(-Inf,-Inf), yend=c(-Inf, Inf), size=1.5)+
  theme(strip.background = element_blank(), strip.text = element_text(color="black"))

Stochasticity_turnover_p2


#=====================================
# turnover group volcano plot at family
#=====================================
RC_family <- OTU_seq_tax_counts_table[, c(7,50:89)]
RC_family <- melt(RC_family)
RC_family <- aggregate(value ~ family + variable, RC_family, sum)
RC_family <- dcast(RC_family, family ~ variable)
rownames(RC_family) <- RC_family$family
RC_family <- RC_family[,-1]

RC_family <- RC_family[rowSums(RC_family==0) < 20,] #remove excess zeros
RC_family
RC_family_dds <- DESeqDataSetFromMatrix(countData = round(RC_family),
                                    colData = rhizo_cycling_SAM,
                                    design = ~ turnover_group) 

RC_family_dds <- DESeq(RC_family_dds)
?results
RC_family_res <- results(RC_family_dds, pAdjustMethod = "fdr")
RC_family_res <- data.frame(RC_family_res[order(RC_family_res$padj),])

write.csv(RC_family_res, "metagenome_analysis_results/family_deseq_Cycling_before_turnover_group.csv")

RC_family_res$centrality <- "Eigencentrality < 0.6"
RC_family_res$centrality[which(rownames(RC_family_res) %in% c(keystonetaxa, "Pseudomonadaceae"))]  <-  "Eigencentrality > 0.6"
RC_family_res$family <- rownames(RC_family_res)
RC_family_res$y <- ifelse(-log10(RC_family_res$padj) > 10, -log10(RC_family_res$padj)-9.5, -log10(RC_family_res$padj))
RC_family_volcano <-
  ggplot(RC_family_res, aes(x = log2FoldChange, y = y, shape=centrality)) +
  geom_point(data=subset(RC_family_res, log2FoldChange > 1 & padj < 0.01), color = "#F8766D", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(RC_family_res, log2FoldChange < - 1 & padj < 0.01), color = "#F8766D", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(RC_family_res, (log2FoldChange > - 1& log2FoldChange < 1) | padj > 0.01), aes(shape=centrality), color="grey50")+
  geom_text_repel(data= RC_family_res[which(rownames(RC_family_res) %in% c(keystonetaxa, "Pseudomonadaceae")),] , aes(label=family))+
  scale_x_continuous(breaks=-5:3*2, labels = -5:3*2, limits = c(-10,4), expand = c(0,0))+
  scale_y_continuous(breaks=c(0:8, 8.55, 8.95, 9.5, 10.5), labels = c(0:8,"","",19,20), limits = c(0,10.5), expand = c(0,0))+
  scale_shape_manual(values = c(1,19))+
  ylab(expression(-log[10](italic(P[adj]))))+
  xlab(expression(log[2](FC)))+
  labs(color="", shape="")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x= 0, xend= 0, y = Inf, yend=-Inf, size = 0.1)+
  annotate(geom="segment", x=-Inf, xend=Inf, y = 2, yend = 2, size = 0.2, color="red")+
  annotate(geom="text", x=-6, y=2, vjust=-0.5, label=expression(italic(P[adj]) == 0.01), color="red")+
  annotate(geom="segment", x=-Inf, xend=Inf, y = -log10(0.07671760), yend = -log10(0.07671760), size = 0.2, color="grey50")+
  annotate(geom="segment", x=c(1, -1), xend = c(1, -1), y = Inf, yend = -Inf, size = 0.2, color="grey30", linetype=3)+
  annotate(geom="text", x=-6, y=-log10(0.07671760), vjust=-0.5, label=expression(P == 0.05), color="grey50")+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=8.55, yend=8.95, color="white", size=2)+
  annotate(geom="segment", x=c(-1,0,1), xend=c(-1,0,1), y=8.55, yend=8.95, color="white", size=2)+
  annotate(geom="curve", x=-5:1*2, xend=-5:1*2+1, curvature = 0.3, y=8.55, yend=8.55, size=0.3)+
  annotate(geom="curve", x=-5:1*2+1, xend=-5:1*2+2, curvature = -0.3, y=8.55, yend=8.55, size=0.3)+
  annotate(geom="curve", x=-5:1*2, xend=-5:1*2+1, curvature = 0.3, y=8.95, yend=8.95, size=0.3)+
  annotate(geom="curve", x=-5:1*2+1, xend=-5:1*2+2, curvature = -0.3, y=8.95, yend=8.95, size=0.3)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks.length=unit(0.2,"cm"),
        axis.ticks = element_line(color = "black", size=0.3))+
  theme(legend.position="top")

RC_family_volcano

EC_family <- OTU_seq_tax_counts_table[, c(7,10:49)]
EC_family <- melt(EC_family)
EC_family <- aggregate(value ~ family + variable, EC_family, sum)
EC_family <- dcast(EC_family, family ~ variable)
rownames(EC_family) <- EC_family$family
EC_family <- EC_family[,-1]

EC_family <- EC_family[rowSums(EC_family==0) < 20,] #remove excess zeros
EC_family
EC_family_dds <- DESeqDataSetFromMatrix(countData = round(EC_family),
                                        colData = endo_cycling_SAM,
                                        design = ~ turnover_group) 

EC_family_dds <- DESeq(EC_family_dds)
EC_family_res <- results(EC_family_dds, pAdjustMethod = "fdr")
EC_family_res <- data.frame(EC_family_res[order(EC_family_res$padj),])
EC_family_res
write.csv(EC_family_res, "metagenome_analysis_results/family_deseq_Cycling_before_turnover_group.csv")

EC_family_res$centrality <- "Eigencentrality < 0.6"
EC_family_res$centrality[which(rownames(EC_family_res) %in% c(keystonetaxa, "Pseudomonadaceae"))]  <-  "Eigencentrality > 0.6"
EC_family_res$family <- rownames(EC_family_res)

EC_family_res

EC_family_volcano <-
  ggplot(EC_family_res, aes(x = log2FoldChange, y = -log10(padj), shape=centrality)) +
  geom_point(data=subset(EC_family_res, log2FoldChange > 1 & padj < 0.01), color = "#F8766D", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(EC_family_res, log2FoldChange < - 1 & padj < 0.01), color = "#F8766D", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(EC_family_res, (log2FoldChange > - 1& log2FoldChange < 1) | padj > 0.01), aes(shape=centrality), color="grey50")+
  geom_text_repel(data= EC_family_res[which(rownames(EC_family_res) %in% c(keystonetaxa, "Pseudomonadaceae")),] , aes(label=family))+
  scale_x_continuous(breaks=-6:2, labels = -6:2, limits = c(-6,2))+
  scale_y_continuous(breaks=0:3, labels = 0:3, limits = c(0,10.5), expand = c(0,0))+
  scale_shape_manual(values = c(1,19))+
  ylab(expression(-log[10](italic(P[adj]))))+
  xlab(expression(log[2](FC)))+
  labs(shape="", color="")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x= 0, xend= 0, y = Inf, yend=-Inf, size = 0.1)+
  annotate(geom="segment", x=-Inf, xend=Inf, y = 2, yend = 2, size = 0.2, color="red")+
  annotate(geom="text", x=1, y=2, vjust=-0.5, label=expression(italic(P[adj]) == 0.01), color="red")+
  annotate(geom="segment", x=-Inf, xend=Inf, y = -log10(0.07671760), yend = -log10(0.07671760), size = 0.2, color="grey50")+
  annotate(geom="segment", x=c(1, -1), xend = c(1, -1), y = Inf, yend = -Inf, size = 0.2, color="grey30", linetype=3)+
  annotate(geom="text", x=1, y=-log10(0.07671760), vjust=-0.5, label=expression(P == 0.05), color="grey50")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks.length=unit(0.2,"cm"),
        axis.ticks = element_line(color = "black", size=0.3))+
  theme(legend.position="top")

EC_family_volcano

ggarrange(EC_family_volcano, RC_family_volcano, labels = c("A","B"))

#=========
# picrust
#=========
#rhizosphere KEGG onthology predction
rhizo_cycling_SAM <- SAM[41:80]
rhizo_cycling_SAM$turnover_group <- ifelse(rhizo_cycling_SAM$Cycling %in% 1:3, "1 ~ 3", "4 ~ 10")
rhizo_cycling_SAM$turnover_group <- factor(rhizo_cycling_SAM$turnover_group, levels= rev(unique(rhizo_cycling_SAM$turnover_group)))

RC_KO <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/rhizosphere_cycling/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"), sep="\t", row.names = 1)

RC_KO <- RC_KO[rowSums(RC_KO==0) < 20,] #remove excess zeros
RC_KO_dds <- DESeqDataSetFromMatrix(countData = round(RC_KO),
                                    colData = rhizo_cycling_SAM,
                                    design = ~ turnover_group) 
RC_KO_dds
RC_KO_dds <- DESeq(RC_KO_dds)
RC_KO_res <- results(RC_KO_dds)
RC_KO_res <- data.frame(RC_KO_res[order(RC_KO_res$padj),])

write.csv(RC_KO_res[which(RC_KO_res$log2FoldChange > 1 & RC_KO_res$padj < 0.01),], "metagenome_analysis_results/picrust2_results/rhizosphere_cycling/picrust_KO_deseq_Cycling_before_turnover_group.csv")
write.csv(RC_KO_res[which(RC_KO_res$log2FoldChange < -1 & RC_KO_res$padj < 0.01),], "metagenome_analysis_results/picrust2_results/rhizosphere_cycling/picrust_KO_deseq_Cycling_after_turover_group.csv")

RC_KO_res$nif_cluster <- "Not nif gene"
RC_KO_res$nif_cluster_gene <- NA
nifs_KO <- c("K02584","K02585","K02586","K02587","K02588","K02589","K02590","K02591","K02592","K02593","K02594","K02595","K02596","K02597")
names(nifs_KO) <- c("nifA", "nifB", "nifD", "nifE", "nifH", "nifHD1", "nifHD2", "nifK", "nifN", "nifT", "nifV", "nifW", "nifX", "nifZ")

RC_KO_res[nifs_KO,]$nif_cluster <- "nif cluster"
RC_KO_res[nifs_KO,]$nif_cluster_gene <-
  c("nifA", "nifB", "nifD", "nifE", "nifH", "nifHD1", "nifHD2", "nifK", "nifN", "nifT", "nifV", "nifW", "nifX", "nifZ")

RC_KO_res$reconfigured_y <- -log10(RC_KO_res$padj)
RC_KO_res[!is.na(RC_KO_res$nif_cluster_gene),]
RC_KO_res


RC_KO_volcano <-
  ggplot(RC_KO_res, aes(x = ifelse(log2FoldChange > 6, log2FoldChange -16, log2FoldChange),
                        y = reconfigured_y,
                        color = nif_cluster, shape = nif_cluster)) +
  geom_point(data=subset(RC_KO_res, log2FoldChange > 1 & padj < 0.01), color = "pink2", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(RC_KO_res, log2FoldChange < - 1 & padj < 0.01), color = "pink2", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(RC_KO_res, nif_cluster != "nif cluster" & (log2FoldChange <= 1 & log2FoldChange >= -1) | padj > 0.01), alpha = 1, stroke = 0.5) +
  labs(color = "Nif gene cluster", shape = "Nif gene cluster")+
  geom_text_repel(data=subset(RC_KO_res, nif_cluster == "nif cluster"), aes(label = nif_cluster_gene), color="black")+
  geom_point(data=subset(RC_KO_res, nif_cluster == "nif cluster")) +
  scale_shape_manual(values=c(19,1)) +
  scale_color_manual(values=c("#F8766D","grey80"))+
  scale_x_continuous(breaks=-5:3*2, labels = -5:3*2, limits = c(-10,5.5))+
  scale_y_continuous(breaks=0:5*5, labels = 0:5*5, limits = c(0,28), expand = c(0,0))+
  ylab(expression(-log[10](italic(P[adj]))))+
  xlab(expression(log[2](FC)))+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x= 0, xend= 0, y = Inf, yend=-Inf, size = 0.1, linetype=3)+
  annotate(geom="segment", x=-Inf, xend=Inf, y = 2, yend = 2, size = 0.2, color="red")+
  annotate(geom="text", x=-6, y=2, vjust=-0.5, label=expression(italic(P[adj]) == 0.01), color="red")+
  annotate(geom="segment", x=-Inf, xend=Inf, y = 0.9325962, yend = 0.9325962, size = 0.2, color="grey50")+
  annotate(geom="segment", x=c(1, -1), xend = c(1, -1), y = Inf, yend = -Inf, size = 0.2, color="grey50")+
  annotate(geom="text", x=-6, y=0.9019, vjust=-0.5, label=expression(P == 0.05), color="grey50")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks.length=unit(0.2,"cm"),
        axis.ticks = element_line(color = "black", size=0.3),
        legend.position = "top")

RC_KO_volcano

RC_KO_res_before_turnover_kegg_mapper_module <- read.csv("metagenome_analysis_results/picrust2_results/rhizosphere_cycling/Cycling_before_turnover_rhizosphere_picrust2_KO_deseq(adjp1_log2foldchage0.5)_kegg_mapper_searchModule_hit.txt", sep="*", header = F)
RC_KO_res_before_turnover_kegg_mapper_module <- data.frame("modul_results" = regmatches(RC_KO_res_before_turnover_kegg_mapper_module$V1, 
                                                                                        regexpr(pattern = "M[0-9]{5,5}.{1,}", RC_KO_res_before_turnover_kegg_mapper_module$V1)))
RC_KO_res_before_turnover_kegg_mapper_module$module_ID<- substr(RC_KO_res_before_turnover_kegg_mapper_module$modul_results,1,6)
RC_KO_res_before_turnover_kegg_mapper_module$before_turnout_hit <- as.numeric(
  substr(
    RC_KO_res_before_turnover_kegg_mapper_module$modul_results, 
    nchar(RC_KO_res_before_turnover_kegg_mapper_module$modul_results)-1,
    nchar(RC_KO_res_before_turnover_kegg_mapper_module$modul_results)-1))
RC_KO_res_before_turnover_kegg_mapper_module$modul_results <- substr(RC_KO_res_before_turnover_kegg_mapper_module$modul_results,8, nchar(RC_KO_res_before_turnover_kegg_mapper_module$modul_results)-4)

RC_KO_res_after_turnover_kegg_mapper_module <- read.csv("metagenome_analysis_results/picrust2_results/rhizosphere_cycling/Cycling_after_turnover_rhizosphere_picrust2_KO_deseq(adjp1_log2foldchage-0.5)_kegg_mapper_searchModule_hit.txt", sep="*", header = F)
RC_KO_res_after_turnover_kegg_mapper_module <- data.frame("modul_results" = regmatches(RC_KO_res_after_turnover_kegg_mapper_module$V1, 
                                                                                       regexpr(pattern = "M[0-9]{5,5}.{1,}", RC_KO_res_after_turnover_kegg_mapper_module$V1)))

RC_KO_res_after_turnover_kegg_mapper_module$module_ID<- substr(RC_KO_res_after_turnover_kegg_mapper_module$modul_results,1,6)
RC_KO_res_after_turnover_kegg_mapper_module$after_turnout_hit <- as.numeric(
  substr(
    RC_KO_res_after_turnover_kegg_mapper_module$modul_results, 
    nchar(RC_KO_res_after_turnover_kegg_mapper_module$modul_results)-1,
    nchar(RC_KO_res_after_turnover_kegg_mapper_module$modul_results)-1))
RC_KO_res_after_turnover_kegg_mapper_module$modul_results <- substr(RC_KO_res_after_turnover_kegg_mapper_module$modul_results,8, nchar(RC_KO_res_after_turnover_kegg_mapper_module$modul_results)-4)

RC_KO_res_kegg_mapper_module <- merge(RC_KO_res_before_turnover_kegg_mapper_module, RC_KO_res_after_turnover_kegg_mapper_module, all=T, key=module_ID)
RC_KO_res_kegg_mapper_module$before_turnout_hit[is.na(RC_KO_res_kegg_mapper_module$before_turnout_hit)] <- 0
RC_KO_res_kegg_mapper_module$after_turnout_hit[is.na(RC_KO_res_kegg_mapper_module$after_turnout_hit)] <- 0
RC_KO_res_kegg_mapper_module <- RC_KO_res_kegg_mapper_module[order(-RC_KO_res_kegg_mapper_module$before_turnout_hit),c(2,1,3,4)]
colnames(RC_KO_res_kegg_mapper_module)[2] <- "modul_discription"
colnames(RC_KO_res_kegg_mapper_module) <- c("KEGG modul ID", "modul discriotion", "before hit", "after hit")

RC_KO_res_kegg_mapper_module$`KEGG modul ID` <- as.factor(RC_KO_res_kegg_mapper_module$`KEGG modul ID`)
RC_KO_res_kegg_mapper_module <- RC_KO_res_kegg_mapper_module[order(RC_KO_res_kegg_mapper_module$`KEGG modul ID`),]
RC_KO_res_kegg_mapper_module <- RC_KO_res_kegg_mapper_module[order(RC_KO_res_kegg_mapper_module$`before hit`, decreasing = T),]

rownames(RC_KO_res_kegg_mapper_module) <-NULL

datatable(RC_KO_res_kegg_mapper_module,
          options = list(
            pageLength = 10
          )) %>% 
  formatStyle(
    'KEGG modul ID',
    target = "row",
    background = styleEqual("M00175","yellow")) %>%
  formatStyle(
    'before hit',
    background = styleColorBar(RC_KO_res_kegg_mapper_module$`before hit`, 'lightblue'),
    backgroundSize = '100% 90%',
    backgroundRepeat = 'no-repeat') %>%
  formatStyle(
    'after hit',
    background = styleColorBar(RC_KO_res_kegg_mapper_module$`after hit`, 'lightblue'),
    backgroundSize = '100% 90%',
    backgroundRepeat = 'no-repeat')

# endosphere KEGG onthology prediction
endo_cycling_SAM <- SAM[1:40,]
endo_cycling_SAM$turnover_group <- ifelse(endo_cycling_SAM$Cycling %in% 1:3, "1 ~ 3", "4 ~ 10")
endo_cycling_SAM$turnover_group <- factor(endo_cycling_SAM$turnover_group, levels= rev(unique(endo_cycling_SAM$turnover_group)))

EC_KO <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/endosphere_cycling/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"), sep="\t", row.names = 1)

EC_KO <- EC_KO[rowSums(EC_KO==0) < 20,] #remove excess zeros

EC_KO_dds <- DESeqDataSetFromMatrix(countData = round(EC_KO),
                                    colData = endo_cycling_SAM,
                                    design = ~ turnover_group) 
EC_KO_dds <- DESeq(EC_KO_dds)
EC_KO_res <- results(EC_KO_dds)
EC_KO_res <- data.frame(EC_KO_res[order(EC_KO_res$padj),])
EC_KO_res
EC_KO_res$nif_cluster <- "Not nif gene"
EC_KO_res$nif_cluster_gene <-NA
EC_KO_res[nifs_KO,]$nif_cluster <- "nif cluster"
EC_KO_res[nifs_KO,]$nif_cluster_gene <-
  c("nifA", "nifB", "nifD", "nifE", "nifH", "nifHD1", "nifHD2", "nifK", "nifN", "nifT", "nifV", "nifW", "nifX", "nifZ")


EC_KO_description <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/endosphere_cycling/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"), sep="\t", row.names = 1)[,1]
names(EC_KO_description) <- rownames(EC_KO)
EC_KO_description
EC_KO_res$gene_description <- EC_KO_description[rownames(EC_KO_res)]
EC_KO_res$gene_name <- gsub(";.{1,}$","",EC_KO_res$gene_description)
EC_KO_res$gene_name <- gsub(",.{1,}$","",EC_KO_res$gene_name)
EC_KO_res$gene_name
EC_KO_res$reconfigured_y <- -log10(EC_KO_res$padj)

EC_KO_res[EC_KO_res$nif_cluster!="Not nif gene",]
EC_KO_res
EC_KO_volcano <-
  ggplot(EC_KO_res, aes(x = log2FoldChange,
                        y = reconfigured_y,
                        color = nif_cluster, shape = nif_cluster)) + 
  geom_point(data=subset(EC_KO_res, log2FoldChange > 1 & padj < 0.01), color = "pink2", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(EC_KO_res, log2FoldChange < - 1 & padj < 0.01), color = "pink2", alpha = 1, stroke = 0.5)+
  geom_point(data=subset(EC_KO_res, nif_cluster != "nif cluster" & (log2FoldChange <= 1 & log2FoldChange >= -1) | padj > 0.01), alpha = 1, stroke = 0.5) +
  labs(color = "Nif gene cluster", shape = "Nif gene cluster")+
  geom_text_repel(data=subset(EC_KO_res, nif_cluster == "nif cluster"), aes(label = nif_cluster_gene), color="black")+
  geom_point(data=subset(EC_KO_res, nif_cluster == "nif cluster")) +
  scale_shape_manual(values=c(19,1)) +
  scale_color_manual(values=c("#F8766D","grey80"))+
  scale_x_continuous(breaks = -5:5, limits=c(-4,4))+
  scale_y_continuous(breaks=0:5*5, labels = 0:5*5, limits = c(0,28), expand = c(0,0))+
  ylab(expression(-log[10](italic(P[adj]))))+
  xlab(expression(log[2](FC)))+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1.5)+
  annotate(geom="segment", x= 0, xend= 0, y = Inf, yend=-Inf, size = 0.1)+
  annotate(geom="segment", x=-Inf, xend=Inf, y = 2, yend = 2, size = 0.2, color="red")+
  annotate(geom="text", x=-2, y=2, vjust=-0.5, label=expression(italic(P[adj]) == 0.01), color="red")+
  annotate(geom="segment", x=-Inf, xend=Inf, y = 0.8153857, yend = 0.8153857, size = 0.2, color="grey50")+
  annotate(geom="text", x=-2, y=0.8153857, vjust=-0.5, label=expression(P == 0.05), color="grey50")+
  annotate(geom="segment", x=c(1, -1), xend = c(1, -1), y = Inf, yend = -Inf, size = 0.2, color="grey30", linetype=3)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks.length=unit(0.2,"cm"),
        axis.ticks = element_line(color = "black", size=0.3),
        legend.position = "top")
EC_KO_volcano

ggarrange(EC_KO_volcano, RC_KO_volcano, nrow =1, labels=c("A","B"))

RC_MetaCYC <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/rhizosphere_cycling/pathways_out/path_abun_unstrat_descrip.tsv.gz"), sep="\t", row.names = 2)

amino_acid_synthesis <- 
  c("superpathway of arginine and polyamine biosynthesis",
    "L-arginine biosynthesis I (via L-ornithine)",
    "L-arginine biosynthesis II (acetyl cycle)",
    "superpathway of L-aspartate and L-asparagine biosynthesis",
    "superpathway of branched amino acid biosynthesis",
    "superpathway of aromatic amino acid biosynthesis",
    "L-lysine biosynthesis I",
    "L-ornithine biosynthesis",
    "L-histidine biosynthesis",
    "L-methionine biosynthesis III",
    "L-isoleucine biosynthesis I (from threonine)",
    "superpathway of L-lysine, L-threonine and L-methionine biosynthesis I",
    "L-lysine biosynthesis II",
    "L-lysine biosynthesis III",
    "superpathway of L-isoleucine biosynthesis I",
    "L-lysine biosynthesis V",
    "L-isoleucine biosynthesis II",
    "L-lysine biosynthesis VI",
    "L-isoleucine biosynthesis III",
    "L-isoleucine biosynthesis IV",
    "L-arginine biosynthesis III (via N-acetyl-L-citrulline)",
    "superpathway of L-methionine biosynthesis (by sulfhydrylation)",
    "superpathway of L-methionine biosynthesis (transsulfuration)",
    "L-glutamate and L-glutamine biosynthesis",
    "superpathway of L-phenylalanine biosynthesis",
    "superpathway of L-tyrosine biosynthesis",
    "L-arginine biosynthesis IV (archaebacteria)",
    "superpathway of L-alanine biosynthesis",
    "superpathway of L-serine and glycine biosynthesis I",
    "superpathway of L-threonine biosynthesis",
    "L-tryptophan biosynthesis",
    "L-valine biosynthesis"
  )

amino_acid_degradation <- 
  c("superpathway of L-arginine, putrescine, and 4-aminobutanoate degradation",
    "L-arginine degradation II (AST pathway)",
    "glucose degradation (oxidative)",
    "L-histidine degradation I",
    "L-leucine degradation I",
    "superpathway of L-arginine and L-ornithine degradation",
    "superpathway of ornithine degradation",
    "L-glutamate degradation V (via hydroxyglutarate)",
    "L-histidine degradation II",
    "L-glutamate degradation VIII (to propanoate)",
    "L-tryptophan degradation to 2-amino-3-carboxymuconate semialdehyde",
    "L-tryptophan degradation IX",
    "L-tryptophan degradation XII (Geobacillus)",
    "L-tyrosine degradation I",
    "L-valine degradation I"
  )

RC_amino_acid_biosynt <- RC_MetaCYC[amino_acid_synthesis,]
dim(RC_amino_acid_biosynt)
RC_amino_acid_biosynt$pathway <- row.names(RC_amino_acid_biosynt)
RC_amino_acid_biosynt <- melt(RC_amino_acid_biosynt, id.vars = 1)
colnames(RC_amino_acid_biosynt)[2] <- "sample"
RC_amino_acid_biosynt$group <- ifelse(RC_amino_acid_biosynt$sample %in% rownames(rhizo_cycling_SAM)[rhizo_cycling_SAM$turnover_group == "1 ~ 3"], "1 ~ 3", "4 ~ 10")
RC_amino_acid_biosynt$group2 <- paste(RC_amino_acid_biosynt$group, RC_amino_acid_biosynt$pathway)

RC_amino_acid_biosynt$pathway <- gsub(" [(]", "\n(", RC_amino_acid_biosynt$pathway)
RC_amino_acid_biosynt$pathway <- gsub(" of", "\nof", RC_amino_acid_biosynt$pathway)
RC_amino_acid_biosynt$pathway <- gsub(" and", "\nand", RC_amino_acid_biosynt$pathway)
RC_amino_acid_biosynt$pathway <- gsub(" biosynthesis", "\nbiosynthesis", RC_amino_acid_biosynt$pathway)


for(i in c(1:64)){
  st.p.value <- shapiro.test(subset(RC_amino_acid_biosynt, group2 == unique(RC_amino_acid_biosynt$group2)[i])$value)[2]
  print(st.p.value)
}

wilcox_tP <- c()
for(i in 1:32){
  wilcox_tP <- c(wilcox_tP, wilcox.test(value ~ group,
                                        subset(RC_amino_acid_biosynt, pathway == unique(RC_amino_acid_biosynt$pathway)[i], alternative = "two.sided"))$p.value)
  
}

wilcox_tP <- as.numeric(unlist(wilcox_tP))
wilcox_tP <- data.frame("pathway" = unique(RC_amino_acid_biosynt$pathway), 
                        "wilcox_testP" = wilcox_tP,
                        "expressionP" = paste0("P = ", sprintf("%.3f", wilcox_tP)))
wilcox_tP$test_ajdP<- p.adjust(wilcox_tP$wilcox_testP, "fdr")

wilcox_tP$expression_adjP <- paste0("italic(P[adj]) == ", sprintf("%.3f", wilcox_tP$test_ajdP))

signif_wilcoxP <- RC_amino_acid_biosynt$pathway %in% wilcox_tP$pathway[wilcox_tP$test_ajdP <0.05]

RC_amino_acid_biosynt$signif <- NA
RC_amino_acid_biosynt[signif_wilcoxP,"signif"] <- "*"

subset(RC_amino_acid_biosynt, signif=="*")

amino_acid_biosynthesis_P <-
  ggplot(RC_amino_acid_biosynt, aes(x=group, y=value, fill = group))+
  geom_rect(data = subset(RC_amino_acid_biosynt, signif=="*"), fill="yellow", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha = 0.004)+
  geom_boxplot()+
  facet_wrap(c("pathway"), scales="free_y", ncol = 4)+
  geom_nodetext(data=wilcox_tP, aes(label = expression_adjP, x=1.5, y=Inf, fill=NULL),vjust = 1.5, size = 3, parse = T)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=1)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1)+
  xlab("")+
  ylab("PICRUSt2 predicted pathway abundances")+
  theme_light()+
  theme(legend.position = "none",
        strip.background = element_rect(fill="grey60"))+
  ggtitle("Amino acid biosynthesis pathway abundance",
          subtitle = "before turnover (cycling 1 ~ 3, n=12), after turnover (cycling 4 ~ 10, n=28), unpaired two-Samples Wilcoxon test")


amino_acid_biosynthesis_P

RC_amino_acid_degrad <- RC_MetaCYC[amino_acid_degradation,]
dim(RC_amino_acid_degrad)
RC_amino_acid_degrad$pathway <- row.names(RC_amino_acid_degrad)
RC_amino_acid_degrad <- melt(RC_amino_acid_degrad, id.vars = 1)
colnames(RC_amino_acid_degrad)[2] <- "sample"
RC_amino_acid_degrad$group <- ifelse(RC_amino_acid_degrad$sample %in% rownames(rhizo_cycling_SAM)[rhizo_cycling_SAM$turnover_group == "1 ~ 3"],"1 ~ 3","4 ~ 10")
RC_amino_acid_degrad$group2 <- paste(RC_amino_acid_degrad$group, RC_amino_acid_degrad$pathway)
RC_amino_acid_degrad$group2

RC_amino_acid_degrad$pathway <- gsub(" [(]", "\n(", RC_amino_acid_degrad$pathway)
RC_amino_acid_degrad$pathway <- gsub(" of", "\nof", RC_amino_acid_degrad$pathway)
RC_amino_acid_degrad$pathway <- gsub(" and", "\nand", RC_amino_acid_degrad$pathway)
RC_amino_acid_degrad$pathway <- gsub(" to", "\nto", RC_amino_acid_degrad$pathway)
RC_amino_acid_degrad$pathway <- gsub(" degradation", "\ndegradation", RC_amino_acid_degrad$pathway)
RC_amino_acid_degrad$pathway <- gsub(" semialde", "\nsemialde", RC_amino_acid_degrad$pathway)

for(i in c(1:30)){
  st.p.value <- shapiro.test(subset(RC_amino_acid_degrad, group2 == unique(RC_amino_acid_degrad$group2)[i])$value)[2]
  print(st.p.value)
}

wilcox_tP2 <- c()
for(i in 1:32){
  wilcox_tP2 <- c(wilcox_tP2, wilcox.test(value ~ group,
                                          subset(RC_amino_acid_degrad, pathway == unique(RC_amino_acid_degrad$pathway)[i], alternative = "two.sided"))$p.value)
}

wilcox_tP2 <- as.numeric(unlist(wilcox_tP2))
wilcox_tP2 <- data.frame("pathway" = unique(RC_amino_acid_degrad$pathway), 
                         "wilcox_testP" = wilcox_tP2,
                         "expressionP" = paste0("P = ", sprintf("%.3f", wilcox_tP2)))

wilcox_tP2
wilcox_tP2$test_adjP <- p.adjust(wilcox_tP2$wilcox_testP, "fdr")
wilcox_tP2$expression_adjP <- paste0("italic(P[adj]) == ", sprintf("%.3f", wilcox_tP2$test_adjP))
signif_wilcoxP2 <- RC_amino_acid_degrad$pathway %in% wilcox_tP2$pathway[wilcox_tP2$test_adjP < 0.05]



RC_amino_acid_degrad$signif <- NA
RC_amino_acid_degrad[signif_wilcoxP2,"signif"] <- "*"

subset(RC_amino_acid_degrad, signif=="*")

amino_acid_degradation_P <- 
  ggplot(RC_amino_acid_degrad, aes(x=group, y=value, fill = group))+
  geom_rect(data = subset(RC_amino_acid_degrad, signif=="*"), fill="yellow", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha = 0.004)+
  geom_boxplot()+
  facet_wrap(c("pathway"), scales="free_y")+
  geom_nodetext(data=wilcox_tP2, aes(label = expression_adjP, x=1.5, y=Inf, fill=NULL), vjust = 1.5, size = 3, parse = T)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size=1)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1)+
  xlab("")+
  ylab("PICRUSt2 predicted pathway abundances")+
  theme_light()+
  theme(legend.position = "none",
        strip.background = element_rect(fill="grey60"))+
  ggtitle("Amino acid degradation pathway abundance",
          subtitle = "before turnover (cycling 1 ~3, n=12), after turnover (cycling 4 ~ 10, n=28), unpaired two-samples Wilcoxon test")

amino_acid_degradation_P


statformat <- function(val,z){
  sub("^(-?)0.", "\\1.", sprintf(paste("%.",z,"f", sep = ""), val))
}

#========
# PCoA
#===========
# endosphere Bray-Curtis
endosphere_OTU_rel <- t(OTU_relative[,1:40]) # relative abundance of endosphere OTUs 

endosphere_bray <- as.matrix(vegdist(endosphere_OTU_rel, method="bray"))

endosphere_pcoa_bray <-cmdscale(endosphere_bray, eig=TRUE)
colnames(endosphere_pcoa_bray$points) <- c("PCoA1", "PCoA2")
endosphere_pcoa_bray$points <- cbind(endosphere_pcoa_bray$points, SAM[1:40,])

endosphere_pcoa_bray_mean <- endosphere_pcoa_bray$points %>%
  group_by(Cycling) %>%
  summarise(mean_PCoA1 = mean(PCoA1, na.rm=T),mean_PCoA2 = mean(PCoA2, na.rm=T))

endosphere_pcoa_bray$points <- merge(endosphere_pcoa_bray$points, endosphere_pcoa_bray_mean, by="Cycling")

Endosphere_PCoA_bray_p <- 
  ggplot(endosphere_pcoa_bray$point, aes(PCoA1, PCoA2, color=as.factor(Cycling)))+
  geom_path(aes(x=mean_PCoA1, y=mean_PCoA2, color=NULL), alpha=0.2, size=5)+
  geom_encircle(aes(fill=as.factor(Cycling)), expand = 0, s_shape=1, alpha =0.1)+
  geom_segment(aes(xend=mean_PCoA1, yend=mean_PCoA2))+
  geom_point(size=3)+
  geom_nodes(aes(x=mean_PCoA1, y=mean_PCoA2), size = 6 ,shape = 21, fill="white", stroke= 1)+
  geom_nodetext(aes(x=mean_PCoA1, y=mean_PCoA2, label=Cycling), color="black")+
  labs(x=paste0("PCoA 1 (",
                round(endosphere_pcoa_bray$eig[1]/sum(endosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       y=paste0("PCoA 2 (",
                round(endosphere_pcoa_bray$eig[2]/sum(endosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       color="Cycling"
  )+
  scale_x_continuous(breaks = -5:5*0.1, limits=c(-0.4,0.5))+
  scale_y_continuous(breaks = -5:5*0.1, limits=c(-0.5,0.4))+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  guides(fill = FALSE)+
  coord_fixed()+
  theme_light()
Endosphere_PCoA_bray_p

# endosphere weighted score
endosphere_pcoa_bray_wascore <- as.data.frame(wascores(endosphere_pcoa_bray$points[,2:3] , t(OTU_relative[,41:80])))

endosphere_pcoa_bray_wascore$keystone <- "Not influencer OTUs"
endosphere_pcoa_bray_wascore$keystone[which(OTU_tax_table$family %in% keystonetaxa)] <- "Influencer OTU"

endosphere_pcoa_bray_wascore$average_abundance <- rowMeans(OTU_relative[,41:80])
endosphere_pcoa_bray_wascore$average_abundance2 <- "Less than 0.1%" 
endosphere_pcoa_bray_wascore$average_abundance2[endosphere_pcoa_bray_wascore$average_abundance > 0.001] <- "More than 0.1%" 

endosphere_pcoa_bray_wascore <- endosphere_pcoa_bray_wascore[!is.nan(endosphere_pcoa_bray_wascore$PCoA1),]

endosphere_weighted_scores_p <- 
  ggplot(endosphere_pcoa_bray_wascore, aes(x=PCoA1, y=PCoA2, color=keystone))+
  geom_point()+
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks = -4:4*0.1, limits=c(-0.4,0.5))+
  scale_y_continuous(breaks = -4:4*0.1, limits=c(-0.5,0.4))+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  theme_light()+
  theme(legend.position="bottom")
ggMarginal(endosphere_weighted_scores_p, groupFill = T)

endosphere_weighted_scores_over0.001_p <- 
  ggplot(subset(endosphere_pcoa_bray_wascore, average_abundance2 == "More than 0.1%"), 
         aes(x=PCoA1, y=PCoA2, color=keystone))+
  geom_point()+
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks = -4:4*0.1, limits=c(-0.4,0.5))+
  scale_y_continuous(breaks = -4:4*0.1, limits=c(-0.5,0.4))+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  labs(x=paste0("PCoA 1 (",
                round(endosphere_pcoa_bray$eig[1]/sum(endosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       y=paste0("PCoA 2 (",
                round(endosphere_pcoa_bray$eig[2]/sum(endosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       color = "OTUs")+
  theme_light()+
  theme(legend.position="bottom")
ggMarginal(endosphere_weighted_scores_over0.001_p, groupFill = T)


# rhziosphere Bray-Curtis
rhizosphere_OTU_rel <- t(OTU_relative[,41:80]) # relative abundance of rhizosphere OTUs 

rhizosphere_bray <- as.matrix(vegdist(rhizosphere_OTU_rel, method="bray"))

rhizosphere_pcoa_bray <-cmdscale(rhizosphere_bray, eig=TRUE)
colnames(rhizosphere_pcoa_bray$points) <- c("PCoA1", "PCoA2")
rhizosphere_pcoa_bray$points <- cbind(rhizosphere_pcoa_bray$points, SAM[41:80,])

rhizosphere_pcoa_bray_mean <- rhizosphere_pcoa_bray$points %>%
  group_by(Cycling) %>%
  summarise(mean_PCoA1 = mean(PCoA1, na.rm=T),mean_PCoA2 = mean(PCoA2, na.rm=T))

rhizosphere_pcoa_bray$points <- merge(rhizosphere_pcoa_bray$points, rhizosphere_pcoa_bray_mean, by="Cycling")

Rhizosphere_PCoA_bray_p <- 
  ggplot(rhizosphere_pcoa_bray$point, aes(PCoA1, PCoA2, color=as.factor(Cycling)))+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  geom_path(aes(x=mean_PCoA1, y=mean_PCoA2, color=NULL), alpha=0.2, size=5)+
  geom_encircle(aes(fill=as.factor(Cycling)), expand = 0, s_shape=1, alpha =0.1)+
  geom_segment(aes(xend=mean_PCoA1, yend=mean_PCoA2))+
  geom_point(size=3)+
  geom_nodes(aes(x=mean_PCoA1, y=mean_PCoA2), size = 6 ,shape = 21, fill="white", stroke= 1)+
  geom_nodetext(aes(x=mean_PCoA1, y=mean_PCoA2, label=Cycling), color="black")+
  labs(x=paste0("PCoA 1 (",
               round(rhizosphere_pcoa_bray$eig[1]/sum(rhizosphere_pcoa_bray$eig), 3)*100,
               "%)"),
       y=paste0("PCoA 2 (",
                round(rhizosphere_pcoa_bray$eig[2]/sum(rhizosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       color="Cycling"
       )+
  guides(fill = FALSE)+
  scale_x_continuous(breaks = -4:4*0.1, limits=c(-0.5,0.4))+
  scale_y_continuous(breaks = -4:4*0.1, limits=c(-0.4,0.5))+
  coord_fixed(ratio=1)+
  theme_light()
Rhizosphere_PCoA_bray_p

#rhizosphere weighted score--------
rhizosphere_pcoa_bray_wascore <- as.data.frame(wascores(rhizosphere_pcoa_bray$points[,2:3] , t(OTU_relative[,41:80])))

rhizosphere_pcoa_bray_wascore$keystone <- "Not influencer OTU"
rhizosphere_pcoa_bray_wascore$keystone[which(OTU_tax_table$family %in% keystonetaxa)] <- "Influencer OTU"

rhizosphere_pcoa_bray_wascore$average_abundance <- rowMeans(OTU_relative[,41:80])
rhizosphere_pcoa_bray_wascore$average_abundance2 <- "Less than 0.1%" 
rhizosphere_pcoa_bray_wascore$average_abundance2[rhizosphere_pcoa_bray_wascore$average_abundance > 0.001] <- "More than 0.1%" 

rhizosphere_pcoa_bray_wascore <- rhizosphere_pcoa_bray_wascore[!is.nan(rhizosphere_pcoa_bray_wascore$PCoA1),]

rhizosphere_weighted_scores_p <- 
  ggplot(rhizosphere_pcoa_bray_wascore, aes(x=PCoA1, y=PCoA2, color=keystone))+
  geom_point()+
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks = -4:4*0.1, limits=c(-0.5,0.4))+
  scale_y_continuous(breaks = -4:4*0.1, limits=c(-0.4,0.5))+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  theme_light()+
  theme(legend.position="bottom")
ggMarginal(rhizosphere_weighted_scores_p, groupFill = T)

rhizosphere_weighted_scores_over0.001_p <- 
  ggplot(subset(rhizosphere_pcoa_bray_wascore, average_abundance2 == "More than 0.1%"), 
         aes(x=PCoA1, y=PCoA2, color=keystone))+
  geom_point()+
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks = -4:4*0.1, limits=c(-0.5,0.4))+
  scale_y_continuous(breaks = -4:4*0.1, limits=c(-0.4,0.5))+
  annotate(geom="segment", y=0, yend = 0, x=Inf, xend=-Inf, size=0.3, alpha=0.5)+
  annotate(geom="segment", x=0, xend = 0, y=Inf, yend=-Inf, size=0.3, alpha=0.5)+
  labs(x=paste0("PCoA 1 (",
                round(rhizosphere_pcoa_bray$eig[1]/sum(rhizosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       y=paste0("PCoA 2 (",
                round(rhizosphere_pcoa_bray$eig[2]/sum(rhizosphere_pcoa_bray$eig), 3)*100,
                "%)"),
       color = "OTUs")+
  theme_light()+
  theme(legend.position="bottom")
ggMarginal(rhizosphere_weighted_scores_over0.001_p, groupFill = T)


ggarrange(Endosphere_PCoA_bray_p+
               theme(legend.position = "bottom"),
             Rhizosphere_PCoA_bray_p+
               theme(legend.position = "bottom"),
             ggMarginal(endosphere_weighted_scores_over0.001_p+
                          annotate(geom="text", x=-Inf, y=Inf, label=" Endosphere,\n weighted score of OTUs\nmore than 0.1% abundance", 
                                   vjust=1.2, hjust=0),
                          groupFill = T
                        ),
             ggMarginal(rhizosphere_weighted_scores_over0.001_p+
                          annotate(geom="text", x=-Inf, y=Inf, label=" Rhizosphere,\n weighted score of OTUs\nmore than 0.1% abundance",
                                   vjust=1.2, hjust=0),
                          ggtitle("D"), 
                        groupFill = T
             ),
             labels = c("A", "B", "C", "D"))



keystonetaxa_abundance <- melt(colSums(OTU_relative[which(OTU_tax_table$family %in% keystonetaxa),1:80]))
colnames(keystonetaxa_abundance) <- "abundance"

keystonetaxa_abundance$cycling <- gsub("[ER]A", "",rownames(keystonetaxa_abundance))
keystonetaxa_abundance$cycling <- as.numeric(gsub("C[0-9]", "", keystonetaxa_abundance$cycling))

keystonetaxa_abundance$group <- ifelse(keystonetaxa_abundance$cycling %in% 1:3, "1 ~ 3", "4 ~ 10")

keystonetaxa_abundance$PotID <- gsub("[ER]A[0-9]{1,}", "", rownames(keystonetaxa_abundance))

keystonetaxa_abundance$Site <- ifelse(rownames(keystonetaxa_abundance) %in% sample_ord[1:40], "Endosphere", "Rhizosphere")

keystonetaxa_abundance$DI <- unlist(SAM[1:80,8])


shapiro.test(subset(keystonetaxa_abundance, group=="1 ~ 3" & Site =="Endosphere")$abundance)
shapiro.test(subset(keystonetaxa_abundance, group=="4 ~ 10" & Site =="Endosphere")$abundance)

shapiro.test(subset(keystonetaxa_abundance, group=="1 ~ 3" & Site =="Rhizosphere")$abundance)
shapiro.test(subset(keystonetaxa_abundance, group=="4 ~ 10" & Site =="Rhizosphere")$abundance)

bartlett.test(abundance ~ group, subset(keystonetaxa_abundance, Site=="Rhizosphere"))
bartlett.test(abundance ~ group, subset(keystonetaxa_abundance, Site=="Endosphere"))

# keystonetaxa abundance, turnover group anova
Rhizosphere_turnover_group_keystone_t.test <- t.test(abundance ~ group, subset(keystonetaxa_abundance, Site=="Rhizosphere"), var.equal = T)
Rhizosphere_turnover_group_keystone_t.test

Endosphere_turnover_group_keystone_t.test <- t.test(abundance ~ group, subset(keystonetaxa_abundance, Site=="Endosphere"), var.equal = T)
Endosphere_turnover_group_keystone_t.test

# DI turnover group DI welch's anova
shapiro.test(subset(keystonetaxa_abundance, group=="1 ~ 3" & Site =="Endosphere")$DI) #정규성 만족
shapiro.test(subset(keystonetaxa_abundance, group=="4 ~ 10" & Site =="Endosphere")$DI) #정규성 만족 못함

bartlett.test(DI ~ group, subset(keystonetaxa_abundance, Site=="Endosphere")) #등분산 만족 못함

turnover_group_DI_wilcox.test <- wilcox.test(DI ~ group, subset(keystonetaxa_abundance, Site=="Endosphere"))
turnover_group_DI_wilcox.test

#
endosphere_keystonetaxa_abundance_turnover_group <- 
  ggplot(subset(keystonetaxa_abundance, Site == "Endosphere"), aes(group, abundance, group=group))+ 
  stat_summary(geom = "bar", position= position_dodge(.8), alpha=0.5, width=0.8, fun.data = mean_sd)+
  stat_summary(geom = "errorbar", position= position_dodge(.8), width=0.4, fun.data = mean_sd)+
  geom_jitter(shape=1, stroke=1.2, alpha=0.5)+
  facet_wrap(Site ~ ., scales = "free")+
  theme_light()+
  annotate(geom= "text", x=-Inf, y=0.9, label=expression(italic(P) == 2.59%*%10^-2), hjust=-0.1)+
  annotate(geom= "segment", x=1, xend=2, y=0.85, yend=0.85, size=0.5)+
  annotate(geom= "segment", x=1, xend=1, y=0.85, yend=0.83, size=0.5)+
  annotate(geom= "segment", x=2, xend=2, y=0.85, yend=0.83, size=0.5)+
  annotate(geom= "text", x=1.5, y=0.85, label= "*", size=7)+
  annotate(geom= "segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate(geom= "segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1)+
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.93), expand=c(0,0))+
  ylab("Relative abundance\nof the influencer groups")+
  xlab("")+
  theme(strip.background = element_rect(fill="grey90"), strip.text = element_text(color="black"))

endosphere_keystonetaxa_abundance_turnover_group

rhizosphere_keystonetaxa_abundance_turnover_group <- 
  ggplot(subset(keystonetaxa_abundance, Site == "Rhizosphere"), aes(group, abundance, group=group))+ 
  stat_summary(geom = "bar", position= position_dodge(.8), alpha=0.5, width=0.8, fun.data = mean_sd)+
  stat_summary(geom = "errorbar", position= position_dodge(.8), width=0.4, fun.data = mean_sd)+
  geom_jitter(shape=1, stroke=1.2, alpha=0.5)+
  facet_wrap(Site ~ ., scales = "free")+
  theme_light()+
  annotate(geom= "text", x=-Inf, y=0.9, label=expression(italic(P) == 8.91%*%10^-8), hjust=-0.1)+
  annotate(geom= "segment", x=1, xend=2, y=0.55, yend=0.55, size=0.5)+
  annotate(geom= "segment", x=1, xend=1, y=0.55, yend=0.53, size=0.5)+
  annotate(geom= "segment", x=2, xend=2, y=0.55, yend=0.53, size=0.5)+
  annotate(geom= "text", x=1.5, y=0.55, label= "***", size=7)+
  annotate(geom= "segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate(geom= "segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1)+
  scale_y_continuous(labels = scales::percent, expand=c(0,0), limits = c(0, 0.93))+
  ylab("Relative abundance\nof the influencer groups")+
  xlab("")+
  theme(strip.background = element_rect(fill="grey90"), strip.text = element_text(color="black"))

rhizosphere_keystonetaxa_abundance_turnover_group

DI_turnover_group <- 
  ggplot(subset(keystonetaxa_abundance, Site =="Endosphere"), aes(group, DI, group=group))+
  stat_summary(geom = "bar", position= position_dodge(.8), alpha=0.5, width=0.8)+
  stat_summary(geom = "errorbar", position= position_dodge(.8), width=0.4)+
  geom_jitter(shape=1, stroke=1.2, alpha=0.5)+
  theme_light()+
  annotate(geom= "text", x=-Inf, y=75, label=expression(italic(P) == 2.25%*%10^-3), hjust=-0.1)+
  annotate(geom= "segment", x=1, xend=2, y=60, yend=60, size=0.5)+
  annotate(geom= "segment", x=1, xend=1, y=60, yend=57, size=0.5)+
  annotate(geom= "segment", x=2, xend=2, y=60, yend=57, size=0.5)+
  annotate(geom= "text", x=1.5, y=60, label= "**", size=7)+
  annotate(geom= "segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate(geom= "segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size=1)+
  ylab("Disease index average")+
  xlab("")+
  scale_y_continuous(limits = c(0,80),expand=c(0,0))
DI_turnover_group

ggarrange(DI_turnover_group, endosphere_keystonetaxa_abundance_turnover_group,
          rhizosphere_keystonetaxa_abundance_turnover_group, nrow =1, labels = c("A", "B" ,"C"))


rhizosphere_kestonegroup_abundance_and_nitrofixation <- cbind(keystonetaxa_abundance[41:80,], t(RC_KO[c("K02584","K02585","K02586","K02587","K02588","K02591","K02592","K02593","K02594","K02595","K02596","K02597"),]))
rhizosphere_kestonegroup_abundance_and_nitrofixation$Nifgenecluster <- rowSums(rhizosphere_kestonegroup_abundance_and_nitrofixation[,7:18])

for_correlation <- rhizosphere_kestonegroup_abundance_and_nitrofixation[,c(1,6,19)]
colnames(for_correlation) <- c("The influencers\nrelative abundance (%)", "Average of disease index", "Nif gene cluster prediction")
for_correlation <- cbind(for_correlation, SAM[41:80,])
for_correlation$`The influencers\nrelative abundance (%)` <- for_correlation$`The influencers\nrelative abundance (%)`*100

key_vs_DI<- ggscatter(for_correlation, x="The influencers\nrelative abundance (%)", 
                      y="Average of disease index",
                      add = "reg.line",
                      add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  #stat_cor(method = "pearson", label.x = 0.1, label.y = 65)
  annotate(geom = "richtext", x=15, y = 65, label="*&rho;* = -0.52, *P* = 0.00057", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)+
  ylab("Disease\nindex")
key_vs_DI

key_vs_nif <- ggscatter(for_correlation, x="The influencers\nrelative abundance (%)", 
                        y="Nif gene cluster prediction",
                        add = "reg.line",
                        add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  #stat_cor(method = "pearson", label.x = 0.1, label.y = 23000)
  annotate(geom = "richtext", x=15, y = 23000, label="*&rho;* = 0.6, *P* = 0.00005", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)
  
key_vs_nif 

DI_vs_nif <- ggscatter(for_correlation, x="average of disease index", 
                       y="Nif gene cluster prediction",
                       add = "reg.line",
                       add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
  xlab("Disease\nindex")+
  #stat_cor(method = "pearson", label.x = 15, label.y = 23000)
  annotate(geom = "richtext", x=20, y = 23000, label="*&rho;* = -0.5, *P* = 0.001", label.padding = grid::unit(rep(0, 4), "pt"), fill = NA, label.color = NA)
DI_vs_nif

phenotype_keystone_spearman <- ggarrange(key_vs_DI, key_vs_nif, DI_vs_nif, labels = c("A","B","C"), nrow = 1)
phenotype_keystone_spearman

keystonetaxa_family_abundance <- data.frame("keystone_family"=NA, "rel_abundance"=NA)
for(i in 1:length(keystonetaxa)){
a <- melt(colSums(OTU_relative[which(OTU_tax_table$family %in% keystonetaxa[i]),41:80]))
a$keystone_family <- keystonetaxa[i]
colnames(a) <- c("rel_abundance","keystone_family") 
keystonetaxa_family_abundance <- rbind(keystonetaxa_family_abundance, a[,2:1])
}

keystonetaxa_family_abundance <- keystonetaxa_family_abundance[-1,]

keystonetaxa_family_abundance$sample <- sample_ord[41:80]
rownames(keystonetaxa_family_abundance) <- NULL
keystonetaxa_family_abundance <- keystonetaxa_family_abundance[,c(3,1,2)]

keystonetaxa_family_abundance$DI <- keystonetaxa_abundance$DI[41:80]
keystonetaxa_family_abundance$rel_abundance <- keystonetaxa_family_abundance$rel_abundance*100
keystone_info
DI_vs_each_key_pearson <- data.frame(keystone_family=c("Beijerinckiaceae", "Comamonadaceae","Devosiaceae","Rhizobiaceae","Sphingobacteriaceae", "Sphingomonadaceae","Xanthomonadaceae"),
                                     rho = c("*&rho;* = -0.33", "*&rho;* = -0.54","*&rho;* = -0.36","*&rho;* = -0.46","*&rho;* = -0.32","*&rho;* = -0.47","*&rho;* = -0.52"),
                                     P=c("*P* = 0.035", "*P* = 0.000", "*P* = 0.021", "*P* = 0.000", "*P* = 0.048", "*P* = 0.000", "*P* = 0.001"))


DI_vs_each_key <- ggscatter(keystonetaxa_family_abundance, x="rel_abundance", 
                           y="DI",
                           add = "reg.line",
                           add.params = list(color = "blue", fill = "lightgray"),conf.int = F)+
#  stat_cor(method = "pearson", label.x = Inf, label.y = 65, hjust=1)+
  geom_richtext(data=DI_vs_each_key_pearson, aes(x=Inf, y=65, label=rho), fill = NA, label.color = NA, label.padding = grid::unit(0, "pt"), hjust=1.1)+
  geom_richtext(data=DI_vs_each_key_pearson, aes(x=Inf, y=57, label=P), fill = NA, label.color = NA, label.padding = grid::unit(0, "pt"), hjust=1.1)+
  facet_wrap(.~keystone_family, scales="free", ncol =4)+
  ylab("Disease index")+
  xlab("Relative abundance (%)")

DI_vs_each_key


#nif test

nif_ab <- read.csv("Nfb/Nfb_medium_ab602.csv")
nif_ab$X602_Absorbance <- nif_ab$X602_Absorbance - mean(nif_ab$X602_Absorbance[1:4])

order <- unique(nif_ab$Treatment)


nif_ab$Treatment <- factor(nif_ab$Treatment, levels = order)

is.data.frame(nif_ab)

require(ggplot2)

for(i in 1:8){
  print(shapiro.test((nif_ab$X602_Absorbance[which(nif_ab$Treatment == unique(nif_ab$Treatment)[i])])))
} 

bartlett.test(X602_Absorbance ~ Treatment, nif_ab) # equal var

kruskal.test(data=nif_ab, X602_Absorbance ~ Treatment)
nif_tukey <- nparcomp(X602_Absorbance ~ Treatment, data=nif_ab, type="Tukey")

nif_tukey_p <- nif_tukey$Analysis

nif_tukey_p <- nif_tukey_p[which(nif_tukey_p$p.Value < 0.05),]
nif_tukey_p$siglabel <- "***"

nif_tukey_p$x1 <- gsub("^p[(] ","", nif_tukey_p$Comparison)
nif_tukey_p$x1 <- gsub(" , [A-Za-z0-9 ]{1,}[)]","", nif_tukey_p$x1)
nif_tukey_p$x1 <- factor(nif_tukey_p$x1, levels=unique(nif_ab$Treatment))


nif_tukey_p$x2 <- gsub("^p[()] [A-Za-z0-9 ]{1,} , ","", nif_tukey_p$Comparison)
nif_tukey_p$x2 <- gsub(" [)]","", nif_tukey_p$x2)
nif_tukey_p$x2 <- factor(nif_tukey_p$x2, levels=unique(nif_ab$Treatment))

nif_tukey_p$y <- 1:18*0.03+0.3

nitrogen_fixation_p <- ggplot(data=nif_ab, mapping=aes(x=Treatment, y=X602_Absorbance, fill=Treatment, color=Treatment))+
  stat_summary(geom="bar", fun.data = mean_se, alpha=0.5)+
  stat_summary(geom="errorbar", fun.data = mean_se, alpha=0.5, width=0.5)+
  geom_jitter()+
  geom_segment(data=nif_tukey_p[1:4,], aes(x=x1, xend=x2, y=y, yend=y, fill=NULL, color=NULL))+
  geom_segment(data=nif_tukey_p[1:4,], aes(x=x1, xend=x1, y=y, yend=y-0.01, fill=NULL, color=NULL))+
  geom_segment(data=nif_tukey_p[1:4,], aes(x=x2, xend=x2, y=y, yend=y-0.01, fill=NULL, color=NULL))+
  geom_text(data=nif_tukey_p[1:4,], aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y+0.002, label=siglabel, fill=NULL, color=NULL))+
  annotate(geom="text", x= 1, y=0.25, label=expression(paste("Kruskal-Wallis test, ", italic(P)==3.713%*%10^4)), hjust=0)+
  ylab("602 nm absorbance")+
  xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none")

photo_path <- "photo_image/nif_test_1_month(none).png"

keystone_nif_test_strain <- data.frame("treatment"=c("Control", "WR72a", "Gsoil 3165", "N7R2", "Gsoil 034", "N10R5", "BO184", "Gsoil 124"))
keystone_nif_test_strain$treatment <- factor(keystone_nif_test_strain$treatment, levels = keystone_nif_test_strain$treatment)

keystone_Nfb_test_label <- ggplot(keystone_nif_test_strain, aes(treatment))+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_x_discrete(expand = c(0.065,0))+
  xlab("")


nif_photo <- ggdraw() +
  draw_image(file.path(getwd(),photo_path))+
  draw_plot(keystone_Nfb_test_label)+
  coord_fixed(ratio=0.33)

nif_photo

ggarrange(nitrogen_fixation_p, nif_photo, ncol=1, heights=c(4,2), labels = c("A", "B"))

keystone_info <- read.csv("phenotype_recode/keystone_info.csv", header = F)
keystone_info

keystone_info_ft <- flextable(keystone_info[-1:-2,])
keystone_info_ft <- delete_part(keystone_info_ft, part = "header")
keystone_info_ft <- hline_top(keystone_info_ft, j = NULL, border = NULL, part = "body")
keystone_info_ft <- add_header_row(keystone_info_ft, top = TRUE, values = keystone_info[2,])
keystone_info_ft <- add_header_row(keystone_info_ft, top = TRUE, values = keystone_info[1,])

keystone_info_ft <- merge_at(keystone_info_ft, i=1, j = 5:6, part="header")
keystone_info_ft <- merge_at(keystone_info_ft, i=1, j = 7:8, part="header")
keystone_info_ft <- merge_at(keystone_info_ft, i=1:2, j = 1, part="header")
keystone_info_ft <- merge_at(keystone_info_ft, i=1:2, j = 2, part="header")
keystone_info_ft <- merge_at(keystone_info_ft, i=1:2, j = 3, part="header")
keystone_info_ft <- merge_at(keystone_info_ft, i=1:2, j = 4, part="header")
keystone_info_ft <- hline(keystone_info_ft, i = 1, j = 5:6, border = NULL, part = "header")
keystone_info_ft <- hline(keystone_info_ft, i = 1, j = 7:8, border = NULL, part = "header")
keystone_info_ft <- hline(keystone_info_ft, i = 2, j = 1:8, border = NULL, part = "header")
keystone_info_ft <- hline_top(keystone_info_ft, j = NULL, border = fp_border(color="grey40", width = 2), part = "header")%>%
  align(i = NULL, j = NULL, align = "center", part = "header") %>%
  align(i = NULL, j = 3, align = "center", part="body") %>%
  align(i = c(3,5), j = 4, align = "center", part="body")

keystone_info_ft%>%autofit()%>%colformat_md(j=1)