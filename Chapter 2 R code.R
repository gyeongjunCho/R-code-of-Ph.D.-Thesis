Packages <- c("dada2", "DECIPHER", "doParallel", "ggplot2", "phyloseq", 
              "iNEXT", "dplyr", "RColorBrewer", "picante", "reshape2", 
              "ggnetwork", "network", "sna", "keras", "tidyr", "tidygraph", 
              "entropart", "gtable", "grid", "RColorBrewer", "picante", 
              "ggExtra", "dplyr", "ggrepel", "DESeq2", "DT", "gridExtra", "agricolae", 
              "progress", "ecodist","ggalt", "ggpubr", "userfriendlyscience", "compiler", "nparcomp")
attachRequiredPackages <-  function() {lapply(Packages, 
                                              FUN = function(Packages){
                                                do.call("require", list(Packages))
                                              })
}
attachRequiredPackages()
registerDoParallel(cores = 64)

require(DunnettTests)

#========
#root DI
#=======
DI <- read.csv("Ginseng DI.csv")
DI$DI_mean <-  (DI$DI0*0 + DI$DI1*1 + DI$DI2*5 + DI$DI3*30 + DI$DI4*75)/rowSums(DI[,-(1:2)])

DI$Treatment <- factor(DI$Treatment, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

#정규성(p > 0.05면 정규성)
for(i in 1:6){
  print(
    shapiro.test(subset(DI, Treatment==c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val")[i])$DI_mean)
  )
}
  
#등분산 검정 (p > 0.05면 등분산)
bartlett.test(DI_mean ~ Treatment, data= DI)

#등분산만족 시 anova
root_DI_aov <- aov(DI_mean ~ Treatment, data= DI)
summary(root_DI_aov)

summary_root_DI_aov <- unlist(summary(root_DI_aov))
root_DI_aov_p <- as.numeric(as.character(summary_root_DI_aov[9]))
root_DI_aov_p <- signif(root_DI_aov_p,3)
root_DI_aov_p

TukeyHSD_root_DI <- TukeyHSD(root_DI_aov)
TukeyHSD_root_DI <- as.data.frame(TukeyHSD_root_DI$Treatment)

Treatment_compare_group <- strsplit(rownames(TukeyHSD_root_DI), split = "-")
Treatment_compare_group <- unlist(Treatment_compare_group)

TukeyHSD_root_DI <- cbind(TukeyHSD_root_DI, compare_group1 = Treatment_compare_group[(1:(length(Treatment_compare_group)/2)*2)-1])
TukeyHSD_root_DI <- cbind(TukeyHSD_root_DI, compare_group2 = Treatment_compare_group[1:(length(Treatment_compare_group)/2)*2])
TukeyHSD_root_DI$compare_group1 <- factor(TukeyHSD_root_DI$compare_group1, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
TukeyHSD_root_DI$compare_group2 <- factor(TukeyHSD_root_DI$compare_group2, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

TukeyHSD_root_DI$p_sig_x <- (as.numeric(TukeyHSD_root_DI$compare_group1) +
                             as.numeric(TukeyHSD_root_DI$compare_group2))/2

TukeyHSD_root_DI$p_sig <- 
  ifelse(TukeyHSD_root_DI$`p adj` < 0.05, no = "ns", 
         yes = ifelse(TukeyHSD_root_DI$`p adj` < 0.01, no = "*",
                      yes = ifelse(TukeyHSD_root_DI$`p adj` < 0.001, no = "**", yes= "***")))

TukeyHSD_root_DI <- subset(TukeyHSD_root_DI, p_sig != "ns")

TukeyHSD_root_DI


plot_root_DI <- ggplot(data = DI, mapping = aes(x = Treatment, y = DI_mean, color = Treatment, fill= Treatment))+
  stat_summary(geom = "bar", fun.data = "mean_se", alpha=0.5)+
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5, size=1)+
  geom_jitter(alpha=0.3, shape=21)+
  annotate(geom="text", label= expression(paste("ANOVA, ", ~italic(P) == 0.013, "; Tukey's HSD")), x=2.5, y= 20)+
  annotate(geom="segment", 
           x = TukeyHSD_root_DI$compare_group1,
           xend = TukeyHSD_root_DI$compare_group2,
           y = 1:length(TukeyHSD_root_DI$compare_group1)+20.5,
           yend = 1:length(TukeyHSD_root_DI$compare_group1)+20.5)+
  annotate(geom = "segment", 
           x = TukeyHSD_root_DI$compare_group1,
           xend = TukeyHSD_root_DI$compare_group1,
           y = 1:length(TukeyHSD_root_DI$compare_group1)+20.5,
           yend = 1:length(TukeyHSD_root_DI$compare_group1)+ 20
  )+
  annotate(geom = "segment", 
           x = TukeyHSD_root_DI$compare_group2,
           xend = TukeyHSD_root_DI$compare_group2,
           y = 1:length(TukeyHSD_root_DI$compare_group1)+20.5,
           yend = 1:length(TukeyHSD_root_DI$compare_group1)+ 20
  )+
  annotate(geom = "text", label = TukeyHSD_root_DI$p_sig,
           x = TukeyHSD_root_DI$p_sig_x,
           y = 1:length(TukeyHSD_root_DI$compare_group1)+20.5)+
#  scale_y_continuous(limits = c(0,3), expand = c(0,0))+
  theme_light()+
  ylab("Disease index")+
  guides(fill="none", color="none")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size= 1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size= 1.5)

plot_root_DI

#===============
# stem phenotype
#===============
# length
#등분산 검정 (p > 0.05면 등분산)
bartlett.test(stem_length ~ Treatment, data= ginseng_stem)

ginseng_stem <- read.csv("phenotype_recode/amino acid treatment stem phenotype.csv")
ginseng_stem$Treatment<- factor(ginseng_stem$Treatment, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))


stem_length_kruskal <- aov(stem_length ~ Treatment, data= ginseng_stem)

summary_stem_length_aov <- unlist(summary(stem_length_aov))
stem_length_aov_p <- as.numeric(as.character(summary_stem_length_aov[9]))
stem_length_aov_p <- signif(stem_length_aov_p,3)

TukeyHSD_stem_length <- TukeyHSD(stem_length_aov)
TukeyHSD_stem_length <- as.data.frame(TukeyHSD_stem_length$Treatment)

Treatment_compare_group <- strsplit(rownames(TukeyHSD_stem_length), split = "-")
Treatment_compare_group <- unlist(Treatment_compare_group)

TukeyHSD_stem_length <- cbind(TukeyHSD_stem_length, compare_group1 = Treatment_compare_group[(1:(length(Treatment_compare_group)/2)*2)-1])
TukeyHSD_stem_length <- cbind(TukeyHSD_stem_length, compare_group2 = Treatment_compare_group[1:(length(Treatment_compare_group)/2)*2])
TukeyHSD_stem_length$compare_group1 <- factor(TukeyHSD_stem_length$compare_group1, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
TukeyHSD_stem_length$compare_group2 <- factor(TukeyHSD_stem_length$compare_group2, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
TukeyHSD_stem_length$p_sig_x <- (as.numeric(TukeyHSD_stem_length$compare_group1) +
                                 as.numeric(TukeyHSD_stem_length$compare_group2))/2

TukeyHSD_stem_length$p_sig <- 
ifelse(TukeyHSD_stem_length$`p adj` < 0.05, no = "ns", 
       yes = ifelse(TukeyHSD_stem_length$`p adj` < 0.01, no = "*",
                    yes = ifelse(TukeyHSD_stem_length$`p adj` < 0.001, no = "**", yes= "***")))

TukeyHSD_stem_length <- subset(TukeyHSD_stem_length, p_sig != "ns")

plot_stem_length <- ggplot(data = ginseng_stem, mapping = aes(x = Treatment, y = stem_length, color = Treatment, fill= Treatment))+
  stat_summary(geom = "bar", fun.data = "mean_se", alpha=0.5)+
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5, size=1)+
  geom_jitter(alpha=0.3, shape=21)+
  annotate(geom="text", label= expression(paste("ANOVA, ", italic(P) == 0, ".000; Tukey's HSD")), x=3, y= 23)+
  annotate(geom="segment", 
           x = TukeyHSD_stem_length$compare_group1,
           xend = TukeyHSD_stem_length$compare_group2,
           y = 1:length(TukeyHSD_stem_length$compare_group1)+15,
           yend = 1:length(TukeyHSD_stem_length$compare_group1)+15)+
  annotate(geom = "segment", 
           x = TukeyHSD_stem_length$compare_group1,
           xend = TukeyHSD_stem_length$compare_group1,
           y = 1:length(TukeyHSD_stem_length$compare_group1)+15,
           yend = 1:length(TukeyHSD_stem_length$compare_group1)+ 14.7
           )+
  annotate(geom = "segment", 
           x = TukeyHSD_stem_length$compare_group2,
           xend = TukeyHSD_stem_length$compare_group2,
           y = 1:length(TukeyHSD_stem_length$compare_group1)+15,
           yend = 1:length(TukeyHSD_stem_length$compare_group1)+ 14.7
           )+
  annotate(geom = "text", label = TukeyHSD_stem_length$p_sig,
           x = TukeyHSD_stem_length$p_sig_x,
           y = 1:length(TukeyHSD_stem_length$compare_group1)+15)+
  scale_y_continuous(limits = c(0,25), expand = c(0,0))+
  theme_light()+
  ylab("Stem length (cm)")+
  guides(fill="none", color="none")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size= 1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size= 1.5)

plot_stem_length


# weight
#등분산 검정 (p > 0.05면 등분산)
bartlett.test(stem_weight ~ Treatment, data= ginseng_stem)

#등분산을 만족하므로 일반적인 anova

stem_weight_aov <- aov(stem_weight ~ Treatment, data= ginseng_stem)

summary_stem_weight_aov <- unlist(summary(stem_weight_aov))
stem_weight_aov_p <- as.numeric(as.character(summary_stem_weight_aov[9]))
stem_weight_aov_p <- signif(stem_weight_aov_p,3)

TukeyHSD_stem_weight <- TukeyHSD(stem_weight_aov)
TukeyHSD_stem_weight <- as.data.frame(TukeyHSD_stem_weight$Treatment)

Treatment_compare_group <- strsplit(rownames(TukeyHSD_stem_weight), split = "-")
Treatment_compare_group <- unlist(Treatment_compare_group)

TukeyHSD_stem_weight <- cbind(TukeyHSD_stem_weight, compare_group1 = Treatment_compare_group[(1:(length(Treatment_compare_group)/2)*2)-1])
TukeyHSD_stem_weight <- cbind(TukeyHSD_stem_weight, compare_group2 = Treatment_compare_group[1:(length(Treatment_compare_group)/2)*2])
TukeyHSD_stem_weight$compare_group1 <- factor(TukeyHSD_stem_weight$compare_group1, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
TukeyHSD_stem_weight$compare_group2 <- factor(TukeyHSD_stem_weight$compare_group2, levels = c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
TukeyHSD_stem_weight$p_sig_x <- (as.numeric(TukeyHSD_stem_weight$compare_group1) +
                                   as.numeric(TukeyHSD_stem_weight$compare_group2))/2

TukeyHSD_stem_weight$p_sig <- 
  ifelse(TukeyHSD_stem_weight$`p adj` < 0.05, no = "ns", 
         yes = ifelse(TukeyHSD_stem_weight$`p adj` < 0.01, no = "*",
                      yes = ifelse(TukeyHSD_stem_weight$`p adj` < 0.001, no = "**", yes= "***")))

TukeyHSD_stem_weight <- subset(TukeyHSD_stem_weight, p_sig != "ns")


plot_stem_weight <- ggplot(data = ginseng_stem, mapping = aes(x = Treatment, y = stem_weight, color = Treatment, fill= Treatment))+
  stat_summary(geom = "bar", fun.data = "mean_se", alpha=0.5)+
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.5, size=1)+
  geom_jitter(alpha=0.3, shape=21)+
  annotate(geom="text", label= expression(paste("ANOVA, ", italic(P) == "0.000; Tukey's HSD" )), x=3, y= 1.8)+
  annotate(geom="segment",
           x = TukeyHSD_stem_weight$compare_group1,
           xend = TukeyHSD_stem_weight$compare_group2,
           y = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.8,
           yend = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.8)+
  annotate(geom = "segment", 
           x = TukeyHSD_stem_weight$compare_group1,
           xend = TukeyHSD_stem_weight$compare_group1,
           y = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.8,
           yend = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.77
  )+
  annotate(geom = "segment", 
           x = TukeyHSD_stem_weight$compare_group2,
           xend = TukeyHSD_stem_weight$compare_group2,
           y = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.8,
           yend = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.77
  )+
  annotate(geom = "text", label = TukeyHSD_stem_weight$p_sig,
           x = TukeyHSD_stem_weight$p_sig_x,
           y = 1:length(TukeyHSD_stem_weight$compare_group1)/8 + 0.8)+
  scale_y_continuous(limits = c(0,2), expand = c(0,0))+
  theme_light()+
  ylab("Stem weight (g)")+
  guides(fill="none", color="none")+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf, yend=-Inf, size= 1.5)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=Inf, yend=-Inf, size= 1.5)

plot_stem_weight

gridExtra::grid.arrange(plot_root_DI, plot_stem_length, plot_stem_weight, ncol =3)



#==================
# fusarium density
#==================

inputsoil <- c(0.5, 0.5, 0.5, 0.5,
               NA,
               0.53, 0.49, 0.5, 0.51, #bulk 
               0.37, 0.31, 0.24, 0.42, 0.23, #Con
               0.35, 0.40, 0.33, 0.44, 0.41, #ammonium
               0.29, 0.37, 0.41, 0.30, 0.35, #Glu
               0.34, 0.36, 0.29, 0.32, 0.36, #Asp
               0.22, 0.31, 0.29, 0.34, 0.53, #Asn
               0.25, 0.23, 0.34,  0.32, 0.29) #Val


qPCR_sample_info <- read.csv("F.solani_qPCR/F.solani -  Quantification Plate View Results_SYBR.csv")[0:7*4+2, c(-2,-14)]
colnames(qPCR_sample_info) <- c("loc",1:11)


qPCR_sample_info <- melt(qPCR_sample_info, id.vars = "loc", variable_name = "well_num")
qPCR_sample_info <- qPCR_sample_info[which(qPCR_sample_info$value!=""),]

rownames(qPCR_sample_info) <- NULL

qPCR_sample_info$Standard <- F 
qPCR_sample_info$Standard[1:8] <- T


qPCR_sample_info$inputsoil <- NA
qPCR_sample_info$inputsoil[1:39*2] <- inputsoil
qPCR_sample_info$inputsoil[1:39*2-1] <- inputsoil

qPCR_sample_info

qPCR_sample_info$value <- gsub("^R","",qPCR_sample_info$value)

qPCR_sample_info$well_num <- as.character(qPCR_sample_info$variable)
qPCR_sample_info$well_num2  <- qPCR_sample_info$well_num

qPCR_sample_info$well_num2[which(nchar(qPCR_sample_info$well_num)==1)] <- paste0("0", qPCR_sample_info$well_num[which(nchar(qPCR_sample_info$well_num)==1)])

qPCR_sample_info$well_name <- paste0(qPCR_sample_info$loc, qPCR_sample_info$well_num)
qPCR_sample_info$well_name2 <- paste0(qPCR_sample_info$loc, qPCR_sample_info$well_num2)


qPCR_sample_info$well_name

melting_curve <- read.csv("F.solani_qPCR/F.solani -  Melt Curve Derivative Results_SYBR.csv")[,-1]
melting_curve <- melt(melting_curve, id.vars = "Temperature", variable_name = "well_name")

colnames(melting_curve)[2] <- "well_name"
melting_curve <- merge(x = melting_curve, y = qPCR_sample_info, by = "well_name")

melting_curve
melting_curve$treatment <- melting_curve$value.y

melting_curve$treatment <- gsub("n[0-9]", "n", melting_curve$treatment)
melting_curve$treatment <- gsub("u[0-9]", "u", melting_curve$treatment)
melting_curve$treatment <- gsub("p[0-9]", "p", melting_curve$treatment)
melting_curve$treatment <- gsub("l[0-9]", "l", melting_curve$treatment)
melting_curve$treatment <- gsub("C[0-9]", "C", melting_curve$treatment)

melting_curve

melting_plot<- 
  ggplot(melting_curve, aes(x=Temperature, y=value.x, group=well_name, color=treatment))+
  geom_line()+
  theme_classic()+
  ylab("-(d/dT)")

Ct <- read.csv("F.solani_qPCR/F.solani -  Quantification Cq Results_0.csv")
Ct <- Ct[,c(2,8)]

colnames(Ct) <- c("well_name2", "Ct")

Ct <- merge(qPCR_sample_info, Ct, by = "well_name2")

Ct_standard <- subset(Ct, Standard==TRUE)
Ct_standard
Ct_standard$`CFU/soilg` <- c(9.8*10^5, 9.8*10^5, 
                             4.9*10^5, 4.9*10^5, 
                             9.8*10^4, 9.8*10^4, 
                             4.9*10^4, 4.9*10^4) 
standard <- 
  ggplot(Ct_standard[-1,], aes(`CFU/soilg`,Ct))+
  stat_summary(geom="point")+
  stat_summary(geom="line")+
  geom_point(shape=1)+
  stat_summary(geom="errorbar", fun.data = "mean_sdl", width=0.5)+
  ylab(expression(Cycle[threshold]))+
  geom_text(aes(x=10^5.8, y=30, label="R² = 0.998"))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(.x)))+
  xlab(expression(paste("Standard ", (log[10](CFU/g)))))+
  theme_classic()

standard

lm_formular <- lm(Ct_standard$Ct[-1] ~ log10(Ct_standard$`CFU/soilg`[-1]))
lm_formular

summary(lm_formular)
a <-lm_formular$coefficients[2]
b <-lm_formular$coefficients[1]
CFU_g<- function(x){10^((x-b)/a)}

Ct$Ct
Ct$CFU_per_soil <- CFU_g(Ct$Ct)*0.5/Ct$inputsoil


Ct_sample <- subset(Ct, Standard == F & value != "H2O")

Ct_sample

Ct_sample$treatment <- Ct_sample$value
Ct_sample$treatment <- gsub("[0-9]", "", Ct_sample$treatment)
Ct_sample$treatment <- gsub("NHCl", "NH₄Cl", Ct_sample$treatment)
Ct_sample$treatment <- gsub("BC", "Bulk soil", Ct_sample$treatment)
Ct_sample$treatment <- gsub("Con", "Control", Ct_sample$treatment)
Ct_sample$treatment <- factor(Ct_sample$treatment, levels = unique(Ct_sample$treatment))

Ct_sample_summary <- 
  Ct_sample%>%
  group_by(value, treatment)%>%
  summarise(sample_CFU_per_soil= mean(CFU_per_soil), SD= sd(CFU_per_soil))

Ct_sample_summary$log10CFU <- log10(Ct_sample_summary$sample_CFU_per_soil)

shapiro.test(subset(Ct_sample_summary, treatment == "Bulk soil")$log10CFU)
shapiro.test(subset(Ct_sample_summary, treatment == "Control")$log10CFU)
shapiro.test(subset(Ct_sample_summary, treatment == "NH₄Cl")$log10CFU)
shapiro.test(subset(Ct_sample_summary, treatment == "Glu")$log10CFU)
shapiro.test(subset(Ct_sample_summary, treatment == "Asp")$log10CFU)
shapiro.test(subset(Ct_sample_summary, treatment == "Asn")$log10CFU)
shapiro.test(subset(Ct_sample_summary, treatment == "Val")$log10CFU)

bartlett.test(log10CFU~treatment, Ct_sample_summary)

CFU_aov <- aov(log10CFU ~ treatment, Ct_sample_summary)

summarise_CFU_aov<- summary(CFU_aov)

CFU_aov_pvalue <- unlist(summarise_CFU_aov)[9]
CFU_aov_pvalue

CFU_post_hoc <- TukeyHSD(CFU_aov)

CFU_post_hoc

CFU_Tukey_sig <- CFU_post_hoc$treatment[which(CFU_post_hoc$treatment[,4] < 0.05),]
CFU_Tukey_sig <- as.data.frame(CFU_Tukey_sig)
CFU_Tukey_sig$sig_label <- ifelse(CFU_Tukey_sig$`p adj` >= 0.001, ifelse(CFU_Tukey_sig$`p adj` >= 0.01,ifelse(CFU_Tukey_sig$`p adj` >= 0.05,"ns","*"),"**"),"***")

CFU_Tukey_sig

CFU_Tukey_sig$xmin <- row.names(CFU_Tukey_sig)
CFU_Tukey_sig$xmin <- gsub("-.{1,}","",CFU_Tukey_sig$xmin)
CFU_Tukey_sig$xmin <- as.numeric(
  factor(CFU_Tukey_sig$xmin, levels = unique(Ct_sample$treatment))
)

CFU_Tukey_sig$xmax <- row.names(CFU_Tukey_sig)
CFU_Tukey_sig$xmax <- gsub(".{1,}-","",CFU_Tukey_sig$xmax)
CFU_Tukey_sig$xmax <- as.numeric(
  factor(CFU_Tukey_sig$xmax, levels = unique(Ct_sample$treatment))
)


CFU_Tukey_sig$y<- c(10^(8+-1:4*0.2))
CFU_Tukey_sig$yend <- c(10^(8+-1:4*0.2-0.05))

Fusarium_density_plot <-        
  Ct_sample_summary %>%
  ggplot(aes(treatment,sample_CFU_per_soil))+
  stat_summary(geom="bar", fun.data = mean_sdl)+
  stat_summary(geom="errorbar", fun.data = mean_sdl, width=0.5)+
  geom_point(position = position_jitter())+
  geom_segment(data=CFU_Tukey_sig, aes(x=xmin, xend=xmax, y=y, yend=y))+
  geom_segment(data=CFU_Tukey_sig, aes(x=xmin, xend=xmin, y=y, yend=yend))+
  geom_segment(data=CFU_Tukey_sig, aes(x=xmax, xend=xmax, y=y, yend=yend))+
  geom_text(data=CFU_Tukey_sig, aes(x=(xmin+xmax)/2, y=y, label=sig_label))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format(.x)))+
  coord_cartesian(ylim=c(10^4.5, 10^9))+
  ylab(expression(paste(italic(F.~solani)~density~(log[10](CFU/soil~g)))))+
  xlab("")+
  theme_classic()

Fusarium_density_plot


FW <- read.csv("F.solani_aa_treatment/F.solani weight.csv")
FW$treatment <- gsub("NH4Cl", "NH₄Cl", FW$treatment)
FW$treatment <- gsub("Con", "Control", FW$treatment)
FW$treatment <- factor(FW$treatment, levels = unique(FW$treatment))

unique(FW$treatment)

for(i in c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val")){
  x <- subset(FW, treatment == i)$F_w
  res <- shapiro.test(x)
  print(res)
}

bartlett.test(F_w ~ treatment, FW)

kruskal.test(F_w ~ treatment, FW)
FW_tukey <- nparcomp(F_w~treatment, FW)  

FW_tukey$Analysis$Comparison <- gsub("[()]", "", FW_tukey$Analysis$Comparison)

FW_tukey$Analysis$Comparison <- gsub("^[p] ", "", FW_tukey$Analysis$Comparison)
FW_tukey_p <- as.tibble(FW_tukey$Analysis[, c(1, 6)])

FW_tukey_p

table_theme <- gridExtra::ttheme_default(
  core = list(padding = unit(c(2.5, 2.5), "mm")))

FW_tukey_p$p.Value <- round(FW_tukey_p$p.Value, 3)

colnames(FW_tukey_p)[2] <- "Tukey's HSD\nP value"

FW_tukey_t <- tableGrob(FW_tukey_p, theme = table_theme, rows = NULL)

FW_tukey_t

Fusarium_dry_w_p <- 
  ggpubr::ggarrange(
  grid.arrange(
    ggplot(data = FW, mapping = aes(x = treatment, y = F_w*1000))+
      stat_summary(fun.data = mean_se, geom = "bar")+
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5)+
      geom_jitter()+
      annotate(geom = "text", x=1, y=23, label=expression(paste(italic(P)==0,".000")))+
      theme(panel.background = element_rect(fill="white"),
            axis.line.x.bottom = element_line(color="black"),
            axis.line.y.left = element_line(color="black"))+
      ylab(expression(italic(Fusarium~solani)~dry~weight~(mg)))+
      xlab("")+
      scale_y_continuous(expand = c(0,0), limits = c(0, 25)),
    FW_tukey_t,
    nrow = 1, widths=c(1,0.5))
  )

Fusarium_dry_w_p
#============
# metagenome
#============

# Input raw files directory
path <- "raw_fastq_files"

# Extract raw files list
list.files(path)

# Foward reads list
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))

# Reverse reads list
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extraction sample names from file list
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Apply sample name
names(fnFs) <- sample.names 
names(fnRs) <- sample.names

sample_ord <- c("BC1", "BC3", "BC5", "BC7",
               "RCon1", "RCon2", "RCon3", "RCon4", "RCon5",
               "RAsn1", "RAsn2", "RAsn3", "RAsn4", "RAsn5",
               "RAsp1", "RAsp2", "RAsp3", "RAsp4", "RAsp5",
               "RGlu1", "RGlu2", "RGlu3", "RGlu4", "RGlu5",
               "RVal1", "RVal2", "RVal3", "RVal4", "RVal5",
               "RNH4Cl1", "RNH4Cl2", "RNH4Cl3", "RNH4Cl4", "RNH4Cl5"
               )

#===================================
# 2. Inspect read quality profiles
#===================================
# Visulaization raw foward qulity
raw_foward_Q <- plotQualityProfile(fnFs[sample_ord])
raw_foward_Q

# Visulaization raw reverse qulity
raw_reverse_Q <- plotQualityProfile(fnRs[sample_ord])
raw_reverse_Q

#=====================
# 3. Filter and trim
#=====================
# Place filtered files in filtered/ subdirectory
filt_path <- "filtered_fastq_files"
filtFs <- file.path(filt_path , paste0(sample.names, "_F_filt.fastq.gz")) # set filtered foward path info
filtRs <- file.path(filt_path , paste0(sample.names, "_R_filt.fastq.gz")) # set filtered reverse path info

# Set sample names for filtered
names(filtFs) <- sample.names 
names(filtRs) <- sample.names  

?filterAndTrim

out <- filterAndTrim(fwd = fnFs, # foward reads
                     filt = filtFs, # filterd foward reads save info
                     rev = fnRs, # reverse reads 
                     filt.rev = filtRs, # filterd reverse reads save info
                     truncLen=c(270,200), # trim to Qulity score 30 or higher
                     trimLeft = c(nchar("GTGYCAGCMGCCGCGGTAA"),
                                  nchar("GGACTACNVGGGTWTCTAAT")), # trim the 515F, 806R primer
                     maxN = 0, maxEE = c(1,1), rm.phix = TRUE, # Other statistical quality check criteria
                     n=1e8, # sampling reads number
                     compress = TRUE, # compress
                     multithread = FALSE) # On Windows set multithread=FALSE; An error has been reported when the operation has an unacceptable amount of RAM. In this case, "FALSE" is recommended.


# Visulaization filtered foward qulity
filt_foward_Q <- plotQualityProfile(filtFs[sample_ord])
filt_foward_Q


# Visulaization filtered reverse qulity
filt_reverse_Q <- plotQualityProfile(filtRs[sample_ord])
filt_reverse_Q

#=========================
# 4. Learn the error rates
#=========================
# This step took about 5 hours when 1700X (16 threads) and 3000Mhz DDR4 32Gib ram.
# This step took about 5 hours when 3970X (64 threads) and 3000Mhz DDR4 64Gib ram.
# If you want to save time, reduce the nbase value.

errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e+10)
errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errF, nominalQ=TRUE)
errF_plot
errR_plot

#=====================================================================
# 5. Sample Inference using Divisive Amplicon Denoising Algorithm (DADA)
#=====================================================================
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#======================
# 6. Merge paired reads
#======================
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#=========================
# 7. Construct sequence table
#=========================
seqtab <- makeSequenceTable(mergers)
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 200:300]
#====================
# 8. Remove chimeras
#====================
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#===============================
# 9. Assign taxonomy with IDTAXA
#===============================
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

#==================================
# 10. no-bacterial OTUs removing
#==================================
taxid_OTUtable <- cbind(data.frame(taxid), as.data.frame(t(seqtab.nochim)[,sample_ord]))
taxid_OTUtable <- subset(taxid_OTUtable, domain == "Bacteria")
taxid_OTUtable <- subset(taxid_OTUtable, order != "Chloroplast")
taxid_OTUtable <- subset(taxid_OTUtable, family != "Mitochondria")

#======================================
# 11. Track reads through the pipeline
#======================================
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

track
dir.create("metagenome_analysis_results")
dir.create("metagenome_analysis_results/DADA2_results")
write.csv(track, "metagenome_analysis_results/DADA2_results/Reads_track.csv")

#================
# 12. output files
#================
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

rhizosphere_fasta <- foreach(i = which(rowSums(OTU_counts_table[,6:35]) != 0), .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], rownames(taxid_OTUtable)[i])
  print(match)
}

rhizosphere_fasta <- as.vector(rhizosphere_fasta)
rhizosphere_fasta <- as.data.frame(gsub("^OTU_", ">OTU_", rhizosphere_fasta))
write.table(rhizosphere_fasta, file = "metagenome_analysis_results/DADA2_results/rhizosphere_OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)


# OTU counts table
OTU_counts_table <- cbind(OTUsID, taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
dim(OTU_counts_table)
row.names(OTU_counts_table) <- NULL
colnames(OTU_counts_table)[1] <- "#OTU ID"
write.table(OTU_counts_table, file = "metagenome_analysis_results/DADA2_results/OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(OTU_counts_table[which(rowSums(OTU_counts_table[, 6:35]) > 0), c(1, 6:35)], file = "metagenome_analysis_results/DADA2_results/rhizosphere_OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

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


#===================
# 1. Phylogenetic tree
#===================
# Load fast file
fasta_path <- "metagenome_analysis_results/DADA2_results/OTUs_sequence.fasta"
eqs <- readDNAStringSet(fasta_path)
seqs <- OrientNucleotides(seqs)
aligned_seqs <- AlignSeqs(seqs)
BrowseSeqs(aligned_seqs, highlight = 0)
dir.create("metagenome_analysis_results/metagenomic sequence pylogenetic tree")
writeXStringSet(aligned_seqs, file="metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs_aligned_sequence.fasta")

rhizosphere_fasta_path <- "metagenome_analysis_results/DADA2_results/rhizosphere_OTUs_sequence.fasta"
rhizosphere_seqs <- readDNAStringSet(rhizosphere_fasta_path)
rhizosphere_seqs <- OrientNucleotides(rhizosphere_seqs)
rhizosphere_aligned_seqs <- AlignSeqs(rhizosphere_seqs)
BrowseSeqs(rhizosphere_aligned_seqs, highlight = 0)
writeXStringSet(rhizosphere_aligned_seqs, file="metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere_OTUs_aligned_sequence.fasta")

# Aligned sequences into a database
#---All
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(aligned_seqs, type = "XStringSet", dbFile = dbConn, "")

Add2DB(myData=data.frame(identifier=OTUsID,
                         stringsAsFactors=FALSE),
       dbConn)

cons <- IdConsensus(dbConn,
                    threshold=0.3,
                    minInformation=0.1)

#---rhizosphere
rhizosphere_dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(rhizosphere_aligned_seqs, type = "XStringSet", dbFile = rhizosphere_dbConn, "")


length(which(rowSums(OTU_counts_table[,-1:-5])>0))

Add2DB(myData=data.frame(identifier=OTUsID[which(rowSums(OTU_counts_table[,-1][,5:34])>0)],
                         stringsAsFactors=FALSE),
       rhizosphere_dbConn)

rhizosphere_cons <- IdConsensus(rhizosphere_dbConn,
                                threshold=0.3,
                                minInformation=0.1)


# calculate a maximum like-hood tree
#---All
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


#---rhizosphere
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


# out new newick file
WriteDendrogram(dend, file = "metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum likelihood tree (JC69+G4).newick",
                quoteLabels = F,
                internalLabels = F)
newickfile <- read.table("metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum likelihood tree (JC69+G4).newick", header = F)$V1
newickfile <- gsub('italic[(][\"]',"",newickfile)
newickfile <-gsub('\")',"",newickfile)
write.table(newickfile, quote = F, col.names = F,row.names = F,"metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum likelihood tree (JC69+G4).newick")

WriteDendrogram(rhizosphere_dend, file = "metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere OTUs maximum likelihood tree (T92+G4).newick",
                quoteLabels = F,
                internalLabels = F)
rhizosphere_newickfile <- read.table("metagenome_analysis_results/metagenomic sequence pylogenetic tree/rhizosphere OTUs maximum likelihood tree (T92+G4).newick", header = F)$V1
rhizosphere_newickfile <- gsub('italic[(][\"]',"",newickfile)
rhizosphere_newickfile <-gsub('\")',"",newickfile)
write.table(rhizosphere_newickfile, quote = F, col.names = F,row.names = F,"metagenome_analysis_results/metagenomic sequence pylogenetic tree/OTUs maximum likelihood tree (T92+G4).newick")



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
  read.csv("phenotype_recode/sample_infomation.csv", row.names = 1)
)
SAM$treatment <- gsub("NH[?]{3,}", "NH₄", SAM$treatment)
SAM$treatment <- factor(SAM$treatment, levels = c("bulk", "control", "NH₄Cl", "Glu",  "Asp", "Asn", "Val"))

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

# rarefaction
rarefaction_curve_data <- ggiNEXT( 
  iNEXT(as.data.frame(
    t(OTU)
  )
  )
)$data

rarefaction_curve_data$treatment <- gsub("[0-9]$", "", rarefaction_curve_data$col)
rarefaction_curve_data$treatment <- gsub("^[R]", "", rarefaction_curve_data$treatment)
rarefaction_curve_data$treatment <- gsub("BC", "Bulk soil", rarefaction_curve_data$treatment)
rarefaction_curve_data$treatment <- gsub("Con", "Control", rarefaction_curve_data$treatment)
rarefaction_curve_data$treatment <- gsub("NH4Cl", "NH₄Cl", rarefaction_curve_data$treatment)
rarefaction_curve_data$treatment <- factor(rarefaction_curve_data$treatment, levels <- c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
rarefaction_curve_max_data <- aggregate(x ~ site + treatment + lty, subset(rarefaction_curve_data, lty=="interpolated"), max)
rarefaction_curve_max_data$y <- aggregate(y ~ site + treatment + lty, subset(rarefaction_curve_data, lty=="interpolated"), max)$y
rarefaction_curve_max_data

rarefaction_curve_plot_c3 <- ggplot(data = rarefaction_curve_data, mapping = aes(x=x , y=y, color=site, linetype=lty, label=site))+
  geom_line()+
  geom_point(data=rarefaction_curve_max_data)+
  theme_light()+
  scale_x_continuous(limits=c(0,110000))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  facet_wrap(treatment ~ ., nrow =2)+
  xlab("Reads")+
  ylab("Inferenced OTU number")+
  labs(linetype="Method")+
  guides(color=F)+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background = element_rect(fill="grey90"),
        strip.text = element_text(color="black"))

rarefaction_curve_plot_c3
#------------------
# Alphyadiversity
#----------------
# Calcuation
alpha_diversity <- plot_richness(phy_obj_reads)$data
alpha_diversity$treatment <- gsub("control", "Control", alpha_diversity$treatment)
alpha_diversity$treatment <- gsub("bulk", "Bulk soil", alpha_diversity$treatment)
alpha_diversity$treatment <- factor(alpha_diversity$treatment, levels=c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

# visualization
alpha_diversity_plot_c3 <- ggplot(data = subset(alpha_diversity, alpha_diversity$variable %in% c("Observed", "Shannon", "Simpson")), mapping = aes(x = treatment, y=value))+
  stat_summary(geom = "point", fun.data = "mean_se")+
  stat_summary(geom = "errorbar",size=0.2, fun.data = "mean_se")+
  geom_point(alpha=0.5, aes(color=treatment))+
  theme_light()+
  facet_wrap(. ~ variable, scales = "free_y", nrow=1)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  ylab("Value of alpha diversity index")+
  xlab("")+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1),
        strip.background = element_rect(fill="grey90"),
        strip.text = element_text(color="black"))
alpha_diversity_plot_c3

#=================================================
# Relative abandunce bar graph (beta diversity)
#=================================================
# sorting relative abandunce data  
relative_abundance <- plot_bar(phy_obj_relative)$data
relative_abundance <- relative_abundance[as.character(1:dim(relative_abundance)[1]),]


# output relative abandance data and NA replacement
write.table(relative_abundance,  file = "metagenome_analysis_results/relative_abundance.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t", na = "Not assigned")
relative_abundance <- read.csv("metagenome_analysis_results/relative_abundance.tsv", sep="\t")

relative_abundance$treatment <- gsub("control", "Control", relative_abundance$treatment)
relative_abundance$treatment <- gsub("bulk", "Bulk soil", relative_abundance$treatment)

# Sums data by taxa levels
phylum_relative_abundance <- aggregate(Abundance ~ phylum + Sample + treatment, relative_abundance, sum)
class_relative_abundance <- aggregate(Abundance ~ class + Sample + treatment, relative_abundance, sum)
order_relative_abundance <- aggregate(Abundance ~ order + Sample + treatment, relative_abundance, sum)
family_relative_abundance <- aggregate(Abundance ~ family + Sample + treatment, relative_abundance, sum)
genus_relative_abundance <- aggregate(Abundance ~ genus + Sample + treatment, relative_abundance, sum)


# top 10 sorting
# phylum_level
sum_phylum<- aggregate(Abundance ~ phylum, relative_abundance, sum)
sum_phylum <- sum_phylum %>% 
  arrange(desc(Abundance))
top_phylum <- as.vector(sum_phylum[1:10,1])

# class_level
sum_class <- aggregate(Abundance ~ class, relative_abundance, sum)
sum_class <- sum_class %>% 
  arrange(desc(Abundance))
top_class <- as.vector(sum_class[1:10,1])


# order_level
sum_order <- aggregate(Abundance ~ order, relative_abundance, sum)
sum_order <- sum_order %>% 
  arrange(desc(Abundance))
top_order <- as.vector(sum_order[1:10,1])


# family_level
sum_family <- aggregate(Abundance ~ family, relative_abundance, sum)
sum_family <- sum_family %>% 
  arrange(desc(Abundance))
top_family <- as.vector(sum_family[1:10,1])


# genus_level
sum_genus <- aggregate(Abundance ~ genus, relative_abundance, sum)
sum_genus <- sum_genus %>% 
  arrange(desc(Abundance))
top_genus <- as.vector(sum_genus[1:10,1])

SAM
# taxa top bar graph
top_bar <- function(df, fill_label="Top taxa", top_toxa_names){

  df$top_taxa <- df[,1]
  df$top_taxa <- factor(df$top_taxa, levels = c(top_toxa_names,"others"))
  df$top_taxa[is.na(df$top_taxa)] <- "others"
  df$treatment <- factor(df$treatment, )
    
  df$Sample <- factor(df$Sample, levels = sample_ord) 

  agg_df <- aggregate(Abundance ~ top_taxa + Sample + treatment, df, sum)
  agg_df$treatment <- factor(agg_df$treatment, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu",  "Asp", "Asn", "Val"))
     
  p <- ggplot(data = agg_df, mapping = aes(x = Sample, y=Abundance, fill=top_taxa))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c(brewer.pal(n=10, name="Paired"),"#999999"))+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank())+
    facet_wrap(treatment~., scales = "free_x", nrow = 1)+
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), expand = c(0,0))+ # fill option
    ylab("Relative abundance")+
    xlab("")+
    labs(fill=fill_label)
}

top_10_phylum_bar <- top_bar(phylum_relative_abundance, fill_label = "Top 10 phylum", top_toxa_names = top_phylum )
top_10_phylum_bar 

top_10_class_bar <- top_bar(class_relative_abundance,  fill_label = "Top 10 class", top_toxa_names = top_class )
top_10_class_bar

top_10_order_bar <- top_bar(order_relative_abundance, fill_label = "Top 10 order", top_toxa_names = top_order)
top_10_order_bar

top_10_family_bar <- top_bar(family_relative_abundance, fill_label = "Top 10 family", top_toxa_names = top_family)
top_10_family_bar

top_10_genus_bar <- top_bar(genus_relative_abundance, fill_label = "Top 10 genus", top_toxa_names = top_genus)
top_10_genus_bar

ggarrange(top_10_phylum_bar, top_10_class_bar, top_10_family_bar, top_10_genus_bar, ncol=1, labels=c("A","B","C","D"))

OTU_relative <- OTU/rowSums(OTU)

OTU_seq_tax_relative_table <- cbind(OTU_seq_tax_counts_table[,1:9], t(OTU_relative))

OTU_seq_tax_relative_ldf <- melt(OTU_seq_tax_relative_table)
colnames(OTU_seq_tax_relative_ldf)[10:11] <- c("sample", "relative") 
OTU_family_relative_ldf <- aggregate(relative ~ family + sample, OTU_seq_tax_relative_ldf, sum)
OTU_family_relative_ldf$sample

keystone <- c("Xanthomonadaceae", "Rhizobiaceae", "Sphingobacteriaceae", "Devosiaceae", "Comamonadaceae", "Sphingomonadaceae", "Beijerinckiaceae")
OTU_keystone_relative_ldf <- subset(OTU_family_relative_ldf, family %in% keystone)
OTU_keystone_relative_ldf$treatment <- OTU_keystone_relative_ldf$sample

OTU_keystone_relative_ldf$treatment <- gsub("[R0-9]", "", OTU_keystone_relative_ldf$treatment)
OTU_keystone_relative_ldf$treatment <- gsub("NHCl", "NH₄Cl", OTU_keystone_relative_ldf$treatment)
OTU_keystone_relative_ldf$treatment <- gsub("BC", "Bulk soil", OTU_keystone_relative_ldf$treatment)
OTU_keystone_relative_ldf$treatment <- gsub("Con", "Control", OTU_keystone_relative_ldf$treatment)
OTU_keystone_relative_ldf$treatment <- factor(OTU_keystone_relative_ldf$treatment, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu",  "Asp", "Asn", "Val"))

keystone_relative1 <- 
ggplot(OTU_keystone_relative_ldf, aes(treatment, relative, fill=family))+
  stat_summary(geom="bar", fun.data = mean_se, position = position_stack())+
  stat_summary(geom="errorbar", fun.data = mean_se, width=0.5)+
  geom_jitter(size=0.1)+
  facet_wrap(. ~ family, scales="free", ncol=4)+
  scale_y_continuous(labels = function(x) paste0(x*100), expand=c(0,NA),)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf,yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=-Inf,yend=Inf, size=2)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(size = 10),
        legend.position = "none",
        panel.background = element_blank())+
  labs(x="", y="Relative abundance (%)")

OTU_keystone_relative_summ_ldf<- 
OTU_keystone_relative_ldf%>%
  group_by(sample, treatment) %>%
  summarise(Abundance = sum(relative))

OTU_keystone_relative_summ_ldf

#정규성(p > 0.05면 정규성)
for(i in 1:6){
  print(
    shapiro.test(subset(OTU_keystone_relative_summ_ldf, treatment==c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val")[i])$Abundance)
  )
} # Val is not normality

kruskal.test(Abundance ~ treatment, OTU_keystone_relative_summ_ldf)

keystone_relative2 <-
ggplot(OTU_keystone_relative_ldf, aes(treatment, relative, fill=family))+
  stat_summary(geom="bar", fun.data = mean_se, position = position_stack(), color="white", size=0.2)+
  geom_jitter(data=OTU_keystone_relative_summ_ldf, aes(y=Abundance, fill=NULL))+
  stat_summary(data=OTU_keystone_relative_summ_ldf, geom="errorbar", fun.data = mean_se, aes(y=Abundance, fill=NULL), width=0.5)+
  scale_y_continuous(labels = function(x) paste0(x*100), expand=c(0,0), limits=c(0, 0.4))+
  annotate(geom = "text", x=0.5, y= 0.37, label=expression(paste("Kruskal-Wallis rank sum test,",~italic(P)==0.087)), hjust=0)+
  annotate(geom = "segment", x=Inf, xend=-Inf, y=-Inf,yend=-Inf, size=2)+
  annotate(geom = "segment", x=-Inf, xend=-Inf, y=-Inf,yend=Inf, size=2)+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        strip.text.x = element_text(size = 12),
        legend.position = "none",
        panel.background = element_blank())+
  labs(x="", y="Relative abundance (%)")


ggarrange(keystone_relative2 ,keystone_relative1, widths = c(1,2), labels=c("A", "B"))



#==========picrust2
KO <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"), sep="\t", row.names = 1)
KO

nifs_KO <- c("K02584","K02585","K02586","K02587","K02588","K02589","K02590","K02591","K02592","K02593","K02594","K02595","K02596","K02597")
nifs_KO

names(nifs_KO) <- c("nifA", "nifB", "nifD", "nifE", "nifH", "nifHD1", "nifHD2", "nifK", "nifN", "nifT", "nifV", "nifW", "nifX", "nifZ")

KO_nif <- KO[nifs_KO,]
KO_nif
rownames(KO_nif) <- names(nifs_KO)
KO_nif$nifs <- names(nifs_KO)
KO_nif_ldf <- melt(KO_nif)

colnames(KO_nif_ldf)[2] <- "samples"

KO_nif_ldf$variable <- gsub("[0-9]", "", KO_nif_ldf$samples)
KO_nif_ldf$variable <- gsub("BC", "Bulk soil", KO_nif_ldf$variable)
KO_nif_ldf$variable <- gsub("R", "", KO_nif_ldf$variable)
KO_nif_ldf$variable <- gsub("Con", "Control", KO_nif_ldf$variable)
KO_nif_ldf$variable <- gsub("NHCl", "NH₄Cl", KO_nif_ldf$variable)

KO_nif_ldf$variable <- factor(KO_nif_ldf$variable, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn",  "Val"))

ggplot(KO_nif_ldf, aes(variable, value))+
  stat_summary(geom="bar")+
  geom_point()+
  facet_wrap(.~nifs, scales = "free")+
  theme_light()

subset(KO_nif_ldf, nifs %in% c("nifH", "nifD", "nifK"))%>%
  group_by(samples, variable)%>%
  summarise(value = sum(value))%>%
  ggplot(aes(variable, value))+
  stat_summary(geom="bar", position=position_dodge())+
  stat_summary(geom="errorbar", position=position_dodge(width=0.9), width = 0.5)+
  geom_jitter(position=position_dodge(width=0.9))

ggplot(subset(KO_nif_ldf, nifs %in% c("nifH", "nifD", "nifK")), aes(variable, value, shape=nifs, fill=nifs))+
  stat_summary(geom="bar", position=position_dodge())+
  stat_summary(geom="errorbar", position=position_dodge(width=0.9), width = 0.5)+
  geom_jitter(position=position_dodge(width=0.9))+
  facet_wrap(.~nifs)+
  annotate(geom="segment", x=Inf, xend=-Inf, y=-Inf,yend=-Inf, size=2)+
  annotate(geom="segment", x=-Inf, xend=-Inf, y=-Inf,yend=Inf, size=2)+
  coord_cartesian(ylim = c(500,2400))+
  theme_light()+
  ylab("PICRUSt2 gene score")+
  theme(panel.grid = element_blank())

KO2 <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/KO_predicted.tsv.gz"), sep="\t", row.names = 1)
nifHDK_OTU <- KO2[, nifs_KO[c("nifH","nifD","nifK")]]
nifHDK_OTU <- nifHDK_OTU[which(rowSums(nifHDK_OTU)>0),]

nifHDK_OTU$G <- rowSums(nifHDK_OTU)

nifHDK_OTU_relative <- OTU_seq_tax_relative_ldf[which(OTU_seq_tax_relative_ldf$OTUsID %in% rownames(nifHDK_OTU)),]

nifHDK_OTU_relative$nifHDK_GOO <- nifHDK_OTU_relative$relative * nifHDK_OTU[nifHDK_OTU_relative$OTUsID, 4]

nifHDK_family_relative <- nifHDK_OTU_relative%>%
  group_by(family, sample)%>%
  summarise(nifHDK_GOO=sum(nifHDK_GOO),
            rel_abundance = sum(relative))

nifHDK_family_relative$treatment <- gsub("[0-9]", "", nifHDK_family_relative$sample)
nifHDK_family_relative$treatment <- gsub("BC", "Bulk soil", nifHDK_family_relative$treatment)
nifHDK_family_relative$treatment <- gsub("^R", "", nifHDK_family_relative$treatment)
nifHDK_family_relative$treatment <- gsub("Con", "Control", nifHDK_family_relative$treatment)
nifHDK_family_relative$treatment <- gsub("NHCl", "NH₄Cl", nifHDK_family_relative$treatment)
nifHDK_family_relative$treatment

nifHDK_family_relative$treatment <- factor(nifHDK_family_relative$treatment, levels =  c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

nifHDK_family_GOO <- nifHDK_family_relative
rm(nifHDK_family_relative)

nifHDK_GOO_ind_plot <- ggplot(nifHDK_family_GOO, aes(x=treatment, y=nifHDK_GOO, fill=family))+
  stat_summary(geom="bar")+
  stat_summary(geom="errorbar", width=0.3, size=0.3)+
  annotate(geom="segment", x=c(-Inf,-Inf), xend = c(Inf, -Inf), y=c(-Inf, -Inf), yend=c(-Inf,Inf))+
  geom_jitter(size=0.2)+
  facet_wrap(family~., scales = "free")+
  ylab(expression(italic(nifHDK)~GOO~index))+
  xlab("")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1),
        panel.background = element_blank())

nifHDK_GOO_ind_plot

nifHDK_family_GOO_sum <- nifHDK_family_GOO_sum %>% group_by(treatment, sample) %>% summarise(nifHDK_GOO = sum(nifHDK_GOO))

for(i in 1:6){
  print(
    shapiro.test(subset(nifHDK_family_GOO_sum, treatment==c("Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val")[i])$nifHDK_GOO)
  )
} # some GOO is not normality

kruskal.test(nifHDK_GOO ~ treatment, nifHDK_family_GOO_sum)
nifHDK_GOO_Tukey<- nparcomp(nifHDK_GOO ~ treatment, nifHDK_family_GOO_sum)

nifHDK_GOO_Tukey_sig<- nifHDK_GOO_Tukey$Analysis[which(nifHDK_GOO_Tukey$Analysis$p.Value <0.05), ]

nifHDK_GOO_Tukey_sig$x1 <- "Bulk soil"
nifHDK_GOO_Tukey_sig$x2 <- c("Glu", "Asp", "Asn", "Val")
nifHDK_GOO_Tukey_sig$y <- 1:4*0.01+0.1

nifHDK_GOO_Tukey_sig$x1 <- factor(nifHDK_GOO_Tukey_sig$x1, levels = c("Bulk soil","Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
nifHDK_GOO_Tukey_sig$x2 <- factor(nifHDK_GOO_Tukey_sig$x2, levels = c("Bulk soil","Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
nifHDK_family_GOO$treatment

nifHDK_GOO_plot <- ggplot(nifHDK_family_GOO, aes(x=treatment, y=nifHDK_GOO, fill=family))+
  stat_summary(geom="bar", position = position_stack(), color="white", size=0.2)+
  stat_summary(data=nifHDK_family_GOO_sum, aes(fill=NULL), geom="errorbar", width=0.5)+
  geom_jitter(data=nifHDK_family_GOO_sum, aes(fill=NULL))+
  annotate(geom="segment", x=c(-Inf,-Inf), xend = c(Inf, -Inf), y=c(-Inf, -Inf), yend=c(-Inf,Inf))+
  annotate(geom="text", x=0.6, y = 0.15, label=expression(paste("Kruskal-Wallis sum rank test, ", italic(P)==0.013)), hjust=0,size=3.5)+
  geom_segment(data=nifHDK_GOO_Tukey_sig, aes(x=x1, xend=x2, y=y, yend=y, fill=NULL))+
  geom_segment(data=nifHDK_GOO_Tukey_sig, aes(x=x1, xend=x1, y=y, yend=y-0.002, fill=NULL))+
  geom_segment(data=nifHDK_GOO_Tukey_sig, aes(x=x2, xend=x2, y=y, yend=y-0.002, fill=NULL))+
  geom_text(data=nifHDK_GOO_Tukey_sig, aes(x=(as.numeric(x2)+as.numeric(x1))/2, y=y, fill=NULL, label="***"), vjust=0)+
  ylab(expression(italic(nifHDK)~GOO~index))+
  scale_y_continuous(expand=c(0,0), limits = c(0,0.16))+
  xlab("")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1),
        panel.background = element_blank())

nifHDK_GOO_plot

ggarrange(nifHDK_GOO_plot, nifHDK_GOO_ind_plot, widths=c(1,2), labels = c("A","B"))

#===========================================picrust2 pathway
PW <-read.csv(gzfile("metagenome_analysis_results/picrust2_results/pathways_out/path_abun_unstrat_descrip.tsv.gz"), sep="\t")
gsub("R","",colnames(PW)[-c(1:2)]) == Ct_sample_summary[c(11:19,1:10,20:24,30:34,25:29),]

PW_amino_rhizo <- PW[which(rowSums(PW[,-c(1:6)] > 0) > 15),-c(3:6)]

Fusarium_dens <- Ct_sample_summary[c(11:19,1:10,20:24,30:34,25:29),]

Fd_DI <- cbind(Fusarium_dens[-c(1:4),], DI[c(1:5,16:20, 6:15, 26:30, 21:25),8])
colnames(Fd_DI)[6] <- "PotDI"

Fd_DI_pearson_res <- cor.test(Fd_DI$log10CFU, Fd_DI$PotDI, alternative = "two.sided")
Fd_DI_pearson_res <- data.frame(versus1="Fusarium density", versus2="DI", path_description="", rho = Fd_DI_pearson_res$estimate, P = Fd_DI_pearson_res$p.value)
Fd_DI_pearson_res

Fd_DI_pearson_pot_c3<- ggplot(Fd_DI, aes(PotDI, log10(sample_CFU_per_soil),  color=treatment, shape=treatment))+
  geom_point(size=3)+
  stat_smooth(aes(color=NULL, shape=NULL), color = "#FC4E07", method = "lm",  level = 0)+
  theme_light()+
  xlab("Disease index")+
  ylab(expression(italic(F.~solani)~log[10](CFU/soil~g)))+
  annotate("text", x=2, y=7, label=expression(paste("Pearson's correlation, ", rho == 0.403,", ", italic(P) == 0.027)), size=4, hjust=0)

#Fusarium vs Pathway cor.test

P <- c()
rho <- c()
pathID <- c()
path_description <- c()
for(i in 1:nrow(PW_amino_rhizo)){
  cor_res <- cor.test(log10(Fusarium_dens$sample_CFU_per_soil[-c(1:4)]), as.numeric(PW_amino_rhizo[i,-c(1:2)]), alternative = "two.sided")
  P <- c(P, cor_res$p.value)
  rho <- c(rho, cor_res$estimate)
  pathID <- c(pathID, PW_amino_rhizo$pathway[i])
  path_description <- c(path_description, PW_amino_rhizo$description[i])
  }

pearson_amino_Fd_vs_Pathway <- data.frame(versus2 = pathID, path_description, rho, P)

#DI vs Pathway cor.test

P <- c()
rho <- c()
pathID <- c()
path_description <- c()
for(i in 1:nrow(PW_amino_rhizo)){
  cor_res <- cor.test(Fd_DI$PotDI, as.numeric(PW_amino_rhizo[i,-c(1:2)]), alternative = "two.sided")
  P <- c(P, cor_res$p.value)
  rho <- c(rho, cor_res$estimate)
  pathID <- c(pathID, PW_amino_rhizo$pathway[i])
  path_description <- c(path_description, PW_amino_rhizo$description[i])
}

pearson_amino_DI_vs_Pathway <- data.frame(versus2 = pathID, path_description, rho, P)

pearson_amino_DI_vs_Pathway <- cbind(versus1 = "DI", pearson_amino_DI_vs_Pathway)
pearson_amino_Fd_vs_Pathway <- cbind(versus1 = "Fusarium density", pearson_amino_Fd_vs_Pathway)

pearson_amino_results <- rbind(Fd_DI_pearson_res, pearson_amino_DI_vs_Pathway, pearson_amino_Fd_vs_Pathway)

pearson_amino_results$FDR_q <- p.adjust(pearson_amino_results$P, method="fdr")

write.csv(subset(pearson_amino_results, FDR_q < 0.01 & rho > 0), "metagenome_analysis_results/picrust2_results/pearson_plus_rho_significant_pathway.csv")
write.csv(subset(pearson_amino_results, FDR_q < 0.01 & rho < 0), "metagenome_analysis_results/picrust2_results/pearson_minus_rho_significant_pathway.csv")

negative_cor_PW <- subset(pearson_amino_results, FDR_q < 0.01 & rho < 0)[,2]
names(negative_cor_PW) <- subset(pearson_amino_results, FDR_q < 0.01 & rho < 0)[,3]
negative_cor_pw_table <- subset(pearson_amino_results, FDR_q < 0.01 & rho < 0)[,-2]

Fd_DI_negative_cor_PW <- data.frame(Fd_DI, t(PW_amino_rhizo[which(PW_amino_rhizo$pathway%in%negative_cor_PW),][,-c(1,2)]))
colnames(Fd_DI_negative_cor_PW)[-c(1:6)] <- PW[which(PW$pathway%in%negative_cor_PW),2]
Fd_DI_negative_cor_PW <- melt(Fd_DI_negative_cor_PW, id.vars = c("value", "treatment","sample_CFU_per_soil","SD","log10CFU","PotDI"))
colnames(Fd_DI_negative_cor_PW)[1] <- "sample_ID"

Fd_DI_negative_cor_PW$negative_to <- ifelse(Fd_DI_negative_cor_PW$variable %in% names(negative_cor_PW)[1:5], "DI negative", "F. solani negative")
Fd_DI_negative_cor_PW$variable <- gsub(" [(]", "\n(", Fd_DI_negative_cor_PW$variable)

negative_cor_pw_table$label <- c("2.268%*%10^-3", "7.296%*%10^-5", "8.767%*%10^-3", "2.368%*%10^-4", "2.692%*%10^-3", "4.409%*%10^-3", "4.119%*%10^-3", "4.292%*%10^-3")

negative_cor_pw_table

colnames(negative_cor_pw_table)[2] <- "variable"
negative_cor_pw_table$variable <- gsub(" [(]", "\n(", negative_cor_pw_table$variable)

rm(deparse2)

DI_negative_PW_p <- 
  ggplot(subset(Fd_DI_negative_cor_PW, negative_to =="DI negative"), aes(x=PotDI, y= value, color=treatment))+
  geom_point()+
  facet_wrap(variable~., ncol=3, scales = "free")+
  stat_smooth(aes(group=NULL, color=NULL), method = "lm", level = 0, show.legend = FALSE)+
  ylab("PICRUSt2 pathway score")+
  xlab("Root rot disease index")+
  geom_text(data = negative_cor_pw_table[1:5,], aes(x=11, y=Inf, label = paste("italic(P[adj]) ==", label), color=NULL), parse = T, vjust=2,  show.legend = FALSE)+
  labs(color="")+
  theme_light()+
  theme(strip.text = element_text(color="black"),
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(color="black"),
        axis.ticks = element_line(color="black"))

DI_negative_PW_p

Fd_DI_negative_cor_PW$variable <- gsub("superpathway of geranyl", "superpathway of\ngeranyl", Fd_DI_negative_cor_PW$variable)
negative_cor_pw_table$variable <- gsub("superpathway of geranyl", "superpathway of\ngeranyl", negative_cor_pw_table$variable)
Fd_negative_PW_p <- 
  ggplot(subset(Fd_DI_negative_cor_PW, negative_to =="F. solani negative"), aes(x=log10(sample_CFU_per_soil), y= value, color=treatment))+
  geom_point()+
  stat_smooth(aes(group=NULL, color=NULL), method = "lm", level = 0, show.legend = FALSE)+
  facet_wrap(variable~., ncol=3, scales = "free")+
  ylab("PICRUSt2 pathway score")+
  xlab(expression(italic(F.~solani)~density~(log[10](CFU/soil~g))))+
  geom_text(data = negative_cor_pw_table[6:8,], aes(x=6.3, y=Inf, label = paste("italic(P[adj]) ==", label), color=NULL), parse = T, vjust=2,  show.legend = FALSE)+
  labs(color="")+
  theme_light()+
  theme(strip.text = element_text(color="black"),
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(color="black"),
        axis.ticks = element_line(color="black"))

Fd_negative_PW_p

ggarrange(DI_negative_PW_p, Fd_negative_PW_p, nrow = 2, heights = c(2, 1.1),common.legend = TRUE, legend="right", labels = c("A","B"))

DI_negative_PW_bar <- 
  ggplot(subset(Fd_DI_negative_cor_PW, negative_to =="DI negative"), aes(x=treatment, y= value, color=treatment))+
  geom_point()+
  facet_nested_wrap(variable~., ncol=3, scales = "free")+
  ylab("PICRUSt2 pathway score")+
  xlab("Root rot disease index")+
  labs(color="")+
  theme_light()+
  theme(strip.text = element_text(color="black"),
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(color="black"),
        axis.ticks = element_line(color="black"))

DI_negative_PW_bar

Fd_negative_PW_bar <- 
  ggplot(subset(Fd_DI_negative_cor_PW, negative_to =="F. solani negative"), aes(x=treatment, y= value, color=treatment))+
  geom_point()+
  facet_wrap(variable~., ncol=3, scales = "free")+
  ylab("PICRUSt2 pathway score")+
  xlab(expression(log[10](italic(F.~solani)~CFU/soil~g)))+
  labs(color="")+
  theme_light()+
  theme(strip.text = element_text(color="black"),
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(color="black"),
        axis.ticks = element_line(color="black"))
Fd_negative_PW_bar

ggarrange(DI_negative_PW_p, Fd_negative_PW_p, nrow = 2, heights = c(2,1.1),common.legend = TRUE, legend="right")


KO_OTU <- read.csv("metagenome_analysis_results/picrust2_results/KO_predicted.tsv.gz", sep = "\t")

# AST-pathway
AST_PWY <- c(
  "K00673", "EC:2.3.1.109",
  "K01484", "EC:3.5.3.23",
  "K00840", "EC:2.6.1.81",
  "K06447", "EC:1.2.1.71",
  "K05526", "EC:3.5.1.96"
)

# ARGORNPROST-PWY

ARGORNPROST_PWY <- c(
  "K01478", "EC:3.5.3.6",
  "K00611", "EC:2.1.3.3",
  "K00926", "EC:2.7.2.2",
  "K21898", "EC:5.1.1.12",
  "K17899", "EC:5.4.3.5",
  "K21672", "EC:1.4.1.12",
  "K21399", "EC:2.3.1.263",
  "K21400", "EC:2.3.1.263",
  "K00819", "EC:2.6.1.13",
  "K00286", "EC:1.5.1.2",
  "K01777", "EC:5.1.1.4",
  "K10793", "EC:1.21.4.1",
  "K10794",
  "K01750", "EC:4.3.1.12"
)

#DHGLUCONATE-PYR-CAT-PWY

DHGLUCONATE_PYR_CAT_PWY <- 
  c("K00117", "EC:1.1.5.2",
    "K01053", "EC:3.1.1.17",
    "K06151", "K06152", "EC:1.1.99.3",
    "K11441", "EC:2.7.1.13",
    "K00032", "EC:1.1.1.43"
    )


# pylymyxin resistance
PWY0_1338 <-
  c("K10011", "EC:1.1.1.305",
    "K07806", "EC:2.6.1.87",
    "K15895", "EC:2.6.1.92",
    "K10011", "EC:2.1.2.13",
    "K10012", "EC:2.4.2.53",
    "K07264", "EC:2.4.2.43"
    )

# pyridoxal 5'-phosphate biosynthesis I
PYRIDOXSYN_PWY <-
  c("K03472", "EC:1.2.1.72",
    "K03473", "EC:1.1.290",
    "K00831", "EC:2.6.1.52",
    "K00097", "EC:1.1.1.262",
    "K01662", "EC:2.2.1.7",
    "K03474", "EC:2.6.99.2",
    "K00275", "EC:1.4.3.5"
  )

# superpathway of geranylgeranyldiphosphate biosynthesis I (via mevalonate)
GG_PW<-
  c("K01662", "EC:2.2.1.7",
    "K00099", "EC:1.1.1.267",
    "K00991", "K12506", "EC:2.7.7.60",
    "K00919", "EC:2.7.1.148",
    "K01770", "K12506", "EC:4.6.1.12",
    "K03526", "EC:1.17.7.3",
    "K03527", "EC:1.17.7.4",
    "K01823", "EC:5.3.3.2",
    "K13787", "EC:2.5.1.1",
    "K13787", "EC:2.5.1.10",
    "K13787", "EC:2.5.1.29"
  )


# isoprene biosynthesis II (engineered)
IP_PW<-
  c("K00626","EC:2.3.1.9",
    "K00632","K07508","K07509","K07513","EC:2.3.1.16",
    "K01641","EC:2.3.3.10",
    "K00021","EC:1.1.1.34",
    "K00869","EC:2.7.1.36",
    "K00938","K13273","EC:2.7.4.2",
    "K01597","EC:4.1.1.33",
    "K01823","EC:5.3.3.2",
    "K12742","EC:4.2.3.27"
  )

# mevalonate pathway I
MEV_PW<-
  c("K00626","EC:2.3.1.9",
    "K00632","K07508","K07509","K07513","EC:2.3.1.16",
    "K01641","EC:2.3.3.10",
    "K00021","EC:1.1.1.34",
    "K00869","EC:2.7.1.36",
    "K00938","K13273","EC:2.7.4.2",
    "K01597","EC:4.1.1.33",
    "K01823","EC:5.3.3.2"
  )

PW_table <- 
rbind(
  data.frame("MetaCyc pathway ID"="AST-PWY", "pathway description" = "L-arginine degradation II (AST pathway)", "KEGG KO ID" = AST_PWY),
  data.frame("MetaCyc pathway ID"="ARGORNPROST-PWY", "pathway description" = "L-arginine degradation (Stickland reaction)", "KEGG KO ID" = ARGORNPROST_PWY),
  data.frame("MetaCyc pathway ID"="DHGLUCONATE-PYR-CAT-PWY", "pathway description" = "Glucose degradation (oxidative)", "KEGG KO ID" = DHGLUCONATE_PYR_CAT_PWY),
  data.frame("MetaCyc pathway ID"="PWY0_1338", "pathway description" = "Polymyxin resistance", "KEGG KO ID" = PWY0_1338),
  data.frame("MetaCyc pathway ID"="PYRIDOXSYN-PWY", "pathway description" = "pyridoxal 5'-phosphate biosynthesis I", "KEGG KO ID" = PYRIDOXSYN_PWY),
  data.frame("MetaCyc pathway ID"="PWY-5910", "pathway description" = "Superpathway of geranylgeranyldiphosphate biosynthesis I (via mevalonate)", "KEGG KO ID" = IP_PW),
  data.frame("MetaCyc pathway ID"="PWY-7391", "pathway description" = "Isoprene biosynthesis II (engineered)", "KEGG KO ID" = GG_PW),
  data.frame("MetaCyc pathway ID"="PWY-922", "pathway description" = "Mevalonate pathway I (eukaryotes and bacteria)", "KEGG KO ID" = MEV_PW)
)

PW_table_KO <- subset(PW_table, substr(KEGG.KO.ID,1,1)=="K")
KO2 <- read.csv(gzfile("metagenome_analysis_results/picrust2_results/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"), sep="\t", row.names = 1)

PW_table_KO$KEGG.KO.ID.description <- KO2[PW_table_KO$KEGG.KO.ID,1]
PW_table_KO$KEGG.KO.ID.description[is.na(PW_table_KO$KEGG.KO.ID.description)] <- "not included PICRUSt2 (version 2.3.0-b)"

PW_OTUrelative <-
  function(PWY){
    PWY_KO <- cbind(Fd_DI, data.frame(t(KO[which(rownames(KO)%in%PWY),-(1:4)])))
    colnames(PWY_KO)[1] <- "sampleID"
    Fd_DI_PWY_KO_ldf <- melt(PWY_KO[,c(1,2,3,6:ncol(PWY_KO))], id.vars = c("sampleID","treatment","sample_CFU_per_soil","PotDI"))
                              
                              plot1 <- 
                                ggplot(Fd_DI_PWY_KO_ldf, aes(PotDI, value, color=treatment))+
                                geom_point()+
                                facet_wrap(.~variable, scales = "free", nrow=1)
                              plot2 <-
                                ggplot(Fd_DI_PWY_KO_ldf, aes(sample_CFU_per_soil, value, color=treatment))+
                                geom_point()+
                                facet_wrap(.~variable, scales = "free", nrow=1)
                              plot3 <-
                                ggplot(Fd_DI_PWY_KO_ldf, aes(treatment, value))+
                                geom_point()+
                                facet_wrap(.~variable, scales = "free", nrow=1)
                 
                              
                              PWY_KO_OTU <- KO_OTU[,colnames(KO_OTU) %in% c("sequence", PWY)]
                              
                              PWY_KO_OTU <- PWY_KO_OTU[
                                rowSums(PWY_KO_OTU[,-1]) > 0 ,
                              ]
                              
                              PWY_KO_OTU <- cbind(PWY_KO_OTU, TAX[PWY_KO_OTU$sequence])
                              PWY_KO_OTU$sum_score <- rowSums(PWY_KO_OTU[,2:(which(colnames(PWY_KO_OTU)=="domain")-1)])
                              rownames(PWY_KO_OTU) <- PWY_KO_OTU$sequence
                              
                              PWY_OTU_rel <- relative_abundance[which(relative_abundance$OTU %in% PWY_KO_OTU$sequence),]
                              
                              PWY_OTU_rel_individual <- PWY_OTU_rel %>%
                                group_by(Sample, treatment, OTU) %>%
                                summarise("Abundance"=sum(Abundance))
                              
                              PWY_OTU_rel_individual$sum_score <- as.character(PWY_KO_OTU[PWY_OTU_rel_individual$OTU,"sum_score"])
                              PWY_OTU_rel_individual$sum_score <- factor(PWY_OTU_rel_individual$sum_score, levels=as.character(1:30))
                              PWY_OTU_rel_individual$treatment <- gsub("bulk","Bulk soil",PWY_OTU_rel_individual$treatment)
                              PWY_OTU_rel_individual$treatment <- gsub("control","Control",PWY_OTU_rel_individual$treatment)
                              PWY_OTU_rel_individual$treatment <- factor(PWY_OTU_rel_individual$treatment, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
                              PWY_OTU_rel_individual$AkGk <- PWY_OTU_rel_individual$Abundance*as.numeric(PWY_OTU_rel_individual$sum_score)
                              
                              PWY_OTU_rel_individual <- PWY_OTU_rel_individual%>%
                                group_by(Sample, treatment, sum_score) %>%
                                summarise("Abundance"=sum(Abundance),
                                          "AkGk"=sum(AkGk))
                              
                              PWY_OTUrel_summary<- PWY_OTU_rel_individual%>%
                                group_by(Sample, treatment) %>%
                                summarise("Abundance"=sum(Abundance),
                                          "SumAkGk"=sum(AkGk))
                              
                              PWY_OTUrel_summary$PotDI <- c(rep(NA,4), DI$DI_mean[c(16:20, 6:10, 1:5, 11:15, 21:30)])
                              PWY_OTUrel_summary$Fusarium_CFU <- c(rep(NA,4), Fd_DI$sample_CFU_per_soil[c(6:10, 11:15, 1:5, 16:20, 26:30, 21:25)])
                              
                              print(
                                list("PWY_OTU_rel_individual"=PWY_OTU_rel_individual, 
                                     "PWY_OTUrel_sum"=PWY_OTUrel_summary, 
                                     "DI_KOvalue_p"=plot1, 
                                     "Fdensity_KOvalue_p"=plot2,
                                     "treatment_KOvalue_p"=plot3,
                                     "OTU"=PWY_KO_OTU[,1])
                              )
  }

AST_PWY_OTU_list <- PW_OTUrelative(AST_PWY)
AGR_PWY_OTU_list <- PW_OTUrelative(ARGORNPROST_PWY)
DG_PWY_OTU_list <- PW_OTUrelative(DHGLUCONATE_PYR_CAT_PWY)
PR_PWY_OTU_list <- PW_OTUrelative(PWY0_1338)
P5P_PWY_OTU_list <- PW_OTUrelative(PYRIDOXSYN_PWY)
GG_PWY_OTU_list <- PW_OTUrelative(GG_PW)
IP_PWY_OTU_list <- PW_OTUrelative(IP_PW)
MEV_PWY_OTU_list <- PW_OTUrelative(MEV_PW)

AST_PWY_OTU_list$PWY_OTUrel_sum$treatment

AST_PWY_bar_P <- 
  ggplot(AST_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=AST_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadance\nof bacteria contributing \n to L-arginine degradation II (%)")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))

for(i in 1:6){
  print(
    shapiro.test(
      subset(AST_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(AST_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

bartlett.test(SumAkGk ~ treatment, AST_PWY_OTU_list$PWY_OTUrel_sum)

kruskal.test(SumAkGk ~ treatment, AST_PWY_OTU_list$PWY_OTUrel_sum)

AST_PWY_OTUrel_sum_tukey <- nparcomp(SumAkGk ~ treatment, AST_PWY_OTU_list$PWY_OTUrel_sum)$Analysis
AST_PWY_OTUrel_sum_tukey <- AST_PWY_OTUrel_sum_tukey[AST_PWY_OTUrel_sum_tukey$p.Value < 0.05,]

AST_PWY_OTUrel_sum_tukey <- 
  cbind(AST_PWY_OTUrel_sum_tukey,
        gsub(" [)]",
             "",
             gsub("p[(] ",
                  "",
                  matrix(
                    unlist(
                      strsplit(AST_PWY_OTUrel_sum_tukey$Comparison, split=" , ")
                    ),
                    ncol=2, 
                    byrow=T
                  )
             )
        )
  )

colnames(AST_PWY_OTUrel_sum_tukey)[7:8] <- c("x1", "x2")

AST_PWY_OTUrel_sum_tukey$x1 <- factor(AST_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
AST_PWY_OTUrel_sum_tukey$x2 <- factor(AST_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
AST_PWY_OTUrel_sum_tukey$y <- 1:9*0.15+2

AST_PWY_bar_P2 <- 
  ggplot(AST_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_segment(data=AST_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=AST_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=AST_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=AST_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label="***"))+
  annotate(geom="text", x=4, y=3.6, label=expression(paste("Kruskal-Wallis rank sum test, ", italic(P)==0.001)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 3.8))

#----------------------
AGR_PWY_bar_P <- 
  ggplot(AGR_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=AGR_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadanceof bacteria contributing \nto arginine, ornithine\nand proline interconversion (%)")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))

for(i in 1:6){
  print(
    shapiro.test(
      subset(AGR_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(AGR_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

kruskal.test(SumAkGk ~ treatment, AGR_PWY_OTU_list$PWY_OTUrel_sum)

AGR_PWY_OTUrel_sum_tukey <- nparcomp(SumAkGk ~ treatment, AGR_PWY_OTU_list$PWY_OTUrel_sum)$Analysis
AGR_PWY_OTUrel_sum_tukey <- AGR_PWY_OTUrel_sum_tukey[AGR_PWY_OTUrel_sum_tukey$p.Value < 0.05,]

AGR_PWY_OTUrel_sum_tukey <- 
  cbind(AGR_PWY_OTUrel_sum_tukey,
        gsub(" [)]",
             "",
             gsub("p[(] ",
                  "",
                  matrix(
                    unlist(
                      strsplit(AGR_PWY_OTUrel_sum_tukey$Comparison, split=" , ")
                    ),
                    ncol=2, 
                    byrow=T
                  )
             )
        )
  )

colnames(AGR_PWY_OTUrel_sum_tukey)[7:8] <- c("x1", "x2")

AGR_PWY_OTUrel_sum_tukey$x1 <- factor(AGR_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
AGR_PWY_OTUrel_sum_tukey$x2 <- factor(AGR_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
AGR_PWY_OTUrel_sum_tukey$y <- 1:3*0.15+5


AGR_PWY_bar_P2 <- 
  ggplot(AGR_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_segment(data=AGR_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=AGR_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=AGR_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=AGR_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label="***"))+
  annotate(geom="text", x=4, y=5.8, label=expression(paste("Kruskal-Wallis rank sum test, ", italic(P)==0.013)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(3.5, 6))

AGR_PWY_bar_P2 
#----------------------
DG_PWY_bar_P <- 
  ggplot(DG_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=DG_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadance\nof bacteria contributing \n to glucose degradation (%)")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))

for(i in 1:6){
  print(
    shapiro.test(
      subset(DG_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(DG_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

kruskal.test(SumAkGk ~ treatment, DG_PWY_OTU_list$PWY_OTUrel_sum)

DG_PWY_OTUrel_sum_tukey <- nparcomp(SumAkGk ~ treatment, DG_PWY_OTU_list$PWY_OTUrel_sum)$Analysis
DG_PWY_OTUrel_sum_tukey <- DG_PWY_OTUrel_sum_tukey[DG_PWY_OTUrel_sum_tukey$p.Value < 0.05,]

DG_PWY_OTUrel_sum_tukey <- 
  cbind(DG_PWY_OTUrel_sum_tukey,
        gsub(" [)]",
             "",
             gsub("p[(] ",
                  "",
                  matrix(
                    unlist(
                      strsplit(DG_PWY_OTUrel_sum_tukey$Comparison, split=" , ")
                    ),
                    ncol=2, 
                    byrow=T
                  )
             )
        )
  )

colnames(DG_PWY_OTUrel_sum_tukey)[7:8] <- c("x1", "x2")

DG_PWY_OTUrel_sum_tukey$x1 <- factor(DG_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
DG_PWY_OTUrel_sum_tukey$x2 <- factor(DG_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
DG_PWY_OTUrel_sum_tukey$y <- 1:7*0.15+2.5

DG_PWY_bar_P2 <- 
  ggplot(DG_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_segment(data=DG_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=DG_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=DG_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=DG_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label="***"))+
  annotate(geom="text", x=4, y=3.8, label=expression(paste("Kruskal-Wallis rank sum test, ", italic(P)==0.002)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0.5, 4))
DG_PWY_bar_P2 
#----------------------
PR_PWY_bar_P <- 
  ggplot(PR_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=PR_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadance\nof bacteria contributing \n to polymyxin resistance (%)")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))

for(i in 1:6){
  print(
    shapiro.test(
      subset(PR_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(PR_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

bartlett.test(SumAkGk ~ treatment, PR_PWY_OTU_list$PWY_OTUrel_sum)

kruskal.test(SumAkGk ~ treatment, PR_PWY_OTU_list$PWY_OTUrel_sum)

PR_PWY_OTUrel_sum_tukey <- nparcomp(SumAkGk ~ treatment, PR_PWY_OTU_list$PWY_OTUrel_sum)$Analysis
PR_PWY_OTUrel_sum_tukey <- PR_PWY_OTUrel_sum_tukey[PR_PWY_OTUrel_sum_tukey$p.Value < 0.05,]

PR_PWY_OTUrel_sum_tukey <- 
  cbind(PR_PWY_OTUrel_sum_tukey,
        gsub(" [)]",
             "",
             gsub("p[(] ",
                  "",
                  matrix(
                    unlist(
                      strsplit(PR_PWY_OTUrel_sum_tukey$Comparison, split=" , ")
                    ),
                    ncol=2, 
                    byrow=T
                  )
             )
        )
  )

colnames(PR_PWY_OTUrel_sum_tukey)[7:8] <- c("x1", "x2")

PR_PWY_OTUrel_sum_tukey$x1 <- factor(PR_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
PR_PWY_OTUrel_sum_tukey$x2 <- factor(PR_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
PR_PWY_OTUrel_sum_tukey$y <- 1:8*0.1+1.3

PR_PWY_bar_P2 <- 
  ggplot(PR_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  labs(fill="number\nof assosicated genes")+
  geom_segment(data=PR_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=PR_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=PR_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=PR_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label="***"))+
  annotate(geom="text", x=4, y=2.3, label=expression(paste("Kruskal-Wallis rank sum test, ", italic(P)==0.001)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 2.5))

#----------------------------------------------
P5P_PWY_bar_P <- 
  ggplot(P5P_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=P5P_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadance\nof bacteria contributing \nto pyridoxal 5'-phosphate biosynthesis I (%)")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))

for(i in 1:6){
  print(
    shapiro.test(
      subset(P5P_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(P5P_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

bartlett.test(SumAkGk ~ treatment, P5P_PWY_OTU_list$PWY_OTUrel_sum)


P5P_PWY_SumAkGk_aov <- aov(SumAkGk ~ treatment, P5P_PWY_OTU_list$PWY_OTUrel_sum)
unlist(summary(P5P_PWY_SumAkGk_aov))

P5P_PWY_OTUrel_sum_tukey <- data.frame(TukeyHSD(P5P_PWY_SumAkGk_aov)$treatment)

P5P_PWY_OTUrel_sum_tukey <- P5P_PWY_OTUrel_sum_tukey[P5P_PWY_OTUrel_sum_tukey$p.adj < 0.05,]

P5P_PWY_OTUrel_sum_tukey$p.adj_mark <- "*"
P5P_PWY_OTUrel_sum_tukey$p.adj_mark[which(P5P_PWY_OTUrel_sum_tukey$p.adj < 0.01)] <- "**"
P5P_PWY_OTUrel_sum_tukey$p.adj_mark[which(P5P_PWY_OTUrel_sum_tukey$p.adj < 0.001)] <- "***"

P5P_PWY_OTUrel_sum_tukey

P5P_PWY_OTUrel_sum_tukey <- 
  cbind(P5P_PWY_OTUrel_sum_tukey,
        matrix(
          unlist(
            strsplit(rownames(P5P_PWY_OTUrel_sum_tukey), split="-")
          ),
          ncol=2, 
          byrow=T
        )
  )

colnames(P5P_PWY_OTUrel_sum_tukey)[6:7] <- c("x1", "x2")

P5P_PWY_OTUrel_sum_tukey$x1 <- factor(P5P_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
P5P_PWY_OTUrel_sum_tukey$x2 <- factor(P5P_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
P5P_PWY_OTUrel_sum_tukey$y <- 1:9*0.12+5.5

P5P_PWY_bar_P2 <- 
  ggplot(P5P_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_segment(data=P5P_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=P5P_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=P5P_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=P5P_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label=p.adj_mark))+
  annotate(geom="text", x=4, y=6.8, label=expression(paste("ANOVA, ", italic(P)==6.542%*%10^-8)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(3.5, 7))


#----------------------------------------------
GG_PWY_bar_P <- 
  ggplot(GG_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=GG_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadance\nof bacteria contributing to superpathway\nof geranylgeranyl diphosphate biosynthesis II(%)")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))


for(i in 1:6){
  print(
    shapiro.test(
      subset(GG_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(GG_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

bartlett.test(SumAkGk ~ treatment, GG_PWY_OTU_list$PWY_OTUrel_sum)

GG_PWY_SumAkGk_aov <- aov(SumAkGk ~ treatment, GG_PWY_OTU_list$PWY_OTUrel_sum)
unlist(summary(GG_PWY_SumAkGk_aov))

GG_PWY_OTUrel_sum_tukey <- data.frame(TukeyHSD(GG_PWY_SumAkGk_aov)$treatment)

GG_PWY_OTUrel_sum_tukey <- GG_PWY_OTUrel_sum_tukey[GG_PWY_OTUrel_sum_tukey$p.adj < 0.05,]

GG_PWY_OTUrel_sum_tukey$p.adj_mark <- "*"
GG_PWY_OTUrel_sum_tukey$p.adj_mark[which(GG_PWY_OTUrel_sum_tukey$p.adj < 0.01)] <- "**"
GG_PWY_OTUrel_sum_tukey$p.adj_mark[which(GG_PWY_OTUrel_sum_tukey$p.adj < 0.001)] <- "***"

GG_PWY_OTUrel_sum_tukey

GG_PWY_OTUrel_sum_tukey <- 
  cbind(GG_PWY_OTUrel_sum_tukey,
        matrix(
          unlist(
            strsplit(rownames(GG_PWY_OTUrel_sum_tukey), split="-")
          ),
          ncol=2, 
          byrow=T
        )
  )

colnames(GG_PWY_OTUrel_sum_tukey)[6:7] <- c("x1", "x2")

GG_PWY_OTUrel_sum_tukey$x1 <- factor(GG_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
GG_PWY_OTUrel_sum_tukey$x2 <- factor(GG_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
GG_PWY_OTUrel_sum_tukey$y <- 6:1*0.12+8.2

GG_PWY_bar_P2 <- 
  ggplot(GG_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_segment(data=GG_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=GG_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=GG_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=GG_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label=p.adj_mark))+
  annotate(geom="text", x=4, y=9, label=expression(paste("ANOVA, ", italic(P)==3.356%*%10^-6)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(6.3, 9.2))

#----------------------------------------------
MEVIP_PWY_bar_P <- 
  ggplot(IP_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=Abundance*100))+
  #  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(data=IP_PWY_OTU_list$PWY_OTU_rel_individual, aes(fill=sum_score), 
               color="white", size=0.1, geom="bar", fun.data = "mean_se", position=position_stack())+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_jitter()+
  ylab("Relative abuncadance\nof bacteria contributing \nIsoprene biosynthesis II\nand mevalonate pathway (%) ")+
  xlab("")+
  theme_classic()+
  labs(fill="number\nof assosicated genes")+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0, 103))

for(i in 1:6){
  print(
    shapiro.test(
      subset(IP_PWY_OTU_list$PWY_OTUrel_sum, treatment == unique(IP_PWY_OTU_list$PWY_OTUrel_sum$treatment)[i])$SumAkGk
    )
  )
}

bartlett.test(SumAkGk ~ treatment, IP_PWY_OTU_list$PWY_OTUrel_sum)


IP_PWY_SumAkGk_aov <- aov(SumAkGk ~ treatment, IP_PWY_OTU_list$PWY_OTUrel_sum)
unlist(summary(IP_PWY_SumAkGk_aov))

IP_PWY_OTUrel_sum_tukey <- data.frame(TukeyHSD(IP_PWY_SumAkGk_aov)$treatment)

IP_PWY_OTUrel_sum_tukey <- IP_PWY_OTUrel_sum_tukey[IP_PWY_OTUrel_sum_tukey$p.adj < 0.05,]

IP_PWY_OTUrel_sum_tukey$p.adj_mark <- "*"
IP_PWY_OTUrel_sum_tukey$p.adj_mark[which(IP_PWY_OTUrel_sum_tukey$p.adj < 0.01)] <- "**"
IP_PWY_OTUrel_sum_tukey$p.adj_mark[which(IP_PWY_OTUrel_sum_tukey$p.adj < 0.001)] <- "***"

IP_PWY_OTUrel_sum_tukey

IP_PWY_OTUrel_sum_tukey <- 
  cbind(IP_PWY_OTUrel_sum_tukey,
        matrix(
          unlist(
            strsplit(rownames(IP_PWY_OTUrel_sum_tukey), split="-")
          ),
          ncol=2, 
          byrow=T
        )
  )

colnames(IP_PWY_OTUrel_sum_tukey)[6:7] <- c("x1", "x2")

IP_PWY_OTUrel_sum_tukey$x1 <- factor(IP_PWY_OTUrel_sum_tukey$x1, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
IP_PWY_OTUrel_sum_tukey$x2 <- factor(IP_PWY_OTUrel_sum_tukey$x2, levels = c("Bulk soil", "Control","NH₄Cl", "Glu", "Asp", "Asn", "Val"))
IP_PWY_OTUrel_sum_tukey$y <- 2:1*0.12+8.2

MEVIP_PWY_bar_P2 <- 
  ggplot(IP_PWY_OTU_list$PWY_OTUrel_sum, aes(x=treatment, y=SumAkGk))+
  stat_summary(geom="bar", fun.data = "mean_se")+
  stat_summary(geom="errorbar", fun.data = "mean_se", width=0.5)+
  geom_segment(data=IP_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x2, y=y, yend=y))+
  geom_segment(data=IP_PWY_OTUrel_sum_tukey, aes(x=x1, xend=x1, y=y, yend=y-0.025))+
  geom_segment(data=IP_PWY_OTUrel_sum_tukey, aes(x=x2, xend=x2, y=y, yend=y-0.025))+
  geom_text(data=IP_PWY_OTUrel_sum_tukey, aes(x=(as.numeric(x1)+as.numeric(x2))/2, y=y, label=p.adj_mark))+
  annotate(geom="text", x=4, y=8.8, label=expression(paste("ANOVA, ", italic(P)==0.013)), size=5)+
  geom_jitter()+
  ylab("GOO index")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(6, 9))

#----------------------------------------------
ggarrange(AST_PWY_bar_P, AST_PWY_bar_P2,
          AGR_PWY_bar_P, AGR_PWY_bar_P2,
          DG_PWY_bar_P, DG_PWY_bar_P2,
          PR_PWY_bar_P, PR_PWY_bar_P2,
          P5P_PWY_bar_P, P5P_PWY_bar_P2, nrow=5, ncol=2, legend = "top", labels=c("A","","B","","C","","D","","E",""))

ggarrange(GG_PWY_bar_P, GG_PWY_bar_P2,
          MEVIP_PWY_bar_P, MEVIP_PWY_bar_P2,
          nrow=2, ncol=2, legend = "top", labels=c("A","","B",""))

#======================
important_PWY_OTU <- 
intersect(
  intersect(
    intersect(AST_PWY_OTU_list$OTU, 
              AGR_PWY_OTU_list$OTU),
    intersect(DG_PWY_OTU_list$OTU, 
              PR_PWY_OTU_list$OTU)
  ),
  intersect(
    intersect(P5P_PWY_OTU_list$OTU, 
              GG_PWY_OTU_list$OTU),
    intersect(IP_PWY_OTU_list$OTU, 
              MEV_PWY_OTU_list$OTU)
  )
)


relative_abundance%>%
  group_by(treatment)%>%
  summarise(sum(Abundance))

datatable(TAX[important_PWY_OTU,], options = list(pageLength = 100))

important_PWY_relative_abundance <- 
  relative_abundance[which(relative_abundance$OTU %in% important_PWY_OTU,),] %>%
  group_by(family, Sample, treatment)%>%
  summarise(Abundance=sum(Abundance))
important_PWY_relative_abundance 

important_PWY_relative_abundance$treatment <- factor(important_PWY_relative_abundance$treatment, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

for(i in 1:6){
  for(j in 1:7){
    shapiro.test(
      subset(important_PWY_relative_abundance, family == unique(important_PWY_relative_abundance$family)[i]& treatment == unique(important_PWY_relative_abundance$treatment)[j])$Abundance
    )
  }
}

imp_family_A <- c()
pvalues_A<- c()
for(i in 1:6){
  fam <- unique(important_PWY_relative_abundance$family)[i]
  p <- kruskal.test(Abundance ~ treatment, subset(important_PWY_relative_abundance, family == fam))$p.value
  imp_family_A <- c(imp_family_A, fam)
  pvalues_A <- c(pvalues_A, p)
}

important_PWY_A_kruskal <- data.frame(imp_family_A, pvalues_A)
important_PWY_A_kruskal$padj <- p.adjust(important_PWY_A_kruskal$pvalues, method="fdr")
important_PWY_A_kruskal$padj_expression <-paste("italic(P)[adj] ==", round(important_PWY_A_kruskal$padj,3))
colnames(important_PWY_A_kruskal)[1] <- "family"


p_df <- data.frame(NA, NA, NA, NA, NA,NA, NA)
colnames(p_df) <- c("Comparison", "Estimator", "Lower", "Upper", "Statistic", "p.Value", "family")

for(i in 1:6){
  fam <- unique(important_PWY_relative_abundance$family)[i]
  p<- data.frame(nparcomp(Abundance ~ treatment, subset(important_PWY_relative_abundance, family == fam))$Analysis)
  p$family <- fam
  p_df<-rbind(p_df, p)
}
p_df <- na.omit(p_df)

Comparison <- gsub("[()]","", matrix(unlist(strsplit(p_df$Comparison, split=" , ")), ncol=2, byrow = T))
Comparison <- gsub("^p","", Comparison)
Comparison <- gsub("^ ","", Comparison)
Comparison <- gsub(" $","", Comparison)
colnames(Comparison) <- c("Comparison1", "Comparison2")

p_df<- cbind(p_df, Comparison)

p_df_sig_for_A <- subset(p_df, p.Value <= 0.05)
p_df_sig_for_A$Comparison1 <- factor(p_df_sig_for_A$Comparison1, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
p_df_sig_for_A$Comparison2 <- factor(p_df_sig_for_A$Comparison2, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

p_df_sig_for_A$y <- NA

p_df_sig_for_A[which(p_df_sig_for_A$family=="Burkholderiaceae"),]$y <- 9:1*0.03+0.5
p_df_sig_for_A[which(p_df_sig_for_A$family=="Enterobacteriaceae"),]$y <- 6:1*0.4+6
p_df_sig_for_A[which(p_df_sig_for_A$family=="Myxococcaceae"),]$y <- 4:1*0.004+0.07
p_df_sig_for_A[which(p_df_sig_for_A$family=="Oxalobacteraceae"),]$y <- 0.17
p_df_sig_for_A[which(p_df_sig_for_A$family=="Pseudomonadaceae"),]$y <- 8:1*1.3+35
p_df_sig_for_A[which(p_df_sig_for_A$family=="Sumerlaeaceae"),]$y <- 6:1 *0.001+0.025

p_df_sig_for_A$y2 <- NA

p_df_sig_for_A[which(p_df_sig_for_A$family=="Burkholderiaceae"),]$y2 <- 0.01
p_df_sig_for_A[which(p_df_sig_for_A$family=="Enterobacteriaceae"),]$y2 <- 0.15
p_df_sig_for_A[which(p_df_sig_for_A$family=="Myxococcaceae"),]$y2 <- 0.0011
p_df_sig_for_A[which(p_df_sig_for_A$family=="Oxalobacteraceae"),]$y2 <- 0.01
p_df_sig_for_A[which(p_df_sig_for_A$family=="Pseudomonadaceae"),]$y2 <- 0.5
p_df_sig_for_A[which(p_df_sig_for_A$family=="Sumerlaeaceae"),]$y2 <- 0.00035



p_df_sig_for_A$x3 <- (as.numeric(p_df_sig_for_A$Comparison1)+as.numeric(p_df_sig_for_A$Comparison2))/2

ggplot(important_PWY_relative_abundance, aes(treatment, Abundance*100, group=NULL))+
  stat_summary(geom="bar")+
  stat_summary(geom="errorbar", width=0.5)+
  geom_jitter()+
  geom_text(data = important_PWY_A_kruskal, aes(x=6, y=Inf, label=padj_expression), parse = T, hjust=0.5, vjust=2)+
  geom_segment(data=p_df_sig_for_A, aes(x=Comparison1, xend=Comparison2, y=y, yend=y))+
  geom_segment(data=p_df_sig_for_A, aes(x=Comparison1, xend=Comparison1, y=y, yend=y-y2/2))+
  geom_segment(data=p_df_sig_for_A, aes(x=Comparison2, xend=Comparison2, y=y, yend=y-y2/2))+
  geom_text(data=p_df_sig_for_A, aes(x=x3, y=y, label="***"))+
  facet_wrap(.~family, scale="free")+
  ylab("Relative abundance (%)")+
  xlab("")+
  theme(panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"))

#-------GOO----------

important_PWY_relative_abundance_OTU <- 
  relative_abundance[which(relative_abundance$OTU %in% important_PWY_OTU,),] %>%
  group_by(family, Sample, treatment, OTU)%>%
  summarise(Abundance=sum(Abundance))

important_PWY <- unique(c(AST_PWY, ARGORNPROST_PWY, DHGLUCONATE_PYR_CAT_PWY, PWY0_1338, PYRIDOXSYN_PWY, GG_PW, IP_PW, MEV_PW))

KO_OTU_important <- KO_OTU[which(KO_OTU$sequence%in%important_PWY_OTU), 
       c(1, which(colnames(KO_OTU)%in%important_PWY))]

important_OTU_Gk <- rowSums(KO_OTU_important[,-1])
names(important_OTU_Gk) <- KO_OTU_important$sequence

important_PWY_relative_abundance_OTU$Gk <- NA

for(k in names(important_OTU_Gk)){
  OTU_k <- which(important_PWY_relative_abundance_OTU$OTU == k)
  important_PWY_relative_abundance_OTU$Gk[OTU_k] <-  important_OTU_Gk[k]
}

important_PWY_relative_abundance_OTU$AkGk <- important_PWY_relative_abundance_OTU$Gk * important_PWY_relative_abundance_OTU$Abundance

important_PWY_relative_SumAkGk <- 
  important_PWY_relative_abundance_OTU %>%
  group_by(family, Sample, treatment) %>%
  summarise(SumAkGk=sum(AkGk))

important_PWY_relative_SumAkGk$treatment <- factor(important_PWY_relative_SumAkGk$treatment, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

for(i in 1:6){
  for(j in 1:7){
    shapiro.test(
    subset(important_PWY_relative_SumAkGk, family == unique(important_PWY_relative_SumAkGk$family)[i]& treatment == unique(important_PWY_relative_SumAkGk$treatment)[j])$SumAkGk
    )
  }
}

imp_family <- c()
pvalues<- c()
for(i in 1:6){
  fam <- unique(important_PWY_relative_SumAkGk$family)[i]
  p <- kruskal.test(SumAkGk ~ treatment, subset(important_PWY_relative_SumAkGk, family == fam))$p.value
  imp_family <- c(imp_family, fam)
  pvalues <- c(pvalues, p)
}

important_PWY_kruskal <- data.frame(imp_family, pvalues)
important_PWY_kruskal$padj <- p.adjust(important_PWY_kruskal$pvalues, method="fdr")
important_PWY_kruskal$padj_expression <-paste("italic(P)[adj] ==", round(important_PWY_kruskal$padj,3))
colnames(important_PWY_kruskal)[1] <- "family"

p_df <- data.frame(NA, NA, NA, NA, NA,NA, NA)
colnames(p_df) <- c("Comparison", "Estimator", "Lower", "Upper", "Statistic", "p.Value", "family")

for(i in 1:6){
  fam <- unique(important_PWY_relative_SumAkGk$family)[i]
  p<- data.frame(nparcomp(SumAkGk ~ treatment, subset(important_PWY_relative_SumAkGk, family == fam))$Analysis)
  p$family <- fam
  p_df<-rbind(p_df, p)
}
p_df <- na.omit(p_df)

Comparison <- gsub("[()]","", matrix(unlist(strsplit(p_df$Comparison, split=" , ")), ncol=2, byrow = T))
Comparison <- gsub("^p","", Comparison)
Comparison <- gsub("^ ","", Comparison)
Comparison <- gsub(" $","", Comparison)
colnames(Comparison) <- c("Comparison1", "Comparison2")

p_df<- cbind(p_df, Comparison)

p_df_sig <- subset(p_df, p.Value <= 0.05)
p_df_sig$Comparison1 <- factor(p_df_sig$Comparison1, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))
p_df_sig$Comparison2 <- factor(p_df_sig$Comparison2, levels = c("Bulk soil", "Control", "NH₄Cl", "Glu", "Asp", "Asn", "Val"))

p_df_sig$y <- NA

p_df_sig[which(p_df_sig$family=="Burkholderiaceae"),]$y <- 9:1*0.02+0.2
p_df_sig[which(p_df_sig$family=="Enterobacteriaceae"),]$y <- 6:1*0.3+2
p_df_sig[which(p_df_sig$family=="Myxococcaceae"),]$y <- 4:1*0.0022+0.021
p_df_sig[which(p_df_sig$family=="Oxalobacteraceae"),]$y <- 0.17
p_df_sig[which(p_df_sig$family=="Pseudomonadaceae"),]$y <- 8:1*1+15
p_df_sig[which(p_df_sig$family=="Sumerlaeaceae"),]$y <- 6:1 *0.0007+0.008

p_df_sig$y2 <- NA

p_df_sig[which(p_df_sig$family=="Burkholderiaceae"),]$y2 <- 0.01
p_df_sig[which(p_df_sig$family=="Enterobacteriaceae"),]$y2 <- 0.15
p_df_sig[which(p_df_sig$family=="Myxococcaceae"),]$y2 <- 0.0011
p_df_sig[which(p_df_sig$family=="Oxalobacteraceae"),]$y2 <- 0.01
p_df_sig[which(p_df_sig$family=="Pseudomonadaceae"),]$y2 <- 0.5
p_df_sig[which(p_df_sig$family=="Sumerlaeaceae"),]$y2 <- 0.00035

p_df_sig$x3 <-(as.numeric(p_df_sig$Comparison1)+as.numeric(p_df_sig$Comparison2))/2

ggplot(important_PWY_relative_SumAkGk, aes(treatment, SumAkGk, group=NULL))+
  stat_summary(geom="bar")+
  stat_summary(geom="errorbar", width=0.5)+
  geom_jitter()+
  geom_text(data = important_PWY_kruskal, aes(x=6, y=Inf, label=padj_expression), parse = T, hjust=0.5, vjust=2)+
  geom_segment(data=p_df_sig, aes(x=Comparison1, xend=Comparison2, y=y, yend=y))+
  geom_segment(data=p_df_sig, aes(x=Comparison1, xend=Comparison1, y=y, yend=y-y2/2))+
  geom_segment(data=p_df_sig, aes(x=Comparison2, xend=Comparison2, y=y, yend=y-y2/2))+
  geom_text(data=p_df_sig, aes(x=x3, y=y, label="***"))+
  facet_wrap(.~family, scale="free")+
  ylab("GOO index of important pathways")+
  xlab("")+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        panel.background = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

