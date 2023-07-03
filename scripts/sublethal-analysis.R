## Import packages
library(survival)
library("survminer")
library(dplyr)
library("ggplot2")
library(ggsignif)
library(ggpubr)
library(rstatix)
library(cowplot)
library("lemon")
library(reprex)
library(AICcmodavg)
library(DescTools)
library("plotrix")
library(data.table)
library("rockchalk")
library(truncnorm)
library("BSDA")
library(ggtext)

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

######################### READ FILES #############################
sublethal.csv = read.csv("sublethality-data.csv", header = T)
sublethal.csv

## Treatment only
sublethal.csv.treat = subset(sublethal.csv, Concentration != 'Control')

## Remove mg/ml in Concentration
sublethal.csv.treat$Concentration <- as.numeric(gsub(" mg/ml", "", 
                                                     sublethal.csv.treat$Concentration))

sublethal.csv.treat$Concentration <- factor(sublethal.csv.treat$Concentration,
                                      levels = c(50, 100, 250, 400, 500),
                                      labels = c("50 mg/ml", "100 mg/ml", "250 mg/ml", "400 mg/ml", "500 mg/ml"))

sublethal.csv.treat$Time <- factor(sublethal.csv.treat$Time, levels = c(0,1,2, 3,4),
                            labels = c("0h", "24h", "48h", "72h", "96h"))

## Initial count is 30
sublethal.csv.treat$Survival =  sublethal.csv.treat$Alive/30.0

sublethal.csv.summary <- sublethal.csv.treat %>% 
  group_by(Treatment, Concentration, Time) %>%
  summarise(Survival.mean = mean(Survival),
            Survival.std = std.error(Survival),
            obs.rawint = n()) %>%
  as.data.frame()

sublethal.csv.summary
sublethal.csv.summary = subset(sublethal.csv.summary, Time != "96h")
write.csv(sublethal.csv.summary, "Sublethal-values.csv")

## Visualize
# color_scheme = c("#92ffc0",  "#5285d1", "#002661")
# color_scheme = c("#62cff4",  "#548ecc", "#ff78d6") # Bubble gum
# color_scheme = c("#13c2fc",  "#7f40ff", "#ff78d6") # Bubble gum

## Itogon
sublethal.csv.summary.itogon = subset(sublethal.csv.summary, Treatment == "Itogon")
q1 = ( ggplot(sublethal.csv.summary.itogon, aes(x=Time, y=Survival.mean, fill=Concentration)) +
         theme_classic()  + labs(x = '', y = 'Survival Rate', title = "Itogon") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ Concentration', scales="free", switch = "x") + 
         geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         # scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               axis.text.x = element_text(size=9),
               axis.text.y = element_text(size=11, 
                                          color = c("black", "black", "black", "black", "#960101", "black")),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q1

png(filename = file.path(getwd(),"/Itogon.png"), width = 7.5, height = 3.0, units = "in", res = 600)
plot(q1)
dev.off()
print(paste("Saved at ", getwd()))

## Kapangan
## Kapangan
sublethal.csv.summary.kp = subset(sublethal.csv.summary, Treatment == "Kapangan")
q2 = ( ggplot(sublethal.csv.summary.kp, aes(x=Time, y=Survival.mean, fill=Concentration)) +
         theme_classic()  + labs(x = '', y = 'Survival Rate', title = "Kapangan") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ Concentration', scales="free", switch = "x") + 
         geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         # scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               axis.text.x = element_text(size=9),
               axis.text.y = element_text(size=11, 
                                          color = c("black", "black", "black", "black", "#960101", "black")),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q2

png(filename = file.path(getwd(),"/Kapangan.png"), width = 7.5, height = 3.0, units = "in", res = 600)
plot(q2)
dev.off()
print(paste("Saved at ", getwd()))

## Sablan
sublethal.csv.summary.sb = subset(sublethal.csv.summary, Treatment == "Sablan")
q3 = ( ggplot(sublethal.csv.summary.sb, aes(x=Time, y=Survival.mean, fill=Concentration)) +
         theme_classic()  + labs(x = '', y = 'Survival Rate', title = "Sablan") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ Concentration', scales="free", switch = "x") + 
         geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         # scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               axis.text.x = element_text(size=9),
               axis.text.y = element_text(size=11, 
                                          color = c("black", "black", "black", "black", "#960101", "black")),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q3

png(filename = file.path(getwd(),"/Sablan.png"), width = 7.5, height = 3.0, units = "in", res = 600)
plot(q3)
dev.off()
print(paste("Saved at ", getwd()))

## La Trinidad
sublethal.csv.summary.lt = subset(sublethal.csv.summary, Treatment == "La Trinidad")
q4 = ( ggplot(sublethal.csv.summary.lt, aes(x=Time, y=Survival.mean, fill=Concentration)) +
         theme_classic()  + labs(x = '', y = 'Survival Rate', title = "La Trinidad") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ Concentration', scales="free", switch = "x") + 
         geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         # scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               axis.text.x = element_text(size=9),
               axis.text.y = element_text(size=11, 
                                          color = c("black", "black", "black", "black", "#960101", "black")),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q4

png(filename = file.path(getwd(),"/La Trinidad.png"), width = 7.5, height = 3.0, units = "in", res = 600)
plot(q4)
dev.off()
print(paste("Saved at ", getwd()))

treatment.summary = sublethal.csv.summary

## Controls
sublethal.csv.ctrl = subset(sublethal.csv, Concentration == 'Control')

sublethal.csv.ctrl$Time <- factor(sublethal.csv.ctrl$Time, levels = c(0,1,2, 3,4),
                                   labels = c("0h", "24h", "48h", "72h", "96h"))

## Initial count is 30
sublethal.csv.ctrl$Survival =  sublethal.csv.ctrl$Alive/30.0

sublethal.csv.summary <- sublethal.csv.ctrl %>% 
  group_by(Treatment, Time) %>%
  summarise(Survival.mean = mean(Survival),
            Survival.std = std.error(Survival),
            obs.rawint = n()) %>%
  as.data.frame()


sublethal.csv.summary$Treatment <- factor(sublethal.csv.summary$Treatment,
                                          levels = c("E. coli", "DMSO", "EGCG"),
                                   labels = c("OP50 <i>E. coli<i>", "0.1% DMSO", "200 uM EGCG"))

sublethal.csv.summary = arrange(sublethal.csv.summary, Treatment)
sublethal.csv.summary = subset(sublethal.csv.summary, Time!="96h")
ctrls.summary = sublethal.csv.summary

## controls
q5 = ( ggplot(sublethal.csv.summary, aes(x=Time, y=Survival.mean, fill=Treatment)) +
         theme_classic()  + labs(x = '', y = 'Survival Rate', title = "Control") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ Treatment', scales="free", switch = "x") + 
         geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         # scale_fill_manual(values=color_scheme) +
         # scale_x_discrete(labels = c("OP50 <i>E. coli</i>",
         #                             "0.1% DMSO",
         #                             "200 uM EGCG")) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               axis.text.x = element_text(size=10),
               axis.text.y = element_text(size=11, 
                                          color = c("black", "black", "black", "black", "#960101", "black")),
               legend.position="none",
               strip.text.x = element_markdown(size=13, face = "bold"),
               strip.background = element_blank()) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q5

png(filename = file.path(getwd(),"/Control.png"), width = 6.0, height = 3.0, units = "in", res = 600)
plot(q5)
dev.off()
print(paste("Saved at ", getwd()))

## MERGE FILES
ctrls.summary$Concentration = ctrls.summary$Treatment
all.summary = rbind(ctrls.summary, subset(treatment.summary, Concentration == "100 mg/ml"))

color_scheme = c("OP50 <i>E. coli<i>" = "#f8766d", "0.1% DMSO"="#00ba38", "200 uM EGCG" = "#619cff", "100 mg/ml"="#a3a500" )
q6 = ( ggplot(all.summary, aes(x=Time, y=Survival.mean, fill=Concentration)) +
         theme_classic()  + labs(x = '', y = 'Survival Rate', title = "Best Concentration") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Survival.mean-Survival.std, ymax=Survival.mean+Survival.std), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ Treatment', scales="free", switch = "x") +
         geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         # scale_x_discrete(labels = c("OP50 <i>E. coli</i>",
         #                             "0.1% DMSO",
         #                             "200 uM EGCG")) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               axis.text.x = element_text(size=10),
               axis.text.y = element_text(size=11, 
                                          color = c("black", "black", "black", "black", "#960101", "black")),
               legend.position="none",
               strip.text.x = element_markdown(size=13, face = "bold"),
               strip.background = element_blank()) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q6

png(filename = file.path(getwd(),"/Best-Concentration.png"), width = 10.50, height = 3.5, units = "in", res = 600)
plot(q6)
dev.off()
print(paste("Saved at ", getwd()))