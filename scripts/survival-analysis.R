## Install packages for Kaplan-Meier plots
# install.packages("survivalAnalysis")
# install.packages("survminer")

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
library("plotrix")     
library(AICcmodavg)
library(DescTools)
library(ggpattern)
library("BSDA")

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Load survival csv file
data = read.csv("lifespan-data.csv", header = T)

## Define a function that will create a survival/event table
## Per row will be a biological sample that either survived (0) or dead (1).
## No censored will be applied here. 

data$group = ifelse(data$Treatment %in% c('EGCG', 'E. coli', 'DMSO'),
                    data$Treatment,
                    paste(data$Treatment,"-", substr(data$Concentration,1,1)))

to_survival_from_table = function(dat) {
  newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','group', 'Trial'), drop = FALSE], row.names=NULL)
  newdat$group = factor(newdat$group, levels=c("E. coli", "DMSO", "EGCG",
                                               "Itogon - L", "Itogon - M", "Itogon - H",
                                               "Kapangan - L", "Kapangan - M", "Kapangan - H",
                                               "La Trinidad - L", "La Trinidad - M", "La Trinidad - H",
                                               "Sablan - L", "Sablan - M", "Sablan - H"),
                        labels=c("E. coli", "DMSO", "EGCG",
                                 "Itogon - L", "Itogon - M", "Itogon - H",
                                 "Kapangan - L", "Kapangan - M", "Kapangan - H",
                                 "La Trinidad - L", "La Trinidad - M", "La Trinidad - H",
                                 "Sablan - L", "Sablan - M", "Sablan - H"))
  newdat$status = 1
  return(newdat)
}

to_survival_from_table_2 = function(dat) {
  newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','Concentration','Treatment', 'Trial'), drop = FALSE], row.names=NULL)
  newdat$Treatment= factor(newdat$Treatment, levels=c("E. coli", "DMSO", "EGCG",
                                                      "Itogon", 
                                                      "Kapangan",
                                                      "La Trinidad", 
                                                      "Sablan"), 
                           labels=c("E. coli", "DMSO", "EGCG",
                                    "Itogon",
                                    "Kapangan",
                                    "La Trinidad",
                                    "Sablan"))
  newdat$Concentration= factor(newdat$Concentration, levels=c("Control", "Low", "Mid", "High"), 
                               labels=c("Control", "Low", "Mid", "High"))
  newdat$status = 1
  return(newdat)
}

newdat = to_survival_from_table(data)
newdat

## Fit survival function
fit <- survfit(Surv(newdat$Time, newdat$status) ~ group, data = newdat)
fit


color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")
g = ggsurvplot(fit, 
               # risk.table = TRUE, # Add risk table
               # risk.table.col = "group", # Change risk table color by groups
               # linetype = "strata", # Change line type by groups
               ggtheme = theme_classic(), # Change ggplot2 theme
               # palette = color_scheme,
               xlab="Time in days",
               size = 0.8,
               legend = c(1.38, 0.50),
               legend.title = element_blank(),
               
               legend.labs = c( "OP50 <i>E. coli</i>", "0.1% DMSO", "200 uM EGCG",
                                "It - L", "It - M", "It - H",
                                "Kp - L", "Kp - M", "Kp - H",
                                "LT - L", "LT - M", "LT- H",
                                "Sa - L", "Sa - M", "Sa - H"),
               xlim = c(0,20),
               font.legend = c(9, "plain", "black"),
               font.x = c(15, 'plain', 'black'),
               font.y = c(15, 'plain', 'black'),
               font.tickslab = c(11, 'plain', 'black'),
               surv.scale = 'percent'
) + guides(colour = guide_legend(nrow = 5, byrow = TRUE))

g$plot = g$plot + theme(plot.margin = unit(c(0, 8.5, 0, 0), "cm"),
                        legend.text = element_markdown())

# ggsave("survival_plot.png", plot = g, dpi = 600)

png(filename = file.path(getwd(),"/Lifespan-km-plot.png"), width = 8, height = 3.5, units = "in", res = 600)
plot(g$plot)
dev.off()
print(paste("Saved at ", getwd()))
print(g$plot)


## Plot mean lifetime
## Prepare dataframe with mean and standard dev
newdat.p = to_survival_from_table_2(data)
newdat.p$lifetime = newdat.p$Time -1
df = newdat.p %>% group_by(Treatment, Concentration) %>% summarise(mean = mean(lifetime), stdv = std.error(lifetime), count = n())
df$upper = df$mean + df$stdv
df$lower = df$mean - df$stdv


## STATISTICS


max_y = (max(df$mean) + 1.5)

control.tox = df[df$Concentration == "Control", ] 
y1 = ggplot(control.tox, aes(Treatment, mean, fill=Concentration)) +
  theme_classic() + labs(x = '', y = 'Mean Lifespan (Days)', title="Control") +
  
  geom_bar(position=position_dodge(), stat="identity", colour='black') +  
  scale_y_continuous(limits = c(0, max_y)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>",
                              "0.1% <br>DMSO",
                              "200 uM <br>EGCG")) +
  scale_fill_manual(values=color_scheme[1]) +
  geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
  # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
  scale_color_manual(values="#960101") +
  theme(legend.position="none", axis.text.x = element_markdown(size=8),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"))
#       plot.title = element_text(hjust = 0.5),
#       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# ) 

y1




control.tox = df[df$Treatment == "Itogon", ] 
y2 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
  theme_classic() + labs(x = '', y = '', title="Itogon") +
  scale_y_continuous(limits = c(0, max_y)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +  
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  scale_fill_manual(values=color_scheme[2]) +
  geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
  # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
  scale_color_manual(values="#960101") +
  theme(legend.position="none", axis.text.x = element_text(size=8),
        # axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"))
#       plot.title = element_text(hjust = 0.5),
#       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# ) 

y2

control.tox = df[df$Treatment == "Kapangan", ] 
y3 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
  theme_classic() + labs(x = '', y = '', title="Kapangan") +
  
  geom_bar(position=position_dodge(), stat="identity", colour='black') +  
  scale_y_continuous(limits = c(0,max_y)) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  scale_fill_manual(values=color_scheme[3]) +
  geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
  # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
  scale_color_manual(values="#960101") +
  theme(legend.position="none", axis.text.x = element_text(size=8),
        # axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"))
y3

control.tox = df[df$Treatment == "La Trinidad", ] 
y4 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
  theme_classic() + labs(x = '', y = '', title="La Trinidad") +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +  
  scale_y_continuous(limits = c(0,max_y)) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  scale_fill_manual(values=color_scheme[4]) +
  geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
  # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
  scale_color_manual(values="#960101") +
  theme(legend.position="none", axis.text.x = element_text(size=8),
        # axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0.0,0,0,0), "cm"))
#       plot.title = element_text(hjust = 0.5),
#       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# ) 

y4

control.tox = df[df$Treatment == "Sablan", ] 
y5 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
  theme_classic() + labs(x = '', y = '', title="Sablan") +
  
  geom_bar(position=position_dodge(), stat="identity", colour='black') +  
  scale_y_continuous(limits = c(0,max_y)) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  scale_fill_manual(values=color_scheme[5]) +
  geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
  # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
  scale_color_manual(values="#960101") +
  theme(legend.position="none", axis.text.x = element_text(size=8),
        plot.title = element_text(hjust = 0.5),
        # axis.text.y = element_blank(),
        plot.margin = unit(c(0.0,0,0,0), "cm"))
#       plot.title = element_text(hjust = 0.5),
#       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
#       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# ) 
y5


newdat$lifetime = newdat$Time -1
res = DunnettTest(x=newdat$lifetime, g=newdat$group)
print(res)
write.csv(as.data.frame(res$`E. coli`), file="Lifespan-stat-results.csv")
write.csv(df, file="Lifespan-values.csv")

y1 = y1 + annotate("text", x = 3, y = 13.65, label = "***", size = 7, color = "#960101") 
y2 = y2 + annotate("text", x = 3, y = 13.0, label = "*", size = 7, color = "#960101")


short.term.x = plot_grid(y1, y2, y3, y4, y5, align = "h", ncol = 5, rel_widths = c(0.25, 0.20, 0.20,0.20,0.20)) +
  theme(plot.margin = unit(c(0.0,0,0,0), "cm")) 
short.term.x

png(filename = file.path(getwd(),"/Lifespan-barplot-sig.png"), width = 8, height = 3.5, units = "in", res = 600)

plot(short.term.x)
dev.off()
print(paste("Saved at ", getwd()))

# 
# 
# ###########################################################################
# ######################### HEAT STRESS #####################################
# ###########################################################################
# ## Load survival csv file
# data = read.csv("heat-stress-data.csv", header = T)
# 
# ## Define a function that will create a survival/event table
# ## Per row will be a biological sample that either survived (0) or dead (1).
# ## No censored will be applied here. 
# 
# data$group = ifelse(data$Treatment %in% c('EGCG', 'E. coli', 'DMSO'),
#                     data$Treatment,
#                     paste(data$Treatment,"-", substr(data$Concentration,1,1)))
# 
# to_survival_from_table = function(dat) {
#   newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','group', 'Trial'), drop = FALSE], row.names=NULL)
#   newdat$group = factor(newdat$group, levels=c("E. coli", "DMSO", "EGCG",
#                                                "Itogon - L", "Itogon - M", "Itogon - H",
#                                                "La Trinidad - L", "La Trinidad - M", "La Trinidad - H"),
#                         labels=c("E. coli", "DMSO", "EGCG",
#                                  "Itogon - L", "Itogon - M", "Itogon - H",
#                                  "La Trinidad - L", "La Trinidad - M", "La Trinidad - H"))
#   newdat$status = 1
#   return(newdat)
# }
# 
# to_survival_from_table_2 = function(dat) {
#   newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','Concentration','Treatment', 'Trial'), drop = FALSE], row.names=NULL)
#   newdat$Treatment= factor(newdat$Treatment, levels=c("E. coli", "DMSO", "EGCG",
#                                                       "Itogon", "La Trinidad"), 
#                            labels=c("E. coli", "DMSO", "EGCG",
#                                     "Itogon",
#                                     "La Trinidad"))
#   newdat$Concentration= factor(newdat$Concentration, levels=c("Control", "Low", "Mid", "High"), 
#                                labels=c("Control", "Low", "Mid", "High"))
#   newdat$status = 1
#   return(newdat)
# }
# 
# newdat = to_survival_from_table(data)
# newdat
# 
# ## Fit survival function
# fit <- survfit(Surv(newdat$Time, newdat$status) ~ group, data = newdat)
# fit
# 
# ## Graph the K-M plot.
# color_scheme = c("#AA8B56",  "#9BA17B", "61764B","#4E6C50", "#395144")
# g = ggsurvplot(fit, 
#                # risk.table = TRUE, # Add risk table
#                # risk.table.col = "group", # Change risk table color by groups
#                # linetype = "strata", # Change line type by groups
#                ggtheme = theme_classic(), # Change ggplot2 theme
#                # palette = color_scheme,
#                xlab="Time in days",
#                size = 0.8,
#                legend = c(0.30, 0.18),
#                legend.title = element_blank(),
#                legend.labs = c("E. coli", "DMSO", "EGCG",
#                                "Itogon - L", "Itogon - M", "Itogon - H",
#                                "La Trinidad - L", "La Trinidad - M", "La Trinidad - H"),
#                xlim = c(0,17),
#                font.legend = c(10, "plain", "black"),
#                font.x = c(15, 'plain', 'black'),
#                font.y = c(15, 'plain', 'black'),
#                font.tickslab = c(11, 'plain', 'black'),
#                surv.scale = 'percent'
# ) + guides(colour = guide_legend(nrow = 3)) 
# g
# 
# 
# ## Plot mean lifetime
# ## Prepare dataframe with mean and standard dev
# newdat.p = to_survival_from_table_2(data)
# newdat.p$lifetime = newdat.p$Time -1
# df = newdat.p %>% group_by(Treatment, Concentration) %>% summarise(mean = mean(lifetime), stdv = std.error(lifetime), count = n())
# df$upper = df$mean + df$stdv
# df$lower = df$mean - df$stdv
# 
# control.tox = df[df$Concentration == "Control", ] 
# y1 = ggplot(control.tox, aes(Treatment, mean, fill=Concentration)) +
#   theme_classic() + labs(x = '', y = 'Mean Lifespan', title="Control") +
#   
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +  
#   scale_y_continuous(limits = c(0,15)) +
#   scale_x_discrete(labels = c("E. coli", "DMSO", "EGCG")) +
#   scale_fill_manual(values=color_scheme[1]) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
#   # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
#   scale_color_manual(values="#960101") +
#   theme(legend.position="none", axis.text.x = element_text(size=8),
#         plot.title = element_text(hjust = 0.5))
# #       plot.title = element_text(hjust = 0.5),
# #       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
# #       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# # ) 
# 
# y1
# 
# control.tox = df[df$Treatment == "Itogon", ] 
# y2 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
#   theme_classic() + labs(x = '', y = '', title="Itogon") +
#   
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +  
#   scale_y_continuous(limits = c(0,15), labels = c('','','')) +
#   scale_x_discrete(labels = c("Low", "Mid", "High")) +
#   scale_fill_manual(values=color_scheme[2]) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
#   # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
#   scale_color_manual(values="#960101") +
#   theme(legend.position="none", axis.text.x = element_text(size=8),
#         plot.title = element_text(hjust = 0.5))
# #       plot.title = element_text(hjust = 0.5),
# #       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
# #       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# # ) 
# 
# y2
# 
# 
# control.tox = df[df$Treatment == "La Trinidad", ] 
# y4 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
#   theme_classic() + labs(x = '', y = '', title="La Trinidad") +
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +  
#   scale_y_continuous(limits = c(0,15), labels = c('','','')) +
#   scale_x_discrete(labels = c("Low", "Mid", "High")) +
#   scale_fill_manual(values=color_scheme[4]) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
#   # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
#   scale_color_manual(values="#960101") +
#   theme(legend.position="none", axis.text.x = element_text(size=8),
#         plot.title = element_text(hjust = 0.5))
# #       plot.title = element_text(hjust = 0.5),
# #       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
# #       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# # ) 
# 
# y4
# 
# 
# newdat$lifetime = newdat$Time -1
# DunnettTest(x=newdat$lifetime, g=newdat$group)
# 
# y1 = y1 + annotate("text", x = 3, y = 12.2, label = "***", size = 7, color = "#960101")
# y2 = y2 + annotate("text", x = 1, y = 11.5, label = "***", size = 7, color = "#960101") +
#   annotate("text", x = 2, y = 11.0, label = "*", size = 7, color = "#960101") +
#   annotate("text", x = 3, y = 11.4, label = "**", size = 7, color = "#960101")
# y4 = y4 + annotate("text", x = 3, y = 11.1, label = "***", size = 7, color = "#960101")
# 
# 
# short.term.x = plot_grid(y1, y2, y4, align = "h", ncol = 3, rel_widths = c(0.33, 0.33, 0.33)) +
#   theme(plot.margin = unit(c(0.7,0,0,0), "cm")) 
# short.term.x
# 
# png(filename = file.path(getwd(),"/Heat-stress-barplot-sig.png"), width = 8, height = 3.5, units = "in", res = 600)
# 
# plot(short.term.x)
# dev.off()
# print(paste("Saved at ", getwd()))
# 
# png(filename = file.path(getwd(),"/Heat-stress-km-plot.png"), width = 8, height = 5, units = "in", res = 600)
# 
# plot(g$plot)
# dev.off()
# print(paste("Saved at ", getwd()))
# 
# ###########################################################################
# ######################### UV STRESS #####################################
# ###########################################################################
# ## Load survival csv file
# data = read.csv("UV-stress-data.csv", header = T)
# 
# ## Define a function that will create a survival/event table
# ## Per row will be a biological sample that either survived (0) or dead (1).
# ## No censored will be applied here. 
# 
# data$group = ifelse(data$Treatment %in% c('EGCG', 'E. coli', 'DMSO'),
#                     data$Treatment,
#                     paste(data$Treatment,"-", substr(data$Concentration,1,1)))
# 
# to_survival_from_table = function(dat) {
#   newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','group', 'Trial'), drop = FALSE], row.names=NULL)
#   newdat$group = factor(newdat$group, levels=c("E. coli", "DMSO", "EGCG",
#                                                "Itogon - L", "Itogon - M", "Itogon - H",
#                                                "La Trinidad - L", "La Trinidad - M", "La Trinidad - H"),
#                         labels=c("E. coli", "DMSO", "EGCG",
#                                  "Itogon - L", "Itogon - M", "Itogon - H",
#                                  "La Trinidad - L", "La Trinidad - M", "La Trinidad - H"))
#   newdat$status = 1
#   return(newdat)
# }
# 
# to_survival_from_table_2 = function(dat) {
#   newdat = data.frame(dat[rep(seq_len(dim(dat)[1]), dat$Dead), c('Time','Dead','Concentration','Treatment', 'Trial'), drop = FALSE], row.names=NULL)
#   newdat$Treatment= factor(newdat$Treatment, levels=c("E. coli", "DMSO", "EGCG",
#                                                       "Itogon", "La Trinidad"), 
#                            labels=c("E. coli", "DMSO", "EGCG",
#                                     "Itogon",
#                                     "La Trinidad"))
#   newdat$Concentration= factor(newdat$Concentration, levels=c("Control", "Low", "Mid", "High"), 
#                                labels=c("Control", "Low", "Mid", "High"))
#   newdat$status = 1
#   return(newdat)
# }
# 
# newdat = to_survival_from_table(data)
# newdat
# 
# ## Fit survival function
# fit <- survfit(Surv(newdat$Time, newdat$status) ~ group, data = newdat)
# fit
# 
# ## Graph the K-M plot.
# color_scheme = c("#AA8B56",  "#9BA17B", "61764B","#4E6C50", "#395144")
# g = ggsurvplot(fit, 
#                # risk.table = TRUE, # Add risk table
#                # risk.table.col = "group", # Change risk table color by groups
#                # linetype = "strata", # Change line type by groups
#                ggtheme = theme_classic(), # Change ggplot2 theme
#                # palette = color_scheme,
#                xlab="Time in days",
#                size = 0.8,
#                legend = c(0.27, 0.18),
#                legend.title = element_blank(),
#                legend.labs = c("E. coli", "DMSO", "EGCG",
#                                "Itogon - L", "Itogon - M", "Itogon - H",
#                                "La Trinidad - L", "La Trinidad - M", "La Trinidad - H"),
#                xlim = c(0,8),
#                font.legend = c(10, "plain", "black"),
#                font.x = c(15, 'plain', 'black'),
#                font.y = c(15, 'plain', 'black'),
#                font.tickslab = c(11, 'plain', 'black'),
#                surv.scale = 'percent'
# ) + guides(colour = guide_legend(nrow = 3)) 
# g
# 
# 
# ## Plot mean lifetime
# ## Prepare dataframe with mean and standard dev
# newdat.p = to_survival_from_table_2(data)
# newdat.p$lifetime = newdat.p$Time -1
# df = newdat.p %>% group_by(Treatment, Concentration) %>% summarise(mean = mean(lifetime), stdv = std.error(lifetime), count = n())
# df$upper = df$mean + df$stdv
# df$lower = df$mean - df$stdv
# 
# control.tox = df[df$Concentration == "Control", ] 
# y1 = ggplot(control.tox, aes(Treatment, mean, fill=Concentration)) +
#   theme_classic() + labs(x = '', y = 'Mean Lifespan', title="Control") +
#   
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +  
#   scale_y_continuous(limits = c(0,4)) +
#   scale_x_discrete(labels = c("E. coli", "DMSO", "EGCG")) +
#   scale_fill_manual(values=color_scheme[1]) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
#   # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
#   scale_color_manual(values="#960101") +
#   theme(legend.position="none", axis.text.x = element_text(size=8),
#         plot.title = element_text(hjust = 0.5))
# #       plot.title = element_text(hjust = 0.5),
# #       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
# #       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# # ) 
# 
# y1
# 
# control.tox = df[df$Treatment == "Itogon", ] 
# y2 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
#   theme_classic() + labs(x = '', y = '', title="Itogon") +
#   
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +  
#   scale_y_continuous(limits = c(0,4), labels = c('','','','', '')) +
#   scale_x_discrete(labels = c("Low", "Mid", "High")) +
#   scale_fill_manual(values=color_scheme[2]) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
#   # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
#   scale_color_manual(values="#960101") +
#   theme(legend.position="none", axis.text.x = element_text(size=8),
#         plot.title = element_text(hjust = 0.5))
# #       plot.title = element_text(hjust = 0.5),
# #       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
# #       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# # ) 
# 
# y2
# 
# 
# control.tox = df[df$Treatment == "La Trinidad", ] 
# y4 = ggplot(control.tox, aes(Concentration, mean, fill=Treatment)) +
#   theme_classic() + labs(x = '', y = '', title="La Trinidad") +
#   
#   geom_bar(position=position_dodge(), stat="identity", colour='black') +  
#   scale_y_continuous(limits = c(0,4), labels = c('','','','','')) +
#   scale_x_discrete(labels = c("Low", "Mid", "High")) +
#   scale_fill_manual(values=color_scheme[4]) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9))+
#   # geom_hline(aes(yintercept = healthy.baseline.mean, color = "U"), linetype='dashed', color='#960101') +
#   scale_color_manual(values="#960101") +
#   theme(legend.position="none", axis.text.x = element_text(size=8),
#         plot.title = element_text(hjust = 0.5))
# #       plot.title = element_text(hjust = 0.5),
# #       axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
# #       axis.text.y = element_text(color = c("black", "black", "black", "black", "#960101", "black"))
# # ) 
# 
# y4
# 
# 
# newdat$lifetime = newdat$Time -1
# DunnettTest(x=newdat$lifetime, g=newdat$group)
# 
# # y1 = y1 + annotate("text", x = 3, y = 12.2, label = "***", size = 7, color = "#960101")
# # y2 = y2 + annotate("text", x = 1, y = 11.5, label = "***", size = 7, color = "#960101") +
# #   annotate("text", x = 2, y = 11.0, label = "*", size = 7, color = "#960101") +
# #   annotate("text", x = 3, y = 11.4, label = "**", size = 7, color = "#960101")
# # y4 = y4 + annotate("text", x = 3, y = 11.1, label = "***", size = 7, color = "#960101")
# # 
# 
# short.term.x = plot_grid(y1, y2, y4, align = "h", ncol = 3, rel_widths = c(0.33, 0.33, 0.33)) +
#   theme(plot.margin = unit(c(0.7,0,0,0), "cm")) 
# short.term.x
# 
# png(filename = file.path(getwd(),"/UV-stress-barplot-sig.png"), width = 8, height = 3.5, units = "in", res = 600)
# 
# plot(short.term.x)
# dev.off()
# print(paste("Saved at ", getwd()))
# 
# png(filename = file.path(getwd(),"/UVs-stress-km-plot.png"), width = 8, height = 5, units = "in", res = 600)
# 
# plot(g$plot)
# dev.off()
# print(paste("Saved at ", getwd()))