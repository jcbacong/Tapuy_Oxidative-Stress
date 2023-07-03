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
library(AICcmodavg)
library(DescTools)


## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Prepare group name
rename_wine_or_lees = function(x) {
  if(substr(x,start=1,stop=2) == "FW"){
    return("Field Wine")
  } else if (substr(x,start=1,stop=2) == "FL") {
    return("Field Lees")
  } else if(substr(x,start=1,stop=2) == "LW") {
    return("Lab Wine")
  } else if(substr(x,start=1,stop=2) == "LL") {
    return("Lab Lees")
  } else {
    return("Control")
  }
}

rename_concentration = function(x) {
  if(substr(x,start=3,stop=3) == "L"){
    return("Low")
  } else if (substr(x,start=3,stop=3) == "M") {
    return("Mid")
  } else if(substr(x,start=3,stop=3) == "H") {
    return("High")
  } else {
    return("Control")
  }
}

##########################################################################################################
########################################### SHORT TERM MEMORY ############################################
##########################################################################################################

## Load survival csv file
learn.csv = read.csv("short-term-learning.csv", header = T)
memory.csv = read.csv("short-term-memory.csv", header = T)

## Learning
learn.csv$CI_index = (learn.csv$N_butanone - learn.csv$N_EtOH) / ((learn.csv$N_total - learn.csv$N_Origin) + 1e-6)
learn.csv$Type = "Learning"

## Memory
memory.csv$CI_index = (memory.csv$N_butanone - memory.csv$N_EtOH) / ((memory.csv$N_total - memory.csv$N_Origin) + 1e-6)
memory.csv$Type = "Short-term Memory"

## Summarize
stm.csv = rbind(learn.csv,memory.csv)
stm.csv$Treatment = factor(stm.csv$Treatment, 
                           levels = c("Control", "Control", "Control","Itogon", "Kapangan", "La Trinidad", "Sablan"),
                           labels = c("Control", "Control", "Control","Itogon", "Kapangan", "La Trinidad", "Sablan"))
stm.csv$Concentration = factor(stm.csv$Concentration, 
                           levels = c("E. coli", "DMSO", "EGCG","Low", "Mid", "High"),
                           labels = c("E. coli", "DMSO", "EGCG","Low", "Mid", "High"))
stm.csv$Type = factor(stm.csv$Type,
                      levels=c("Learning", "Short-term Memory"),
                      labels=c("Learning", "Short-term Memory"))


stm.df <- stm.csv %>% 
  group_by(Treatment,Concentration, Type) %>%
  summarise(mean.CI = mean(CI_index) + 1e-6,
            std.CI = std.error(CI_index),
            count = n()) %>%
  as.data.frame()
stm.df
stm.df[stm.df < 1e-6] = 5e-3
write.csv(stm.df, "Short-term-values.csv")

## Statistics
## Learning
dat = subset(stm.csv, Type == "Learning")
dat$name = paste(dat$Treatment, "-", dat$Concentration)
dat$name = factor(dat$name, 
                  levels = c("Control - E. coli", "Control - DMSO", "Control - EGCG",
                             "Itogon - Low", "Itogon - Mid", "Itogon - High",
                             "Kapangan - Low","Kapangan - Mid", "Kapangan - High",
                             "La Trinidad - Low", "La Trinidad - Mid", "La Trinidad - High",
                             "Sablan - Low", "Sablan - Mid", "Sablan - High"))
res = DunnettTest(x=dat$CI_index, g=dat$name)
write.csv(res$`Control - E. coli`, "Learning-short-term-stat-results.csv")
res


## Memory
dat = subset(stm.csv, Type == "Short-term Memory")
dat$name = paste(dat$Treatment, "-", dat$Concentration)
dat$name = factor(dat$name, 
                  levels = c("Control - E. coli", "Control - DMSO", "Control - EGCG",
                             "Itogon - Low", "Itogon - Mid", "Itogon - High",
                             "Kapangan - Low","Kapangan - Mid", "Kapangan - High",
                             "La Trinidad - Low", "La Trinidad - Mid", "La Trinidad - High",
                             "Sablan - Low", "Sablan - Mid", "Sablan - High"))
res = DunnettTest(x=dat$CI_index, g=dat$name)
write.csv(res$`Control - E. coli`, "Memory-short-term-stat-results.csv")
res



## Visualize boxplots
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

## Controls - No Toxin
dat = stm.df[stm.df$Treatment == "Control", ]
x1 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = 'Chemotaxis Index (C.I)', title="Control") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" )) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "0.1% <br>DMSO", "200uM <br>EGCG")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.y = element_text(size=15))
x1


## Field Wine - No Toxin
dat = stm.df[stm.df$Treatment == "Itogon", ]
x2 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Itogon") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x2

## Field Wine - No Toxin
dat = stm.df[stm.df$Treatment == "Kapangan", ]
x3 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Kapangan") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x3

## Field Wine - No Toxin
dat = stm.df[stm.df$Treatment == "La Trinidad", ]
# dat[dat < 1e-6] = 5e-3
x4 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="La Trinidad") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x4

## Field Wine - No Toxin
dat = stm.df[stm.df$Treatment == "Sablan", ]
# dat[dat < 1e-6] = 5e-3
x5 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Sablan") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  # guides(fill = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x5

# x5 = x5 +
#   annotate("text", x = 0.78, y = 0.7, label = '*', size = 5, color = "#960101") 
# x5

## Field Wine - No Toxin
dat = stm.df[stm.df$Treatment == "Sablan", ]
# dat[dat < 1e-6] = 5e-3
x6 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Sablan", fill=NULL) +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  # guides(fill = guide_legend(title = " ")) 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x6

legend <- get_legend(x6)


# short.term.x = ggarrange(plotlist = list(x1, x2, x3, x4, x5), ncol = 5, nrow = 1)
short.term.x = plot_grid(x1, x2, x3, x4, x5, align = "h", ncol = 5, 
                         rel_widths = c(0.24,0.19,0.19,0.19,0.19))

p <- plot_grid(short.term.x, legend, ncol = 2, rel_widths = c(1.0, .2))
plot(p)

png(filename = file.path(getwd(),"/Short-term-all.png"), width = 11, height = 3.5, units = "in", res = 600)

plot(p)
dev.off()
print(paste("Saved at ", getwd()))


##########################################################################################################
########################################### LONG TERM MEMORY ############################################
##########################################################################################################

## Load survival csv file
learn.csv = read.csv("Long-term Learning.csv", header = T)
memory.csv = read.csv("Long-term Memory.csv", header = T)

## Learning
learn.csv$CI_index = (learn.csv$N_butanone - learn.csv$N_EtOH) / ((learn.csv$N_total - learn.csv$N_Origin) + 1e-6)
learn.csv$Type = "Learning"

## Memory
memory.csv$CI_index = (1e-6)
memory.csv$Type = "Long-term Memory"

## Summarize
ltm.csv = rbind(learn.csv,memory.csv)
ltm.csv$Treatment = factor(ltm.csv$Treatment, 
                           levels = c("Control", "Control", "Control","Itogon", "Kapangan", "La Trinidad", "Sablan"),
                           labels = c("Control", "Control", "Control","Itogon", "Kapangan", "La Trinidad", "Sablan"))
ltm.csv$Concentration = factor(ltm.csv$Concentration, 
                               levels = c("E. coli", "DMSO", "EGCG","Low", "Mid", "High"),
                               labels = c("E. coli", "DMSO", "EGCG","Low", "Mid", "High"))
ltm.csv$Type = factor(ltm.csv$Type,
                      levels=c("Learning", "Long-term Memory"),
                      labels=c("Learning", "Long-term Memory"))


ltm.df <- ltm.csv %>% 
  group_by(Treatment,Concentration, Type) %>%
  summarise(mean.CI = mean(CI_index) + 1e-6,
            std.CI = std.error(CI_index),
            count = n()) %>%
  as.data.frame()
ltm.df
ltm.df[ltm.df < 1e-6] = 5e-3
write.csv(ltm.df, "Long-term-values.csv")


##############################################################################
########################## CREATE PLOT #######################################
##############################################################################

ltm.df.2 = subset(ltm.df, Type != "Learning")
three.df = rbind(stm.df, ltm.df.2)
three.df

three.df$Treatment = factor(three.df$Treatment, 
                           levels = c("Control", "Control", "Control","Itogon", "Kapangan", "La Trinidad", "Sablan"),
                           labels = c("Control", "Control", "Control","Itogon", "Kapangan", "La Trinidad", "Sablan"))
three.df$Concentration = factor(three.df$Concentration, 
                               levels = c("E. coli", "DMSO", "EGCG","Low", "Mid", "High"),
                               labels = c("E. coli", "DMSO", "EGCG","Low", "Mid", "High"))
three.df$Type = factor(three.df$Type,
                      levels=c("Learning", "Short-term Memory","Long-term Memory"),
                      labels=c("Learning", "Short-term Memory","Long-term Memory"))


## Statistics
## Learning
# dat = subset(ltm.csv, Type == "Learning")
# dat$name = paste(dat$Treatment, "-", dat$Concentration)
# dat$name = factor(dat$name, 
#                   levels = c("Control - E. coli", "Control - DMSO", "Control - EGCG",
#                              "Itogon - Low", "Itogon - Mid", "Itogon - High",
#                              "Kapangan - Low","Kapangan - Mid", "Kapangan - High",
#                              "La Trinidad - Low", "La Trinidad - Mid", "La Trinidad - High",
#                              "Sablan - Low", "Sablan - Mid", "Sablan - High"))
# res = DunnettTest(x=dat$CI_index, g=dat$name)
# write.csv(res$`Control - E. coli`, "Learning-long-term-stat-results.csv")
# res


# ## Memory
# dat = subset(stm.csv, Type == "Long-term Memory")
# dat$name = paste(dat$Treatment, "-", dat$Concentration)
# dat$name = factor(dat$name, 
#                   levels = c("Control - E. coli", "Control - DMSO", "Control - EGCG",
#                              "Itogon - Low", "Itogon - Mid", "Itogon - High",
#                              "Kapangan - Low","Kapangan - Mid", "Kapangan - High",
#                              "La Trinidad - Low", "La Trinidad - Mid", "La Trinidad - High",
#                              "Sablan - Low", "Sablan - Mid", "Sablan - High"))
# res = DunnettTest(x=dat$CI_index, g=dat$name)
# write.csv(res$`Control - E. coli`, "Memory-long-term-stat-results.csv")
# res



## Visualize boxplots
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

## Controls - No Toxin
dat = three.df[three.df$Treatment == "Control", ]
x1 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = 'Chemotaxis Index (C.I)', title="Control") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" , "#f25c78")) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>", "0.1% <br>DMSO", "200uM <br>EGCG")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.y = element_text(size=15))
x1




## Field Wine - No Toxin
dat = three.df[three.df$Treatment == "Itogon", ]
x2 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Itogon") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" , "#f25c78")) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x2

# x2 = x2 +
#   annotate("text", x = 0.78, y = 0.85, label = '*', size = 5, color = "#960101")
# x2

## Field Wine - No Toxin
dat = three.df[three.df$Treatment == "Kapangan", ]
x3 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Kapangan") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" , "#f25c78" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x3

# x3 = x3 +
#   annotate("text", x = 1.78, y = 0.90, label = '*', size = 5, color = "#960101")
# x3

## Field Wine - No Toxin
dat = three.df[three.df$Treatment == "La Trinidad", ]
# dat[dat < 1e-6] = 5e-3
x4 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="La Trinidad") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" , "#f25c78" )) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  guides(fill = guide_legend(title = " ")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x4

## Field Wine - No Toxin
dat = three.df[three.df$Treatment == "Sablan", ]
# dat[dat < 1e-6] = 5e-3
x5 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Sablan") +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" , "#f25c78")) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  # guides(fill = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x5

# x5 = x5 +
#   annotate("text", x = 0.78, y = 0.7, label = '*', size = 5, color = "#960101") 
# x5

## Field Wine - No Toxin
dat = three.df[three.df$Treatment == "Sablan", ]
# dat[dat < 1e-6] = 5e-3
x6 =ggplot(dat, aes(x = Concentration,y = mean.CI, fill = Type)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  theme_classic() +labs(x = '', y = ' ', title="Sablan", fill=NULL) +
  geom_errorbar(aes(ymin=mean.CI, ymax=mean.CI+std.CI), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("#ffa600","#b64cdf" , "#f25c78")) +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
  # guides(fill = guide_legend(title = " ")) 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "", "","" ),
                     limits = c(0,1)) +
  theme(axis.text.x = element_markdown(size=10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
x6

legend <- get_legend(x6)


# short.term.x = ggarrange(plotlist = list(x1, x2, x3, x4, x5), ncol = 5, nrow = 1)
short.term.x = plot_grid(x1, x2, x3, x4, x5, align = "h", ncol = 5, 
                         rel_widths = c(0.24,0.19,0.19,0.19,0.19))

p <- plot_grid(short.term.x, legend, ncol = 2, rel_widths = c(1.0, .2))
plot(p)

png(filename = file.path(getwd(),"/All.png"), width = 13, height = 3.5, units = "in", res = 600)

plot(p)
dev.off()
print(paste("Saved at ", getwd()))








