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
pos = read.csv("meron H2O2 Results.csv", header = T)
pos$type = "Treated"
neg = read.csv("wala H2O2 Results.csv", header = T)
neg$type = "Untreated"

ndata = rbind(neg, pos)


## Pattern the table from the previous KM plots.
## Drop the missing and consider only those alive and dead
## Arrange the tables as necessary
# data = pos
# ndata = subset(data, select = -c(X))
# # ndata$control = as.factor(ndata$control)
ndata$control <- ordered(ndata$control,levels=c("OP50", "DMSO", "EGCG", 
                                                "It L", "It M", "It H", 
                                                "Lt L", "Lt M", "Lt H",
                                                "Kp L", "Kp M", "Kp H",
                                                "Sb L", "Sb M", "Sb H"))
ndata = ndata %>% arrange(factor(control, levels=c("OP50", "DMSO", "EGCG", 
                                                            "It L", "It M", "It H", 
                                                            "Lt L", "Lt M", "Lt H",
                                                            "Kp L", "Kp M", "Kp H",
                                                            "Sb L", "Sb M", "Sb H")))
ndata$type <- ordered(ndata$type,levels=c("Untreated", "Treated"),
                      labels = c("- H2O2", "+ H2O2"))

ndata.summ = ndata %>% group_by(control, type) %>% summarise(mean = mean(Mean), sem = std.error(Mean), count = n())
ndata.summ = ndata.summ %>% arrange(factor(control, levels=c("OP50", "DMSO", "EGCG", 
                                                             "It L", "It M", "It H", 
                                                             "Lt L", "Lt M", "Lt H",
                                                             "Kp L", "Kp M", "Kp H",
                                                             "Sb L", "Sb M", "Sb H")))
ndata.summ
write.csv(ndata.summ, "Oxidative-Stress-values.csv")

n.plot = ggplot(ndata.summ, aes(control, mean, fill=type)) +
  theme_classic() + labs(x = '', y = 'Meean GFP Intensity (A.U)', title="") +
  geom_bar(position=position_dodge(), stat="identity", colour="black") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sem), width=.2,position=position_dodge(.9)) +
  scale_x_discrete(labels = c("OP50 <br><i>E. coli</i>",
                              "0.1% <br>DMSO",
                              "200 uM <br>EGCG",
                              "It L", "It M", "It H", 
                              "Lt L", "Lt M", "Lt H",
                              "Kp L", "Kp M", "Kp H",
                              "Sb L", "Sb M", "Sb H")) +
  # facet_grid('. ~ control', scales="free", switch = "x") +
  scale_fill_manual(values=c("#c7f5a6", "#068f44")) +
  guides(fill = guide_legend(title = "Oxidative Stress")) +
  theme(axis.text.x = element_markdown(size=10))
n.plot

# png(filename = file.path(getwd(),"/temp.png"), width = 8, height = 3, units = "in", res = 600)
# 
# plot(n.plot)
# dev.off()
# print(paste("Saved at ", getwd()))

### STAT TEST
## OP50 alternative hypothesis: -H202 < +H2O2
summary_df <- data.frame(control = character(),
                         t.stat = numeric(),
                         p.value = numeric(),
                         stringsAsFactors = FALSE)

## If not negative control: we want -H2O2 > +H2O2 (because protective na siya from the stress)
for (ctrl in unique(ndata.summ$control)) {
  dat.test = subset(ndata.summ, control==ctrl)
  if (ctrl == "OP50") {
    res = tsum.test(mean.x=dat.test[dat.test$type=="- H2O2", "mean"]$mean,
                    s.x=dat.test[dat.test$type=="- H2O2", "sem"]$sem,
                    n.x=dat.test[dat.test$type=="- H2O2", "count"]$count,
                    
                    mean.y=dat.test[dat.test$type=="+ H2O2", "mean"]$mean,
                    s.y=dat.test[dat.test$type=="+ H2O2", "sem"]$sem,
                    n.y=dat.test[dat.test$type=="- H2O2", "count"]$count,
                    alternative = "less")
  }
  if (ctrl != "OP50"){
    
    res = tsum.test(mean.x=dat.test[dat.test$type=="- H2O2", "mean"]$mean,
                    s.x=dat.test[dat.test$type=="- H2O2", "sem"]$sem,
                    n.x=dat.test[dat.test$type=="- H2O2", "count"]$count,
                    
                    mean.y=dat.test[dat.test$type=="+ H2O2", "mean"]$mean,
                    s.y=dat.test[dat.test$type=="+ H2O2", "sem"]$sem,
                    n.y=dat.test[dat.test$type=="+ H2O2", "count"]$count,
                    alternative = "greater")
  }
  summary_row <- data.frame(control = ctrl,
                            t.stat = res$statistic,
                            p.value = res$p.value,
                            stringsAsFactors = FALSE)
  summary_df <- rbind(summary_df, summary_row)
}
print(summary_df)
write.csv(summary_df, file="stat_summary.csv")

## List of significant (***)
print(subset(summary_df, p.value < 0.05)$control)

## PLOT WITH STAT
op50 <- tibble(
  x = c(0.78, 0.78, 1.22, 1.22),
  y = c(65, 100, 100, 85)
) ## p = 1.321677e-74

dmso <- tibble(
  x = c(1.78, 1.78, 2.22, 2.22),
  y = c(155, 160, 160, 133)
) ## p = 0

itl <- tibble(
  x = c(3.78, 3.78, 4.22, 4.22),
  y = c(120, 125, 125, 110)
) ## p = 0

itm <- tibble(
  x = c(4.78, 4.78, 5.22, 5.22),
  y = c(193, 202, 202, 80)
) ## p = 0

ith <- tibble(
  x = c(5.78, 5.78, 6.22, 6.22),
  y = c(160, 170, 170, 143)
) ## p = 0

ltl <- tibble(
  x = c(6.78, 6.78, 7.22, 7.22),
  y = c(158, 168, 168, 135)
) ## p = 0

ltm <- tibble(
  x = c(7.78, 7.78, 8.22, 8.22),
  y = c(145, 150, 150, 139)
) ## p = 0

lth <- tibble(
  x = c(8.78, 8.78, 9.22, 9.22),
  y = c(142, 150, 150, 139)
) ## p = 0

kpl <- tibble(
  x = c(9.78, 9.78, 10.22, 10.22),
  y = c(146, 155, 155, 140)
) ## p = 0

kpm <- tibble(
  x = c(10.78, 10.78, 11.22, 11.22),
  y = c(155, 167, 167, 150)
) ## p = 0

kph <- tibble(
  x = c(11.78, 11.78, 12.22, 12.22),
  y = c(145, 153, 153, 122)
) ## p = 0

sbl <- tibble(
  x = c(12.78, 12.78, 13.22, 13.22),
  y = c(190, 200, 200, 152)
) ## p = 0

n.plot.1 = n.plot +
  geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 1, y = 105, label = '***', size = 5, color = "#960101") +
  geom_line(data = dmso, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 2, y = 165, label = '***', size = 5, color = "#960101") +
  geom_line(data = itl, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 4, y = 127, label = '***', size = 5, color = "#960101") +
  geom_line(data = itm, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 5, y = 204, label = '***', size = 5, color = "#960101") +
  geom_line(data = ith, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 6, y = 173, label = '***', size = 5, color = "#960101") +
  geom_line(data = ltl, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 7, y = 171, label = '***', size = 5, color = "#960101") +
  geom_line(data = ltm, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 8, y = 153, label = '***', size = 5, color = "#960101") +
  geom_line(data = lth, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 9, y = 153, label = '***', size = 5, color = "#960101") +
  geom_line(data = kpl, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 10, y = 157, label = '***', size = 5, color = "#960101") +
  geom_line(data = kpm, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 11, y = 169, label = '***', size = 5, color = "#960101") +
  geom_line(data = kph, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 12, y = 155, label = '***', size = 5, color = "#960101") +
  geom_line(data = sbl, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 13, y = 202, label = '***', size = 5, color = "#960101")
  

n.plot.1


png(filename = file.path(getwd(),"/stat_ASAP.png"), width = 9, height = 3, units = "in", res = 600)

plot(n.plot.1)
dev.off()
print(paste("Saved at ", getwd()))
