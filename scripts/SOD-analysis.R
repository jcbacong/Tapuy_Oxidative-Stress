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
df.csv= read.csv("Final SOD data.csv", header = T)
# df.csv.summary = subset(df.csv.summary, Type!="Lab Wine")
# df.csv.summary$Sample[1] = "Ascorbic Acid"

df.csv$Sample = factor(df.csv$Sample,
                           levels = c("OP50", "DMSO", "EGCG", "Itogon", "Kapangan", "La Trinidad", "Sablan"),
                           labels = c("OP50 <br><i>E. coli</i>", "0.1% <br>DMSO", "200uM <br>EGCG",
                                      "Itogon", "Kapangan", "La Trinidad", "Sablan"))

res = DunnettTest(x=df.csv$SOD_activity, g=df.csv$Sample)
write.csv(res$`OP50 <br><i>E. coli</i>`, "SOD-stat-results.csv")
res

df.summary = df.csv %>% 
  group_by(Sample) %>%
  summarise(Mean = mean(SOD_activity),
            Std = sd(SOD_activity),
            count = n()) %>%
  as.data.frame()
write.csv(df.summary, "SOD-values.csv")


color_scheme = c("#84d94f","#84d94f","#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

## Itogon
# ascorbic.acid = subset(df.csv.summary, Treatment == "Ascorbic acid")
q1 = ( ggplot(df.summary, aes(x=Sample, y=Mean, fill=Sample)) +
         theme_classic()  + labs(x = '', y = 'SOD Activity (U/mL)') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Mean-Std, ymax=Mean+Std), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ Sample', scales="free", switch = "x", space='free') +
         # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=13, face="bold"),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) + 

         scale_y_continuous(limits = c(0,2))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q1

q1 = q1 +
  annotate("text", x = 3, y = 0.97, label = "*", size = 7, color = "#960101") +
  annotate("text", x = 4, y = 1.63, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 5, y = 1.85, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 6, y = 1.70, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 7, y = 1.75, label = "***", size = 7, color = "#960101")
q1

## Statistics
png(filename = file.path(getwd(),"/SOD-Parallel.png"), width = 8, height = 3.5, units = "in", res = 600)

plot(q1)
dev.off()
print(paste("Saved at ", getwd()))

