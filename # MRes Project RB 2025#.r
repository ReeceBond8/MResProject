# MRes Project RB 2025#
#### Library ####
install.packages("binom")
library("ggplot2")
library("tidyr")
library("dplyr")
library("patchwork")
library("tidyverse")
library("binom")

#### Data ####
d1 = read.csv("Page11R.csv") #Clinically proven variants from Week 4
d1.1 = d1 %>% filter(Complex == 1) #Complex I variants only
d1.3 = d1 %>% filter(Complex == 3) #Complex III variants only
d1.4 = d1 %>% filter(Complex == 4) #Complex IV variants only
d1.5 = d1 %>% filter(Complex == 5) #Complex V variants only
d2 = read.csv("MRCUKCategorical0.27.csv") #% of variants in simple table for UK ME Cases/Controls
d2.1 = d2 %>% #longer table which groups columns for easier structure
  pivot_longer(cols = c(Cases, Controls),
               names_to = "Group",
               values_to = "Count")
d2.2 = read.csv("SAPBACategorical0.27.csv")
d2.3 = d2.2 %>% #same as 2.1 but for SAPBA cohort
  pivot_longer(cols = c(Cases, Controls),
               names_to = "Group",
               values_to = "Count")
d3 = read.csv("SAPBACategoricalCountsForChiSq.csv") #Counts from 0.27 ME SAPBA cohort
d3.1 = d3 %>% select(Cases, Controls,)
d3.2 = read.csv("UKMRCMEcountsforChiSquared.csv") #Counts from 0.27 ME MRC UK cohort
d3.3 = d3.2 %>% select(Cases, Controls,)
d4 = read.csv("ComplexPredictorData.csv")
d4.1 = pivot_longer(d4, cols = -complex, names_to = "Tool", values_to = "percentage")
#### MLC vs APOGEE 2 vs MutPred ####

#APOGEE 2 Pathogenic Predicting By Complex

p3 = ggplot(d1.1, aes(x = 1, y = APOGEE2.Score)) +
  #sets the original plot backing with blank x variable as only analysing the APOGEE score
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  #jitter spreads the data points out so they are not in a straight line
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  #sets the boxplot and alpha dictates the transparency of the box
  ggtitle("Complex I APOGEE2 Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  #removes the axis numbers and title as unneeded in this case
  theme(axis.title.x = element_blank())

p4 = ggplot(d1.3, aes(x = 1, y = APOGEE2.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  ggtitle("Complex III APOGEE2 Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p5 = ggplot(d1.4, aes(x = 1, y = APOGEE2.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  ggtitle("Complex IV APOGEE2 Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p6 = ggplot(d1.5, aes(x = 1, y = APOGEE2.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  ggtitle("Complex V APOGEE2 Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p3 + p4 + p5 + p6

#MLC by Complex#

p3.1 = ggplot(d1.1, aes(x = 1, y = MLC.Score)) +
  #sets the original plot backing with blank x variable as only analysing the APOGEE score
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  #jitter spreads the data points out so they are not in a straight line
  geom_boxplot(alpha = 0.5, fill = "salmon", colour="black") +
  #sets the boxplot and alpha dictates the transparency of the box
  ggtitle("Complex I MLC Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  #removes the axis numbers and title as unneeded in this case
  theme(axis.title.x = element_blank())

p4.1 = ggplot(d1.3, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  geom_boxplot(alpha = 0.5, fill = "salmon", colour="black") +
  ggtitle("Complex III MLC Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p5.1 = ggplot(d1.4, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  geom_boxplot(alpha = 0.5, fill = "salmon", colour="black") +
  ggtitle("Complex IV MLC Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p6.1 = ggplot(d1.5, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  geom_boxplot(alpha = 0.5, fill = "salmon", colour="black") +
  ggtitle("Complex V MLC Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p3.1 + p4.1 + p5.1 + p6.1

#MutPred by Complex#

p3.2 = ggplot(d1.1, aes(x = 1, y = MutPred.Score)) +
  #sets the original plot backing with blank x variable as only analysing the APOGEE score
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="green") +
  #jitter spreads the data points out so they are not in a straight line
  geom_boxplot(alpha = 0.5, fill = "chartreuse1", colour="black") +
  #sets the boxplot and alpha dictates the transparency of the box
  ggtitle("Complex I MutPred Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  #removes the axis numbers and title as unneeded in this case
  theme(axis.title.x = element_blank())

p4.2 = ggplot(d1.3, aes(x = 1, y = MutPred.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="green") +
  geom_boxplot(alpha = 0.5, fill = "chartreuse1", colour="black") +
  ggtitle("Complex III MutPred Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p5.2 = ggplot(d1.4, aes(x = 1, y = MutPred.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="green") +
  geom_boxplot(alpha = 0.5, fill = "chartreuse1", colour="black") +
  ggtitle("Complex IV MutPred Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p6.2 = ggplot(d1.5, aes(x = 1, y = MutPred.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="green") +
  geom_boxplot(alpha = 0.5, fill = "chartreuse1", colour="black") +
  ggtitle("Complex V MutPred Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

p3.2 + p4.2 + p5.2 + p6.2

#Proven pathogenic variants: MLC vs APOGEE 2 vs MutPred - ALL COMPLEXES#

p3.3 = ggplot(d1, aes(x = 1, y = APOGEE2.Score)) +
  #sets the original plot backing with blank x variable as only analysing the APOGEE score
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  #jitter spreads the data points out so they are not in a straight line
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0,1)) + 
  #sets the boxplot and alpha dictates the transparency of the box
  ggtitle("APOGEE2 Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  #removes the axis numbers and title as unneeded in this case
  theme(axis.title.x = element_blank())
p4.3 = ggplot(d1, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  geom_boxplot(width=0.3, alpha = 0.5, fill = "salmon", colour="black") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0,1)) + 
  ggtitle("MLC Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) 
#adds the plots either side of each other 
p5.3 = ggplot(d1, aes(x = 1, y = MutPred.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour = "green") +
  geom_boxplot(alpha = 0.5, fill = "chartreuse1", colour = "black") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0,1)) + 
  ggtitle("MutPred Scores") +
  ylab("Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())
p3.3 + p4.3 + p5.3

#Predictor by Complex Graph# 

p4= ggplot(d4.1, aes(x = factor(complex), y = percentage, 
                          shape = Tool, fill = Tool)) +  # Using fill aesthetic for colors
  geom_point(size = 7.5, position = position_dodge(width = 0.5),
             color = "black",  # Border color
             stroke = 0.5) +   # Border thickness
  
  # Y-axis settings
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10),
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  # X-axis with Roman numerals
  scale_x_discrete(
    labels = c("1" = "I", "3" = "III", "4" = "IV", "5" = "V")
  ) +
  
  # Shapes that can be filled (21-25 are fillable shapes)
  scale_shape_manual(values = c(21, 24, 22)) +  # Circle, Triangle, Square
  
  # Fill colors for tools
  scale_fill_manual(
    values = c("APOGEE2" = "lightblue", 
               "MLC" = "salmon", 
               "MutPred" = "chartreuse1"),
    name = "Tool"  # Legend title
  ) +
  
  # Labels and title
  labs(
    x = "Complex", 
    y = "Variants Classified as Pathogenic(%)",
    title = "Proven Pathogenic Variant Prediction Accuracy per Tool by Complex"
  ) +
  
  # Theme settings
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.key = element_rect(fill = "white")  # White background for legend
  )
pdf("Predictor.pdf", width = 8, height = 6)
p4
dev.off()

#### ME in SABPA and MRC UK cohorts ####

#MRC UK Bar Chart for % of deleterious variants harboured at APOGEE 2 threshold 0.27#

# Calculating CI (Wilsons CIs) first

figMRCUKME = ggplot(d2.1, aes(x = ï..Variants, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", linewidth = 0.5) +
  labs(title = "Number of deleterious variants (%) in UK ME cases and controls",
       x = "Number of deleterious variants",
       y = "% of Individuals",
       fill = "Group") +
  scale_fill_manual(values = c("Cases" = "steelblue", "Controls" = "salmon")) +
  theme_bw()
pdf("figMRCUKME", width =6, height = 6)
figMRCUKME
dev.off()

#SAPBA Bar Chart for % of deleterious variants harboured at APOGEE 2 threshold 0.27

figSAPBAME = ggplot(d2.3, aes(x = ï..Variants, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", linewidth = 0.5) +
  labs(title = "Number of deleterious variants (%) in SAPBA ME cases and controls",
       x = "Number of deleterious variants",
       y = "% of Individuals",
       fill = "Group") +
  scale_fill_manual(values = c("Cases" = "steelblue", "Controls" = "salmon")) +
  theme_bw()
pdf("figSAPBAME", width = 6, height = 6)
figSAPBAME
dev.off()

# Chi Squared testing of categorical data #

#SAPBA cohort 0.27 APOGEE 2 

head(d3)
colnames(d3)[colnames(d3) ==  "ï..Variants"] = "No of Variants"
colnames(d3)
contigency_table = table(d3.1$Cases, d3.1$Controls)
print(contigency_table)
r1 = chisq.test(d3.1)
print(r1)

# UK MRC cohort 0.27 APOGEE 2 #

head(d3.2)
colnames(d3.2)[colnames(d3.2) ==  "ï..Variants"] = "No of Variants"
colnames(d3.2)
contigency_table = table(d3.3$Cases, d3.3$Controls)
print(contigency_table)
r2 = chisq.test(d3.3)
print(r2)
