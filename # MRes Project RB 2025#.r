# MRes Project RB 2025#
#### Library ####
install.packages("binom")
install.packages("RColorBrewer")
install.packages("ggpubr")
library("ggplot2")
library("tidyr")
library("dplyr")
library("patchwork")
library("tidyverse")
library("binom")
library("RColorBrewer")
library("Forecats")
library("pROC")

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
d5 = read.csv("tRNAPathogenicMutationsFD.csv")
d6 = read.csv("Distribition-tRNAs-In-56.csv")
d7 = read.csv("56MostFrequent-tRNA.csv")
d8 = read.csv("CommonVariants-Classified.csv") %>%
  mutate(Classification = factor(Classification, levels = c("Benign", "VUS", "Pathogenic")))
d9 = read.csv("PathogenicVariantsClassified.csv") %>%
  mutate(Classification = factor(Classification, levels = c("Benign", "VUS", "Pathogenic")))
d10 = read.csv("CommonVariants-Classified-0.457.csv") %>%
  mutate(Classification = factor(Classification, levels = c("Benign", "VUS", "Pathogenic")))
d11 = read.csv("PathogenicVariantsClassified-0.457.csv") %>%
  mutate(Classification = factor(Classification, levels = c("Benign", "VUS", "Pathogenic")))
d12 = read.csv("ROC2.csv")
d13 = read.csv("LogReg.csv")
d14 = read.csv("CFS_Samples_UK_with_mtVLM+FIS+FS36.csv")
FIS = d14$FIS
Score = d14$mtVLM..0.16..tRNA.score
SF = d14$SF.36
d15 = read.csv("tRNA+PE+Combined+MEUK.csv")
FIS2 = d15$FIS
SF2 = d15$SF36
Score2 = d15$Variant.Load.Score

#### Colours ####

point_colors <- brewer.pal(3, "Dark2")  
fill_colors  <- brewer.pal(3, "Pastel2")
palette_set3 = brewer.pal(3, "Set1")  
palette_pastel3 = brewer.pal(3, "Pastel1")

#### MLC vs APOGEE 2 vs MutPred for Protein Encoding  ####

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

#### tRNA APOGEE (0.716) vs MLC vs MitoTIP #####

#Basic Graphs Comparing Scores of Pathogenic Proven Variants# 

#See Week 12 Folder for use of RColorBrewerSet#

h1 = ggplot(d5, aes(x = 1, y = Unbiased.APOGEE.score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  ggtitle("APOGEE tRNA Pathogenicity Scores for Cfrm Pathogenic Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  #removes the axis numbers and title as unneeded in this case
  theme(axis.title.x = element_blank())

h2 = ggplot(d5, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  geom_boxplot(alpha = 0.5, fill = "salmon", colour="black") +
  ggtitle("MLC tRNA Pathogenicity Scores for Cfrm Pathogenic Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

h3 = ggplot(d5, aes(x = 1, y = MitoTIP.Raw.Scores)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour = "green") +
  geom_boxplot(alpha = 0.5, fill = "Chartreuse1", colour="black") +
  scale_fill_brewer(palette = "set2") +
  ggtitle("MitoTIP Raw tRNA Pathogenicity Scores for Cfrm Pathogenic Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())
h1 + h2 + h3

#Pathogenic Variants Against Pathogenic Thresholds#

h1.1 = sum(d5$APOGEE.Score >=0.716, na.rm = TRUE)
h2.1 = sum(d5$MLC.Score >=0.5, na.rm = TRUE)
h3.1 = sum(d5$MitoTIP.Raw.Scores > 12.66, na.rm = TRUE)

#Pathogenic Variants Against Benign Thresholds#

h1.2 = sum(d5$APOGEE.Score <= 0.396, na.rm = TRUE) # Use of the VUS-/LB/B threshold here from protein encoding
h2.2 = sum(d5$MLC.Score < 0.5, na.rm = TRUE) 
h3.2 = sum(d5$MitoTIP.Raw.Scores < 8.44, na.rm = TRUE)

#Illustration of the above classifications#

palette_set3 = brewer.pal(3, "Set1")  
palette_pastel3 = brewer.pal(3, "Pastel1")

ggplot(d9, aes(fill = Classification, y = ï..Count, x = Tool.Used)) +
  geom_col(linewidth = 0.5, colour = "black", width = 0.6, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name = NULL,
                    breaks = c("Pathogenic", "VUS", "Benign"), 
                    values = c(brewer.pal(3, "Set2"), "grey")) +
  scale_x_discrete(breaks = c("APOGEE","MLC","MitoTIP"),
                   labels = c("APOGEE","MLC","MitoTIP")) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60),
                     labels = c("0", "10", "20", "30", "40", "50", "60"),
                     limits = c(0,60)) +
  labs(x = "Pathogenicity Tool",
       y = "tRNA confirmed Pathogenic Variant Count",
       title = "Classification of 57 confirmed pathogenic variants by pathogenicity predictor tool") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"))

#Pathogenic Proven Variants using RColourBrewer for tRNA#

h1 = ggplot(d5, aes(x = 1, y = Unbiased.APOGEE.score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="blue") +
  geom_boxplot(alpha = 0.5, fill = "lightblue", colour="black") +
  ggtitle("APOGEE tRNA Pathogenicity Scores for Cfrm Pathogenic Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  #removes the axis numbers and title as unneeded in this case
  theme(axis.title.x = element_blank())

h2 = ggplot(d5, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour="red") +
  geom_boxplot(alpha = 0.5, fill = "salmon", colour="black") +
  ggtitle("MLC tRNA Pathogenicity Scores for Cfrm Pathogenic Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

h3 = ggplot(d5, aes(x = 1, y = MitoTIP.Raw.Scores)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour = "green") +
  geom_boxplot(alpha = 0.5, fill = "Chartreuse1", colour="black") +
  scale_fill_brewer(palette = "set2") +
  ggtitle("MitoTIP Raw tRNA Pathogenicity Scores for Cfrm Pathogenic Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())
h1 + h2 + h3

# Testing the benign variants - the distribution of tRNA in 56 commmon variants from GenBank #

tRNA_Benign_Palette = colorRampPalette(brewer.pal(9, "Set1"))(22)
ggplot(d6, aes(x = reorder(x = ï..tRNA, -Count), y = Count, fill = ï..tRNA)) +
  geom_bar(stat = "identity", width = 0.7, colour = "Black", linewidth = 0.5) +
  scale_fill_manual(values = tRNA_Benign_Palette) +
  labs(title = "tRNA loci distribution amongst 57 common variants",
       x = "tRNA loci",
       y = "Count within 57 most frequent variants (GenBank)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1.5),
    panel.grid.major.y = element_line(colour = "grey80"),
    legend.position = "none"
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)))

# APOGEE vs MLC vs MitoTIP in 'Benign common population variants # 

palette_set2 = brewer.pal(3, "Accent")  
palette_pastel2 = brewer.pal(3, "Accent")  #Opposing colours to that used in Pathogenic variants for easy comparison

h11 = ggplot(d7, aes(x = 1, y = Unbiased.APOGEE.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour = palette_set2[2]) + 
  geom_boxplot(alpha = 0.5, fill = palette_pastel2[2], colour = "black") +
  ggtitle("APOGEE tRNA Pathogenicity Scores for 'Common' population Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

h12 = ggplot(d7, aes(x = 1, y = MLC.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour = palette_set2[1]) + 
  geom_boxplot(alpha = 0.5, fill = palette_pastel2[1], colour = "black") +
  ggtitle("MLC tRNA Pathogenicity Scores for 'Common' population Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

h13 = ggplot(d7, aes(x = 1, y = MitoTIP.Score)) +
  geom_jitter(width = 0.2, height = 0, alpha = 1, colour = palette_set2[3]) + 
  geom_boxplot(alpha = 0.5, fill = palette_pastel2[3], colour = "black") +
  ggtitle("MitoTIP Raw tRNA Pathogenicity Scores for 'Common' population Variants") +
  ylab("Pathogenicity Scores") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

h11 + h12 + h13

# Common population variant counts against 'Pathogenic' Thresholds#
h11.1 = sum(d7$Unbiased.APOGEE.Score >=0.716, na.rm = TRUE)
h12.1 = sum(d7$MLC.Score >=0.5, na.rm = TRUE)
h13.1 = sum(d7$MitoTIP.Score > 12.66, na.rm = TRUE)

#Common population variant counts against 'Benign' Thresholds#

h11.2 = sum(d7$Unbiased.APOGEE.Score <= 0.396, na.rm = TRUE) # Use of the VUS-/LB/B threshold here from protein encoding
h12.2 = sum(d7$MLC.Score < 0.5, na.rm = TRUE) 
h13.2 = sum(d7$MitoTIP.Score < 8.44, na.rm = TRUE) #this is the likely benign threshold

#Illustrated by a Stacked Bar Chart#

ggplot(d8, aes(fill = Classification, y = ï..Count, x = Tool.Used)) +
  geom_col(linewidth = 0.5, colour = "black", width = 0.6, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(name = NULL,
                    breaks = c("Pathogenic", "VUS", "Benign"), 
                    values = c(brewer.pal(3, "Set2"), "grey")) +
  scale_x_discrete(breaks = c("APOGEE","MLC","MitoTIP"),
                   labels = c("APOGEE","MLC","MitoTIP")) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60),
                     labels = c("0", "10", "20", "30", "40", "50", "60"),
                     limits = c(0,60)) +
  labs(x = "Pathogenicity Tool",
       y = "tRNA Common Variant Count",
       title = "Classification of 57 Common mt-tRNA population variants by pathogenicity predictor tool") +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    
    # Plots with all tools common and pathogenic variants # 
    
    #APOGEE Scores for common and pathogenic
    
    h14 = ggplot() +
      geom_boxplot(data = d7, aes(x = 1, y = Unbiased.APOGEE.Score), alpha = 0.5, fill = palette_pastel2[2], colour = "black") +
      geom_jitter(data = d7, aes(x = 1, y = Unbiased.APOGEE.Score), width = 0.2, height = 0, alpha = 1, colour = palette_set2[2]) + 
      geom_boxplot(data = d5, aes(x = 1, y = Unbiased.APOGEE.Score), alpha = 0.5, fill = fill_colors[1], colour = "black") +
      geom_jitter(data = d5, aes(x = 1, y = Unbiased.APOGEE.Score), width = 0.2, height = 0, alpha = 1, colour = point_colors[1]) +  
      ggtitle("APOGEE tRNA Pathogenicity Scores for Common and Pathogenic  Variants") +
      ylab("Pathogenicity Scores") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    #MLC Scores for common and pathogenic
    
    h15 = ggplot() +
      geom_jitter(data = d7, aes(x = 1, y = MLC.Score), width = 0.2, height = 0, alpha = 1, colour = palette_set2[1]) + 
      geom_boxplot(data = d7, aes(x = 1, y = MLC.Score), alpha = 0.5, fill = palette_pastel2[1], colour = "black") +
      geom_jitter(data = d5, aes(x = 1, y = MLC.Score), width = 0.2, height = 0, alpha = 1, colour = point_colors[2]) +  # Dark orange
      geom_boxplot(data = d5, aes(x = 1, y = MLC.Score),alpha = 0.5, fill = fill_colors[2], colour = "black") +
      ggtitle("MLC tRNA Pathogenicity Scores for Common and Pathogenic Variants") +
      ylab("Pathogenicity Scores") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    #MitoTIP scores for common and pathogenic
    
    h16 = ggplot(data = d7, aes(x = 1, y = MitoTIP.Score)) +
      geom_jitter(data = d7, aes(x = 1, y = MitoTIP.Score), width = 0.2, height = 0, alpha = 1, colour = palette_set2[3]) + 
      geom_boxplot(data = d7, aes(x = 1, y = MitoTIP.Score), alpha = 0.5, fill = palette_pastel2[3], colour = "black") +
      geom_jitter(data = d5, aes(x = 1, y = MitoTIP.Score), width = 0.2, height = 0, alpha = 1, colour = point_colors[3]) +  # Dark lavender
      geom_boxplot(data = d5, aes(x = 1, y = MitoTIP.Score), alpha = 0.5, fill = fill_colors[3], colour = "black") +        # Light lavender
      ggtitle("MitoTIP Raw tRNA Pathogenicity Scores for Common and Pathogenic Variants") +
      ylab("Pathogenicity Scores") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    h14 + h15 + h16
#### tRNA APOGEE (0.457) vs MLC vs MitoTIP #####
    
    #Calculating the optimum threshold for tRNA in APOGEE via ROC and Youden's J stat#
    
    
stopifnot(all(d12$Status %in% c(0,1)))
    # Calculate ROC curve
roc_obj = roc(response = d12$Status, 
                  predictor = d12$ï..APOGEE.Score,
                  direction = "<") # Higher scores = more pathogenic
    # Get optimal threshold (Youden's J statistic)
optimal = coords(roc_obj, x = "best", 
                     ret = c("threshold", "specificity", "sensitivity"))
roc_plot = ggroc(roc_obj, legacy.axes = TRUE, size = 1) + 
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   color = "grey50", linetype = "dashed") +
      geom_point(aes(x = 1 - optimal$specificity, 
                     y = optimal$sensitivity),
                 color = "blue", size = 5) +
      annotate("text", x = 0.7, y = 0.2, 
               label = paste("Optimal threshold:", round(optimal$threshold, 3)),
               color = "blue", size = 5, face = "bold") +
      labs(title = "ROC Curve for APOGEE Score Pathogenicity Prediction",
           x = "False Positivity Rate (1 - Specificity)",
           y = "True Positivity Rate (Sensitivity)") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold")
      ) +
      scale_x_reverse(limits = c(1,0)) # Invert x-axis' so beginning at 0 (1-spec)
    print(roc_plot)
cat("\nOptimal Threshold Results:\n")
cat("-------------------------\n")
cat("Threshold (APOGEE score):", optimal$threshold, "\n")
cat("Sensitivity (TP rate):   ", optimal$sensitivity, "\n")
cat("Specificity (TN rate):   ", optimal$specificity, "\n")
cat("FPR (1 - Specificity):   ", 1 - optimal$specificity, "\n") #This outputs False and True rates and optimum thresholds 
    
    #Log Regression Analysis of MitoTIP and APOGEE values#
    
    
logistic = glm(Status ~ ï..APOGEE.Score, data = d13, family="binomial")
summary(logistic)
logreggraph = data.frame(pathprobability = logistic$fitted.values, path=d13$Status)
logreggraph = logreggraph[
      order(logreggraph$pathprobability, decreasing=FALSE),]
logreggraph$rank = 1:nrow(logreggraph) 
reg1 = ggplot(logreggraph, aes(x = rank, y = pathprobability)) +
      geom_point(aes(colour=path), alpha=1, shape=4, stroke=2) +
      xlab("mt-tRNA variant rank by APOGEE score") +
      ylab("Predicted Probability of mt-tRNA variant classified as Pathogenic by APOGEE score") +
      theme_bw()
    
    
logistic2 = glm(Status ~ MitoTIP.Score, data = d13, family="binomial")
summary(logistic2)
logreggraph2 = data.frame(pathprobability = logistic2$fitted.values, path=d13$Status)
logreggraph2 = logreggraph2[
      order(logreggraph2$pathprobability, decreasing=FALSE),]
logreggraph2$rank = 1:nrow(logreggraph2) 
reg2 = ggplot(logreggraph2, aes(x = rank, y = pathprobability)) +
      geom_point(aes(colour=path), alpha=1, shape=4, stroke=2) +
      xlab("mt-tRNA variant rank by MitoTIP score") +
      ylab("Predicted Probability of mt-tRNA variant classified as Pathogenic by MitoTIP score") +
      theme_bw()
    
reg1 + reg2
    
    #Pathogenic Variants Against 0.457 Pathogenic Thresholds#
    
    h21.1 = sum(d5$APOGEE.Score >=0.457, na.rm = TRUE)
    h22.1 = sum(d5$MLC.Score >=0.5, na.rm = TRUE)
    h23.1 = sum(d5$MitoTIP.Raw.Scores > 12.66, na.rm = TRUE)
    
    #Pathogenic Variants Against 0.457 Benign Thresholds#
    
    h21.2 = sum(d5$APOGEE.Score < 0.457, na.rm = TRUE) # Use of the VUS-/LB/B threshold here from protein encoding
    h22.2 = sum(d5$MLC.Score < 0.5, na.rm = TRUE) 
    h23.2 = sum(d5$MitoTIP.Raw.Scores < 8.44, na.rm = TRUE)
    
    #Illustration of the above classifications#
    
    ggplot(d10, aes(fill = Classification, y = ï..Count, x = Tool.Used)) +
      geom_col(linewidth = 0.5, colour = "black", width = 0.6, position = position_stack(reverse = TRUE)) +
      scale_fill_manual(name = NULL,
                        breaks = c("Pathogenic", "VUS", "Benign"), 
                        values = c(brewer.pal(3, "Set2"), "grey")) +
      scale_x_discrete(breaks = c("APOGEE","MLC","MitoTIP"),
                       labels = c("APOGEE(0.457)","MLC","MitoTIP")) +
      scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60),
                         labels = c("0", "10", "20", "30", "40", "50", "60"),
                         limits = c(0,60)) +
      labs(x = "Pathogenicity Tool",
           y = "tRNA confirmed Pathogenic Variant Count",
           title = "Classification of 57 confirmed pathogenic variants by pathogenicity predictor tool") +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Common population variant counts against 'Pathogenic' Thresholds#
    h31.1 = sum(d7$Unbiased.APOGEE.Score >=0.457, na.rm = TRUE)
    h32.1 = sum(d7$MLC.Score >=0.5, na.rm = TRUE)
    h33.1 = sum(d7$MitoTIP.Score > 12.66, na.rm = TRUE)
    
    #Common population variant counts against 'Benign' Thresholds#
    
    h31.2 = sum(d7$Unbiased.APOGEE.Score < 0.457, na.rm = TRUE) # Use of the VUS-/LB/B threshold here from protein encoding
    h32.2 = sum(d7$MLC.Score < 0.5, na.rm = TRUE) 
    h33.2 = sum(d7$MitoTIP.Score < 8.44, na.rm = TRUE) #this is the likely benign threshold
    
    #Illustrated by a Stacked Bar Chart#
    
    z1 = ggplot(d11, aes(fill = Classification, y = ï..Count, x = Tool.Used)) +
      geom_col(linewidth = 0.5, colour = "black", width = 0.6, position = position_stack(reverse = TRUE)) +
      scale_fill_manual(name = NULL,
                        breaks = c("Pathogenic", "VUS", "Benign"), 
                        values = c(brewer.pal(3, "Set2"), "grey")) +
      scale_x_discrete(breaks = c("APOGEE","MLC","MitoTIP"),
                       labels = c("APOGEE (0.457)","MLC","MitoTIP")) +
      scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60),
                         labels = c("0", "10", "20", "30", "40", "50", "60"),
                         limits = c(0,60)) +
      labs(x = "Pathogenicity Tool",
           y = "tRNA Common Variant Count",
           title = "Classification of 57 cfrm pathogenic variants by pathogenicity predictor tool") +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))
    
    z2 = ggplot(d10, aes(fill = Classification, y = ï..Count, x = Tool.Used)) +
      geom_col(linewidth = 0.5, colour = "black", width = 0.6, position = position_stack(reverse = TRUE)) +
      scale_fill_manual(name = NULL,
                        breaks = c("Pathogenic", "VUS", "Benign"), 
                        values = c(brewer.pal(3, "Set2"), "grey")) +
      scale_x_discrete(breaks = c("APOGEE","MLC","MitoTIP"),
                       labels = c("APOGEE (0.457)","MLC","MitoTIP")) +
      scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60),
                         labels = c("0", "10", "20", "30", "40", "50", "60"),
                         limits = c(0,60)) +
      labs(x = "Pathogenicity Tool",
           y = "tRNA confirmed Pathogenic Variant Count",
           title = "Classification of 57 common population variants by pathogenicity predictor tool") +
      theme_bw() +
      theme(
        text = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))
    z1 + z2
    
    
    # Plots with all tools common and pathogenic variants # 
    
    #APOGEE Scores for common and pathogenic
    
    h14 = ggplot() +
      geom_boxplot(data = d7, aes(x = 1, y = Unbiased.APOGEE.Score), alpha = 0.5, fill = palette_pastel2[2], colour = "black") +
      geom_jitter(data = d7, aes(x = 1, y = Unbiased.APOGEE.Score), width = 0.2, height = 0, alpha = 1, colour = palette_set2[2]) + 
      geom_boxplot(data = d5, aes(x = 1, y = Unbiased.APOGEE.Score), alpha = 0.5, fill = fill_colors[1], colour = "black") +
      geom_jitter(data = d5, aes(x = 1, y = Unbiased.APOGEE.Score), width = 0.2, height = 0, alpha = 1, colour = point_colors[1]) +  
      ggtitle("APOGEE tRNA Pathogenicity Scores for Common and Pathogenic  Variants") +
      ylab("Pathogenicity Scores") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    #MLC Scores for common and pathogenic
    
    h15 = ggplot() +
      geom_jitter(data = d7, aes(x = 1, y = MLC.Score), width = 0.2, height = 0, alpha = 1, colour = palette_set2[1]) + 
      geom_boxplot(data = d7, aes(x = 1, y = MLC.Score), alpha = 0.5, fill = palette_pastel2[1], colour = "black") +
      geom_jitter(data = d5, aes(x = 1, y = MLC.Score), width = 0.2, height = 0, alpha = 1, colour = point_colors[2]) +  # Dark orange
      geom_boxplot(data = d5, aes(x = 1, y = MLC.Score),alpha = 0.5, fill = fill_colors[2], colour = "black") +
      ggtitle("MLC tRNA Pathogenicity Scores for Common and Pathogenic Variants") +
      ylab("Pathogenicity Scores") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    #MitoTIP scores for common and pathogenic
    
    h16 = ggplot(data = d7, aes(x = 1, y = MitoTIP.Score)) +
      geom_jitter(data = d7, aes(x = 1, y = MitoTIP.Score), width = 0.2, height = 0, alpha = 1, colour = palette_set2[3]) + 
      geom_boxplot(data = d7, aes(x = 1, y = MitoTIP.Score), alpha = 0.5, fill = palette_pastel2[3], colour = "black") +
      geom_jitter(data = d5, aes(x = 1, y = MitoTIP.Score), width = 0.2, height = 0, alpha = 1, colour = point_colors[3]) +  # Dark lavender
      geom_boxplot(data = d5, aes(x = 1, y = MitoTIP.Score), alpha = 0.5, fill = fill_colors[3], colour = "black") +        # Light lavender
      ggtitle("MitoTIP Raw tRNA Pathogenicity Scores for Common and Pathogenic Variants") +
      ylab("Pathogenicity Scores") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    h14 + h15 + h16

#### UK ME cohort assessment (severity vs tRNA mtVLM score) ####
    
# Testing Linearity of tRNA mtVLM score vs FIS # 
    
PearsonsCor = cor.test(FIS, Score, method = "pearson", use = "complete.obs") #Complete obs used to ignore NA values 
PearsonsCor #Prints Pearsons Cor to test for linearity (rem: 0 indicates no linearity)
    
# Plot showing r'ship between mtVLM tRNA score and FIS with the regression line
model1=lm(FIS ~ Score, data = d14)
r_squared1=summary(model1)$r.squared
r_squared_label1=paste0("R² = ", round(r_squared1, 4))

sh1 = ggplot(d14, aes(x = Score, y = FIS)) +
      geom_point(size = 3, alpha = 0.6, color = "steelblue") +  # Larger, semi-transparent points
      geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1.5) +  # Thicker regression line
      annotate("text", x = Inf, y = Inf, label = r_squared_label1,  # Add R²
               hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
      labs(
        title = "Effect of tRNA mtVLM Score on FIS",
        x = "mtVLM Score (tRNA only) ",
        y = "FIS"
      ) +
      theme_bw(base_size = 14) +  # Clean theme with larger text
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Bold, centered title
        axis.title = element_text(face = "bold"),  # Bold axis labels
        panel.grid.major = element_line(color = "grey90"),  # Lighter grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
      )

#Above repeated for the SF-36 Score

PearsonsCor2 = cor.test(SF, Score, method = "pearson", use = "complete.obs") 
PearsonsCor2 

model2=lm(SF ~ Score, data = d14)
r_squared2=summary(model2)$r.squared
r_squared_label2=paste0("R² = ", round(r_squared2, 4))

sh2 = ggplot(d14, aes(x = Score, y = SF)) +
  geom_point(size = 3, alpha = 0.6, color = "steelblue") +  
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1.5) +  
  annotate("text", x = Inf, y = Inf, label = r_squared_label2,  
           hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
  labs(
    title = "Effect of tRNA mtVLM Score on SF-36",
    x = "mtVLM Score (tRNA only) ",
    y = "SF-36"
  ) +
  theme_bw(base_size = 14) +  
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    axis.title = element_text(face = "bold"), 
    panel.grid.major = element_line(color = "grey90"),  
    panel.grid.minor = element_blank()
    )

sh1 + sh2
#### UK ME cohort PE + tRNA mtVLM scores vs FIS/SF36 scores
#### UK ME tRNA + PE scores vs FIS/SF36 scores ####

PearsonsCor3 = cor.test(FIS2, Score2, method = "pearson", use = "complete.obs") 
PearsonsCor3

# Plot showing r'ship between mtVLM score (tRNA + PE) and FIS with the regression line

model3=lm(FIS2 ~ Score2, data = d15)
r_squared3=summary(model3)$r.squared
r_squared_label3=paste0("R² = ", round(r_squared3, 4))

p1 = ggplot(d15, aes(x = Score2, y = FIS2)) +
  geom_point(size = 3, alpha = 0.6, color = "steelblue") +  
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1.5) +
  annotate("text", x = Inf, y = Inf, label = r_squared_label3, 
           hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5),
                     labels = c("0", "0.5", "1", "1.5", "2", "2.5"),
                     limits = c(0,2.5)) + #To fix axis from large scores with FIS NA values
  labs(
    title = "Effect of mtVLM (tRNA + PE) Score on FIS",
    x = "mtVLM Score (tRNA + PE) ",
    y = "FIS"
  ) +
  theme_bw(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    axis.title = element_text(face = "bold"),  
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank() 
  )

#Above repeated for the SF-36 Score

PearsonsCor4 = cor.test(SF2, Score2, method = "pearson", use = "complete.obs") 
PearsonsCor4 

model4=lm(SF2 ~ Score2, data = d15)
r_squared4=summary(model4)$r.squared
r_squared_label4=paste0("R² = ", round(r_squared4, 4))

p2 = ggplot(d15, aes(x = Score2, y = SF2)) +
  geom_point(size = 3, alpha = 0.6, color = "steelblue") +  
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1.5) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5),
                     labels = c("0", "0.5", "1", "1.5", "2", "2.5"),
                     limits = c(0,2.5)) +
  annotate("text", x = Inf, y = Inf, label = r_squared_label4,  
           hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
  labs(
    title = "Effect of mtVLM (tRNA + PE) Score on SF-36",
    x = "mtVLM Score (tRNA + PE) ",
    y = "SF-36"
  ) +
  theme_bw(base_size = 14) +  
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    axis.title = element_text(face = "bold"), 
    panel.grid.major = element_line(color = "grey90"),  
    panel.grid.minor = element_blank()
  )

p1 + p2
