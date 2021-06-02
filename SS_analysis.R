library(tidyverse)
library(cowplot)
library(extrafont)
library(pROC)
library(lme4)
library(gt)

##############Multi-Level Logistic Regression Suspended Sediment Dose-Response Model############

#Code written by Greg Courtice, last updated June 2, 2021

#Database consists of Adult and Juvenile Salmonids, in addition to some non-salmonid observations
#Non-salmonid, egg/larvae, and habitat-related observations were removed to assess stress-responses in salmonids to ambient SS exposures
#Predictors: Selected model for Figure 3 uses ln(dose). Alternate predictors: ln(SSC), ln(DoE)
#Life stage also considered to account for unique sensitivities between adult and juveniles
#Unconstrained slopes among successive regressions (i.e. not proportional-odds assumption)

#Effect categories used in study with respect to SEV score (Newcombe and Jensen 1996):
###EOI 1: Behavioural (0<SEV<=4)
###EOI 2: Minor Physiological (4<SEV<=7, no observation @ SEV=7)
###EOI 3: Major Physiological and Lethal (7<SEV<=14)


#import figure font and add it to general theme function
loadfonts(device = "win")
theme_greg <- function () { 
  theme_nothing() %+replace% 
    theme(
      text = element_text(size = 12,
                          family = "Arial")
      #add stuff here
    )
}

#import data

data <- read_csv("data.csv") %>%
  filter(Group %in% c("AS", "JS")) %>%
  filter(Habitat_Related == "No")


# Define Model effect categories and predictor variables
# Study_Year distinguishes studies published in different years by same authors
# NStudy determines how many observations published in each study

# Log_Conc and Log_Dur are natural log of the respective predictor
# LS is life stage (adult or juvenile salmonid)
# Species_NA groups all observations where species is not identified into a single group
# Species_group separates juvenile and adults into distinct groups (we don't expect juveniles and adults to respond similarly to SS exposures)
# Effect is the Effect of Interest (EOI)
# Effect1 is the categorization for EOI<=1 (EOI of 1 = 0, EOI of 2 or 3 = 1)
# Effect2 is the categorization for EOI<=2 (EOI of 1 or 2 = 0, EOI of 3 = 1)

data <- mutate(data, Study_Year = ifelse(Group == "AS", paste0(Study, Year), paste0(Study, Year, "J"))) %>%
  group_by(Study_Year) %>%
  mutate(NStudy = n()) %>%
  ungroup()
data <- mutate(data,
               LS = Group,
               Log_Conc = log(Concentration),
               Log_Dur = log(Duration),
               Dose = Log_Conc + Log_Dur,
               Effect = ifelse(SEV < 5, 1, ifelse(SEV > 4 & SEV < 8, 2, 3)),
               Effect1 = ifelse(SEV < 5, 0, 1),
               Effect2 = ifelse(SEV < 8, 0, 1),
               Species_NA = ifelse(is.na(Species), "Not Available", Species),
               Species_group = ifelse(Group == "AS", Species_NA, paste0(Species_NA, "J"))
)


#Create axis tick labels
dur_labs <- c("1.2 minutes",
              "10 minutes", 
              "1 hour", 
              "6 hours", 
              "1 day",
              "1 week",
              "1 month",
              "6 months",
              "2 years")

conc_labs <- c("0.7",
               "5",
               "10",
               "25",
               "50",
               "100",
               "250",
               "1 000",
               "5 000",
               "50 000",
               "207 000")

#create color blind friendly color palette
palette  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


####### Figure 1: Raw data visualization with univariate distribution projections

#main panel
p_ct <- 
  ggplot(data)+
  theme_greg()+
  geom_point(aes(x = Duration,y = Concentration,
                 shape = as.factor(Effect),
                 color = as.factor(Effect)),
             size = 4.5)+
  xlab("Exposure Duration") +
  ylab(expression(atop("SSC",(mgL^-1))))+
  scale_x_continuous(trans = "log",
                     breaks = c(1.2/60,10/60,1,6,24,168,720,24*30*6,17520),
                     labels = dur_labs) +
  scale_y_continuous(trans = "log",
                     breaks = c(0.66,5,10,25,50,100,250,1000,5000,50000,207000),
                     labels = conc_labs) +
  scale_shape_manual(labels = c(
    "Behavioural   ",
    "Minor Physiological   ",
    "Major Physiological and Lethal   "
  ),
  values = c(1,2,5))+
  scale_color_manual(values = c(palette[2],palette[3],palette[1]))+
  coord_cartesian(
    clip = "off")+
  geom_segment(x = log(min(data$Duration)),
               xend = log(quantile(data$Duration)[2]),
               y = log(0),
               yend = log(0),
               size = 1)+
  geom_segment(x = log(quantile(data$Duration)[4]),
               xend = log(max(data$Duration)),
               y = log(0),
               yend = log(0),
               size = 1)+
  geom_point(x = log(median(data$Duration)),
             y = log(0)+0.1,
             size = 2)+
  geom_segment(y = log(min(data$Concentration)),
               yend = log(quantile(data$Concentration)[2]),
               x = log(0),
               xend = log(0),
               size = 1)+
  geom_segment(y = log(quantile(data$Concentration)[4]),
               yend = log(max(data$Concentration)),
               x = log(0),
               xend = log(0),
               size = 1)+
  geom_point(y = log(median(data$Concentration)),
             x = log(0),
             size = 2)+
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    # axis.line = element_line(colour = "black"),
    legend.position = "bottom",
    axis.text.y = element_text(size = 16,
                               hjust = 1,
                               margin = margin(r = unit(10, "cm"))),
    axis.text.x = element_text(size = 16,
                               margin = margin(t = unit(5, "cm")),
                               angle = -30,
                               hjust = 0),
    axis.title.x = element_text(margin = margin(t = unit(10,"cm"))),
    axis.title.y = element_text(size = 16,
                                angle = 0,
                                vjust = 0.5),
    axis.title = element_text(size = 16),
    axis.ticks.length = unit(0.2,"cm"),
    axis.ticks = element_line(color = "black"))

#duration projection above main panel

set.seed(3)
p_te <- ggplot(data)+
  theme_greg()+
  geom_jitter(aes(x = Duration, y = as.factor(Effect),
                  shape = as.factor(Effect),
                  color = as.factor(Effect)),
              size = 3,
              height = 0.2,
              width = 0)+
  scale_x_continuous(trans = "log"
                     
  )+
  scale_shape_manual(labels = c(
    "Behavioural",
    "Minor Physiological",
    "Major Physiological and Lethal"
  ),
  values = c(1,2,5))+
  scale_color_manual(values = c(palette[2],palette[3],palette[1]))

#concentration projection to the right of main panel
set.seed(1)
p_ce <- ggplot(data)+
  theme_nothing()+
  geom_jitter(aes(x = as.factor(Effect), y = Concentration,
                  shape = as.factor(Effect),
                  color = as.factor(Effect)),
              size = 3,
              height = 0,
              width = 0.25)+
  scale_y_continuous(trans = "log") +
  scale_shape_manual(labels = c(
    "Behavioural",
    "Minor Physiological",
    "Major Physiological and Lethal"
  ),
  values = c(1,2,5))+
  scale_color_manual(values = c(palette[2],palette[3],palette[1]))


#create life stage rugplots between main panel and univariate panels

p_rug_dur <- 
  ggplot(data)+
  theme_greg()+
  geom_segment(aes(x = Duration, 
                   xend = Duration, 
                   y = ifelse(Group == "AS",0.53,0), 
                   yend = ifelse(Group == "AS",1,0.47)),
               colour = ifelse(data$Effect == 1, palette[2],ifelse(data$Effect == 2,palette[3],palette[1])),
               size = 0.3)+
  scale_x_continuous(trans = "log")

p_rug_conc <- 
  ggplot(data)+
  theme_greg()+
  geom_segment(aes(y = Concentration, 
                   yend = Concentration, 
                   x = ifelse(Group == "AS",0.53,0), 
                   xend = ifelse(Group == "AS",1,0.47)),
               colour = ifelse(data$Effect == 1, palette[2],ifelse(data$Effect == 2,palette[3],palette[1])),
               size = 0.3)+
  scale_y_continuous(trans = "log")


#Create sub-plots for appropriate panel spacing
p.top <- plot_grid(p_te,
                   NULL,
                   p_rug_dur,
                   ncol = 1,
                   rel_heights = c(1,0.3,1),
                   align = "v")

p.right <- plot_grid(p_rug_conc,
                     NULL,
                     p_ce,
                     ncol = 3,
                     rel_widths = c(1,0.3,1),
                     align = "h")

#Create annotation panel to serve as legend
p.annot <- ggplot()+
  theme_greg()+
  annotate("rect", xmin = 0, xmax = 5, ymin = 0, ymax =5,
           alpha = 0)+
  annotate("text", x = 0.6, y = 0.4, label = "Adult", hjust = 1,family = "Arial")+
  annotate("text", x = 1.9, y = 1.5, label = "Juvenile", hjust = 1,family = "Arial")+
  annotate("text", x = 3.3, y = 3.3, label = "Behavioural", hjust = 1,family = "Arial")+
  annotate("text", x = 4.1, y = 4.1, label = "Minor Phys.", hjust = 1,family = "Arial")+
  annotate("text", x = 4.8, y = 4.8, label = "Major Phys. and Lethal", hjust = 1,family = "Arial")+
  geom_segment(aes(x = 3.3,
                   xend = 3.3,
                   y = 3,
                   yend = 0),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 1.2,
                   xend = 0,
                   y = 3.3,
                   yend = 3.3),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 4.05,
                   xend = 4.05,
                   y = 3.85,
                   yend = 0),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 2,
                   xend = 0,
                   y = 4.1,
                   yend = 4.1),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 4.8,
                   xend = 4.8,
                   y = 4.5,
                   yend = 0),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 1,
                   xend = 0,
                   y = 4.8,
                   yend = 4.8),
               linetype = "dotted",
               size = 0.5)

#######Plot Figure 1#############
p.f1 <- plot_grid(
  p.top,
  NULL,
  p.annot,
  NULL,NULL,NULL, #row to close up gap between plots
  p_ct+theme(legend.position = "none",
             plot.margin = margin(r = -5)),
  NULL,
  p.right+theme(legend.position = "none",
                plot.margin = margin(l = -5)),
  NULL,NULL,NULL,
  ncol = 3,
  rel_widths = c(1,-0.2,0.5),
  rel_heights = c(0.5,-0.2,1,0.01),
  align = "vh"
)

#Plot univariate dose distribution


dose_vals = 20 #value used for height of dose value identifiers

data_f2 <- data%>%mutate(f2_cat = case_when(
  
  #behavioural
  SEV == 1 ~ 0.6,
  SEV == 3 ~ 0.8,
  SEV == 4 ~ 1,
  
  #minor phys
  SEV == 5 ~ 1.7,
  SEV == 6 ~ 1.9,
  
  #major phys and lethal
  SEV == 8 ~ 2.6,
  SEV == 9 ~ 2.8,
  SEV == 10 ~ 3,
  SEV == 11 ~ 3.2,
  SEV == 12 ~ 3.4,
  SEV == 13 ~ 3.6,
  SEV == 14 ~ 3.8
  
  
  
))

pdose0 <- ggplot(data_f2)+
  theme_greg()+
  geom_point(aes(x = Dose, y = f2_cat,
                 shape = as.factor(Effect),
                 color = as.factor(Effect)),
             size = 4.5)+
  scale_shape_manual(labels = c(
    "Behavioural",
    "Minor Physiological",
    "Major Physiological and Lethal"
  ),
  values = c(1,2,5))+
  xlab(expression(paste(ln(dose, mg%.%h%.%L^-1))))+
  coord_cartesian(xlim =  c(min(data$Dose),max(data$Dose)+1))+
  scale_x_continuous(
    breaks = c(min(data$Dose),5,10,15,max(data$Dose)),
    labels = c(round(min(data$Dose), digits = 1),5,10,15,round(max(data$Dose), digits = 1)),
    sec.axis = dup_axis(name = "Typical Exposure Values in Literature",
                        # labels = c("3 min\n20mg/L",
                        #            "6 min\n100 mg/L",
                        #            "1 hour\n25 mg/L",
                        #            "2.5 hours\n60 mg/L",
                        #            "24 hours\n25 mg/L",
                        #            "9 hours\n500 mg/L",
                        #            "1 week\n300 mg/L",
                        #            "3 days\n185 mg/L",
                        #            "1 month\n300 mg/L",
                        #            "2 weeks\n4 700 mg/L",
                        #            "1 year\n1 040 mg/L"
                        #            
                        # ),
                        labels = c("3 min",
                                   expression(20~mgL^-1),
                                   "1 hour",
                                   expression(25~mgL^-1),
                                   "24 hours",
                                   expression(25~mgL^-1),
                                   "9 hours",
                                   expression(500~mgL^-1),
                                   "1 week",
                                   expression(300~mgL^-1),
                                   "1 month",
                                   expression(300~mgL^-1),
                                   "4 days",
                                   expression(13000~mgL^-1),
                                   "1 year",
                                   expression(1040~mgL^-1)
                                   
                        ),
                        
                        breaks = 
                          c(log(.05*20)+0.1,
                            log(.05*20)-0.6,
                            log(1*25)+0.2,
                            log(1*25)-0.5,
                            log(24*25)+0.1,
                            log(24*25)-0.6,
                            log(9*500)+0.4,
                            log(9*500)-0.3,
                            log(24*7*300)+0.1,
                            log(24*7*300)-0.6,
                            log(720*300)+0.5,
                            log(720*300)-0.2,
                            log(96*13000)+0.5,
                            log(96*13000)-0.2,
                            log(8760*1040)+0.7,
                            log(8760*1040)-0
                          )
    )
  )+
  scale_color_manual(values = c(palette[2],palette[3],palette[1]))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.x.bottom = element_text(size = 20,
                                          margin = margin(t=-15)),
        axis.text.x.top = element_text(size = 14,
                                       angle = 60,
                                       hjust = 0,
                                       margin = margin(b=-25)),
        axis.title.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.ticks.length.x.bottom = unit(0.2,"cm"),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20,
                                   hjust = 1,
                                   family = "Arial"))+
  scale_y_continuous(
    labels = c("Behavioural",
               "Minor Phys.",
               "Lethal\nMajor Phys."),
    breaks = c(0.8,1.8,3.3)
  )+
  geom_segment(x = log(.05*20),
               xend = log(.05*20),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x = log(1*25),
               xend = log(1*25),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x =log(24*25),
               xend = log(24*25),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x = log(9*500),
               xend = log(9*500),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x = log(24*7*300),
               xend = log(24*7*300),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x = log(720*300),
               xend = log(720*300),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x = log(96*13000),
               xend = log(96*13000),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  geom_segment(x = log(8760*1040),
               xend = log(8760*1040),
               y = 0,
               yend = dose_vals,
               linetype = "dashed",
               colour = "lightgrey",
               size = 0.3)+
  coord_cartesian(ylim = c(0,5),
                  xlim = c(min(data$Dose), 23))+
  geom_segment(x = min(data$Dose),
               xend = quantile(data$Dose)[2],
               y = 0,
               yend = 0,
               size = 1)+
  geom_segment(x = quantile(data$Dose)[4],
               xend =max(data$Dose),
               y = 0,
               yend = 0,
               size = 1)+
  geom_point(x = median(data$Dose),
             y = 0,
             size = 2)

pdose <- pdose0 +
  annotate("text", x = 19, y = 0.6, 
           label = "Alarm Reaction",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 0.8, 
           label = "Avoidance Response", 
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 1, 
           label = "Reduction in Feeding", 
           hjust = 0,
           family = "Arial")+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 0.6,
                   yend = 0.6),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 0.8,
                   yend = 0.8),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x =17.8,
                   xend = 18.7,
                   y = 1,
                   yend = 1),
               linetype = "dotted",
               size = 0.5)+
  annotate("text", x = 19, y = 1.55,
           label = "Sublethal Stress and\nIncreased Coughing",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 2.13,
           label = "Long-Term Survivability\nIssues and\nHeightened Blood Glucose",
           hjust = 0,
           family = "Arial")+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 1.7,
                   yend = 1.7),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x =17.8,
                   xend = 18.7,
                   y = 1.9,
                   yend = 1.9),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 18.7,
                   xend = 18.7,
                   y = 1.7,
                   yend = 1.4),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x =18.7,
                   xend = 18.7,
                   y = 1.9,
                   yend = 2.4),
               linetype = "dotted",
               size = 0.5)+
  annotate("text", x = 19, y = 2.6,
           label = "Gill Tissue Damage",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 2.8,
           label = "Growth Rate Reduced",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 3,
           label = "Mortality 0-20%",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 3.2,
           label = "Mortality 20-40%",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 3.4,
           label = "Mortality 40-60%",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 3.6,
           label = "Mortality 60-80%",
           hjust = 0,
           family = "Arial")+
  annotate("text", x = 19, y = 3.8,
           label = "Mortality 80-100%",
           hjust = 0,
           family = "Arial")+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 2.6,
                   yend = 2.6),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 2.8,
                   yend = 2.8),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 3,
                   yend = 3),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 3.2,
                   yend = 3.2),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 3.4,
                   yend = 3.4),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x = 17.8,
                   xend = 18.7,
                   y = 3.6,
                   yend = 3.6),
               linetype = "dotted",
               size = 0.5)+
  geom_segment(aes(x =18.2,
                   xend = 18.7,
                   y = 3.8,
                   yend = 3.8),
               linetype = "dotted",
               size = 0.5)


p_rug_dose <- 
  ggplot(data)+
  theme_greg()+
  geom_segment(aes(x = Dose, 
                   xend = Dose, 
                   y = ifelse(Group == "JS",0.53,0), 
                   yend = ifelse(Group == "JS",1,0.47)),
               colour = ifelse(data$Effect == 1, palette[2],ifelse(data$Effect == 2,palette[3],palette[1])),
               size = 0.3)+
  scale_y_continuous(breaks = c(0.25,0.75),
                     labels = c("Adult","Juvenile"))+
  coord_cartesian(xlim = c(min(data$Dose),23))+
  theme(axis.text.y = element_text(size = 20,
                                   hjust = 1,
                                   family = "Arial"))

######Plot Figure 2###############
p.f2 <- plot_grid(NULL,NULL,
          NULL,p_rug_dose,
          NULL,NULL,
          NULL,pdose,
          NULL,NULL,
          ncol = 2,
          align = "v",
          rel_widths = c(0.01,1),
          # rel_heights = c(0.1,0.08,-0.1,0.5,0.02)
          rel_heights = c(0.11,0.07,-0.14,0.5,0.02)
)








#####Fit logistic regression model ###########

#fit all alternative models

#P(Effect>Behavioural) models:

fit1 <- glmer(as.factor(Effect1) ~ Dose + (1|Study_Year),
              family = "binomial",
              data = filter(data, Effect != 2))

alpha1 = fixef(fit1)[1]
beta1 = fixef(fit1)[2]

sigma1 = as.data.frame(VarCorr(fit1))[5]

fit1.ls <- glmer(as.factor(Effect1) ~ Dose + Group + (1|Study_Year),
                 family = "binomial",
                 data = filter(data, Effect != 2))

alpha1.ls = fixef(fit1.ls)[1]
beta1.ls = c(fixef(fit1.ls)[2], 
             fixef(fit1.ls)[3])

sigma1.ls = as.data.frame(VarCorr(fit1.ls))[5]

fit1.fe <- glm(as.factor(Effect1)~Dose,
               family = "binomial",
               data = filter(data, Effect != 2))
alpha1.fe <- coefficients(fit1.fe)[1]
beta1.fe <- coefficients(fit1.fe)[2]

fit1.fe.ls <- glm(as.factor(Effect1)~Dose+Group,
                  family = "binomial",
                  data = filter(data, Effect != 2))
alpha1.fe.ls <- coefficients(fit1.fe.ls)[1]
beta1.fe.ls <- c(coefficients(fit1.fe.ls)[2],
                 coefficients(fit1.fe.ls)[3])

fit1.conc <- glmer(as.factor(Effect1)~Log_Conc + (1|Study_Year),
                   family = "binomial",
                   data = filter(data, Effect != 2))

alpha1.conc = fixef(fit1.conc)[1]
beta1.conc = fixef(fit1.conc)[2]


sigma1.conc = as.data.frame(VarCorr(fit1.conc))[5]


fit1.dur <- glmer(as.factor(Effect1)~Log_Dur + (1|Study_Year),
                  family = "binomial",
                  data = filter(data, Effect != 2))


alpha1.dur = fixef(fit1.dur)[1]
beta1.dur = fixef(fit1.dur)[2]


sigma1.dur = as.data.frame(VarCorr(fit1.dur))[5]


#Estimate exceedance probabilities, given model fits
m.probs <- mutate(filter(data, Effect !=2),
                  p1 = 1/(1+exp(-(alpha1 + beta1 * Dose))),
                  p1.ls = ifelse(Group == "AS",
                                 1/(1+exp(-(alpha1.ls + beta1.ls[1] * Dose))),
                                 1/(1+exp(-(alpha1.ls + beta1.ls[1] * Dose + beta1.ls[2])))
                  ),
                  p1.fe = 1/(1+exp(-(alpha1.fe + beta1.fe * Dose))),
                  p1.fe.ls = ifelse(Group == "AS",
                                    1/(1+exp(-(alpha1.fe.ls + beta1.fe.ls[1] * Dose))),
                                    1/(1+exp(-(alpha1.fe.ls + beta1.fe.ls[1] * Dose + beta1.fe.ls[2])))
                  ),
                  p1.conc = 1/(1+exp(-(alpha1.conc + beta1.conc * Log_Conc))),
                  p1.dur = 1/(1+exp(-(alpha1.dur + beta1.dur * Log_Dur))))


#generate roc curves
m.pred1 <- roc(m.probs,Effect1,p1)
m.pred1.ls <- roc(m.probs,Effect1,p1.ls)
m.pred1.fe <- roc(m.probs,Effect1,p1.fe)
m.pred1.fe.ls <- roc(m.probs,Effect1,p1.fe.ls)
m.pred1.conc <- roc(m.probs,Effect1,p1.conc)
m.pred1.dur <- roc(m.probs,Effect1,p1.dur)

#model selection
m.summary <- tibble(Model = 1:6,
                    Predictors = c("ln(dose)",
                                   "ln(dose)",
                                   "ln(dose), Life Stage",
                                   "ln(dose), Life Stage",
                                   "ln(SSC)",
                                   "ln(DoE)"
                    ),
                    `Study Level` = c("Yes",
                                      "No",
                                      "Yes",
                                      "No",
                                      "Yes",
                                      "Yes"
                    ),
                    AIC = c(round(AIC(fit1), digits = 1),
                            round(AIC(fit1.fe), digits = 1),
                            round(AIC(fit1.ls), digits = 1),
                            round(AIC(fit1.fe.ls), digits = 1),
                            round(AIC(fit1.conc), digits = 1),
                            round(AIC(fit1.dur), digits = 1)),
                    AUC = c(round(auc(m.pred1), digits = 3),
                            round(auc(m.pred1.fe), digits = 3),
                            round(auc(m.pred1.ls), digits = 3),
                            round(auc(m.pred1.fe.ls), digits = 3),
                            round(auc(m.pred1.conc), digits = 3),
                            round(auc(m.pred1.dur), digits = 3)),
                    Alpha = c(round(alpha1,digits = 1),
                              round(alpha1.fe, digits = 1),
                              round(alpha1.ls, digits = 1),
                              round(alpha1.fe.ls, digits = 1),
                              round(alpha1.conc, digits = 1),
                              round(alpha1.dur, digits = 1)),
                    Beta = c(round(beta1,digits = 1),
                             round(beta1.fe, digits = 1),
                             round(beta1.ls[1], digits = 1),
                             round(beta1.fe.ls[1], digits = 1),
                             round(beta1.conc, digits = 1),
                             round(beta1.dur, digits = 1)),
                    ls1 = c(NA,
                            NA,
                            round(beta1.ls[2], digits = 1),
                            round(beta1.fe.ls[2], digits = 1),
                            NA,
                            NA),
                    Sigma = c(round(sigma1,digits = 1),
                              NA,
                              round(sigma1.ls, digits = 1),
                              NA,
                              round(sigma1.conc, digits = 1),
                              round(sigma1.dur, digits = 1))
                    
                    
                    
)




#######Table 1##########

m.tbl <- gt(m.summary)%>%
  cols_align("center")%>%
  cols_label(
    ls1 = "Life Stage",
    Sigma = "Study SD"
  )%>%
  tab_footnote(
    footnote = "Akaike Information Criterion (AIC) value represents relative model fit performance with a penalty for overfitting. A smaller AIC value indicates a more parsimonious model fit when compared to alternatives. Best fitting value is bolded.",
    locations = cells_column_labels(
      columns = vars(AIC))
  )%>%
  tab_footnote(
    footnote = "Area under the receiver operating characteristic curve (AUC) represents model classification performance. A value of 1.0 indicates a perfect classification while a value of 0.5 indicates a random classification. Best classification values are bolded.",
    locations = cells_column_labels(
      columns = vars(AUC)
    )
  )%>%
  tab_footnote(
    footnote = "Standard deviation of regression error for study level intercepts",
    locations = cells_column_labels(
      columns = vars(Sigma)
    )
  )%>%
  tab_footnote(
    footnote = "Model selected to present in Figure 3.",
    locations = cells_body(
      columns = vars(Model),
      rows = AIC == min(AIC))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = vars(AIC),
      rows = AIC == min(AIC)
    )
  )%>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = vars(AUC),
      rows = AUC == max(AUC)
    )
  )



########Figure 3 ##########




#create data for probability curve and calculate approximate confidence intervals
#CI's are only approximate due to presence of study variability
#use 2*SE instead of 1.96 to emphasize approximate nature of these estimates
#code for confidence intervals adapted from: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#lme4

newdata <- data.frame(Dose = seq(from = 0, to = 20, by = 0.01),
                      Effect1 = 0)

mm <- model.matrix(terms(fit1), newdata)

newdata$Effect1 <- mm %*% fixef(fit1)

pvar1 <- diag(mm %*% tcrossprod(vcov(fit1),mm))

newdata <- data.frame(
  newdata,
  plo = newdata$Effect1 - 1.96*sqrt(pvar1),
  phi = newdata$Effect1 + 1.96*sqrt(pvar1)
)




alpha1 = fixef(fit1)[1]
beta1 = fixef(fit1)[2]


set.seed(1327)
p_f3 <- ggplot(filter(data, Effect %in% c(1,3)))+
  theme_greg()+
  stat_function(fun = function(x) 1/(1+exp(-(alpha1 + beta1*x))),
                mapping = aes(x = Dose, linetype = "P(Effect>Behavioural)"),
                size = 1,
                data = newdata)+
  geom_ribbon(data = newdata,
              mapping = aes(x = Dose,
                            ymin = 1/(1+exp(-plo)), 
                            ymax = 1/(1+exp(-phi))
  ),
  alpha = 0.1)+
  geom_jitter(data = filter(data, Effect == 1),
              aes(x = Dose, y = -0.05),
              shape = 1,
              color = palette[2],
              size = 4,
              height = 0.02,
              width = 0)+
  geom_jitter(data = filter(data, Effect == 3),
              aes(x = Dose, y = 1.05),
              shape = 5,
              color = palette[1],
              size = 4,
              height = 0.02,
              width = 0)+
  geom_jitter(data = filter(data, Effect == 2),
              aes(x = Dose, y = 1.25),
              shape = 2,
              color = palette[3],
              size = 4,
              height = 0.02,
              width = 0)+
  xlab(expression(paste(ln(dose, mg%.%h%.%L^-1))))+
  ylab("Probability")+
  scale_y_continuous(
    breaks = c(0.02,0.1,0.5,0.9),
    labels = c(expression(paste("P(",pi[6.4],") = 0.02-")),"0.1-","0.5-","0.9-"),
    sec.axis = dup_axis(name = "Binary Effect\nClassification",
                        breaks = c(-0.05,1.05,1.25),
                        labels = c("Behavioural\n(BHV)",
                                   "Major Phys.\nand Lethal\n(MaPL)",
                                   "Excluded Data\n(Minor Phys.)")
                        
                        
    ))+
  scale_x_continuous(
    breaks = c(0,5,10,15,20),
    sec.axis = dup_axis(
      
      breaks = c(log(25*24)
      ),
      labels = c("Example Exposure Limit\nln(dose) = 6.4")
    )
  )+
  theme(axis.title.y.left = element_text(size = 20),
        axis.text.y.left = element_text(hjust = 1,
                                        margin = margin(l=-40),
                                        size = 20),
        axis.text.y.right = element_text(hjust = 0,
                                         margin = margin(l = -15)),
        axis.text.x.top = element_text(size = 12,
                                       margin = margin(t = 20,
                                                       b = 5)),
        axis.title.x.bottom = element_text(size = 20),
        axis.title.x.top = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x.bottom = unit(0.2,"cm"),
        axis.ticks.x.bottom = element_line(color = "black"),
        axis.text.x.bottom = element_text(size = 20,
                                          family = "Arial"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank())


tp = nrow(filter(data,Dose >= log(600) & Effect ==3))  #true positive
tn = nrow(filter(data,Dose < log(600) & Effect == 1)) #true negative
fn = nrow(filter(data,Dose < log(600) & Effect ==3))  #false negative
fp = nrow(filter(data,Dose >= log(600) & Effect == 1)) #false positive

tpr = tp/(tp+fn) #true positive rate
tnr = tn/(tn+fp) #true negative rate

tpr = round(tpr, digits = 2)

roc1 <- ggroc(list(
  m.pred1.ls,
  m.pred1.fe,
  m.pred1.fe.ls,
  m.pred1.conc,
  m.pred1.dur,
  m.pred1),
  legacy.axes = TRUE
) +
  scale_color_manual(values = c("grey",
                                "grey",
                                "grey",
                                "grey",
                                "grey",
                                "black")
  )+
  theme_bw()+
  geom_segment(x = 0,
               xend = 1,
               y = 0,
               yend = 1,
               linetype = "dashed",
               size = 0.3)+
  geom_segment(x = 0,
               xend = 0,
               y = 0,
               yend = 1,
               linetype = "dashed",
               size = 0.3)+
  geom_segment(x = 0,
               xend = 1,
               y = 1,
               yend = 1,
               linetype = "dashed",
               size = 0.3)+
  geom_point(aes(x = 1-tnr, y = tpr),
             size = 5,
             color = "black")+
  annotate("text", x = 1-0.3, y = 0.35,
           label = expression(pi[6.4]),
           family = "Arial",
           size = 10)+
  geom_segment(x = 1-tnr,
               xend = 1-0.4,
               y = tpr,
               yend = 0.5,
               size = 0.05)+
  scale_y_continuous(breaks = c(tpr),
                     label = c(paste0(tpr," -")))+
  scale_x_continuous(breaks = c(0,1),
                     labels = c(1,0),
                     sec.axis = dup_axis())+
  ylab("TMR")+
  xlab("1-TBR")+
  theme(axis.title.x.top = element_text(size = 12),
        axis.title.y = element_text(size = 12,
                                    angle = 0,
                                    vjust = 0.5),
        axis.text.y.left = element_text(size = 12,
                                        margin = margin(l = -20)),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12,
                            family = "Arial"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

p.f3 <- p_f3+
  annotation_custom(ggplotGrob(roc1), xmin = 12.5, xmax = 21.5, 
                    ymin = 0.4, ymax = 0.95)+
  geom_segment(x = log(25*24),
               xend = log(25*24),
               y = 1.5,
               yend = -0.1,
               linetype = "dotted")+
  geom_segment(x = log(25*24),
               xend = 10,
               y = 0.4,
               yend = 0.4,
               arrow = arrow(length = unit(0.015, "npc")))+
  geom_segment(x = log(25*24),
               xend = 2*log(25*24)-10,
               y = 0.5,
               yend = 0.5,
               arrow = arrow(length = unit(0.015, "npc")))+
  geom_point(aes(x = log(25*24),
                 y = 1/(1+exp(-(alpha1+beta1*log(25*24))))),
             size = 5)+
  annotate("text",
           x = log(25*24)-3,
           y = 0.25,
           label = expression(pi[6.4]),
           family = "Arial",
           size = 10)+
  annotate("text",
           x = 8.75,
           y = 0.35,
           label = "Classified as MaPL",
           family = "Arial",
           hjust = 0)+
  annotate("text",
           x = 2*log(25*24)-8.5,
           y = 0.56,
           label = "Classified as BHV",
           family = "Arial",
           hjust = 1)+
  annotate("text",
           x = 0.5*(3.16+4.85)-2,
           y = 1.14,
           label = "False BHV (FB)",
           family = "Arial")+
  annotate("text",
           x = 0.5*(3.16+4.85)-2,
           y = 1.08,
           label = "[No Observations]",
           family = "Arial")+
  #False MaPL
  annotate("segment",
           x = 6.6,
           xend = 13.8,
           y = -0.12,
           yend = -0.12,
           size = 0.1)+
  annotate("segment",
           x = 6.6,
           xend = 6.6,
           y = -0.12,
           yend = -0.09,
           size = 0.1)+
  annotate("segment",
           x = 13.8,
           xend = 13.8,
           y = -0.12,
           yend = -0.09,
           size = 0.1)+
  annotate("text",
           x = 0.5*(6.6+13.8),
           y = -0.16,
           label = "False MaPL (FM)",
           family = "Arial")+
  #True BHV
  annotate("segment",
           x = 2*log(25*24)-6.6,
           xend = -0.4,
           y = -0.12,
           yend = -0.12,
           size = 0.1)+
  annotate("segment",
           x = 2*log(25*24)-6.6,
           xend = 2*log(25*24)-6.6,
           y = -0.12,
           yend = -0.09,
           size = 0.1)+
  annotate("segment",
           x = -0.4,
           xend = -0.4,
           y = -0.12,
           yend = -0.09,
           size = 0.1)+
  annotate("text",
           x = 0.5*(2*log(25*24)-6.6-0.4),
           y = -0.16,
           label = "True BHV (TB)",
           family = "Arial")+
  #True MaPL
  annotate("segment",
           x = 8.15,
           xend = max(data$Dose),
           y = 1.1,
           yend = 1.1,
           size = 0.1)+
  annotate("segment",
           x = 8.15,
           xend = 8.15,
           y = 1.07,
           yend = 1.1,
           size = 0.1)+
  annotate("segment",
           x = max(data$Dose),
           xend = max(data$Dose),
           y = 1.07,
           yend = 1.1,
           size = 0.1)+
  annotate("text",
           x = 0.5*(8.15+max(data$Dose)),
           y = 1.14,
           label = "True MaPL (TM)",
           family = "Arial")+
  annotate("text",
           x = 13.5,
           y = 0.3,
           label = "True MaPL Rate (TMR):",
           family = "Arial",
           hjust = 0
  )+
  annotate("text",
           x = 19.5,
           y = 0.3,
           label = expression(frac(TM,TM+FB)),
           family = "Arial",
           hjust = 0.5
  )+
  annotate("text",
           x = 19.5,
           y = 0.15,
           label = expression(frac(TB,TB+FM)),
           family = "Arial",
           hjust = 0.5
  )+
  annotate("text",
           x = 13.5,
           y = 0.15,
           label = "True BHV Rate (TBR):",
           family = "Arial",
           hjust = 0
  )+
  theme(plot.margin = margin(15,15,15,15),
        axis.title.y.left = element_text(margin = margin(r = -50),
                                         vjust = 0.6)
  )+
  geom_segment(x = 20,
               xend = 18.3,
               y = 1.05,
               yend = 1.05,
               linetype = "dotted",
               size = 0.5)+
  geom_segment(x = 20,
               xend = 14.5,
               y = -0.05,
               yend = -0.05,
               linetype = "dotted",
               size = 0.5)+
  geom_segment(x = 20,
               xend = 17.3,
               y = 1.25,
               yend = 1.25,
               linetype = "dotted",
               size = 0.5)+
  geom_segment(x = log(25*24) - 0.1,
               y = 1/(1+exp(-(alpha1+beta1*log(25*24))))+0.01,
               xend = log(25*24)-2.5,
               yend = 0.2,
               size = 0.05)

####### Plot Figure 3 #######

p.f3
