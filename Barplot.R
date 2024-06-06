#### Load packages ----

library(tidyverse)
library(ggbreak) #  set breakpoints for both x and y axes

#### Import data ----

MI_casestudy <- read.csv("MI_casestudy.csv")
MI_casestudy$Group <- factor(MI_casestudy$Group, levels = c("Mean", "RSE"))
MI_casestudy$Percentage <- factor(MI_casestudy$Percentage, levels = c("0", "30", "50", "80"))
MI_casestudy$Covariate <- factor(MI_casestudy$Covariate, levels = c("WT_CL", "WT_Vd"))

#### Barplot ----

# First, divide RSE by 100 to match the unit of Mean

# MI_casestudy <- MI_casestudy %>% 
  #mutate(Value = ifelse(Group == "RSE", Value / 100, Value))
# Should illustrate RSE in the % scale - 010524

# Calculate the difference between mean of 3 approaches with mean of complete data
MI_casestudy <- MI_casestudy %>%
  group_by(Covariate) %>%
  mutate(Rel_bias = ifelse(Group == "Mean", 
                           100 * (Value - Value[Group == "Mean" & Percentage == "0"]) / Value[Group == "Mean" & Percentage == "0"], NA))

# Compare accurary - mean parameter
Mean_diff <- MI_casestudy %>%
  dplyr::filter(Group == "Mean" & !Percentage == 0) %>%
  ggplot(aes(x = Percentage, y = Rel_bias, fill = Technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Covariate, labeller = labeller(Covariate = 
                                         c("WT_CL" = "Body weight on CL",
                                           "WT_Vd" = "Body weight on Vd")
  )) +
  theme_bw() +
  scale_fill_viridis_d() +
  xlab("Missing body weight (%)") +
  ylab("Relative bias (%)") +
  scale_y_continuous(limits = c(-110, 50), 
                     breaks = seq(-110, 50, by = 10), 
                     expand = c(0,0),
                     labels = scales::comma) +
  theme(strip.text.x = element_text(hjust = 0.5, size = 8, face = "bold", family = "sans"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        strip.text.y = element_text(vjust = 0.5, size = 8, face = "bold", family = "sans"), 
        axis.title = element_text(size = 7, family = "sans"),
        axis.text = element_text(size = 7, family = "sans"),
        legend.title = element_text(size = 6, family = "sans"),
        legend.text = element_text(size = 6, family = "sans"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = c(-10, 0, 10), linetype = "dashed", color = "gray50", size = 0.6)
Mean_diff

ggsave("Mean_diff.svg", 
       Mean_diff, 
       dpi = 300, 
       width = 19, 
       height = 9,
       unit = "cm")

# Compare precision - RSE
RSE <- MI_casestudy %>%
  dplyr::filter(Group == "RSE") %>%
  ggplot(aes(x = Percentage, y = Value, fill = Technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Covariate, labeller = labeller(Covariate = 
                                               c("WT_CL" = "Body weight on CL",
                                                 "WT_Vd" = "Body weight on Vd")
  )) +
  theme_bw() +
  scale_fill_viridis_d() +
  xlab("Missing body weight (%)") +
  ylab("Relative standard error (%)") +
  scale_y_continuous(limits = c(0, 850), 
                     breaks = seq(0, 850, by = 50), 
                     expand = c(0,0)) +
  scale_y_break(c(200, 750)) +
  theme(strip.text.x = element_text(hjust = 0.5, size = 8, face = "bold", family = "sans"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        strip.text.y = element_text(vjust = 0.5, size = 8, face = "bold", family = "sans"), 
        axis.title = element_text(size = 7, family = "sans"),
        axis.text = element_text(size = 7, family = "sans"),
        legend.title = element_text(size = 6, family = "sans"),
        legend.text = element_text(size = 6, family = "sans"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", size = 0.6)
RSE

ggsave("RSE.svg", 
       RSE, 
       dpi = 300, 
       width = 19, 
       height = 9,
       unit = "cm")