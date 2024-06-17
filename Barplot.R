# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggbreak) #  set breakpoints for both x and y axes

# Import data -----------------------------------------------------------

MI_casestudy <- read.csv("MI_casestudy.csv")
MI_casestudy$Group <- factor(MI_casestudy$Group, levels = c("Mean", "RSE"))
MI_casestudy$Percentage <- factor(MI_casestudy$Percentage, levels = c("0", "30", "50", "80"))
MI_casestudy$Covariate <- factor(MI_casestudy$Covariate, levels = c("WT_CL", "WT_Vd"))

# Barplot -----------------------------------------------------------

# Should illustrate RSE in the % scale - 010524

# Calculate the difference between mean of 3 approaches with mean of complete data
MI_casestudy <- MI_casestudy |> 
  group_by(Covariate) |>
  mutate(Rel_bias = ifelse(Group == "Mean", 
                           100 * (Value - Value[Group == "Mean" & Percentage == "0"]) / Value[Group == "Mean" & Percentage == "0"], NA))

# Compare accurary - mean parameter
Mean_diff <- MI_casestudy |>
  dplyr::filter(Group == "Mean" & !Percentage == 0) |>
  ggplot(aes(x = Percentage, y = Rel_bias, fill = Technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Covariate, labeller = labeller(Covariate = 
                                         c("WT_CL" = "Body weight on CL",
                                           "WT_Vd" = "Body weight on Vd")
  )) +
  theme_bw() +
  scale_fill_manual(values = c("#00407A", "#1E8DB0", "#F27405")) +
  xlab("Missing body weight (%)") +
  ylab("Relative bias (%)") +
  scale_y_continuous(limits = c(-110, 50), 
                     breaks = seq(-110, 50, by = 10), 
                     expand = c(0,0),
                     labels = scales::comma) +
  theme(strip.text.x = element_text(hjust = 0.5, size = 8, face = "bold", family = "sans", color = "#00407A"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        strip.text.y = element_text(vjust = 0.5, size = 8, face = "bold", family = "sans", color = "#00407A"),
        strip.background = element_rect(fill = "#F27405"),
        axis.title = element_text(size = 7, family = "sans", color = "#00407A"),
        axis.text = element_text(size = 7, family = "sans", color = "#00407A"),
        legend.title = element_text(size = 6, family = "sans", color = "#00407A"),
        legend.text = element_text(size = 6, family = "sans", color = "#00407A"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "#00407A"),
        #legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.2, 'cm'),
        legend.box.margin = unit(0.1, 'cm')) +
  geom_hline(yintercept = c(-10, 0, 10), linetype = "dashed", color = "#00407A", size = 0.6)

Mean_diff

# export the plot
setwd("./Plots/PAGE_2024")
ggsave("Mean_diff.svg", 
       Mean_diff, 
       dpi = 300, 
       width = 19, 
       height = 9,
       unit = "cm")
# return the working directory
Path = getwd()
setwd(dirname(dirname(Path)))

# Compare precision - RSE
RSE <- MI_casestudy |>
  dplyr::filter(Group == "RSE") |>
  ggplot(aes(x = Percentage, y = Value, fill = Technique)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Covariate, labeller = labeller(Covariate = 
                                               c("WT_CL" = "Body weight on CL",
                                                 "WT_Vd" = "Body weight on Vd")
  )) +
  theme_bw() +
  scale_fill_manual(values = c("#00407A", "#1E8DB0", "#F27405")) +
  xlab("Missing body weight (%)") +
  ylab("Relative standard error (%)") +
  scale_y_continuous(limits = c(0, 850), 
                     breaks = seq(0, 850, by = 50), 
                     expand = c(0,0)) +
  scale_y_break(breaks = c(200, 750), expand = FALSE, space = 0.2, scales = "fixed") +
  #scale_y_cut(breaks = c(200, 750)) +
  #scale_y_break(breaks = c(200, 750)) +
  theme(strip.text.x = element_text(hjust = 0.5, size = 8, face = "bold", family = "sans", color = "#00407A"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        strip.text.y = element_text(vjust = 0.5, size = 8, face = "bold", family = "sans", color = "#00407A"),
        strip.background = element_rect(fill = "#F27405"),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title = element_text(size = 7, family = "sans", color = "#00407A"),
        axis.text = element_text(size = 7, family = "sans", color = "#00407A"),
        legend.title = element_text(size = 6, family = "sans", color = "#00407A"),
        legend.text = element_text(size = 6, family = "sans", color = "#00407A"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "#00407A"),
        #legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.2, 'cm'),
        legend.box.margin = unit(0.1, 'cm')) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#00407A", size = 0.6)
RSE

# export the plot
setwd("./Plots/PAGE_2024")
ggsave("RSE.svg", 
       RSE, 
       dpi = 300, 
       width = 19, 
       height = 9,
       unit = "cm")
# return the working directory
Path = getwd()
setwd(dirname(dirname(Path)))
