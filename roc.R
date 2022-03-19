rm(list=ls())

library(tidyverse)
library(forcats)
library(yardstick)
library(broom)
library(plotly)

data_roc <- as_tibble(read.csv("data for roc.csv"))

CURB_fac <- c("A", "B", "C")
death_fac <- c()
data_roc <- data_roc %>%
  dplyr::mutate(CURB_grade = case_when(
    CURB_65 %in% c(0, 1) ~ "A",
    CURB_65 == 2 ~ "B",
    CURB_65 %in% c(3, 4, 5) ~ "C")
    ) %>%
  dplyr::mutate(CURB_factor = factor(CURB_grade, levels = CURB_fac)) %>%
  dplyr::select(-c(CURB_grade))

##plot by outcome
p <- ggplot(data = data_roc, aes(x = as.factor(death), y = APACHEII))
p + geom_boxplot()

glm_out1 <- glm(
  formula = death ~ CURB_factor,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "CURB_65")

glm_out2 <- glm(
  formula = death ~ APACHEII,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "APACHEII")

glm_out3 <- glm(
  formula = death ~ PCT,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "PCT")

glm_out4 <- glm(
  formula = death ~ PCT + CURB_factor,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "PCT + CURB_65")

glm_out5 <- glm(
  formula = death ~ PCT + APACHEII,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "PCT + APACHEII")

glm_out6 <- glm(
  formula = death ~ CURB_factor + APACHEII,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "CURB_65 + APACHEII")

glm_out7 <- glm(
  formula = death ~ PCT + CURB_factor + APACHEII,
  family = binomial,
  data = data_roc
) %>%
  augment() %>%
  dplyr::mutate(model = "PCT + CURB_65 + APACHEII")

glm_out_raw <- bind_rows(glm_out1, glm_out2, glm_out3,
                         glm_out4, glm_out5, glm_out6,
                         glm_out7) %>%
  dplyr::mutate(death_factor = as.factor(death))

glm_out <- bind_rows(glm_out1, glm_out2, glm_out3,
                     glm_out4, glm_out5, glm_out6,
                     glm_out7) %>%
  dplyr::mutate(death_factor = as.factor(death)) %>%
  dplyr::group_by(model) %>%
  roc_curve(event_level = 'second', truth = death_factor, .fitted) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(specificity_comp = 1 - specificity) %>%
  dplyr::arrange(model, specificity_comp, sensitivity)

## differentiate groups by linetype
p.out <- ggplot(
  aes(x = 1 - specificity, y = sensitivity, linetype = model),
  data = glm_out) + # plot with 2 ROC curves for each model
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4, alpha = 0.5)

## differentiate groups by point shape
p.out <- ggplot(
  aes(x = 1 - specificity, y = sensitivity, group = model),
  data = glm_out) +
  geom_line(size = 1) +
  geom_point(aes(shape = model), size = 2)
ggplotly(p.out)
# test.df <- glm_out %>%
#   dplyr::group_by(model) %>% # group to get individual ROC curve for each model
#   roc_curve(event_level = 'second', truth = death_factor, .fitted) %>%
#   dplyr::filter(model == "PCT") %>%
#   dplyr::mutate(specificity_comp = 1 - specificity) %>%
#   dplyr::arrange()
# 
# test.df_sorted <- test.df %>%
#   dplyr::arrange(specificity_comp, sensitivity)
# p1 <- ggplot(aes(x = specificity_comp,
#                  y = sensitivity,
#                  color = model),
#              data = test.df_sorted) + # plot with 2 ROC curves for each model
#   geom_line(size = 1.1) +
#   geom_abline(slope = 1, intercept = 0, size = 0.4)
# ggplotly(p.ggplot)

# plot ROC curves
glm_out %>%
  group_by(model) %>% # group to get individual ROC curve for each model
  roc_curve(event_level = 'second', truth = death_factor, .fitted) %>% # get values to plot an ROC curve
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) #+
  #scale_color_manual(values = c("#48466D", "#3D84A8")) +
  #coord_fixed()



## Compute AUC
# function to obtain variance of AUC
v_AUC <- function(AUC, n_pos, n_neg){
  Q1 <- AUC/(2-AUC)
  Q2 <- (2*AUC^2)/(1+AUC)
  out <- (1/(n_pos*n_neg))*(AUC*(1-AUC) + (n_pos-1)*(Q1-AUC^2) + (n_neg-1)*(Q2-AUC^2))
  return(out)
}
n_d <- data_roc %>% filter(death == 1) %>% nrow()
n_nd <- data_roc %>% filter(death == 0) %>% nrow()
# alpha for CI
CI.alpha <- 0.05
res_AUC <- glm_out_raw %>%
  group_by(model) %>% # group to get individual AUC value for each model
  roc_auc(event_level = 'second', truth = death_factor, .fitted) %>%
  mutate(Var.AUC = v_AUC(AUC = .estimate, n_pos = n_d, n_neg = n_nd)) %>%
  mutate(CI.lower = .estimate + qnorm(p = CI.alpha/2, mean = 0, sd = 1, lower.tail = TRUE)*sqrt(Var.AUC),
         CI.upper = .estimate + qnorm(p = 1 - CI.alpha/2, mean = 0, sd = 1, lower.tail = TRUE)*sqrt(Var.AUC)) %>%
  mutate(CI.upper = ifelse(CI.upper > 1, 1, CI.upper)) %>%
  rename(AUC = .estimate) %>%
  select(-c(.metric, .estimator))
write.csv(res_AUC, "res_AUC.csv", row.names = FALSE)


## Comparison between death group vs nondeath group
res_count <- data_roc %>%
  group_by(death) %>%
  summarise(count = n())

res_summary <- data_roc %>%
  group_by(death) %>%
  select(-subject_id) %>%
  summarise_all(list(mean = mean, sd = sd))

res_compare <- left_join(res_count, res_summary, by = "death")

res_p.value <- data_roc %>%
  summarise_each(funs(t.test(.[death == 1], .[death == 0])$p.value),
                 vars = CURB_65:age)
colnames(res_p.value) <- colnames(data_roc)[-c(1, 2)]

write.csv(res_compare, "res_compare.csv", row.names = FALSE)
write.csv(res_p.value, "res_p.value.csv", row.names = FALSE)

