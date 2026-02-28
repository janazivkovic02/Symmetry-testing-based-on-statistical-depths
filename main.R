library(MASS)
library(ddalpha)
library(mvtnorm)
library(sn)
library(ggplot2)
library(dplyr)
library(tidyr)
source("R/data_generation.R")
source("R/data_generation2.R")
source("R/rn_test_statistic.R")

#### simulacije figure 1 ####

run_experiment <- function(settings, j_vals, sampler,
                           M = 3000, n = 100, alpha = 0.05,
                           setting_labels, panel_order = settings,
                           ylab = "Azzalini skewing") {
  
  monte_carlo_runs <- function(setting, j, M, n, alpha, sampler) {
    reject <- matrix(FALSE, nrow = M, ncol = 5)
    colnames(reject) <- c("Simplicial","L2","Mahalanobis","Projection","Zonoid")
    
    crit <- qnorm(1 - alpha/2) * sqrt(11/3)
    
    for (m in 1:M) {
      x <- sampler(setting, j, n)
      
      r <- c(
        rn_sim(x),
        rn_l2(x),
        rn_mah(x),
        rn_proj(x),
        rn_zonoid(x)
      )
      
      stat <- (4*r - n - 2) / sqrt(n)
      reject[m, ] <- abs(stat) > crit
    }
    
    colMeans(reject)
  }
  
  rezultati <- list()
  for (s in settings) {
    for (j in j_vals) {
      cat("Setting:", s, "j =", j, "\n")
      
      res <- monte_carlo_runs(s, j, M, n, alpha, sampler)
      
      rezultati[[paste(s, "j", j, sep = "_")]] <-
        data.frame(
          Setting = s,
          j = j,
          Test = names(res),
          Rejection = as.numeric(res)
        )
    }
  }
  
  final_table <- do.call(rbind, rezultati)
  
  final_table$Setting <- factor(final_table$Setting, levels = panel_order)
  final_table$SettingLabel <- factor(setting_labels[as.character(final_table$Setting)],
                                     levels = setting_labels[panel_order])
  final_table$Test <- factor(final_table$Test,
                             levels = c("Simplicial","L2","Mahalanobis","Projection","Zonoid"))
  
  p <- ggplot(final_table, aes(x = j, y = Rejection, linetype = Test, color = Test, group = Test)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~SettingLabel, ncol = 2) +
    scale_x_continuous(breaks = j_vals) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    labs(x = "j", y = ylab) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10)
    )
  
  list(final_table = final_table, plot = p)
}

#### figura 1 ####

settings1 <- c(
  "normal_spherical",
  "normal_elliptical",
  "normal_cone",
  "cauchy_spherical",
  "cauchy_elliptical",
  "cauchy_cone"
)

j_vals <- 0:3

setting_labels1 <- c(
  normal_spherical   = "Sferno normalno jezgro",
  cauchy_spherical   = "Sferno Košijevo jezgro",
  normal_elliptical  = "Eliptično normalno jezgro",
  cauchy_elliptical  = "Eliptično Košijevo jezgro",
  normal_cone        = "Centralno simetrično normalno jezgro",
  cauchy_cone        = "Centralno simetrično Košijevo jezgro"
)

panel_order1 <- c(
  "normal_spherical", "cauchy_spherical",
  "normal_elliptical", "cauchy_elliptical",
  "normal_cone", "cauchy_cone"
)

out1 <- run_experiment(
  settings = settings1,
  j_vals = j_vals,
  sampler = generate_sample,   # <--- KLJUČNO
  M = 3000,
  n = 100,
  setting_labels = setting_labels1,
  panel_order = panel_order1,
  ylab = "Azzalini skewing"
)

print(out1$plot)

final_table1 <- out1$final_table

#### figura 2 ####


settings2 <- c("normal_spherical", "cauchy_spherical")
j_vals <- 0:3

setting_labels2 <- c(
  normal_spherical = "Sferno normalno jezgro",
  cauchy_spherical = "Sferno Košijevo jezgro"
)

out2 <- run_experiment(
  settings = settings2,
  j_vals = j_vals,
  sampler = generate_sample_contamination,
  M = 3000, n = 100,
  setting_labels = setting_labels2,
  panel_order = c("normal_spherical", "cauchy_spherical"),
  ylab = "Azzalini skewing sa kontaminacijom"
)

print(out2$plot)
final_table2 <- out2$final_table