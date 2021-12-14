# Script to create functions
# Em Lim
# Updated Winter 2023

# Load packages
library(tidyverse)
library(ggplot2)

# Rename factor levels -----
rename_factors <- function(df){
  new_df <- df %>%
    mutate(
      # experimental design
      nat_num = natural,
      natural = as.factor(natural),
      nat_bin = case_when(
        natural == 1 ~ "Field: unstocked",
        natural == 2 ~ "Field: PP stocked",
        natural == 3 ~ "Field: animals stocked",
        natural == 4 ~ "Field: both stocked",
        natural == 5 ~ "Mesocosm"),
      nat_bin = factor(
        nat_bin,
        levels = c(
          "Field: unstocked",
          "Field: PP stocked",
          "Field: animals stocked",
          "Field: both stocked",
          "Mesocosm")),
      nat_bin2 = as.factor(case_when(natural == 1 ~ "Unstocked",
                           natural == 2 ~ "Unstocked",
                           natural == 3 ~ "Stocked",
                           natural == 4 ~ "Stocked",
                           natural == 5 ~ "Mesocosm")),
      nat_bin2 = factor(nat_bin2, levels = c( "Unstocked",
                                              "Stocked",
                                              "Mesocosm")),
      nat_bin_num = as.factor(case_when(natural == 1 ~ "un",
                                     natural == 2 ~ "un",
                                     natural == 3 ~ "st",
                                     natural == 4 ~ "st",
                                     natural == 5 ~ "me")),
      # consumer groups
      c_group = case_when(
        c_taxa == "Bi" ~ "Birds",
        c_taxa == "F" ~ "Fishes",
        c_taxa == "E" ~ "Sea cucumbers",
        c_taxa == "G" ~ "Molluscs",
        c_taxa == "B" ~ "Molluscs",
        c_taxa == "SP" ~ "Sponges",
        c_taxa == "SI" ~ "Jellies",# this is the 1 obelia
        c_taxa == "MI" ~ "Jellies"),# these are the jellies
      c_group = factor(
        c_group,
        levels = c(
          "Birds",
          "Fishes",
          "Sea cucumbers",
          "Molluscs",
          "Jellies",
          "Sponges")),
      c_activity = case_when(c_taxa == "Bi" | c_taxa == "F" ~ "Vertebrates",
                             c_taxa == "E" | c_taxa == "G" |  c_taxa == "MI" ~ "Mobile Inverts",
                             c_taxa == "SP" | c_taxa == "B" | c_taxa == "SI" ~ "Sessile Inverts"),
      # pp groups
      pp_taxa_old = pp_taxa,
      pp_taxa = case_when(pp_taxa == "CA" ~ "A",
                          pp_taxa == "FA" ~ "A",
                          TRUE ~ as.character(pp_taxa)),
      pp_taxa = as.factor(pp_taxa),
      pp_group = as.factor(
        case_when(
          pp_taxa == "A" ~ 'Algae',
          pp_taxa == "En" | pp_taxa == "Co" ~ 'Coral',
          pp_taxa == "P" ~ 'Phytoplankton',
          pp_taxa == "S" ~ "Plants",
          pp_phyla == "Tracheophyta" ~ "Plants",
          pp_taxa == "O" &
            pp_phyla == "Cnidaria" ~ "Coral")),# the one anemone
      pp_num = case_when(
        pp_group == "Algae" ~ 1,
        pp_group == "Coral" ~ 2,
        pp_group == "Phytoplankton" ~ 3,
        pp_group == "Plants" ~ 4),
      pp = as.factor(
        case_when(
          pp_taxa == "A"  ~ 'Algae',
          pp_taxa == "En" ~ "Endosymbiont",
          pp_taxa == "Co" ~ 'Coral',
          pp_taxa == "P" ~ 'Phytoplankton',
          pp_taxa == "S" ~ "Seagrass",
          pp_taxa == "O" & pp_phyla == "Tracheophyta" ~ "Plants",
          pp_taxa == "O" & pp_phyla == "Cnidaria" ~ "Coral")),
      # responses
      response = case_when(response == "mass" ~ "biomass",
                           response == "dry weight" ~ "biomass",
                           response == "buoyant weight" ~ "biomass",
                           response == "density" ~ "count",
                           TRUE ~ as.character(response)),
      res = as.factor(response),
      measured = ifelse(measured == "blade" & pp_taxa == "S", "leaf", measured),
      grow_time_end = case_when(growth_rate == "Y" ~ "growth_rate_data",
                                growth_rate == "B" ~ "growth_rate_data",
                                growth_inst == "timepoint" ~ "timepoint",
                                growth_inst == "grow" ~ "endpoint"),
      # provisioned nutrients
      prov = case_when(treatment_levels == "biodeposits" ~ "poop",
                          treatment_levels == "faeces" ~ "poop",
                          treatment_levels == "pseudofeces" ~ "pseudofeces",
                          treatment_levels == "effluent" ~ "pee",
                          treatment_levels == "excretion" ~ "pee",
                          treatment_levels == "animals" ~ "animals"),
      # nutrients
      # redo background nut categories with just papers where they described nut levels
      bn_cat = factor(
        case_when(nutrient_levels == "deplete" ~ "Low nut",
                  nutrient_levels == "extremely low" ~ "Low nut",
                  nutrient_levels == "low" ~ "Low nut",
                  nutrient_levels == "nutrient-poor" ~ "Low nut",
                  nutrient_levels == "oligotrophic" ~ "Low nut",
                  nutrient_levels == "nutrient-limited" ~ "Low nut",
                  nutrient_levels == "eutrophic" ~ "High nut",
                  nutrient_levels == "high" ~ "High nut",
                  nutrient_levels == "nutrient-replete" ~ "High nut",
                  TRUE ~ as.character("NA")),
        levels = c("High nut", "Low nut")
      ),
      bn_cat_old = factor(
        ifelse(background_nut_umolL > 4, "High nut", "Low nut"),
        levels = c("High nut", "Low nut")), # paper 46 had been high but they don't actually say it's eutrophic
      # latitude/climate
      abs_lat = abs(lat_dec), # creating an absolute value lat
      climate = factor(
        case_when(
          abs_lat <= 23.5 ~ "Tropical",
          abs_lat <= 35 ~ "Subtropical",
          TRUE ~ as.character("Temperate")),
        # abs_lat <= 66.6 ~ "Temperate",
        #  abs_lat > 66.6 ~ "Polar"),
        levels = c("Tropical", "Subtropical", "Temperate")),
      tt = as.factor(ifelse(abs_lat >= 30, "Temperate", "Tropical")),
      tt2 = ifelse(abs_lat <= 23.5, "Tropical",
                   ifelse(abs_lat >= 30 & abs_lat < 66.6, "Temperate", NA)),
      # temperature
      temp_cat = ifelse(temp >22, "hot", "cold"),
      sst_cat = ifelse(temp >20, "Hot", "Cold")
    )
}

# scale variables
scale_vars <- function(df) {
  new_df <- df %>%
    mutate(lat_stand = c(scale(abs_lat)),
           lat_center = c(scale(abs_lat, scale = FALSE)),
           bn_stand = c(scale(background_nut_umolL)),
           bn_center = c(scale(background_nut_umolL, scale = FALSE)),
           po4_stand = c(scale(po4)),
           po4_stand = c(scale(po4)),
           po4_center = c(scale(po4, scale = FALSE)),
           temp_stand = c(scale(temp)),
           temp_center = c(scale(temp, scale = FALSE)),
           sst_stand = c(scale(sst)),
           sst_center = c(scale(sst, scale = FALSE)),
           duration_stand = c(scale(duration_days)),
           duration_center = c(scale(duration_days, scale = FALSE)),
           duration_log = log(duration_days + 1),
           added_N_stand = c(scale(added_N)),
           added_N_center = c(scale(added_N, scale = FALSE)),
        #   added_N_log = log(added_N + 1),
           indiv = as.factor(indiv),
           paper_num = as.factor(paper_num))
}

# Plotting themes ---------------
## Theme black ----
theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "white"),  
      axis.text.x = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*2.5, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*2.5, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "black",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.5, color = "white"),  
      legend.title = element_text(size = base_size*2.5, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "black"),  
      panel.grid.minor = element_line(color = "black"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*2.5, color = "white"),  
      strip.text.y = element_text(size = base_size*2.5, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

## Same as theme black but white ------
theme_white = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", linewidth  =  0.2),  
      axis.title.x = element_text(size = base_size*2.5, color = "black", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*2.5, color = "black", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*1.5, color = "black"),  
      legend.title = element_text(size = base_size*2.5, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_line(color = "white"),  
      panel.grid.minor = element_line(color = "white"),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey80", color = "grey30"),  
      strip.text.x = element_text(size = base_size*2.5, color = "black"),  
      strip.text.y = element_text(size = base_size*2.5, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "white", fill = "white"),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
}

## Theme for pubs -----
pub_theme <- function() {
  theme_classic() +
    theme(
          # axis
          axis.text = element_text(size = 10, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.line = element_line(colour = "black", linewidth = 0.1),
          axis.ticks = element_line(color = "black", linewidth = 0.2),
          axis.ticks.length = unit(0.1, "lines"),
          # legend
          legend.title = element_text(size = 11, colour = "black"),
          legend.text = element_text(size = 11, colour = "black"),
          # panel border
          panel.background = element_rect(fill = "white", color  =  NA),  
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.spacing = unit(0.1, "lines"),
          # facet
          strip.text.x = element_text(size = 10, color = "black")
          
    )
}


# Pull slope coeff from model ------------------------------------------------
# Slope
slope_coeff_fun <- function(model, slope) {
  int_effect_sum <- coef(summary(model)) %>%
    slice(-(1)) %>%
    slice(-(2:5)) %>%
    mutate(int = factor(c("Algae", "Coral", "Endosymbiont", "Phytoplankton", "Seagrass")),
           effect = ifelse(int != "Algae", estimate + {{slope}},
                           {{slope}}), # add intercept to the estimated coefficient
           ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
           ci.up_adjust = effect + (1.96*se),
           dummy = 1:5)
}

slope_coeff_fun_add <- function(model, slope) {
  int_effect_sum <- coef(summary(model)) %>%
    slice(-(1)) %>%
    slice(-(2:6)) %>%
    mutate(int = factor(c("Algae", "Coral", "Endosymbiont", "Phytoplankton", "Seagrass")),
           effect = ifelse(int != "Algae", estimate + {{slope}},
                           {{slope}}), # add intercept to the estimated coefficient
           ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
           ci.up_adjust = effect + (1.96*se),
           dummy = 1:5)
}

slope_coeff_fun_sge <- function(model, slope) {
  int_effect_sum <- coef(summary(model)) %>%
    slice(-(1)) %>%
    slice(-(2)) %>%
    mutate(int = factor(c("Epiphytes", "Seagrass")),
           effect = ifelse(int != "Epiphytes", estimate + {{slope}},
                           {{slope}}), # add intercept to the estimated coefficient
           ci.lb_adjust = effect - (1.96*se), #adjust the CI based on the new estimate
           ci.up_adjust = effect + (1.96*se),
           dummy = 1:2)
}



# Generate new df ------------------------------------------------------------
new_df <- function(variable, df) {
  # figure out the background nutrient ranges for each pp group
  mm <- df %>%
    group_by(pp) %>%
    summarise(min = min({{variable}}),
              max = max({{variable}})) %>%
    mutate(x = c(1,2,3,4,5)) %>%
    as.data.frame()
  
  # Generate the values
  nut_values <- data.frame(
    var = c(
      seq(from = mm[1,2], to = mm[1,3], length.out = 100),
      seq(from = mm[2,2], to = mm[2,3], length.out = 100),
      seq(from = mm[3,2], to = mm[3,3], length.out = 100),
      seq(from = mm[4,2], to = mm[4,3], length.out = 100),
      seq(from = mm[5,2], to = mm[5,3], length.out = 100)),
    ppCoral = c(rep(0, each = 100), 
                rep(1, each = 100), 
                rep(0, each = 100), 
                rep(0, each = 100),
                rep(0, each = 100)),
    ppEndosymbiont = c(rep(0, each = 100), 
                       rep(0, each = 100), 
                       rep(1, each = 100), 
                       rep(0, each = 100),
                       rep(0, each = 100)),
    ppPhytoplankton = c(rep(0, each = 100), 
                        rep(0, each = 100),
                        rep(0, each = 100), 
                        rep(1, each = 100),
                        rep(0, each = 100)),
    ppSeagrass = c(rep(0, each = 100), 
                   rep(0, each = 100), 
                   rep(0, each = 100), 
                   rep(0, each = 100),
                   rep(1, each = 100)))
}

# another for sge
new_df_sge <- function(variable, df) {
  # figure out the background nutrient ranges for each pp group
  mm <- df %>%
    group_by(pp) %>%
    summarise(min = min({{variable}}),
              max = max({{variable}})) %>%
    mutate(x = c(1,2)) %>%
    as.data.frame()
  
  # Generate the values
  nut_values <- data.frame(
    var = c(
      seq(from = mm[1,2], to = mm[1,3], length.out = 100),
      seq(from = mm[2,2], to = mm[2,3], length.out = 100)),
    ppSeagrass = c(rep(0, each = 100), 
                   rep(1, each = 100)))
}


# Prep for emmeans -----
emmeans_prep <- function(model){
  rownames(model$beta)[1] = "(Intercept)"
  prep = qdrg(object = model) 
}


# Plotting functions -----------------------------------------------------------
## Cat effect plot ---------------------------------------------------------------
effect_plot <- function(sum_data, y_var, facet) {
  
  pal5 <- viridis::viridis(5)
  pal <- pal5[1]
  
  plot <- ggplot(sum_data, aes(emmean, {{y_var}}, colour = {{y_var}})) +
    geom_point(colour = pal) +
    geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, linewidth = 0.5, colour = pal) +
    pub_theme() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, lty = "dashed") +
    ylab("") +
    xlab("Effect size (LnRR)") +
    #scale_color_manual(values = (pal)) +
    # Add text labels inside the plots
    geom_text(aes(x = emmean, label = n), # change to 1 if i add tt
              color = "black", size = 3, vjust = -0.5, hjust = 0)+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

if(facet == 2) {
  plot <- plot + facet_wrap(~ tt + pp, ncol = 5)
} else if(facet == 1) {
  plot <- plot + facet_wrap(~ pp, ncol = 5)
}

return(plot)
}
  


## Slope only plot -------------------------------------------------------------
int_effect_plot <- function(data, theme = "white", pal){
  
  if(theme == "white"){
    theme_var <- pub_theme()
    features <- "black"
    size_var = 2.2
    line_var = 0.5
    hline_var = 0.2
  } else if(theme == "black") {
    theme_var <- theme_black()
    features <- "white"
    size_var = 5
    line_var = 1
    hline_var = 0.5
  }
  
  lab1 <- c("Algae", "Coral", "Endo", "Phyto", "Seagrass")
  ggplot(data, aes(x = effect, y = int, xmin = ci.lb_adjust, xmax = ci.up_adjust, colour = int)) +
    geom_point(size = size_var) +
    geom_errorbar(width = 0, linewidth = line_var) +
    geom_vline(xintercept=0, color = features, linetype="dashed", linewidth = hline_var) +
    labs(x = "", y = " ") +
    theme_var +
    theme(plot.margin = unit(c(0.1,0.15,0,0), "cm"),
          legend.position = "none") +
    scale_y_discrete(labels = lab1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    scale_colour_manual(values = pal)+
    xlab("Slope")
}
# plot.margin = top, right, bottom, left

int_effect_plot_sge <- function(data, pal) {
  
  lab1 <- c("Epiphytes", "Seagrass")
  ggplot(data, aes(x = effect, y = int, xmin = ci.lb_adjust, xmax = ci.up_adjust, colour = int)) +
    geom_point(size = 2.7) +
    geom_errorbar(width = 0, linewidth = 0.5) +
    geom_vline(xintercept=0, color="black", linetype="dashed") +
    labs(x = "", y = " ") +
    pub_theme() +
    theme(plot.margin = unit(c(0.1,0.15,0,0), "cm"),
          legend.position = "none") +
    scale_y_discrete(labels = lab1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    scale_colour_manual(values = rev(pal)) 
}

## Scatter plot function ------------------------------------------------------
scatter_fun <- function(data_sum, data_raw, x, pal, theme = "white"){
  
  if(theme == "white"){
    theme_var <- pub_theme()
    features <- "black"
    size_var = 0.5
    line_var = 0.5
    hline_var = 0.2
    alpha_point = 0.5
    ribbon_alpha = 0.1
  } else if(theme == "black") {
    theme_var <- theme_black()
    features <- "white"
    size_var = 3.5
    line_var = 1.5
    hline_var = 0.7
    alpha_point = 0.8
    ribbon_alpha = 0.2
  }
  
  ggplot() +
    geom_jitter(data = data_raw, 
                aes({{x}}, yi, colour = pp, fill = pp), 
                size = size_var, alpha = alpha_point) +
    geom_line(data = data_sum, aes({{x}}, pred, colour = pp),  linewidth = line_var) +
    geom_ribbon(data = data_sum, aes({{x}}, pred, ymin = ci.lb, ymax = ci.ub, fill = pp), 
                alpha = ribbon_alpha) +
    labs(y = "Effect (LnRR)", colour = "Primary Producer", fill = "Primary Producer") +
    geom_hline(yintercept = 0, color = features, linetype = "dashed", linewidth = hline_var) +
    scale_colour_manual(values = {{pal}},
                        labels=c("Algae" = "Algae",
                                 "Coral" = "Coral", 
                                 "Endosymbiont" = "Endo",
                                 "Phytoplankton" = "Phyto",
                                 "Seagrass" = "Seagrass"),
                        drop = FALSE) +
    scale_fill_manual(values = {{pal}},
                      labels=c("Algae" = "Algae",
                               "Coral" = "Coral", 
                               "Endosymbiont" = "Endo",
                               "Phytoplankton" = "Phyto",
                               "Seagrass" = "Seagrass"),
                      drop = FALSE) +
    theme_var +
    theme(legend.position = "none",
          plot.margin = unit(c(0.1,0,0,0.15), "cm")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
    labs(colour = "", fill = "") 
}



# Save model summary as csv ---------------------------------------------
save_mod_output <- function(model){
  mod_table <-coef(summary(model)) %>%
    rownames_to_column() %>%
    rowwise() %>%
    mutate(rowname = ifelse(rowname == "intrcpt", "Intercept", rowname),
           estimate = ifelse(abs(estimate) < 0.005, sprintf("%.3f",estimate), sprintf("%.2f",estimate)),
           se = ifelse(abs(se) < 0.005, sprintf("%.3f",se), sprintf("%.2f",se)),
           tval = ifelse(abs(tval) < 0.005, sprintf("%.3f",tval), sprintf("%.2f",tval)),
           pval = ifelse(pval < 0.001, sprintf("< 0.001"), 
                            ifelse(pval < 0.005, sprintf("%.3f",pval), sprintf("%.2f",pval)))) %>%
    rename(` ` = rowname,
           SE = se,
           `t-value` = tval,
           `p-value` = pval) %>%
    select(-c(ci.lb, ci.ub))
}

# Save emtrends output as csv ---------
save_slopes_fun <- function(rg, var){
  
  emt <- stats::update(emmeans::contrast(rg, "consec", simple = var), by = "pp")
  
  slope_output <- summary(emt, infer=TRUE) %>%
    as_tibble() %>%
    transmute(pp = pp,
              estimate = ifelse(abs(estimate) < 0.005, sprintf("%.3f",estimate), sprintf("%.2f",estimate)),
              SE = ifelse(abs(SE) < 0.005, sprintf("%.3f",SE), sprintf("%.2f",SE)),
              df = df,
              `t-ratio` = ifelse(abs(t.ratio) < 0.005, sprintf("%.3f",t.ratio), sprintf("%.2f",t.ratio)),
              `p-value` = ifelse(p.value < 0.001, sprintf("< 0.001"), 
                                 ifelse(p.value < 0.005, sprintf("%.3f",p.value), sprintf("%.2f",p.value))))
  
  print(slope_output)
}

# save glm summary as csv  ---------------------------------------------
save_glm_output <- function(mod, filepath){
  
  df <- summary(mod)$coefficients$cond
  
  mod_table <- df %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(estimate = Estimate,
           se = `Std. Error`,
           zval = `z value`,
           pval = `Pr(>|z|)`) %>%
    rowwise() %>%
    mutate(rowname = ifelse(rowname == "(Intercept)", "Intercept", rowname),
           estimate = ifelse(abs(estimate) < 0.005, sprintf("%.3f",estimate), sprintf("%.2f",estimate)),
           se = ifelse(abs(se) < 0.005, sprintf("%.3f",se), sprintf("%.2f",se)),
           zval = ifelse(abs(zval) < 0.005, sprintf("%.3f",zval), sprintf("%.2f",zval)),
           pval = ifelse(pval < 0.001, sprintf("< 0.001"), 
                         ifelse(pval < 0.005, sprintf("%.3f",pval), sprintf("%.2f",pval)))) %>%
    rename(` ` = rowname,
           SE = se,
           `z-value` = zval,
           `p-value` = pval) 
  write_csv(mod_table, filepath)
  
  sprintf("saved to %s", filepath)
}

# save pairwise output --------------------------------------------------
tidy_pairwise_output <- function(df){
table <- df %>%
    mutate(estimate = ifelse(abs(estimate) < 0.005, sprintf("%.3f",estimate), sprintf("%.2f",estimate)),
           SE = ifelse(abs(SE) < 0.005, sprintf("%.3f",SE), sprintf("%.2f",SE)),
           t.ratio = ifelse(abs(t.ratio) < 0.005, sprintf("%.3f",t.ratio), sprintf("%.2f",t.ratio)),
           p.value = ifelse(p.value < 0.001, sprintf("< 0.001"), 
                            ifelse(p.value < 0.005, sprintf("%.3f",p.value), sprintf("%.2f",p.value)))) %>%
    rename(`t-ratio` = t.ratio,
           `p-value` = p.value) %>%
  mutate(pp = case_when(pp == "Phytoplankton" ~ "Phyto",
                        pp == "Endosymbiont" ~ "Endo",
                        TRUE ~ as.factor(pp)))
}

# save pairwise glm output --------------------------------------------------
tidy_glm_pairwise_output <- function(df){
  table <- df %>%
    mutate(estimate = ifelse(abs(estimate) < 0.005, sprintf("%.3f",estimate), sprintf("%.2f",estimate)),
           SE = ifelse(abs(SE) < 0.005, sprintf("%.3f",SE), sprintf("%.2f",SE)),
           z.ratio = ifelse(abs(z.ratio) < 0.005, sprintf("%.3f",z.ratio), sprintf("%.2f",z.ratio)),
           p.value = ifelse(p.value < 0.001, sprintf("< 0.001"), 
                            ifelse(p.value < 0.005, sprintf("%.3f",p.value), sprintf("%.2f",p.value)))) %>%
    rename(`z-ratio` = z.ratio,
           `p-value` = p.value) %>%
    mutate(pp = case_when(pp == "Phytoplankton" ~ "Phyto",
                          pp == "Endosymbiont" ~ "Endo",
                          TRUE ~ as.factor(pp)))
}


# Calculate I2 ------

calc_I2 <- function(mod){
  W <- diag(1/mod$vi)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  100 * sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
}


# AIC comparisson function -----------------------------------------------------
meta_aic <- function(moderator, covar, df){
  
  #print(moderator)
  #print(covar)
  #print(df)
  
  if(moderator == "bn_stand"){
    df_drop <- df %>% 
      drop_na(!!sym(covar)) %>%
      filter(bn_proxy == "measure")
  }
  else if (moderator == "bn_cat") {
    df_drop <- df %>% 
      drop_na(!!sym(covar)) %>% 
      drop_na(bn_cat)
  } else{
    df_drop <- df %>% drop_na(!!sym(covar)) # drop NAs
  }
  
  # create formula 1 for reduced mod
  formula_str1 <- paste("~", moderator, "* pp")
  formula1 <- as.formula(formula_str1)
  
  # create formula 2 for mod with covar
  formula_str2 <- paste("~", moderator, "* pp +", covar)
  formula2 <- as.formula(formula_str2)
  
  mod <- rma.mv(yi = yi,
                V = vi,
                slab = author,
                data = df_drop,
                method = "REML",
                mods = formula1,
                random = ~ 1 | paper_num/indiv, 
                test = "t")
  
  mod2 <- rma.mv(yi = yi,
                 V = vi,
                 slab = author,
                 data = df_drop,
                 method = "REML",
                 mods = formula2,
                 random = ~ 1 | paper_num/indiv, 
                 test = "t")
  
  AIC_tab <- AIC(mod, mod2) %>%
    mutate(model = c(deparse(formula1), deparse(formula2)),
           delta = AIC - min(AIC),
           likelihood = exp( -0.5*delta),
           sum = sum(likelihood),
           AICw = likelihood/sum) %>%
    transmute(model = model, 
              K = df,
              AIC = round(AIC, 2),
              delta = round(delta, 2),
              `AIC Weight` = round(AICw, 2)) %>%
    remove_rownames()
}


# SST functions -------

# original function
download_sst <- function(row) {
  message("Downloading request_id: ", row$request_id,
          " | Lat: ", row$lat_dec,
          " | Lon: ", row$long_dec)
  
  dat <- tryCatch({
    rerddap::griddap(
      datasetx = "ncdcOisst21Agg_LonPM180",
      url = "https://coastwatch.pfeg.noaa.gov/erddap/",
      time = c(format(row$start_date, "%Y-%m-%d"),
               format(row$end_date, "%Y-%m-%d")),
      zlev = c(0, 0),
      latitude = c(row$lat_dec, row$lat_dec),
      longitude = c(row$long_dec, row$long_dec),
      fields = "sst"
    )$data
  }, error = function(e) {
    message("Error for request_id ", row$request_id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  
  dat %>%
    mutate(
      t = as.Date(str_remove(time, "T12:00:00Z")),
      lon = longitude,
      lat = latitude,
      lat_dec = row$lat_dec,
      long_dec = row$long_dec,
      start_date = row$start_date,
      end_date = row$end_date,
      request_id = row$request_id,
      temp = sst
    ) %>%
    select(request_id, lat_dec, long_dec, start_date, end_date, lon, lat, t, temp) %>%
    na.omit()
}

# function for nearest point

# another function
download_sst_nearest <- function(row, grid_size = 0.25, buffer = 0.25) {
  
  # Snap coordinates to nearest OISST 0.25° grid
  lat_snap <- round(row$lat_dec / grid_size) * grid_size
  lon_snap <- round(row$long_dec / grid_size) * grid_size
  
  message("Downloading request_id: ", row$request_id,
          " | Lat: ", lat_snap,
          " | Lon: ", lon_snap)
  
  # First attempt: single grid point
  dat <- tryCatch({
    rerddap::griddap(
      datasetx = "ncdcOisst21Agg_LonPM180",
      url = "https://coastwatch.pfeg.noaa.gov/erddap/",
      time = c(format(row$start_date, "%Y-%m-%d"),
               format(row$end_date, "%Y-%m-%d")),
      zlev = c(0, 0),
      latitude = c(lat_snap, lat_snap),
      longitude = c(lon_snap, lon_snap),
      fields = "sst"
    )$data
  }, error = function(e) NULL)
  
  # If no data, expand to small neighborhood
  if (is.null(dat) || nrow(dat) == 0) {
    message("No data at exact point, expanding to small neighborhood")
    
    lat_min <- lat_snap - buffer
    lat_max <- lat_snap + buffer
    lon_min <- lon_snap - buffer
    lon_max <- lon_snap + buffer
    
    dat <- tryCatch({
      rerddap::griddap(
        datasetx = "ncdcOisst21Agg_LonPM180",
        url = "https://coastwatch.pfeg.noaa.gov/erddap/",
        time = c(format(row$start_date, "%Y-%m-%d"),
                 format(row$end_date, "%Y-%m-%d")),
        zlev = c(0, 0),
        latitude = c(lat_min, lat_max),
        longitude = c(lon_min, lon_max),
        fields = "sst"
      )$data
    }, error = function(e) NULL)
    
    if (is.null(dat) || nrow(dat) == 0) {
      message("No valid data found even after expansion")
      return(NULL)
    }
    
    # Pick the nearest valid SST to original lat/lon
    dat <- dat %>%
      filter(!is.na(sst)) %>%
      mutate(dist = sqrt((latitude - row$lat_dec)^2 + (longitude - row$long_dec)^2)) %>%
      slice_min(dist, n = 1)
  }
  
  # Clean and return
  dat %>%
    mutate(
      t = as.Date(str_remove(time, "T12:00:00Z")),
      lon = longitude,
      lat = latitude,
      lat_dec = row$lat_dec,
      long_dec = row$long_dec,
      start_date = row$start_date,
      end_date = row$end_date,
      request_id = row$request_id,
      temp = sst
    ) %>%
    select(request_id, lat_dec, long_dec, start_date, end_date, lon, lat, t, temp) %>%
    na.omit()
}


# Build an ellipse in Mollweide coordinates -----------------------------------
make_ellipse <- function(center = c(0, 0), a = 18000000, b = 8900000, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  x <- center[1] + a * cos(theta)
  y <- center[2] + b * sin(theta)
  # Ensure the polygon is closed
  ellipse_coords <- cbind(c(x, x[1]), c(y, y[1]))  # Close the loop
  st_polygon(list(ellipse_coords)) %>% st_sfc(crs = "+proj=moll")
}



# Super slick plot fun -------
cont_plot_god_fun <- function(data_frame, model, variable){
  slope <- coef(summary({{model}}))[2,1] # extract the slope
  
  # Calculate slope for each pp
  int_effect_sum <- slope_coeff_fun_sge({{model}}, slope)
  
  # Predictions for scatter plot
  # Generate values
  values <- new_df_sge({{variable}}, {{data_frame}}) %>% # Function to create df with dummy vars
    mutate("var:ppSeagrass" = var*ppSeagrass) %>%
    as.matrix()
  
  # new model with non standardized background nutrients
  model_lat_untran_sge <- rma.mv(yi = yi,
                                 V = vi,
                                 slab = author,
                                 data = meta_sge,
                                 method = "REML",
                                 # will need to specify formula here
                                 mods = ~ {{variable}} + pp + {{variable}}:pp,
                                 random = ~ 1 | paper_num/indiv, 
                                 test = "t")
  
  # use the generated values to plot background nutrients effect for each pp_group, nat_group held to 1
  predict <- as.data.frame(predict(object = model_lat_untran_sge,
                                           newmods = values)) %>%
    cbind(values) %>%
    mutate(pp = c(rep("Algae", each = 100), 
                  rep("Seagrass", each = 100)))
  
  # slope plot
  slope_plot <- int_effect_plot_sge(int_effect_sum) +
    scale_colour_manual(values = rev(pal5)) +
    xlab("Slope")
  
  # scatter plot
  scatter_plot <- scatter_fun(predict, {{data_frame}}, {{variable}}, rev(pal5)) + 
    labs(x = "Absolute latitude") +
    theme(axis.title.y = element_text(vjust = -15))
  
  p <- slope_plot + scatter_plot
}
  