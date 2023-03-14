library(tidyverse)
library(ggplot2)

# Theme black ----
theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "white"),  
      axis.text.x = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
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

# Same as theme black but white ------
theme_white = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_line(colour = "black"),  
      axis.text.x = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*2, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", size  =  0.2),  
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

# Plotting effect sizes -----
effect_plot <- function(data, xmin, xmax, colour, pal, lab1, lab2) {
  ggplot(data, aes(x = effect, y = dummy, xmin = {{xmin}}, xmax = {{xmax}}, colour = {{colour}})) +
    geom_point(size = 6) +
    geom_errorbar(width = .5, size = 2) +
    geom_vline(xintercept=0, color="black", linetype="dashed") +
    labs(x = "", y = " ")  +
    theme_white()+ 
    theme(legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"),
          plot.background = element_rect(colour = "black", fill=NA, size=1)) +
    xlim(c(-0.5, 2)) +
    scale_colour_manual(values = {{pal}}) +
    scale_y_continuous(breaks = 1:length({{lab1}}),
                       labels = {{lab1}},
                       sec.axis = sec_axis(~.,
                                           breaks = 1:length({{lab2}}),
                                           labels = {{lab2}}))
}