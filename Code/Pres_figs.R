## Presentation figures ------
pal <- viridis::viridis(12)
pie(rep(1,12), col = pal)
#pal_b <- c(pal[12], pal[11], pal[10], pal[8], pal[6])
pal_b <- c(pal[12], pal[10], pal[8], pal[7], pal[5])

# slope plot latitude
slope_l_b <- int_effect_plot(int_effect_sum_l, theme = "black", pal_b) +
  facet_wrap(~dummy)

# scatter plot latitude
scatter_l_b <- scatter_fun(lat_predict, 
                           meta_df, 
                           abs_lat, (pal_b), theme = "black") + 
  labs(x = "Absolute latitude") +
  ylim(-1.4, 3.9) +
  xlim(5.8, 74.3)
#  theme(axis.title.y = element_text(vjust = -15))

# slope plot nutrients
slope_n_b <- int_effect_plot(int_effect_sum, theme = "black", pal_b) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  facet_wrap(~dummy)

# scatter plot nutrients
scatter_n_b <- scatter_fun(nut_predict, meta_df %>% filter(bn_proxy == "measure"),
                           background_nut_umolL, (pal_b), theme = "black") + 
  labs(x = expression("Background N"~(mu*M))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  ylim(-1.4, 3.9) +
  xlim(0.03, 17.65)

# slope plot sst
slope_sst_b <- int_effect_plot(int_effect_sum_sst, theme = "black", pal_b) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  facet_wrap(~dummy)

# scatter plot sst
scatter_sst_b <- scatter_fun(sst_predict, meta_df, sst, (pal_b), theme = "black") + 
  labs(x = expression("Annual SST " ( degree*C))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position = "right") +
  ylim(-1.4, 3.9) +
  xlim(-0.8, 28.8)

### Fig.3 plot black background
(slope_l_b + slope_n_b + slope_sst_b +
    plot_layout(ncol = 3)) / 
  (scatter_l_b + scatter_n_b + scatter_sst_b +
     plot_layout(ncol = 3)) +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3_black.png", device = "png", height = 9, width = 14, dpi = 400)

# just scatter panels
scatter_l_b + scatter_n_b + scatter_sst_b +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3def_black.png", device = "png", height = 9, width = 16, dpi = 400)

# just slope panels
slope_l_b + slope_n_b + slope_sst_b +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3abc_black.png", device = "png", height = 7.5, width = 16, dpi = 400)


### JUST ALGAE -------
meta_df_a <- meta_df %>% filter(pp == "Algae")

# coeff plot
pal_ba <- c("#FDE725FF", "black", "black","black","black")

slope_l_b_a <- int_effect_plot(int_effect_sum_l, 
                               theme = "black", pal_ba) +
  facet_wrap(~dummy)

# bn slope
slope_n_b_a <- int_effect_plot(int_effect_sum, theme = "black", pal_ba) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  facet_wrap(~dummy)

# slope sst
slope_sst_b_a <- int_effect_plot(int_effect_sum_sst, theme = "black", pal_ba) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  facet_wrap(~dummy)

# sst with algae and seagrass
pal_bas <- c("#FDE725FF", "black", "black","black", "#2D708EFF")
slope_sst_b_as <- int_effect_plot(int_effect_sum_sst, theme = "black", pal_bas) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  facet_wrap(~dummy)

# just algae
slope_l_b_a + slope_n_b_a + slope_sst_b_a +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3abc_black_a.png", device = "png", height = 7.5, width = 16, dpi = 400)

# all + algae + algae
slope_l_b + slope_n_b_a + slope_sst_b_a +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3abc_black_a2.png", device = "png", height = 7.5, width = 16, dpi = 400)

# all + all + algae
slope_l_b + slope_n_b + slope_sst_b_a +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3abc_black_a3.png", device = "png", height = 7.5, width = 16, dpi = 400)

# all + all + algae/seagrass
slope_l_b + slope_n_b + slope_sst_b_as +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))
# ggsave("Figures/Pres_figs/Fig3abc_black_as.png", device = "png", height = 7.5, width = 16, dpi = 400)

# Scatter plots! 
# scatter plot latitude
scatter_l_b_a <- scatter_fun(lat_predict %>% filter(pp == "Algae"), 
                             meta_df_a,
                             abs_lat, pal_b, theme = "black") + 
  labs(x = "Absolute latitude") +
  ylim(-1.4, 3.9) +
  xlim(5.8, 74.3)

# scatter plot nutrients
scatter_n_b_a <- scatter_fun(nut_predict %>% filter(pp == "Algae"), 
                             meta_df_a %>% filter(bn_proxy == "measure"),
                             background_nut_umolL, pal_b, theme = "black") + 
  labs(x = expression("Background N"~(mu*M))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  ylim(-1.4, 3.9) +
  xlim(0.03, 17.65)

# scatter plot sst
scatter_sst_b_a <- scatter_fun(sst_predict %>% filter(pp == "Algae"), 
                               meta_df_a, 
                               sst, pal_b, theme = "black") + 
  labs(x = expression("Annual SST " ( degree*C))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position = "right") +
  ylim(-1.4, 3.9) +
  xlim(-0.8, 28.8)

# just scatter panels
scatter_l_b_a + scatter_n_b_a + scatter_sst_b_a+
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))

# ggsave("Figures/Pres_figs/Fig3def_black_algae.png", device = "png", height = 9, width = 16, dpi = 400)

### Algae and seagrass -------
meta_df_as <- meta_df %>% filter(pp == "Algae" | pp == "Seagrass")

# scatter plot latitude
scatter_l_b_as <- scatter_fun(lat_predict %>% 
                                filter(pp == "Algae" | pp == "Seagrass"), 
                              meta_df_as, 
                              abs_lat, pal_b, theme = "black") + 
  labs(x = "Absolute latitude")  +
  ylim(-1.4, 3.9) +
  xlim(5.8, 74.3)

# scatter plot nutrients
scatter_n_b_as <- scatter_fun(nut_predict %>% filter(pp == "Algae" | pp == "Seagrass"), 
                              meta_df_as %>% filter(bn_proxy == "measure"),
                              background_nut_umolL, pal_b, theme = "black") + 
  labs(x = expression("Background N"~(mu*M))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  ylim(-1.4, 3.9) +
  xlim(0.03, 17.65)

# scatter plot sst
scatter_sst_b_as <- scatter_fun(sst_predict %>% filter(pp == "Algae" | pp == "Seagrass"), 
                                meta_df_as, 
                                sst, pal_b, theme = "black") + 
  labs(x = expression("Annual SST " ( degree*C))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position = "right")  +
  ylim(-1.4, 3.9) +
  xlim(-0.8, 28.8)

# just scatter panels
scatter_l_b_as + scatter_n_b_as + scatter_sst_b_as +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))

# ggsave("Figures/Pres_figs/Fig3def_black_as.png", device = "png", height = 9, width = 16, dpi = 400)


### Algae and seagrass and phyto -------
meta_df_as <- meta_df %>% filter(pp == "Algae" | pp == "Seagrass" | pp == "Phytoplankton")

# scatter plot latitude
scatter_l_b_asp <- scatter_fun(lat_predict %>% 
                                 filter(pp == "Algae" | pp == "Seagrass" | pp == "Phytoplankton"), 
                               meta_df_as, 
                               abs_lat, pal_b, theme = "black") + 
  labs(x = "Absolute latitude") +
  ylim(-1.4, 3.9)+
  xlim(5.8, 74.3)

# scatter plot nutrients
scatter_n_b_asp <- scatter_fun(nut_predict %>% 
                                 filter(pp == "Algae" | pp == "Seagrass" | pp == "Phytoplankton"), 
                               meta_df_as %>% filter(bn_proxy == "measure"),
                               background_nut_umolL, pal_b, theme = "black") + 
  labs(x = expression("Background N"~(mu*M))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  ylim(-1.4, 3.9) +
  xlim(0.03, 17.65)

# scatter plot sst
scatter_sst_b_asp <- scatter_fun(sst_predict %>% 
                                   filter(pp == "Algae" | pp == "Seagrass" | pp == "Phytoplankton"), 
                                 meta_df_as, 
                                 sst, pal_b, theme = "black") + 
  labs(x = expression("Annual SST " ( degree*C))) +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position = "right")  +
  ylim(-1.4, 3.9) +
  xlim(-0.8, 28.8)

# just scatter panels
scatter_l_b_asp + scatter_n_b_asp + scatter_sst_b_asp +
  plot_annotation(theme = theme(
    plot.background = element_rect(fill = "black", colour = "black"),
    panel.background = element_rect(fill = "black", colour = "black")))

# ggsave("Figures/Pres_figs/Fig3def_black_asp.png", device = "png", height = 9, width = 16, dpi = 400)
