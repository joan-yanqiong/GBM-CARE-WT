#################################
## Title: Figure 3 panel d in Nomura et al - 3-axis per sample
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: show each sample on a 3-axis coordinate system
#################################

####################################################################################################################################
# Upper panel
####################################################################################################################################

ggplot(feature_scores_per_sample %>%
         mutate(Label = ifelse(ID %in% c("P64T2", "P66T2", "P51T1", "P47T1", "P76T2", "P37T2"), ID, "")),
       aes(x = LineagePlot, y = C5)) +
  geom_point(aes(fill = axisColorCode), color = "black", shape = 22, size = 5) +
  geom_text(aes(label = Label), size = 8, hjust = "left", vjust = "outward") +
  geom_smooth(method = "loess", color = "black", size = 1, se = F) +
  scale_fill_gradientn(name = "", colours = axis_color_vec, guide = "colorsteps", rescaler = function(x, ...) { return(x) },
                       breaks = seq(0, 1, .05), labels = NULL,
                       values = c(0, .1, .125, .15, .175, .2, .225, .25, .3,
                                  .35, .45, .475, .5, .525, .55, .575, .6, .65,
                                  .7, .8, .825, .85, .875, .9, .925, .95, 1)) +
  scale_x_continuous(limits = c(-3, 3), oob = squish, breaks = seq(-3, 3, 1)) +
  scale_y_continuous(limits = c(-3, 2), oob = squish) +
  xlab("ECM <-> Neuronal") +
  ylab("-> Glial") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(panel.border = element_blank()) +
  theme(panel.grid.major = element_line())

####################################################################################################################################
# Lower panel
####################################################################################################################################

d <- mdata

d <- d %>%
  ungroup() %>%
  left_join(scp_scores %>%
              select(CellID, MaxAxis, axisColorCode),
            by = "CellID")

ggplot(data = d %>%
         filter(ID %in% c("P64T2", "P66T2", "P51T1", "P47T1", "P76T2", "P37T2")) %>%
         mutate(ID_factor = factor(ID, c("P64T2", "P66T2", "P51T1", "P47T1", "P76T2", "P37T2"))),
       mapping = aes(x = Dx, y = Dy)) +
  facet_grid(cols = vars(ID_factor)) +
  geom_point(aes(fill = axisColorCode), color = "black", shape = 21, size = 3) +
  scale_fill_gradientn(name = "", colours = axis_color_vec, guide = "colorsteps", rescaler = function(x, ...) { return(x) },
                       breaks = seq(0, 1, .05), labels = NULL,
                       values = c(0, .1, .125, .15, .175, .2, .225, .25, .3,
                                  .35, .475, .5, .525, .55, .575, .6, .65,
                                  .7, .8, .825, .85, .875, .9, .925, .95, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.5, 1.5)) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colorsteps(ticks.colour = "black", frame.colour = "black", ticks.linewidth = 1, frame.linewidth = 1)) +
  theme_gbm_pvsr(legend.position = "none")
