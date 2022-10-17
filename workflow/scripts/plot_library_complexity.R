log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

#------------------------------------------------------------------------------------------
# Load in Libraries + custom themes
#------------------------------------------------------------------------------------------

# Lib
library("ggplot2")
library("dplyr")
library("tibble")
library("purrr")
library("tidyr")
library("RColorBrewer")
library("scales")

# Load in custom functions
source("workflow/scripts/custom_functions.R")

# theme
theme_custom <- theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        legend.key.size  = unit(0.4, units = "cm"),
  ) +
  theme(plot.title = element_text(hjust = 0.5)
  ) +
  theme(
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  theme(
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(fill="white")
  )
#------------------------------------------------------------------------------------------
# Read in complexity curve data and estimated future yields
#------------------------------------------------------------------------------------------

c_curve <- snakemake@input$c_curve %>%
	purrr::map(read.delim, header = TRUE) %>%
	purrr::map(setNames, c("total_reads", "distinct_reads"))


names(c_curve) <- as.character(snakemake@params[["sample_names"]])


lc_extrap <- snakemake@input$lc_extrap %>%
	purrr::map(read.delim, header = TRUE) %>%
	purrr::map(setNames, c("TOTAL_READS", "EXPECTED_DISTINCT", "LOWER_0.95CI", "UPPER_0.95CI"))

names(lc_extrap) <- as.character(snakemake@params[["sample_names"]])

#------------------------------------------------------------------------------------------
# Wrangling data
#------------------------------------------------------------------------------------------

# c_curve
comp_curve <- c_curve %>% 
  map2_df(names(c_curve), ~ mutate(.x, Experiment = .y))

# substitute these ones with the package (scales) and the functions commented in the plot
#comp_curve[,1:2] <- comp_curve[, 1:2] /1e6


# lc_extrap
expct_curve <- lc_extrap %>% 
  map2_df(names(lc_extrap), ~ mutate(.x, Experiment = .y)) %>% 
  select(1:2, 5)

# substitute these ones with the package (scales) and the functions commented in the plot
#expct_curve[,1:2] <- expct_curve[,1:2] /1e6


#df_full

final_df <- left_join(comp_curve, expct_curve, c("total_reads" = "TOTAL_READS", "Experiment")) %>% 
  pivot_longer(cols = c("distinct_reads", "EXPECTED_DISTINCT"), names_to = "Distinct_Reads") %>% 
  mutate(Type = case_when(Distinct_Reads %in% "distinct_reads" ~ "Observed",
                          TRUE ~ "Expected"))


#------------------------------------------------------------------------------------------
# Potting data
#------------------------------------------------------------------------------------------

# Deciding how many colors i need
color_count <- length(unique(final_df$Experiment))
get_palette <- colorRampPalette(brewer.pal(9, "Set1"))
 


# plot observed complexity curve
c_curve_plot <- comp_curve %>% 
  ggplot(aes(x = total_reads, y = distinct_reads, color = Experiment)) +
  geom_line() +
  theme_custom +
  #theme_tech(theme = "airbnb") +
  #scale_color_brewer(type = "qual", palette = "Set1") +
  scale_color_manual(values = get_palette(color_count)) +
  geom_abline(linetype = "dashed") +
  scale_y_continuous(labels = unit_format(suffix = "M", scale = 1e-6)) +
  scale_x_continuous(labels = unit_format(suffix = "M", scale = 1e-6)) +
  labs(x = "Total Reads (M)", y = "Distinct Reads (M)", title = "Observed Rarefraction Curve")


lc_extrap_plot <- expct_curve %>% 
  ggplot(aes(x = TOTAL_READS, y = EXPECTED_DISTINCT, color = Experiment)) +
  geom_line() +
  theme_custom +
  #theme_tech(theme = "airbnb") +
  scale_color_manual(values = get_palette(color_count)) +
  #scale_color_brewer(type = "qual", palette = "Set1") +
  geom_vline(xintercept = 50e6, linetype = "dashed") +
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  labs(x = "Total Reads (M)", y = "Expected Distinct Reads (M)", title = "Expected Rarefraction Curve")


obs_exp_plot <- final_df %>% 
  ggplot(aes(x = total_reads, y = value, color = Experiment, linetype = Type)) +
  geom_line() +
  theme_custom +
  #scale_color_brewer(type = "qual", palette = "Set1") +
  scale_color_manual(values = get_palette(color_count)) +
  facet_wrap(~Type) +
  theme(strip.background =element_rect(fill="darkslategray1"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  scale_x_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  labs(x = "Total Reads (M)", y = "Obs/Exp Distinct Reads (M)", title = "Rarefraction Curve of \nObserved and Expected reads")


#------------------------------------------------------------------------------------------
# Save output
#------------------------------------------------------------------------------------------
ggsave(filename = snakemake@output[["c_curve_plot"]], c_curve_plot, width = 14, height = 8)
ggsave(filename = snakemake@output[["lc_extrap_plot"]], lc_extrap_plot, width = 14, height = 8)
ggsave(filename = snakemake@output[["obs_exp_plot"]], obs_exp_plot, width = 14, height = 8)
