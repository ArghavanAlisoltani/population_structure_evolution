# Load the required libraries
library(tidyverse)
library(patchwork)

# Define a function to read, process, and plot the FST data
create_manhattan_plot <- function(filepath, title) {
  
  # Read the tab-separated data file
  fst_data <- read_tsv(filepath, col_names = TRUE, show_col_types = FALSE) %>%
    rename(CHROM = 1, POS = 2, FST = 3) # Rename columns for clarity
  
  # --- Data Cleaning and Preparation ---
  fst_data_cleaned <- fst_data %>%
    mutate(
      FST = as.numeric(FST),
      # Convert negative FST values and NaNs to 0
      FST = if_else(is.na(FST) | FST < 0, 0, FST),
      # Create a numeric version of the scaffold ID for proper sorting
      CHROM_NUM = as.numeric(gsub("[^0-9]", "", CHROM))
    ) %>%
    # Arrange the data by the numeric scaffold and then by position
    arrange(CHROM_NUM, POS)
  
  # --- Cumulative Position Calculation for a Continuous X-axis ---
  data_cumulative <- fst_data_cleaned %>%
    group_by(CHROM_NUM) %>%
    summarise(max_pos = max(as.numeric(POS)), .groups = 'drop') %>%
    mutate(bp_add = lag(cumsum(max_pos), default = 0))
  
  # Join the cumulative information back to the main dataframe
  fst_data_plot <- fst_data_cleaned %>%
    left_join(data_cumulative, by = "CHROM_NUM") %>%
    mutate(bp_cum = as.numeric(POS) + bp_add)
  
  # --- Generate the Plot ---
  plot <- ggplot(fst_data_plot, aes(x = bp_cum, y = FST)) +
    
    # Draw points, coloring by even/odd scaffold number to get alternating colors
    geom_point(alpha = 0.8, size = 0.7, aes(color = as.factor(CHROM_NUM %% 2))) +
    
    # Format the X-axis to display large numbers in Gigabases (Gb)
    scale_x_continuous(name = "Genomic Position", labels = scales::label_number(scale = 1e-9, suffix = " Gb"), expand = c(0.01, 0.01)) +
    
    scale_y_continuous(name = expression(F[ST]), limits = c(0, 0.9), expand = c(0, 0)) +
    
    # Manually set the alternating colors for even and odd scaffolds
    scale_color_manual(values = c("0" = "grey30", "1" = "steelblue")) +
    
    labs(title = title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      
      # --- FONT SIZE CONTROLS ---
      # Increase plot title font, make it bold, and center it
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      
      # Increase the font size for axis titles (e.g., "Fst", "Genomic Position")
      axis.title = element_text(size = 18),
      
      # Increase the font size for the text on the axes (the numbers/labels)
      axis.text = element_text(size = 16),
      
      # Keep the x-axis text angled for readability
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
    )
    
  
  return(plot)
}

# --- Create and Combine the Plots ---

# Specify the file paths for your three data files
file1 <- "PROC3_vs_PROC4.weir.fst"
file2 <- "PROC3_vs_PROC5.weir.fst"
file3 <- "PROC4_vs_PROC5.weir.fst"

# Generate a plot for each file
plot1 <- create_manhattan_plot(file1, "PROC3 vs. PROC4")
plot2 <- create_manhattan_plot(file2, "PROC3 vs. PROC5")
plot3 <- create_manhattan_plot(file3, "PROC4 vs. PROC5")

# Arrange the three plots into a single vertical column
combined_plot <- plot1 + plot2 + plot3 + plot_layout(ncol = 1)

# Display the combined plot
print(combined_plot)

# Save the final image to a file
ggsave("fst_manhattan_grid_vertical.png", combined_plot,
       width = 17, height = 10, dpi = 300)

