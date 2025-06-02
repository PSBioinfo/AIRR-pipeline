# ---------------------------------------------
# clone_analysis.R
# Minimal BCR clone size analysis and visualization
# ---------------------------------------------

# Load required libraries
library(dplyr)
library(ggplot2)

# Define file paths
input_file <- "data/example_db.tsv"
output_file <- "output/clone_histogram.png"

# Read the Change-O formatted dataset
data <- read.delim(input_file, stringsAsFactors = FALSE)

# Group sequences by CLONE and count how many in each
clone_sizes <- data %>%
  group_by(CLONE) %>%
  summarize(size = n()) %>%
  filter(!is.na(CLONE))

# Plot clone size distribution as a histogram
plot <- ggplot(clone_sizes, aes(x = size)) +
  geom_histogram(binwidth = 1, fill = "#377eb8", color = "black") +
  theme_minimal() +
  labs(
    title = "Clone Size Distribution",
    x = "Clone Size (number of sequences)",
    y = "Number of Clones"
  )

# Save the plot to output directory
ggsave(output_file, plot, width = 7, height = 5)

cat("✅ Analysis complete. Plot saved to:", output_file, "\n")

# ---------------------------------------------
# Calculate Shannon Diversity Index
# ---------------------------------------------

# Calculate proportion of sequences in each clone
total_sequences <- sum(clone_sizes$size)
clone_sizes <- clone_sizes %>%
  mutate(p_i = size / total_sequences)

# Calculate Shannon index
shannon_index <- -sum(clone_sizes$p_i * log(clone_sizes$p_i))

# Print the result
cat("🔢 Shannon Diversity Index:", round(shannon_index, 3), "\n")


write(paste("Shannon Diversity Index:", round(shannon_index, 3)),
      file = "output/shannon_index.txt")

# --------------------------------------------------------
# Calculate Shannon Index PER SAMPLE and Plot
# --------------------------------------------------------

library(tibble)

# Function to calculate Shannon index from a vector of clone sizes
calc_shannon <- function(counts) {
  proportions <- counts / sum(counts)
  -sum(proportions * log(proportions))
}

# Calculate Shannon index for each sample
shannon_by_sample <- data %>%
  group_by(SAMPLE, CLONE) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(SAMPLE) %>%
  summarise(shannon = calc_shannon(count))

# Save to file
write.table(shannon_by_sample, "output/shannon_by_sample.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot Shannon index per sample
shannon_plot <- ggplot(shannon_by_sample, aes(x = SAMPLE, y = shannon)) +
  geom_col(fill = "#4daf4a") +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Sample",
       x = "Sample",
       y = "Shannon Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/shannon_by_sample.png", shannon_plot, width = 8, height = 5)

# --------------------------------------------------------
# Calculate Simpson's Reciprocal Index per Sample
# --------------------------------------------------------

# Function to calculate Simpson's reciprocal index
calc_simpson <- function(counts) {
  proportions <- counts / sum(counts)
  1 / sum(proportions^2)
}

# Calculate for each sample
simpson_by_sample <- data %>%
  group_by(SAMPLE, CLONE) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(SAMPLE) %>%
  summarise(simpson = calc_simpson(count))

# Save to file
write.table(simpson_by_sample, "output/simpson_by_sample.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot Simpson's index per sample
simpson_plot <- ggplot(simpson_by_sample, aes(x = SAMPLE, y = simpson)) +
  geom_col(fill = "#984ea3") +
  theme_minimal() +
  labs(title = "Simpson's Reciprocal Diversity Index by Sample",
       x = "Sample",
       y = "Simpson Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/simpson_by_sample.png", simpson_plot, width = 8, height = 5)
