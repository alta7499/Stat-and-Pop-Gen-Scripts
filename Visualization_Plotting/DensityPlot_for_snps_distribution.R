library(ggplot2)
library(scales)

# Read the data
# Replace 'your_file.txt' with the path to your file
snp_data <- read.table("1.5m_snps_pos.txt", header = FALSE, col.names = c("Chromosome", "Position"))

# Convert Chromosome to a factor to maintain order
snp_data$Chromosome <- factor(snp_data$Chromosome, levels = sort(unique(snp_data$Chromosome)))

# Create the density plot
ggplot(snp_data, aes(x = Position, fill = Chromosome)) +
  geom_density(alpha = 0.5) +  # Density plot with transparency
  facet_wrap(~ Chromosome, scales = "free_y", ncol = 4) +  # One plot per chromosome
  labs(
    title = "Density Plot of SNP Physical Positions of 1.5 Million SNPs (imputed)",
    x = "Physical Position (bp)",
    y = "Density"
  ) +
  scale_y_continuous(labels = function(x) x / 1000,  # Convert y-axis to thousands
                     name = "Density (x 1,000)") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    strip.text = element_text(size = 10)  # Adjust facet label size
  )

ggplot(snp_data, aes(x = Chromosome, y = Position)) +
  geom_violin(fill = "skyblue") +
  labs(
    title = "Violin Plot of SNP Positions by Chromosome",
    x = "Chromosome",
    y = "Physical Position (bp)"
  ) +
  theme_minimal()

ggplot(snp_data, aes(x = Position)) +
  geom_rug() +
  labs(
    title = "Rug Plot of SNP Positions",
    x = "Physical Position (bp)",
    y = ""
  ) +
  theme_minimal()
