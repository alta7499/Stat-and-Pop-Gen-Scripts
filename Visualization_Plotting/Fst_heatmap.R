## Creates heatmap from FST values generated from PLINK. 
## Usage: run R script as it is, change paths and files as necessary.
## works for almost all types of FST file so long as the file contains 3 columns, Pop1 Pop2 FST.

# Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Load the data
data <- read.table("PHASE2/Global_fst.SEA_global.fst.summary", header=TRUE)  # Modify filename if needed

# Convert the dataset into a symmetric matrix
fst_matrix <- data %>%
  reshape2::acast(POP1 ~ POP2, value.var = "HUDSON_FST")  # Convert to wide format

# Ensure the matrix is symmetric by copying values
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]

# Convert matrix to long format for ggplot2
fst_long <- melt(fst_matrix, na.rm=TRUE)

# Create the heatmap
ggplot(fst_long, aes(Var1, Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="yellow", high="red", name="Fst") +  # Adjust color scale as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=8),
        axis.text.y = element_text(size=8)) +
  labs(x="Population", y="Population", title="Fst Heatmap of Global Populations") +
  coord_fixed()
