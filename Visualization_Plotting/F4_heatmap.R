library(ggplot2)

# load F4 table, tab delimited. from ADMIXTOOLS2 F4 output, add in one more column of absolute_z.
# Population and F4 configuration must be pre-filtered before use.
df <- read.table("F4_OA_NB.txt", header = TRUE, stringsAsFactors = FALSE)

# Ensure numeric columns
df$est <- as.numeric(df$est)
df$se  <- as.numeric(df$se)

# Order NB populations by mean Z
df$pop1 <- factor(
  df$pop1,
  levels = names(sort(tapply(df$z, df$pop1, mean), decreasing = TRUE))
)

# Order OA populations by mean Z
df$pop3 <- factor(
  df$pop3,
  levels = names(sort(tapply(df$z, df$pop3, mean), decreasing = TRUE))
)

ggplot(df, aes(x = pop3, y = pop1, fill = z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Z score"
  ) +
  labs(
    x = "Orang Asli population (pop3)",
    y = "North Borneo population (pop1)",
    title = "F4 Z-score heatmap: F4(NB, TaiwanAbo; OA, YRI)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )