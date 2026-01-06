library(ggplot2)

# load F4 table, tab delimited. from ADMIXTOOLS2 F4 output, add in one more column of absolute_z.
# Population and F4 configuration must be pre-filtered before use.
df <- read.table("F4_OA_NB.txt", header = TRUE, stringsAsFactors = FALSE)

# Ensure numeric columns
df$est <- as.numeric(df$est)
df$se  <- as.numeric(df$se)

# ordering by direction of gene flow
df$direction <- ifelse(df$z > 0, "NB closer", "TaiwanAbo closer")
df$sig <- abs(df$z) >= 3
df$pop1 <- factor(
  df$pop1,
  levels = names(sort(tapply(df$z, df$pop1, mean), decreasing = TRUE))
)

ggplot(df, aes(x = z, y = pop1,
               color = direction,
               shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 3) +
  facet_wrap(~ pop3, scales = "free_y") +
  scale_color_manual(
    values = c("NB closer" = "#D55E00",        # red
               "TaiwanAbo closer" = "#0072B2") # blue
  ) +
  scale_shape_manual(
    values = c("FALSE" = 1, "TRUE" = 16)
  ) +
  labs(
    x = "Z score  (← TaiwanAbo closer | NB closer →)",
    y = "North Borneo population",
    color = "Relative affinity",
    shape = "Significance",
    title = "Sliding-dot plot of F4(NB, TaiwanAbo; OA, YRI)",
    subtitle = "Color = direction of allele sharing; filled points indicate |Z| ≥ 3"
  ) +
  theme_bw(base_size = 12)
