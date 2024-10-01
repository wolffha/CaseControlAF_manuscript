calc_b <- function(N_case, N_control, AF_total, OR) {
  N_total <- N_case + N_control
  
  return(OR*(1-(N_total/N_case*AF_total))+1/N_case*(N_control+N_total*AF_total))
}
calc_a <- function(N_case, N_control, OR) {
  return((N_control/N_case)*(OR - 1))
}


N_case <- rep(c(5000, 5000, 5000, 5000, 9000, 9000, 9000, 9000, 1000, 1000, 1000, 1000), 6)
N_control <- rep(c(5000, 5000, 5000, 5000, 1000, 1000, 1000, 1000, 9000, 9000, 9000, 9000), 6)
AF <- c(rep(0, 12), rep(0.1, 12), rep(0.25, 12), rep(0.5, 12), rep(0.75, 12), rep(1, 12))

OR <- rep(rep(c(0.1, 0.4, 0.6, 0.8), 3), 6) # for negative a
df1 <- data.frame(OR, AF, N_case, N_control)

OR <- rep(rep(c(1.1, 1.2, 3, 5), 3), 6) # for positive a
df2 <- data.frame(OR, AF, N_case, N_control)

df <- rbind(df1, df2)
#colnames(df) <- c("OR", "AF", "N_case", "N_control")
df$b <- apply(df, 1, function(x)
  calc_b(x[3], x[4], x[2], x[1]))
df$a <- apply(df, 1, function(x)
  calc_a(x[3], x[4], x[1]))

df$sym <- (-df$b)/(2*df$a)
summary(df)

# solve for real roots of a quadratic using the quadratic formula
quad_roots<-function(a,b,c){
  c(((-b-sqrt(b^2-4*a*c))/(2*a)),((-b+sqrt(b^2-4*a*c))/(2*a)))
}

#calculate total sample size
df$N_total <- df$N_control+df$N_case

#set a, b, c of quadratic equation derived in manuscript
df$a <- (df$N_control/df$N_case)*(df$OR-1)
df$b <- (df$OR-((df$N_total/df$N_case)*df$AF*df$OR))+((df$N_control/df$N_case)+(df$N_total*df$AF/df$N_case))
df$c <- -(df$N_total/df$N_case)*df$AF

df$root1 <- 0
df$root2 <- 0
for(i in 1:nrow(df)) {
  roots <- quad_roots(df[i,]$a, df[i,]$b, df[i,]$c)
  df[i,]$root1 <- roots[1]
  df[i,]$root2 <- roots[2]
}
df$root1_valid <- apply(df, 1, function(x)
  ifelse(x[10] >= 0 & round(x[10], 4) <= 1, 1, 0))
df$root2_valid <- apply(df, 1, function(x)
  ifelse(x[11] >= 0 & round(x[11], 4) <= 1, 1, 0))

table(df$root1_valid + df$root2_valid)

summary(df)

out <- df[df$sym < 1 & df$sym > 0, ]
out

pa <- ggplot(df, aes(x = OR, y = b, color = as.factor(AF))) + geom_point()
pb <- ggplot(df, aes(x = AF, y = b, color = as.factor(AF))) + geom_point()
pc <- ggplot(df, aes(x = N_case, y = b, color = as.factor(AF))) + geom_point()
pd <- ggplot(df, aes(x = N_control, y = b, color = as.factor(AF))) + geom_point()
pe <- ggplot(df, aes(x = (N_case + N_control), y = b, color = as.factor(AF))) + geom_point()

plot_grid(pa, pb, pc, pd, pe)

df$simulation <- as.numeric(row.names(df))

# scenario 1 is a+ | b+
# scenario 2: a+ | b-
# scenario 3: a- | b+
# scenario 4: a- | b-

df$a_neg <- ifelse(df$a > 0, 0, 1)
df$b_neg <- ifelse(df$b > 0, 0, 1)

df$scenario <- case_when(
  df$a_neg == 0 & df$b_neg == 0 ~ 1,
  df$a_neg == 0 & df$b_neg == 1 ~ 2,
  df$a_neg == 1 & df$b_neg == 0 ~ 3,
  df$a_neg == 1 & df$b_neg == 1 ~ 4
)

df$a_neg <- factor(df$a_neg, levels = c(0, 1),
                   labels = c("Positive", "Negative"))
df$b_neg <- factor(df$b_neg, levels = c(0, 1),
                   labels = c("Positive", "Negative"))

longroots <- pivot_longer(df, cols = c(root1, root2))
ggplot(longroots, aes(x = name, y = value, group = as.factor(simulation))) +
  annotate(geom = "rect", xmin = 1, xmax = 2, ymin = 0, ymax = 1, fill = "red",
           color = "white", alpha = 0.2) +
  geom_line(linewidth = 0.1, color = "grey30") +
  geom_point(aes(color = as.factor(b_neg)), size = 5, 
             position = position_jitter(h = 0.0, w = 0.05)) +
  theme_bw() +
  ylab("Value (AF_control)") +
  ggtitle("Two Root Solutions by Scenario") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  labs(color = "Sign of b") +
  ggforce::facet_zoom(ylim = c(-.05, 1.05))
ggsave(filename = "RootSolutionSims_20240926_bNegColor.png", width = 14, height = 8, units = "in")

df[df$a < 0 & AF == 1,]

tosave <- df %>% select(OR, AF, N_case, N_control, N_total, a, b, c,
                        root1, root2, simulation)
write.csv(tosave, file = "rootSimulations.csv", row.names = F, quote = F)
