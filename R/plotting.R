# plot the output from simulation
scaleFUN <- function(x) sprintf("%.2f", x)

ap <- sim_res %>%
  group_by(Umax) %>%
  mutate(Dbar2 = 1 - round(mean(Dbar), 2)) %>%
  filter(value == "ar") %>%
  filter(scenario == "depleted") %>%
  ggplot(aes(y = median, x = as.factor(Dbar2))) +
  geom_hline(yintercept = c(ar), lty = 3, color = "black", lwd = 1) +
  geom_violin(width = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  ggtitle(expression(ln ~ alpha)) +
  ylab("Posterior medians") +
  xlab("Mean depletion level") +
  theme_qfc() +
  scale_y_continuous(
    labels = scaleFUN, limits = c(0.5, 2),
    breaks = c(0.50, 1.00, 1.50, 2.00)
  ) +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
ap

bp <- sim_res %>%
  group_by(Umax) %>%
  mutate(Dbar2 = 1 - round(mean(Dbar), 2)) %>%
  filter(value == "b") %>%
  filter(scenario == "depleted") %>%
  ggplot(aes(y = median, x = as.factor(Dbar2))) +
  geom_hline(yintercept = b, lty = 3, color = "black", lwd = 1) +
  geom_violin(width = 1.5) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  ggtitle(expression(beta)) +
  ylab("Posterior medians") +
  xlab("Mean depletion level") +
  theme_qfc() +
  scale_y_continuous(
    labels = scaleFUN, limits = c(0, 13),
    breaks = c(1.00, 5.00, 9.00)
  ) +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
bp

sigo <- sim_res %>%
  group_by(Umax) %>%
  mutate(Dbar2 = 1 - round(mean(Dbar), 2)) %>%
  filter(value == "sdo") %>%
  filter(scenario == "depleted") %>%
  ggplot(aes(y = median, x = as.factor(Dbar2))) +
  geom_hline(yintercept = sdo, lty = 3, color = "black", lwd = 1) +
  geom_violin(width = 0.25) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  ggtitle(expression(sigma[observation])) +
  ylab("Posterior medians") +
  xlab("Mean depletion level") +
  theme_qfc() +
  scale_y_continuous(labels = scaleFUN) +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
sigo

sigp <- sim_res %>%
  group_by(Umax) %>%
  mutate(Dbar2 = 1 - round(mean(Dbar), 2)) %>%
  filter(value == "sdp") %>%
  filter(scenario == "depleted") %>%
  ggplot(aes(y = median, x = as.factor(Dbar2))) +
  geom_hline(yintercept = sdp, lty = 3, color = "black", lwd = 1) +
  geom_violin(width = 0.25) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  ggtitle(expression(sigma[process])) +
  ylab("Posterior medians") +
  xlab("Mean depletion level") +
  theme_qfc() +
  scale_y_continuous(labels = scaleFUN) +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
sigp

s_msy <- sim_res %>%
  group_by(Umax) %>%
  mutate(Dbar2 = 1 - round(mean(Dbar), 2)) %>%
  filter(value == "smsy") %>%
  filter(scenario == "depleted") %>%
  filter(median < 1e2) %>%
  ggplot(aes(y = median, x = as.factor(Dbar2))) +
  geom_hline(yintercept = smsy, lty = 3, color = "black", lwd = 1) +
  geom_violin(width = 0.12) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  ggtitle(expression(S[MSY])) +
  ylab("Posterior medians") +
  xlab("Mean depletion level") +
  theme_qfc() +
  scale_y_continuous(labels = scaleFUN) +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
s_msy

h_msy <- sim_res %>%
  group_by(Umax) %>%
  mutate(Dbar2 = 1 - round(mean(Dbar), 2)) %>%
  filter(value == "hmsy") %>%
  filter(scenario == "depleted") %>%
  ggplot(aes(y = median, x = as.factor(Dbar2))) +
  geom_hline(yintercept = hmsy, lty = 3, color = "black", lwd = 1) +
  geom_violin(width = 0.12) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  ggtitle(expression(h[MSY])) +
  ylab("Posterior medians") +
  xlab("Mean depletion level") +
  theme_qfc() +
  scale_y_continuous(labels = scaleFUN) +
  stat_summary(fun = median, geom = "point", size = 3, col = "darkorange3")
h_msy

p <- plot_grid(ap, bp, sigo, sigp, s_msy, h_msy, nrow = 6, scale = 0.98)
p

ggsave("plots/ts-bias-depleted.pdf", width = 8, height = 11)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# call the function and make some more plots
set.seed(3)
dat <- sr_model()

p1 <- dat %>%
  ggplot(aes(x = year, y = S)) +
  geom_point() +
  geom_line() +
  theme_qfc() +
  ylim(0, 1.0) +
  ggtitle("True S vs. time") +
  geom_hline(yintercept = 0.1, linetype = 2)

p2 <- dat %>%
  ggplot(aes(x = year, y = Ut)) +
  geom_point() +
  geom_line() +
  theme_qfc() +
  ylab("Exploitation rate U(t)") +
  ggtitle("Ut vs. time")

p3 <- dat %>%
  ggplot(aes(x = S, y = ln_RS)) +
  geom_point() +
  ylab("ln(R/S)") +
  xlab("S") +
  ggtitle("True ln(R/S) vs. S relationship") +
  theme_qfc()

p4 <- dat %>%
  ggplot(aes(x = year, y = E)) +
  geom_point() +
  xlab("year") +
  ylab("Escapement") +
  theme_qfc() +
  geom_hline(yintercept = 0.1, linetype = 2) +
  ggtitle("Observed Escapement vs. time")

p5 <- dat %>%
  ggplot(aes(x = year, y = wt)) +
  geom_line() +
  xlab("year") +
  ylab(expression(sigma[process])) +
  theme_qfc() +
  ggtitle(expression(sigma[process] ~ vs. ~ time))

p6 <- dat %>%
  ggplot(aes(x = year, y = vt)) +
  geom_line() +
  xlab("year") +
  ylab(expression(sigma[obs])) +
  theme_qfc() +
  ggtitle(expression(sigma[obs] ~ vs. ~ time))

p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
p

ggsave("plots/ts-simulation-demonstration.pdf", width = 11, height = 8)



# next chunk requires a single fit called fit:

R2 <- fit %>% spread_draws(R[year])
R2 <- R2 %>%
  filter(year > k) %>%
  mutate(year = year - k)
S2 <- fit %>% spread_draws(S[year])
S2 <- S2 %>% filter(year <= (n_year - k))

dat <- left_join(R2, S2)
dat$ln_RS <- log(dat$R / dat$S)

# pluck out the best estimated ar, b
ests <- fit %>%
  gather_draws(ar, b) %>%
  median_qi()

ar_est <- ests$.value[which(ests$.variable == "ar")]
b_est <- ests$.value[which(ests$.variable == "b")]

devs <- fit %>%
  spread_draws(ar, b) %>%
  sample_n(size = 300)

p4 <- dat %>%
  median_qi() %>%
  ggplot(aes(x = S, y = ln_RS, ymin = ln_RS.lower, ymax = ln_RS.upper)) +
  geom_pointinterval(alpha = 0.25) +
  ylab("ln(R/S)") +
  ggqfc::theme_qfc()

p4 <- p4 + geom_abline(
  intercept = devs$ar, slope = -devs$b, color = "black",
  linetype = 1, size = 0.5, alpha = .05
)

p4 <- p4 + geom_abline(
  intercept = ar_est, slope = -b_est, color = "black",
  linetype = 2, size = 1.5
)

p4 <- p4 + geom_abline(
  intercept = ar, slope = -b, color = "steelblue",
  linetype = 2, size = 2
)

p4 <- p4 + ggtitle("Estimated (black/gray) ln(R/S) vs. S  vs. true relationship (blue)")

set.seed(3)
dat2 <- sr_model()
dat2$color <- c(rep("not included", n_year / 2), rep("included", n_year / 2))
p5 <- dat2 %>%
  ggplot(aes(x = S, y = ln_RS, color = color)) +
  geom_point() +
  scale_color_manual(values = c("black", "steelblue")) +
  ylab("ln(R/S)") +
  xlab("S") +
  ggtitle("True ln(R/S) vs. S relationship") +
  theme_qfc()

p5

p6 <- dat2 %>%
  ggplot(aes(y = S, x = year, color = color)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("black", "steelblue")) +
  ylab("Stock Size (S)") +
  xlab("Year") +
  ggtitle("Visualizing what the assessment can and cannot 'see' in time series bias scenario") +
  theme_qfc()
p6

p <- plot_grid(p6, p5, p4, ncol = 1)

ggsave("plots/ts-bias-one-fit.pdf", width = 8, height = 11)

#-------------------------------------------------------------------------------
# even more plotting
#-------------------------------------------------------------------------------
#
# p1 <- fit %>%
#   spread_draws(R[year]) %>%
#   median_qi() %>%
#   ggplot(aes(x = year, y = R, ymin = .lower, ymax = .upper)) +
#   geom_pointinterval()
# p1
# dat <- data.frame(year = 1:(n - k), Rtrue = R[(n_year - (n - (k + 1))):n_year])
# p1 <- p1 + geom_point(
#   data = dat, aes(
#     x = year, y = Rtrue, ymin = Rtrue,
#     ymax = Rtrue
#   ), shape = 16, color = "red",
#   size = 2
# ) +
#   ggqfc::theme_qfc()
#
# p1
#
# p2 <- fit %>%
#   spread_draws(S[year]) %>%
#   median_qi() %>%
#   ggplot(aes(x = year, y = S, ymin = .lower, ymax = .upper)) +
#   geom_pointinterval()
# p2
# dat <- data.frame(year = 1:n, Strue = S[(n_year - (n - 1)):n_year])
# p2 <- p2 + geom_point(
#   data = dat, aes(
#     x = year, y = Strue, ymin = Strue,
#     ymax = Strue
#   ), shape = 16, color = "red",
#   size = 2
# ) +
#   ggqfc::theme_qfc()
# p2
#
# p3 <- fit %>%
#   gather_draws(ar, b, sdp, sdo) %>%
#   median_qi() %>%
#   ggplot(aes(x = .variable, ymin = .lower, ymax = .upper, y = .value)) +
#   geom_pointinterval() +
#   xlab("Parameter") +
#   ylab("Value") +
#   ggqfc::theme_qfc()
# dat <- data.frame(
#   .variable = c("ar", "b", "sdo", "sdp"),
#   .value = c(ar, b, sdo, sdp)
# )
# p3 <- p3 + geom_point(
#   data = dat, aes(
#     x = .variable, y = .value, ymin = .value,
#     ymax = .value
#   ),
#   shape = 16, color = "red",
#   size = 2
# )
# p3
#
# p <- plot_grid(p3, p1, p2, ncol = 1)
# p
# ggsave("plots/self-test-su-peterson.pdf", width = 8, height = 11)

# create a few posteriors for stock-recruit, i.e., ln(R/S) vs. S
