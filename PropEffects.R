# Load required R packages ----
library(tidyverse)
library(nlme)
library(parallel)

# Cross-sectional simulations ----
N <- 100
S <- 10000
BETA <- -c(seq(0, 1, by=0.5),100)
DELTA <- seq(-1, 1, by=0.1)

results <-
  mclapply(BETA, mc.cores = 4, function(Beta){
    set.seed(20240731)
    lapply(DELTA, function(Delta){
      lapply(1:S, function(sim){
        dd <- tibble(
          Tx = c(rep(0, N/2), rep(1,N/2)),
          Y = c(rnorm(N/2, Beta), rnorm(N/2, Beta+Delta)))
        tt <- lm(Y~Tx, data = dd)
        pt <- nls(Y ~ beta*(1-Tx*theta),
          data = dd,
          start = list(beta=coef(tt)[[1]], theta=0))
        ci.pt <- tryCatch(confint(pt), error = function(msg){
          return(tibble(`2.5%`=NA))})

        bind_rows(
          summary(tt)$coef %>%
            as.data.frame() %>%
            rownames_to_column('term') %>%
            bind_cols(confint(tt)) %>%
            filter(term == 'Tx') %>%
            mutate(Approach = 't-test') %>%
            rename(`2.5%`=`2.5 %`, `97.5%`=`97.5 %`),
          summary(pt)$coef %>%
            as.data.frame() %>%
            rownames_to_column('term') %>%
            bind_cols(ci.pt) %>%
            filter(term == 'theta') %>%
            mutate(Approach = 'proportional test')) %>%
          mutate(Simulation = sim, Beta = Beta, Delta = Delta)
        }) %>%
      bind_rows()
    }) %>%
      bind_rows()
  }) %>%
      bind_rows()

## Figure 1. Power for cross-sectional model ----
# Simulated power (and type I error when $\\delta$ = 0) for the proportional 
# model and two-sample t-test based on 10,000 simulations. The control group 
# mean is $\\beta_C$ and $\\delta$ is the mean group difference. 
# With $\\beta_C = 0$, the proportional effect is not identifiable, and 
# simulated power with the proportional test is about 5% regardless of $\\delta$. 
# With $\\beta_C$ near zero ($\\beta_C=-0.5$ or $-1$), the proportional test 
# appears to have better power with $\\delta>0$, and worse power with $\\delta<0$. 
# With $\\beta_C$ far from zero ($\\beta_C=-100$), power for the two methods is 
# identical and the lines are overlapping. Residual variance was simulated to be 
# one, sample size was 50 per group, and the proportional model was fit with R's 
# nonlinear least squares (\\texttt{nls}) function.[@bates1992]
sum.results <- results %>%
  group_by(Beta, Delta, Approach) %>%
  summarise(
    `N sims` = length(Estimate),
    Mean = mean(Estimate),
    `Power (%)` = sum(`Pr(>|t|)`<0.05)/length(`Pr(>|t|)`)*100
  )

sum.results %>%
  select(Beta, Delta, Approach, `Power (%)`) %>%
  mutate(
    Beta = paste0('\"Control group mean, \" * beta[C]==', Beta) %>%
      factor(levels = paste0('\"Control group mean, \" * beta[C]==', 
        c(0, -0.5, -1, -100)))) %>%
  ggplot(aes(x=Delta, y=`Power (%)`)) +
    geom_line(aes(color=Approach)) +
    geom_point(aes(color=Approach, shape=Approach)) +
    scale_shape_manual(values = c(1,3)) +
    facet_wrap(vars(Beta),labeller = label_parsed) +
    xlab(expression("Difference between groups ("*delta*")")) +
    geom_hline(yintercept=5, linetype = 'dashed') +
    geom_vline(xintercept=0, linetype = 'dashed') +
    theme_minimal() +
    theme(legend.position='inside', legend.position.inside = c(0.5, 0.15)) +
    scale_color_manual(values = ggsci::pal_nejm()(n=5)[4:5])


## Figure 2. Zipper plots for cross-sectional model
# Zipper plots [@morris2019using] showing estimates and 95% confidence intervals 
# for cross-sectional simulations with control group mean $\\beta_C=-0.5$. 
# The top row has a treatment effect of $\\delta=0$, and demonstrates Type I error. 
# The bottom row has a treatment effect of $\\delta=0.3$, or equivalently, a 
# proportional effect $\\theta = 0.6$. The intervals are sorted so that those 
# associated with the largest standardized bias are toward the top of each panel, 
# and only 25% of simulations with the largest bias are shown. While t-test 
# estimates are symmetric about the true effect (left), the plots for the 
# proportional test (right) reveal bias and asymmetric confidence interval widths. 
# Vertical dashed lines are at the true value. Horizontal dashed lines are at 
# the target rejection rate of 5\\%. Note that due to the numeric instability of 
# the proportional model, p-values are often inconsistent with confidence interval 
# coverage as can be seen with red intervals (p<0.05) which cover zero in the 
# upper right panel
ttest_t1 <- results %>%
  filter(Beta == -0.5 & Delta == 0 & Approach == 't-test') %>%
  arrange(`Pr(>|t|)`) %>%
  mutate(
    rank = row_number()/length(`Pr(>|t|)`)*100,
    Significance = factor(`Pr(>|t|)` < 0.05,
      levels = c(TRUE, FALSE), labels = c('p<0.05', 'p>0.05'))) %>%
  filter(rank<=25) %>%
ggplot(aes(x=Estimate, y=rank, color = Significance)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`), alpha = 0.1) +
  scale_y_reverse() +
  ylab("P-value percentile") +
  xlab("Estimated mean difference") +
  theme_minimal() +
  theme(legend.position='none') +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ggtitle('t-test Type I error')
    
ptest_t1 <- results %>%
  filter(Beta == -0.5 & Delta == 0 & Approach == 'proportional test') %>%
  arrange(`Pr(>|t|)`) %>%
  mutate(
    rank = row_number()/length(`Pr(>|t|)`)*100,
    Significance = factor(`Pr(>|t|)` < 0.05,
      levels = c(TRUE, FALSE), labels = c('p<0.05', 'p>0.05'))) %>%
  filter(rank<=25) %>%
ggplot(aes(x=Estimate, y=rank, color = Significance)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`), alpha = 0.1) +
  scale_y_reverse() +
  ylab("P-value percentile") +
  xlab("Estimated proportional effect") +
  theme_minimal() +
  theme(legend.position='none') +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ggtitle('Proportional test Type I error') +
  coord_cartesian(xlim = c(-5,3))
  
ttest_pow <- results %>%
  mutate(ddelta = as.factor(Delta)) %>%
  filter(Beta == -0.5 & ddelta == 0.3 & Approach == 't-test') %>%
  mutate(absbias = abs(Estimate - Delta)/`Std. Error`) %>%
  arrange(desc(absbias)) %>%
  mutate(
    rank = row_number()/length(`Pr(>|t|)`)*100,
    Coverage = factor(`2.5%`>Delta | `97.5%`<Delta,
      levels = c(TRUE, FALSE), labels = c('Not covered', 'Covered'))) %>%
  filter(rank<=25) %>%
ggplot(aes(x=Estimate, y=rank, color = Coverage)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`), alpha = 0.1) +
  scale_y_reverse() +
  ylab("Standardized bias percentile") +
  xlab("Estimated mean difference") +
  theme_minimal() +
  theme(legend.position='none') +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  geom_vline(aes(xintercept = Delta), linetype = 'dashed') +
  ggtitle('t-test power')  +
  geom_segment(
    data = tibble(
      arrowx1 = -0.5, arrowx2 = 1, 
      arrowy1 = 36, arrowy2 = 36),
    color = 'black',
    aes(x = arrowx1, y = arrowy1, xend = arrowx2, yend = arrowy2),
    arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(
    data = tibble(
      x = 0.25, y = 38, text = "Favors active group"),
    aes(y=y, x=x, label = text), color = 'black') +
  coord_cartesian(ylim = c(25, 0), clip = "off") +
  theme(plot.margin = unit(c(0,0,2,0), "lines"))

ptest_pow <- results %>%
  mutate(ddelta = as.factor(Delta)) %>%
  filter(Beta == -0.5  & ddelta == 0.3 & Approach == 'proportional test') %>%
  mutate(absbias = abs(Estimate + Delta/Beta)/`Std. Error`) %>%
  arrange(desc(absbias)) %>%
  mutate(
    rank = row_number()/length(`Pr(>|t|)`)*100,
    Coverage = factor(`2.5%`> -Delta/Beta | `97.5%`< -Delta/Beta,
      levels = c(TRUE, FALSE), labels = c('Not covered', 'Covered'))) %>%
  filter(rank<=25) %>%
ggplot(aes(x=Estimate, y=rank, color = Coverage)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`), alpha = 0.1) +
  scale_y_reverse() +
  ylab("Standardized bias percentile") +
  xlab("Estimated proportional effect") +
  theme_minimal() +
  theme(legend.position='none') +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  geom_vline(aes(xintercept = -Delta/Beta), linetype = 'dashed') +
  ggtitle('Proportional test power')  +
  geom_segment(
    data = tibble(
      arrowx1 = -20, arrowx2 = 10, 
      arrowy1 = 36, arrowy2 = 36),
    color = 'black',
    aes(x = arrowx1, y = arrowy1, xend = arrowx2, yend = arrowy2),
    arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(
    data = tibble(
      x = -5, y = 38, text = "Favors active group"),
    aes(y=y, x=x, label = text), color = 'black') +
  coord_cartesian(ylim = c(25, 0), clip = "off") +
  theme(plot.margin = unit(c(0,0,2,0), "lines"))

gridExtra::grid.arrange(
  ttest_t1, ptest_t1, 
  ttest_pow, ptest_pow, ncol=2,
  heights = c(1,1.1))


# Longitudinal simulations ----
## Define group-level means for four scenarios ----
long.means <- expand_grid(
  month = seq(0,18,by=3),
  tx = c(0,1),
  scenario = 1:4) %>%
  mutate(
    # control group means
    Y0 = case_when(
      scenario == 1 ~ -1.5/18 * month,
      scenario == 2 ~ -1.5/2/18 * month,
      scenario == 3 ~ -0.5/18 * month,
      scenario == 4 ~ -0.25/18 * month) %>% 
        as.numeric(),
    # proportional treatment effects
    Y = case_when(
      scenario == 1 ~ Y0 * (1 - (1/3) * tx),
      scenario == 2 ~ Y0 * (1 - (2/3) * tx),
      scenario == 3 ~ Y0 * (1 - (1) * tx),
      scenario == 4 ~ Y0 * (1 - (2) * tx)) %>% 
        as.numeric(),
    Tx = factor(tx, levels = 1:0, labels = c('Active', 'Control')),
    Scenario = factor(scenario, levels = 1:4,
      labels = paste("Scenario", c("A", "B", "C", "D"))))


## Figure 3. Longitudinal mean trend scenarios ----
# Four longitudinal mean trend scenarios all with different proportional effects
# (ranging from 1/3 to 2) but equivalent linear effects (a difference between 
# slopes of 0.5 outcome units per 18 months)
slope_effect <- - (-1.5 - (-1.5 * (1 - (1/3))))/18

ggplot(long.means, aes(x=month, y=Y, color=Tx)) +
  geom_line() +
  geom_point(aes(shape = Tx)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(vars(Scenario)) +
  scale_x_continuous(breaks = seq(0, 18, by = 3)) +
  ylab("Group mean over time") +
  theme_minimal() +
  theme(legend.position='inside', legend.position.inside = c(0.1, 0.1)) +
  ggsci::scale_color_nejm()


## Expand group-level data to individual-level data ----
ind.long.means <- long.means %>%
  cross_join(tibble(id = 1:200)) %>%
  mutate(id = tx*200 + id)
N <- max(ind.long.means$id)

## Define proportional cLDA model ----
# covariates: month, tx (0 control, 1 active)
# parameters: Intercept, beta.month, theta
PcLDA <- function(month, tx, Intercept, beta.month, theta){
  Intercept + month*beta.month*(1-theta*tx) 
}

## Define function to fit models to simulated data ----
# this function is used for both power and type I error simulations
return.model.results <- function(data){
  lapply(1:4, function(scen){
    # fit linear model
    linear.model <- lme(Ysim ~ month + month:tx,
      random = ~ 1 | id, data = data %>% filter(scenario == scen))

    # extract estimate and p-value for difference between slopes
    linear.result <- summary(linear.model)$tTable %>%
      as.data.frame() %>%
      rownames_to_column('term') %>%
      bind_cols(intervals(linear.model)$fixed) %>%
      filter(term == 'month:tx') %>%
      mutate(Model = 'Linear') %>%
      select(-Value)
    
    # use linear model fit for starting values for nonlinear model
    start1 <- fixef(linear.model)[c("(Intercept)", "month")]
    names(start1) <- c("Intercept", "beta.month")
    start1 <- c(start1, theta = 0)
    
    # fit non-linear proportional effect model
    nonlin.model <- nlme(Ysim ~ PcLDA(month, tx, Intercept, 
      beta.month, theta) + b0,
      data = data %>% filter(scenario == scen),
      fixed = Intercept + beta.month + theta ~ 1,
      random = list(b0 ~ 1),
      groups = ~ id,
      start = start1,
      control = nlmeControl(returnObject=TRUE))
    
    # extract proportional effect and p-value
    nonlin.result <- summary(nonlin.model)$tTable %>%
      as.data.frame() %>%
      rownames_to_column('term') %>%
      bind_cols(intervals(nonlin.model)$fixed) %>%
      filter(term == 'theta') %>%
      mutate(Model = 'Proportional effect') %>%
      select(-Value)
    
    # bundle and return results
    linear.result %>%
      bind_rows(nonlin.result) %>%
      mutate(scenario = scen)
  }) %>% bind_rows()
}

## Run longitudinal simulations ----
# fix seeds for each of 10,000 simulated trials
set.seed(20241101)
SEEDS <- sample(1:2e9, size = 10000, replace = FALSE)

long.results <- mclapply(1:10000, mc.cores=10, function(sim){
  # set seed and simulate residuals and random intercepts
  set.seed(SEEDS[sim])
  dd <- tibble(
    id = 1:N,
    rand.int = rnorm(n=N, sd=2)) %>%
    right_join(ind.long.means, by = 'id') %>%
    mutate(
      residual = rnorm(n=nrow(ind.long.means), sd=1.5),
      Y.power = Y + rand.int + residual,
      Y.null = Y0 + rand.int + residual
    )

  # fit models
  power.results <- return.model.results(dd %>% rename(Ysim = Y.power)) %>%
    mutate(Effect = 'Power')
  null.results <- return.model.results(dd %>% rename(Ysim = Y.null)) %>%
    mutate(Effect = 'Type I error')
  bind_rows(power.results, null.results) %>%
    mutate(simulation = sim)
}) %>%
  bind_rows()

## Summarize results ----
long.results <- long.results %>%
  mutate(
    Significance = factor(`p-value` < 0.05,
      levels = c(TRUE, FALSE), labels = c('p<0.05', 'p>0.05')),
    Scenario = factor(scenario, levels = 1:4,
      labels = paste("Scenario", c("A", "B", "C", "D"))))


## Table 1. Simulated power, Type I error, and proportion of rejections ----
# Simulated power, Type I error, and proportion of rejections under the null 
# that favor active treatment from 10,000 simulated trials under the mean trends
# shown in Figure 3
reject.tx <- long.results %>%
  filter(Effect == 'Type I error') %>%
  filter(!is.na(`p-value`) & `p-value`<0.05) %>%
  group_by(Scenario, Model) %>%
  summarise(
    `Proportion of rejections in favor of treatment under the null (%)` = sum(lower > 0)/length(`p-value`)*100)

long.results %>%
  group_by(Scenario, Model, Effect) %>%
  filter(!is.na(`p-value`)) %>%
  summarise(
    Percent = sum(`p-value`<0.05)/length(`p-value`)*100,
    `N sims` = length(`p-value`)) %>%
  pivot_wider(values_from = Percent, names_from = 'Effect') %>%
  left_join(reject.tx, by = c("Scenario", "Model")) %>%
  rename(`Power (%)` = Power, `Type I error (%)` = `Type I error`) %>%
  select(-`N sims`) %>%
  print(comment = FALSE, include.rownames=FALSE)

## Figure 4. Zipper plots for longitudinal power ----
# Zipper plots [@morris2019using] showing estimates and 95% confidence intervals
# for longitudinal power simulations. The intervals are sorted so that those 
# associated with the largest standardized bias are toward the top of each panel, 
# and only the 25% of simulations with the largest bias are shown. While linear 
# model estimates are symmetric about the true effect (left), the plots for the 
# proportional effect model (right) reveal bias and asymetrical confidence 
# interval widths. Confidence intervals also appear to be narrower when 
# proportional effect estimates favor treatment. Vertical dashed lines are at 
# the true value. Horizontal dashed lines are at the target rejection rate of 5\\%
pd <- long.results %>%
  filter(Effect == 'Power') %>%
  group_by(Scenario, Model) %>%
  arrange(`p-value`) %>%
  mutate(Model = paste(Model, "model"))

p1 <- pd %>% filter(Model == 'Linear model') %>%
  mutate(`Standardized bias` = abs(`est.` - slope_effect)/`Std.Error`) %>%
  group_by(Scenario) %>%
  arrange(`Standardized bias`) %>%
  mutate(
    rank = row_number()/length(`p-value`)*100,
    Coverage = case_when(
      lower > slope_effect | upper < slope_effect ~ 'Not covered',
      TRUE ~ 'Covered') %>%
      factor(levels = c('Not covered', 'Covered'))) %>%
  filter(rank>=75) %>%
ggplot(aes(x=`est.`, y=rank, color = Coverage)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=lower, xmax=upper), alpha = 0.1) +
  facet_wrap(vars(Scenario), scales = 'free_x', ncol = 1) +
  ylab("Standardized bias percentile") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position='inside',
    legend.position.inside = c(0.95, 0.06),
    legend.title = element_blank()) +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 95, linetype = 'dashed') +
  geom_vline(xintercept = slope_effect, linetype = 'dashed') +
  ggtitle('Linear model')

maxplot <- 20
p2 <- pd %>% filter(Model == 'Proportional effect model') %>%
  mutate(
    true_theta = case_when(
      scenario == 1 ~ 1/3,
      scenario == 2 ~ 2/3,
      scenario == 3 ~ 1,
      scenario == 4 ~ 2),
    `Standardized bias` = abs(`est.` - true_theta)/`Std.Error`,
    `est.` = case_when(abs(`est.`) > maxplot ~ NA, TRUE ~ `est.`),
    upper = case_when(upper > maxplot ~ maxplot, TRUE ~ upper),
    lower = case_when(lower < -maxplot ~ -maxplot, TRUE ~ lower)) %>%
  group_by(Scenario) %>%
  arrange(`Standardized bias`) %>%
  mutate(
    rank = row_number()/length(`p-value`)*100,
    Coverage = case_when(
      lower > true_theta | upper < true_theta ~ 'Covered',
      TRUE ~ 'Not covered'
    )) %>%
  filter(rank>=75) %>%
  ggplot(aes(x=`est.`, y=rank, color = Coverage)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=lower, xmax=upper), alpha = 0.1) +
  facet_wrap(vars(Scenario), scales = 'free_x', ncol = 1) +
  ylab("Standardized bias percentile") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position='none') +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 95, linetype = 'dashed') +
  geom_vline(aes(xintercept = true_theta), linetype = 'dashed') +
  ggtitle('Proportional effect model')

gridExtra::grid.arrange(p1, p2, ncol=2)

## Figure 5. Zipper plots for for longitudinal Type I error ----
# Zipper plots [@morris2019using] showing estimates and 95% confidence intervals
# for longitudinal Type I error simulations. The intervals are sorted so that 
# those associated with the smallest p-values are toward the top of each panel, 
# and only the smallest 25% of p-values are shown. While linear model estimates 
# are symmetric about the zero (left), the plots for the proportional effect 
# model (right) reveal bias in favor of treatment with a disproportionate number
# of positive dots. Confidence intervals also appear to be narrower when 
# proportional effect estimates favor treatment. Vertical dashed lines are at 
# the true value, zero. Horizontal dashed lines are at the target rejection rate
# of 5\\%
pd <- long.results %>%
  filter(Effect == 'Type I error') %>%
  group_by(scenario, Model) %>%
  arrange(`p-value`) %>%
  mutate(
    rank = row_number()/length(`p-value`)*100, 
    Model = paste(Model, "model")) %>%
  filter(rank<=25)

p1 <- pd %>% filter(Model == 'Linear model') %>%
ggplot(aes(x=`est.`, y=rank, color = Significance)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=lower, xmax=upper), alpha = 0.1) +
  facet_wrap(vars(Scenario), scales = 'free_x', ncol = 1) +
  scale_y_reverse() +
  ylab("P-value percentile") +
  xlab("Estimated difference in slopes") +
  theme_minimal() +
  theme(legend.position='none', 
    plot.margin =  unit(c(0.1,0.1,1,0.1), "in")) +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ggtitle('Linear model')  +
  geom_segment(
    data = tibble(
      arrowx1 = -0.03, arrowx2 = 0.03, 
      arrowy1 = 43, arrowy2 = 43, Scenario = 'Scenario D'),
    color = 'black',
    aes(x = arrowx1, y = arrowy1, xend = arrowx2, yend = arrowy2),
    arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(
    data = tibble(
      x = 0, y = 49, Scenario = 'Scenario D', text = "Favors active group"),
    aes(y=y, x=x, label = text), color = 'black') +
  coord_cartesian(ylim = c(25, 0), clip = "off") +
  theme(plot.margin = unit(c(0,0,2,0), "lines"))

maxplot <- 20
p2 <- pd %>% filter(Model == 'Proportional effect model') %>%
  mutate(
    `est.` = case_when(abs(`est.`) > maxplot ~ NA, TRUE ~ `est.`),
    upper = case_when(upper > maxplot ~ maxplot, TRUE ~ upper),
    lower = case_when(lower < -maxplot ~ -maxplot, TRUE ~ lower)) %>%
  ggplot(aes(x=`est.`, y=rank, color = Significance)) +
  geom_point(alpha=0.1) +
  geom_errorbar(aes(xmin=lower, xmax=upper), alpha = 0.1) +
  scale_y_reverse() +
  facet_wrap(vars(Scenario), scales = 'free_x', ncol = 1) +
  ylab("") +
  xlab("Proportional effect estimate") +
  theme_minimal() +
  theme(legend.position='none') +
  ggsci::scale_color_nejm() +
  geom_hline(yintercept = 5, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ggtitle('Proportional effect model')  +
  geom_segment(
    data = tibble(
      arrowx1 = -5, arrowx2 = 15, 
      arrowy1 = 43, arrowy2 = 43, Scenario = 'Scenario D'),
    color = 'black',
    aes(x = arrowx1, y = arrowy1, xend = arrowx2, yend = arrowy2),
    arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(
    data = tibble(
      x = 5, y = 49, Scenario = 'Scenario D', text = "Favors active group"),
    aes(y=y, x=x, label = text), color = 'black') +
  coord_cartesian(ylim = c(25, 0), clip = "off") +
  theme(plot.margin = unit(c(0,0,2,0), "lines"))

gridExtra::grid.arrange(p1, p2, ncol=2)


