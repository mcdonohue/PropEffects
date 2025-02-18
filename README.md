# Proportional effect models for continuous outcomes

Longitudinal models that assume Gaussian residuals with proportional treatment group means provide direct estimates of the proportional treatment effect in randomized clinical trials. However, these models often make a strong assumption about a fixed proportional treatment effect over time, which can lead to bias and Type I error inflation. We demonstrate that these models are biased and sensitive to the labeling of treatment groups, even when the proportional hazards assumption holds true. Typically, this bias favors the active group and inflates Type I error.

This repository includes R code for the simulations described in [Donohue et al., 2025](https://doi.org/10.48550/arXiv.2502.00214)