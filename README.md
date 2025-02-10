# Survival Analysis Model: From Unidimensional Clustering to Multidimensional Enhancements and Theoretical Integration

This project documents the development of an innovative survival analysis model that evolves from a simple unidimensional clustering approach to a sophisticated multidimensional framework. By integrating theoretical concepts such as MACMAL, dynamic programming for segmentation, and advanced optimization techniques, the model aims to better capture real-world data dynamics and uncertainty in survival prediction.

## Table of Contents

- [Introduction](#introduction)
- [Methodology](#methodology)
  - [Model Fitting and Selection](#model-fitting-and-selection)
  - [Unidimensional Clustering (V0 & V0.5)](#unidimensional-clustering-v0--v05)
- [Model Versions](#model-versions)
  - [V0 – The Original Unidimensional Model](#v0--the-original-unidimensional-model)
  - [V0.5 – Enhanced Unidimensional Model](#v05--enhanced-unidimensional-model)
  - [V1.0 – Multidimensional Clustering](#v10--multidimensional-clustering)
- [Multidimensional Enhancements](#multidimensional-enhancements)
  - [Dynamic Programming for Optimal Segmentation](#dynamic-programming-for-optimal-segmentation)
  - [Incorporating Multiple Covariates](#incorporating-multiple-covariates)
- [Model Averaging and Uncertainty Quantification](#model-averaging-and-uncertainty-quantification)
- [Current Results](#current-results)
- [Challenges and Future Work](#challenges-and-future-work)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Introduction

This project was inspired by an in-depth conversation with Professor Townsend, who highlighted the need for a two-dimensional graphical model that better reflects the asymmetric and dynamic nature of real-world survival data. The model challenges the common assumptions of symmetric or normally distributed data and explores the potential of borrowing methodologies from gene analysis to improve survival predictions.

---

## Methodology

### Model Fitting and Selection

- **Maximum Likelihood Estimation (MLE):** Employed for robust parameter estimation.
- **Information Criteria:** Utilizes AIC, BIC, and AICc to select the best model.
- **Divide-and-Conquer Strategy:** Recursively segments the observation period into distinct time intervals.
- **Model Averaging:** Integrates model uncertainty by weighting candidate models to produce reliable death rate estimates and corresponding confidence intervals at each time point.

### Unidimensional Clustering (V0 & V0.5)

- **V0:** The initial approach applies a single multidimensional Logistic regression over the entire observation period.
  - *Limitation:* Assumes a constant risk across time, which can be unrealistic.
  
- **V0.5 Enhancements:**
  - **Integration of MACMAL Theory:** Embeds theoretical insights to manage complex data distributions.
  - **Support for Multiple Distributions:** Extends from Binomial to Poisson and Negative Binomial, adapting to various data types.
  - **Enhanced Model Averaging:** Provides more accurate weight calculations and confidence interval estimation.
  - **Parallel Computing:** Incorporates OpenMP to accelerate model search and parameter estimation.

---

## Model Versions

### V0 – The Original Unidimensional Model

- **Approach:** Uses a logistic regression model to estimate event probabilities over a homogeneous time interval.
- **Shortcomings:**
  - Fails to capture dynamic risk changes over time.
  - Provides fixed risk estimates without uncertainty quantification.
  - Potentially overfits small or noisy datasets.

### V0.5 – Enhanced Unidimensional Model

- **Improvements:**
  - Integrates MACMAL principles for complex covariate relationships.
  - Supports multiple statistical distributions.
  - Optimizes computation with parallel processing.
  - Introduces confidence interval calculations alongside model averaging.
- **Outcome:** While overcoming several limitations of V0, the unidimensional structure still limits the model’s adaptability.

### V1.0 – Multidimensional Clustering

- **Key Features:**
  - **Dynamic Segmentation:** Divides the observation period into multiple intervals with unique risk profiles.
  - **Multidimensional Covariates:** Considers several influencing factors (e.g., patient age, treatment protocol) simultaneously.
  - **Optimization via Gradient Descent:** Uses enhanced optimization strategies with convergence criteria and adaptive learning rate adjustments.
  - **Minimum Segment Length Constraint:** Ensures stability by avoiding excessively short segments.

---

## Multidimensional Enhancements

### Dynamic Programming for Optimal Segmentation

To efficiently identify segmentation points as the number of segments grows, a dynamic programming (DP) approach is used. A DP table (e.g., `dp[k][s]`) is constructed where:
- **dp[k][s]:** Represents the optimal AIC/BIC value and corresponding log-likelihood when segmenting data starting at time point `s` into `k` segments.
- **Recurrence Relation:**  
  The cost function `Cost(s, s')` evaluates the log-likelihood or AIC/BIC over the segment between time points `s` and `s'`, enabling an iterative build-up of the optimal segmentation scheme.

### Incorporating Multiple Covariates

- **Extended Logistic Regression:**  
  The model is extended to include multiple covariates (e.g., age, gender, treatment protocols), enhancing interpretability and predictive performance.
- **Gradient Descent Optimization:**  
  Implements convergence criteria and learning rate adjustments to efficiently estimate regression coefficients within each segment.

---

## Model Averaging and Uncertainty Quantification

- **Akaike Weights:**  
  Weights for each candidate segmentation model are computed based on AIC/BIC to reflect their relative support.
- **Weighted Averaging:**  
  Final event probability estimates are produced by taking a weighted average across candidate models. This helps mitigate fluctuations and improves robustness, especially near segmentation boundaries.

---

## Current Results

- **V0 (Unidimensional, Origin):**  
  Generates initial graphics and analyses based on a single risk assumption.
- **V0.5 (Enhanced Unidimensional):**  
  Offers improved visualizations and more reliable analysis by incorporating multiple distributions and parallel computing.
- **V1.0 (Multidimensional):**  
  Provides detailed segmentation and risk profiling, capturing dynamic changes in survival data.

*Note: Some results are still under testing and refinement. Future updates will include further predictions and accuracy assessments.*

---

## Challenges and Future Work

### Current Challenges

- **Handling Missing/Censored Data:**  
  The current model assumes complete data, which is a limitation in survival analysis.
- **Dependence on Specific Distribution Assumptions:**  
  Restricting the model to certain distributions may limit flexibility.
- **Visualization and Interpretability:**  
  Existing tools could be improved to provide more intuitive insights.
- **Predictive Accuracy:**  
  Extensive testing across diverse datasets is necessary to validate robustness.

### Potential Improvements

- **Bayesian and Semi-Parametric Integration:**  
  Combining segmentation clustering with Bayesian methods for a more adaptive framework.
- **Nonlinear and Flexible Structures:**  
  Incorporating non-parametric relationships to better capture complex data patterns.
- **Censored Data Support:**  
  Enhancing the model to handle incomplete data.
- **High-Dimensional Feature Selection:**  
  Using advanced techniques to manage and select relevant variables.
- **Spatio-Temporal Extensions:**  
  Expanding to incorporate multidimensional interactions in both space and time.
- **Advanced Optimization:**  
  Further refining algorithms and visualization tools to improve usability and efficiency.
- **Deep Learning Integration:**  
  Exploring deep learning methods for high-dimensional and complex nonlinear relationships.

---
