---
output:
  pdf_document: default
  html_document: default
---
GPCW: Gaussian Process for Estimating Crtical Windows of Susceptibility

##Statistical Model
$$y_i|\boldsymbol{\beta}, \boldsymbol{\theta} \stackrel{\text{iid}}{\sim} \text{Bernoulli}\left(p_i\left(\boldsymbol{\beta}, \boldsymbol{\theta}\right)\right),\ i=1,...,n;$$

$$\log\left(\frac{p_i\left(\boldsymbol{\beta}, \boldsymbol{\theta}\right)}{1 - p_i\left(\boldsymbol{\beta}, \boldsymbol{\theta}\right)}\right) = \textbf{x}_i^{\text{T}} \boldsymbol{\beta} + \sum_{j=1}^{m_i} \text{z}_{ij} \theta\left(j\right);$$

$$\boldsymbol{\theta}=\left(\theta\left(1\right), ..., \theta\left(m\right)\right)^{\text{T}}| \sigma^2_{\theta}, \phi \sim \text{MVN}\left(\boldsymbol{0}_m, \sigma^2_{\theta}\Sigma\left(\phi\right)\right)$$

* $p$: Length of $\textbf{x}_i$ vector (same for all $i$)

* $m = \max\left\{m_i: i=1,...,n\right\}$

##Prior Information
$\beta_j \stackrel{\text{iid}}{\sim}\text{N}\left(0, \sigma^2_{\beta}\right),\ j=1,...,p$

* Default setting: $\sigma^2_{\beta} = 10,000$

$\sigma^2_{\theta} \sim \text{Inverse Gamma}\left(\alpha, \beta\right)$

* Default setting: $\alpha = 3$, $\beta = 2$

$\phi \sim \text{Uniform}\left(a_{\phi}, b_{\phi}\right)$

* Default setting: $a_{\phi} = d$, $b_{\phi} = d$

##Default Initial Values
$\beta_j = 0$ for all $j$

$\theta_j = 0$ for all $j$

$\sigma^2_{\theta} = 0.50$

$\phi = 0.50$

