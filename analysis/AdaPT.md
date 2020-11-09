# AdaPT

Given $n$ hypothesis,  If $R$ is the number of total rejections and $V$ is the number of false rejections (rejections of true null hypotheses), the FDR is 
$$
FDR = E(\frac{V}{R \vee 1}).
$$

## Overview

Suppose we have side information, $x_i$,  of each hypothesis. The proposed method iteratively proposes a threshold $s_t(x)$ and computes $\widehat{FDP}_t$, until its smaller than $\alpha$. Every $H_i$ for which $p_i\leq s_t(x_i)$ are rejected. 
$$
\text{for } t = 1,2,...
\\
\text{if } \widehat{FDP}_t \leq \alpha, \text{reject} \{H_i:p_i\leq s_t(x_i)\};
\\
\text{else } s_{t+1} = update\{x_i,\tilde p_{t,i},A_t,R_t,s_t\}.
$$


The estimator $\widehat{FDP}_t$ is 
$$
\widehat{FDP}_t = \frac{1+A_t}{R_t\vee 1},
$$


where $A_t = |\{i: p_i\geq 1-s_t(x_i)\}|$ and $R_t$ is the number of rejections. 

[Lei and Fithian (2016)](http://proceedings.mlr.press/v48/lei16.pdf) and [Arias-Castro and Chen (2017)](https://projecteuclid.org/euclid.ejs/1494900120) also used this estimator. Might worth a look. 

The threshold $s_t(x)$ should satisfies two constraints: 1. $s_{t+1}\preceq s_t$ for all $x$; 2. large and small p-values are partially masked. The available information when choosing $s_{t+1}$ is $A_t,R_t$ and $(x_i,\tilde p_{t,i})$ where 
$$
\begin{equation}
\tilde{p}_{t, i}=\left\{\begin{array}{ll}
p_i & s_{t}\left(x_{i}\right)<p_{i}<1-s_{t}\left(x_{i}\right) \\
\left\{p_{i}, 1-p_{i}\right\} & \text { otherwise }
\end{array}\right.
\end{equation}
$$
An ideal choice for $s_t(x)$ is a level surface of the local FDR in a two-groups empirical Bayes model.

## Assumptions

Let $H_0$ denote all the null hypotheses. 

1. Assume $(p_i)_{i\in H_0}$ are independent, and independent of $(p_i)_{i\notin H_0}$. 

2. Null p values are mirror conservative. A p-value is conservative if 
   $$
   P(p\in[a_1,a_2])\leq P(o\in [1-a_1,1-a_2]) , 0\leq a_1\leq a_2\leq 0.5.
   $$

## FDR control

Let $V$ be the number of false rejections at $t$ and $U_t = |\{i:i\in H_0,p_i\geq 1-s_t(x_i\}|$ . 

If the null p-values are uniform then $V_t\approx U_t$. 

At the final step $T$ when iterations stop, $FDP_T = \frac{V_T}{R_T\vee 1} = \frac{1+U_T}{R_T\vee 1}\frac{V_T}{1+U_T}\leq \alpha\frac{V_T}{1+U_T}$.  Only need to show $E(\frac{V_T}{1+U_T})\leq 1$. 

## Choose threshold



Two groups model 
$$
\begin{equation}
\begin{array}{c}
H_{i} \mid x_{i} \sim \text { Bernoulli }\left\{\pi_{1}\left(x_{i}\right)\right\} \\
p_{i} \mid H_{i}, x_{i} \sim\left\{\begin{array}{ll}
f_{0}\left(p \mid x_{i}\right) & \text { if } H_{i}=0 \\
f_{1}\left(p \mid x_{i}\right) & \text { if } H_{i}=1
\end{array}\right.
\end{array}
\end{equation}
$$


Assume $f_0(p|x)=1$ then $f(p|x) = 1-\pi_1(x)+\pi_1(x)f_1(p|x)$., and  the conditional lfdr is 
$$
fdr(p|x) = p(H_i \text{ is NULL}|x_i=x,p_i=p) = \frac{1-\pi_1(x)}{f(p|x)}.
$$
Note that $f_1(p|x)$  is continuously non-increasing  in $p$.  The conservative identifying assumption is 
$$
1-\pi_1(x) = f(1|x).
$$
Thus, a conservative estimate $\widehat{fdr}(p|x) = \frac{\hat f(1|x)}{\hat f(p|x)}$. 

The thresholding rule is considered optimal if it's the solution of 
$$
max_s \text{Power}(s) \hspace 1cm \text{subject to FDR}(s)\leq \alpha.    
$$
It can be shown that $fdr(s^*(x)|x) = \frac{1+\lambda\alpha}{1+\lambda}$ where $\lambda$ is Lagrangian multiplier. "The optimal thresholding rules are level surfaces of local FDR". 

### Estimate f(p|x)

Assume $f_1(p|x)$ is from exponential family, denoted as 
$$
h(p;\mu) = \exp\{\eta(\mu)g(p)-A(\mu)\}. 
$$
 The two-group model is 
$$
H_i|x_i\sim Bernoulli(\pi_{1i}), logit(\pi_{1i}) = \theta^T\phi_\pi(x_i),
$$


and  if $H_i=1$, 
$$
p_x|x_i,H_i\sim h(p;\mu_i), \zeta(\mu_i) = \beta^T\phi_\mu(x_i),
$$
where $\zeta$ is the canonical link function and $\phi$ is feature transformation. 

In the log-likelihood $l(\beta,\theta; p,H,x)$, some of the p-values(due to masking) and all H are unknown. AdaPT procedure maximizes the likelihood of the data $D_t=(x_i,\tilde p_{t,i})$. 

The EM algorithm for estimating $\beta,\theta$ proceed as 

1. update expectations $\hat H_i = E_{\theta,\beta}(H_i|D_i)$ and $\hat y_i = E_{\theta,\beta}(y_iH_i|D,H_i = 1)/\hat H_i$. 

2. maximize $l$ over $\beta,\theta$ : two GLM problems.  

   

How to deal with the masked p-values?

Let $p'_{t,i}$ denote the minimum element of $\tilde p_{t,i}$.  For masked p-values, 
$$
p(p_i = p'_{t,i}|\tilde p_{t,i},H_i = 1) = \frac{h(p'_{t,i};\mu_i)}{h(1-p'_{t,i};\mu_i) + h(p'_{t,i};\mu_i)}.
$$


And a similar one holds for $H_i = 0$. So $p(H_i = 1|\tilde p_{t,i}) = \frac{p(p_i=\tilde p_{t,i}|H_i = 1)p(H_i = 1)}{p(p_i = \tilde p_{t,i})}$. 



## Mask z scores

For masking p-values, we can introduce a latent binary variable $z_i$ and let the masked p-value be 
$$
\tilde p_i = p_i*z_i + (1-z_i)*(1-p_i). 
$$


Consider z scores $\hat\beta_i\sim N(\beta_i,1)$ , and $\beta_i\sim g(\cdot)$. 

If we want to mask z-scores, we can introduce a latent variable  $z_i$ and define 
$$
\tilde{\hat\beta}_i = \hat\beta_i*z_i + (1-z_i)*\frac{1}{\hat\beta_i}.
$$
The argument is similar to p-value case that if a z-score is very large(p is very small), we'd like it's counter part to be closer to null, i.e. a z-score closer to 0(a p closer to 1), since a very large z-score(small p value) is more likely to be a non-null( and we want to keep it a null). 

Another question is if we want to mask all z-scores. In AdaPT paper, authors mention that loss of information that is caused by partial masking is small target FDR drops to the ‘practical’ regime. But if we mask all z-scores, the loss of information would be much more. 



























