
**GOAL OF PAPER:**

investigate the ’properness’ of methods with diagnostics for multivariate model parameters.**

**Secondary (optional): When to use which method?** 

# Project Plan:

Why investigate properness? => predictions do not equal right inference. Can't draw conclusions about relationships from imputed data. 

**1) Methods:**
- Mice: PMM, CART 
- MiceRanger: Forest .

**2) Setup:**
- Algorithms (especially Mice)
- Differences (especially Miceranger)

**3) Theory:**
- Why properness? 
- Why methods unproper? 

**4) DGP:**
- MAR (Missing at Random)
- Simulation 

**5) Running Imputation**
- Convergence 
- Model Checks 

**6) Diagnostics: (properness)**
- Between Chain Variance = B 
- If B ~ 0 -> improper 
- Coverage Rate => If good coverage = good (95%) pooling (proper, if lower = Bad (indicator)
- B = Main Judge of Properness
- B => variability across plausible completions of the missing data

Secondary Importance: 
- True Bias 
- FoMi (Fraction of Missing) RF and CART often give low rates 

**7) Conclusion:**
- Why focus on properness? 
- When are CART and RF still better choices? 

**Why?:**
A multiple imputation method is proper *(Rubin, 1987)* if:
After pooling, the MI estimator has the same repeated-sampling distribution as the estimator based on the observed-data posterior.

Simulates the DGP! 

# Question:

Tree-based Multiple Imputation methods:
MI methods based on machine learning algorithms (see e.g. Burgette &
Reiter 2010) have become increasingly popular, although they are not
considered proper MI methods
Objectives:

(a) Describe and discuss the miceranger algorithm. Compare methods
from the miceranger package with methods ’pmm’ and ’cart’ in mice.
Conduct a simulation study and investigate the ’properness’ of these
methods by focusing on diagnostics for multivariate model parame-
ters.
[Packages and vignettes available. Requires implementation of an MC
study.]
