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

**GOAL OF PAPER: investigate the ’properness’ of methods with diagnostics for multivariate model parameters.**

**Secondary (optional): When to use which method?** 

Create outline for showing properness in MC simulation study.

**Project Plan:**


Why investigate properness? => predictions do not equal right inference. Can't draw conclusions a out relationships from imputed data. 

Methods: 
- Mice: PMM, CART 
- MiceRanger: Forest .

Setup: 
- Algorithms (especially Mice)
- Differences (especially Miceranger)

Theory: 
- Why properness? 
- Why methods unproper? 


DGP: 
- MAR (Missing at Random)
- Simulation 

Diagnostics: (properness)
- Between Chain Variance = B 
- If B ~ 0 -> improper 
- Coverage Rate => If good coverage = good (95%) pooling (proper, if lower = Bad 
- True Bias 
- FoMi (Fraction of Missing) RF and CART often give low rates 

Why?: 
A multiple imputation method is proper (Rubin, 1987) if:
After pooling, the MI estimator has the same repeated-sampling distribution as the estimator based on the observed-data posterior.

Simulates the DGP! 
