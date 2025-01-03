# Linear-Birth-Death-Process-with-Catastrophic-Extinction

Topic: Linear Birth-Death Process with Catastrophic Extinction: Analytical Derivations, Parameter Estimation, and Efficient Simulation

# Abstract: 

This study analytically examines the linear birth-death process with catastrophic extinction (B-D-C process), a stochastic model that integrates catastrophic events leading to sudden population extinction into the conventional birth-death framework. Despite its simplicity, deriving the exact transition function and performing robust parameter estimation pose significant challenges, particularly for discretely observed processes. For the first time, we provide the exact derivation and numerical validation of the analytical transition probability function and theoretical moments for the linear B-D-C process. We propose and assess three competing parameter estimation methods: likelihood- and moment-based parameter estimations as well as a stochastic parameter estimation approach while examining the trade-off between computational speed and accuracy under different *in silico* simulation experiments using the exact stochastic simulation algorithm (SSA) for the B-D-C process. Additionally, we present a computationally efficient Monte Carlo simulation technique, leveraging a hybrid tau-leaping method, to accelerate the simulation of the B-D-C process compared to the exact SSA. The methods and findings have wide-ranging applications, particularly in modelling host-parasite systems, where sudden parasite population collapse can occur due to host immune responses, host mortality, or other extinction mechanisms.

**The R codes are attached for reproducibility of the study's results.**

# Below are the labels of the main programming codes/scripts of the study:

1. `BDC-process-Numerical-Proof.r`: R Codes for Numerical Validation of the Derived Transition Probability Function
2. `MLE_BDC.r`: R Codes for implementing Maximum Likelihood Estimation (MLE) of the B-D-C process
3. `BDC_Loglik.jl`: Julia function for evaluating the log-likelihood function for MLE of the B-D-C process
4. `BDCfit.jl`: Julia function for fitting the log-likelihood function (to speed up computational time)
5. `MLE-LogLik-BDC_JuliaFunc.r`: R coupled with Julia codes to speed up MLE of the B-D-C process 
6. `GMM-Estimation-BDC-process.R`: R Codes for implementing Generalised Method of Moments (GMM) Estimation of the B-D-C process
7. `BDC-estimation-GW-process.html`: R codes (from a Jupyter Notebook HTML file) for implementing embedded Galton-Watson (GW) Estimation of the B-D-C process
8. `BDC-LeapSize_StateX_Relationship.r`: R codes to explore the leap size for the two proposed Tau-leaping algorithms of the study
9. `Hybrid-Tauleaping-Algorithms-BDC.r`: Main R codes for running the exact SSA and the two Tau-leaping algorithms 
