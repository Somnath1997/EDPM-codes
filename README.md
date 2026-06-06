# EDPM-codes

R code for implementing Variational Bayes to get the truncation and calculate the cluster values, and then using the cluster values, fitting an EDPM model. Includes methods proposed in the paper "Variational Bayes and Truncation Approximation for Enriched Dirichlet process mixtures: by Somnath Bhadra and Michael J. Daniels. Please see the arxiv version of the paper at https://arxiv.org/abs/2603.12427 to use the codes.


For the initial prior distributions, one can change the parameters (mean and sd) to one's likings.

To run the codes, first run the VB code to calculate the Cluster values for the EDPM truncation. Then use those cluster values as K and L accordingly to plug in and  run the EDPM code.
