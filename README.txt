###########################################################
## Input to CVPAT function
###########################################################
- MV: A dataframe containing only the Manifest Variables (indicators)
      used for the PLS-SEM estimation. The variable names must be the
      same as specified in Model1 and Model2 (see below).

- CVFOLDS: The number of cross-validation folds

Model1: The model specifications as required by the matrixpls
	package. That is, Model1 is a list of three components.
	The first component is a matrix of the inner model,
	the second component a matrix of the reflective measurement model
	and the last component is a matrix of the formative measurement model.

Model2: The structure is the same as for Model1, but specifications is
	for Model2.

hypothesis: The hypothesis to test. If hypothesis = "M1_better_out_of_sample_than_M2"
	    then CVPAT will test if M1 has a significantly better out-of-sample
	    performance than M2. 
	    If hypothesis = "M1!=M2" then CVPAT will test if the out-of-sample
	    performance between M1 and M2 is significantly different from 
	    each other.		

BootSamp: The number of bootstrap resamples

boot.Di: If TRUE, bootstrapping will be conducted on the losses
	 from a single cross-validation. 
	 If FALSE, bootstrapping is on the manifest variables
	 and consequently PLS estimation and cross-validation will be
	 conducted for each bootstrapping run (this options is
	 considerably slower than boot.Di=TRUE)

seed: If FALSE, the random number generator will not be seeded inside
      the CVPAT function. If seed = x (where x is a number, e.g. seed = 42),
      the random number generator will be seeded with x inside the CVPAT
      function. The latter option can be used to replicate results.
      If seed = x, then the seed in the global environment will
      return to its original value after the CVPAT function has run.	



###########################################################
## Output from CVPAT function
###########################################################
The output from the CVPAT function is a list. The list has following components:
- boot.p.values
- losses
- t.stat
- p.value
- conf.int
- conv.fail

----------
The content of each list entry is:
----------

boot.p.values:
	- p.value.perc.t: p-value calculated from the bootstrapped percentiles 
			  of the t-statistics.
	- p.value.var.ttest: p-value calculated from the original t-statistic
		             but replacing the variance of D_bar with the 
			     boostrapped variance of D_bar
	- p.value.perc.D: p-value calculated from the bootstrapped percentiles
			  of D_bar.

losses$case_wise_loss:
	- LossM1: The loss for each case for model 1
	- LossM2: The loss for each case for model 2
	- LossM1_sepLV: The loss for each case on each endogenous construct for
			model 1 
	- LossM2_sepLV: The loss for each case on each endogenous construct for
			model 2

losses$avg_losses:
	- avg_losses_M1: The average loss overall for M1 and the average loss
			 for each endogenous construct in M1
	- avg_losses_M2: The average loss overall for M2 and the average loss
			 for each endogenous construct in M2

t.stat: The non-bootstrapped t-statistic

p.val: The non-bootstrapped p-value

conf.int: The non-bootstrapped confidence interval

conv.fail: Only relevant when boot.Di = FALSE. The proportion of bootstrap
	   runs where the PLS-SEM algorithm had convergence failure in one
	   or several of the cross-validation runs. 
