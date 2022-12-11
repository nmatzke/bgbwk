
#######################################################
# Example simulation and SSEsim
#######################################################
# Load the package (after installation, see above).
library(optimx)         # You need to have some version of optimx available
                        # as it is a BioGeoBEARS dependency; however, if you
                        # don't want to use optimx, and use optim() (from R core) 
                        # you can set:
                        # BioGeoBEARS_run_object$use_optimx = FALSE
                        # ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)


#######################################################
# Set up the simulation parameters
#######################################################
wd = "/drives/GDrive/__GDrive_projects/2020-12-06_Transmitting_Science/lab03/sim/"
setwd(wd)

# Set up simulation parameters
# dvals = c(0, 0.03, 0.15)
# evals = c(0, 0.03, 0.15)
# jvals = c(0, 0.02, 0.1, 0.3)
# brates_sim = 0.3
# drates_sim = c(0, 0.1, 0.3)
# b_exp = c(0, 1)
# d_exp = c(0, -1)

# Set up simulation parameters
dvals = c(0, 0.03)
evals = c(0, 0.03)
jvals = c(0, 0.02, 0.1)
brates_sim = 0.3
drates_sim = c(0, 0.1, 0.3)
b_exp = c(0, 1)
d_exp = c(0, -1)


dej_params = expand.grid(dvals, evals, jvals, brates_sim, drates_sim, b_exp, d_exp)
names(dej_params) = c("d", "e", "j", "brate", "drate", "b_exp", "d_exp")

keepTF = dej_params$d >= dej_params$e
dej_params = dej_params[keepTF, ]

# Remove all zeros
keepTF = (dej_params$d + dej_params$e + dej_params$j) > 0
dej_params = dej_params[keepTF, ]
dej_params

tail(dej_params)

# Remove all zeros
tmpadd = ((dej_params$b_exp == 0) + (dej_params$d_exp == -1))
tail(tmpadd)

remove_TF = tmpadd == 2
keepTF = remove_TF == FALSE
dej_params = dej_params[keepTF, ]
dim(dej_params)


tmpadd = (dej_params$b_exp == 1) + (dej_params$d_exp == 0)
remove_TF = tmpadd == 2
keepTF = remove_TF == FALSE
dej_params = dej_params[keepTF, ]

# Set up the seed
seedval = 54321
seedres = set.seed(seedval)

dim(dej_params)

save(dej_params, file="dej_params.Rdata")




#######################################################
# Set up the simulation inputs
#######################################################
# 100 sims per parameter combination
numsims = 3

simspath = "/drives/GDrive/__GDrive_projects/2020-12-06_Transmitting_Science/lab03/sim/sims/"
param_iter = 1
simnum = 1
simcount = 0

num_cores_to_use = 1
if (.Platform$GUI != "AQUA" && num_cores_to_use > 1)
	{
	cluster_already_open = makeCluster(rep("localhost",num_cores_to_use), type = "SOCK")
	cat("Started cluster with ", num_cores_to_use, " cores.\n\n", sep="")
	} else {
	cluster_already_open = NULL
	}
# This sets the recursion limit
options("expressions") # default 5000
options(expressions=15000) # default 5000



#######################################################
# BEGIN THE SIMULATION FOR-LOOP
#######################################################
# There are only 6 parameter combinations here
# in this _basic_example
# for (param_iter in 6:6)

for (param_iter in 1:1)
#for (param_iter in 1:nrow(dej_params))
	{
	for (simnum in 1:numsims)
		{
		simcount = ((param_iter-1) * numsims) + simnum
		simcount

		#######################################################
		# Set up simulation directory
		#######################################################
		# make e.g. 0001
		ps_count_txt = sprintf("%04.0f", param_iter)
		simcount_txt = sprintf("%03.0f", simnum)

		simdir = paste(simspath, "ps", ps_count_txt, "_sim", simcount_txt, sep="")
		simdir


		txt = paste("\nSimulation #", simcount, " in ", simdir, "/", "\n", sep="")
		cat(txt)


		dir.create(path=simdir, showWarnings=FALSE)
		setwd(simdir)
		list.files(simdir)


		#######################################################
		# Put parameters into simulation
		#######################################################
		# Starter -- don't use
		SSEsim_inputs = SSEsim_setup_inputs()

		# Set up the BD / SSE parameters
		SSEsim_inputs$SSEmodel$brate = dej_params$brate[param_iter]
		SSEsim_inputs$SSEmodel$drate = dej_params$drate[param_iter]
		SSEsim_inputs$SSEmodel$rangesize_b_exponent = dej_params$b_exp[param_iter]
		SSEsim_inputs$SSEmodel$rangesize_d_exponent = dej_params$d_exp[param_iter]


		print(param_iter)
		print(dej_params[param_iter,])


		# Set up the DEC+J parameters
		SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dej_params$d[param_iter]
		SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dej_params$d[param_iter]

		SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = dej_params$e[param_iter]
		SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = dej_params$e[param_iter]

		SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = dej_params$j[param_iter]
		SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = dej_params$j[param_iter]

		# You gotta re-calculate with the new params
		SSEsim_inputs = SSEsim_setup_inputs(SSEmodel=SSEsim_inputs$SSEmodel, BioGeoBEARS_run_object=SSEsim_inputs$BioGeoBEARS_run_object)
		names(SSEsim_inputs)


		# Do a simulation
		#testwd = "/drives/SkyDrive/_________thesis/_doc2/ch2_submission/2014-02-11_reviews/testsim/"
		SSEsim_results = SSEsim_run(SSEsim_inputs, time_stop=100, taxa_stop=50, seed=simcount, printlevel=0, testwd=simdir)
		SSEsim_results

		#names(SSEsim_results)
		#SSEsim_results$trynum
		#SSEsim_results$simtr
		#dej_params[param_iter,]
		#plot(SSEsim_results$simtr)

		# Process the results
		SSEsim_results_processed = SSEsim_to_files(SSEsim_results, simdir=simdir, fossils_older_than=0.001, printlevel=0)
		SSEsim_results_processed
		names(SSEsim_results_processed)
		SSEsim_results_processed$tipranges_observed

		save(SSEsim_results_processed, file="SSEsim_results_processed.Rdata")


		dim(SSEsim_results_processed$table_of_cladogenetic_events_translated)
		dim(SSEsim_results_processed$table_of_cladogenetic_events_observed)
		dim(SSEsim_results_processed$table_of_range_change_events_translated)


		# Plot the results
		# plot(SSEsim_results$simtr, label.offset=0.2)
		# axisPhylo()
		# title("SSE sim, with species labels")
		# nodelabels()
		# tiplabels()






		#######################################################
		# Inference
		#######################################################
		max_range_size = 4

		#######################################################
		# Run DEC
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = FALSE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=4
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx=TRUE
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		runslow = TRUE
		resfn = "DEC_inf.Rdata"
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)
			resDEC = res
			} else {
			# Loads to "res"
			load(resfn)
			resDEC = res
			}

		#######################################################
		# Run DECj
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = FALSE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=4
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx=TRUE
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE

		# Set up DEC+J model
		# Add j as a free parameter
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01

		# Crash fix
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13


		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		resfn = "DECJ_inf.Rdata"
		runslow = TRUE
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)

			resDECj = res
			} else {
			# Loads to "res"
			load(resfn)
			resDECj = res
			}


		# Print the results
		txt = paste("\nSimulation finished. DEC LnL=", resDEC$total_loglikelihood, ", DEC+J LnL=", resDECj$total_loglikelihood, "\n", sep="")
		cat(txt)
		cat("Model params:\n")
		print(dej_params[param_iter,])



		} # end 100 simulations
	} # end dej_params


