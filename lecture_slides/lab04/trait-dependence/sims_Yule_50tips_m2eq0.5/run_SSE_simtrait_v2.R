
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

library(TeachingDemos) # for char2seed()


########################################################
# TO GET THE OPTIMX/OPTIM FIX, AND THE UPPASS FIX, 
# SOURCE THE REVISED FUNCTIONS WITH THESE COMMANDS
#
# CRUCIAL CRUCIAL CRUCIAL: 
# YOU HAVE TO RUN THE SOURCE COMMANDS AFTER 
# *EVERY TIME* YOU DO library(BioGeoBEARS). THE CHANGES ARE NOT "PERMANENT", 
# THEY HAVE TO BE MADE EACH TIME.  IF YOU ARE GOING TO BE OFFLINE, 
# YOU CAN DOWNLOAD EACH .R FILE TO YOUR HARD DRIVE AND REFER THE source()
# COMMANDS TO THE FULL PATH AND FILENAME OF EACH FILE ON YOUR
# LOCAL SYSTEM INSTEAD.
########################################################
# library(BioGeoBEARS)
# source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
# source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
# calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
# calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
    # slight speedup hopefully

#######################################################
# Local source()-ing method -- uses BioGeoBEARS sourceall() function 
# on a directory of .R files, so you don't have to type them out.
# The directories here are on my machine, you would have to make a 
# directory, save the .R files there, and refer to them.
#
# NOTE: it's best to source the "cladoRcpp.R" update first, to avoid warnings like this:
##
## Note: possible error in 'rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs = tmpca_1, ': 
##         unused arguments (m = m, m_null_range = include_null_range, jts_matrix = jts_matrix) 
##
#
# TO USE: Delete or comment out the 'source("http://...")' commands above, and un-comment
#              the below...
########################################################################
# Un-comment (and fix directory paths) to use:
library(BioGeoBEARS)
source("/drives/Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp.R")
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
########################################################################







#######################################################
# Kristina's tree
#######################################################
library(ape)
trfn = "/drives/GDrive/__GDrive_projects/2016-12-07_Kristina_Klaus_Podocarpaceae/2017-05-05_Klaus_all/9Areas_5.0_Schneemann_M0_geog_traits_v2_m2eq0/Podocarpaceae_6markerMCC_158.newick"
tr = read.tree(trfn)

# Infer ML birth/death rates
BD =  birthdeath(tr)
BD

# Estimation of Speciation and Extinction Rates
#             with Birth-Death Models
# 
#      Phylogenetic tree: tr 
#         Number of tips: 158 
#               Deviance: -39.11031 
#         Log-likelihood: 19.55516 
#    Parameter estimates:
#       d / b = 0.8365634   StdErr = 0.05061399 
#       b - d = 0.01492895   StdErr = 0.003992572 
#    (b: speciation rate, d: extinction rate)
#    Profile likelihood 95% confidence intervals:
#       d / b: [0.7883772, 0.8736507]
#       b - d: [0.01193282, 0.0185437]

names(BD)
# [1] "tree" "N"    "dev"  "para" "se"   "CI"  
# Calculate the birthRate and deathRate from the outputs
x1 = unname(BD$para["d/b"])
x2 = unname(BD$para["b-d"])
deathRate = (x2*x1) / (1-x1)
birthRate = deathRate+x2
c(birthRate, deathRate)
# [1] 0.09134398 0.07641503

#plot(tr)
#axisPhylo()






#######################################################
# Set up the simulation parameters
#######################################################
#wd = "/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/examples_new/Psychotria_M0_unconstrained/BGB/"
master_wd = "/drives/GDrive/__GDrive_projects/2016-12-07_Kristina_Klaus_Podocarpaceae/_sims/sims_Yule_50tips_m2eq0.5/"
setwd(master_wd)

# Set up simulation directory
simsdir = slashslash(paste(master_wd, "/sims/", sep="", collapse=""))
dir.create(simsdir, showWarnings=FALSE)

# Set up simulation parameters
# dvals = c(0, 0.03, 0.15)
# evals = c(0, 0.03, 0.15)
# jvals = c(0, 0.02, 0.1, 0.3)
# brates_sim = 0.3
# drates_sim = c(0, 0.1, 0.3)
# b_exp = c(0, 1)
# d_exp = c(0, -1)

# Set up simulation parameters
dvals = c(0.01)
evals = c(0, 0.05)
jvals = c(0.1, 0)
brates_sim = c(0.04491782) #c(0.09)
drates_sim = c(0, 0.01) #c(0.07)
b_exp = c(0)
d_exp = c(0)

# Set up trait simulation parameters
m1 = c(1)
m2 = c(0.5, 0.1, 1)
#m2 = 1
t12 = c(0.005)
t21 = c(0.010)


# Without traits
#dej_params = expand.grid(dvals, evals, jvals, brates_sim, drates_sim, b_exp, d_exp)
#names(dej_params) = c("d", "e", "j", "brate", "drate", "b_exp", "d_exp")

# With traits
dej_params = expand.grid(dvals, evals, jvals, brates_sim, drates_sim, b_exp, d_exp, m1, m2, t12, t21)
names(dej_params) = c("d", "e", "j", "brate", "drate", "b_exp", "d_exp", "m1", "m2", "t12", "t21")

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

# Set up the seed, based on directory name
words = strsplit(master_wd, split="/")[[1]]
lastword = words[length(words)]
seedval = -543 + TeachingDemos::char2seed(lastword, set=FALSE)
seedres = set.seed(seedval)

dim(dej_params)

save(dej_params, file="dej_params.Rdata")




#######################################################
# Set up the simulation inputs
#######################################################
# 100 sims per parameter combination
numsims = 100


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

#trait_fn = NULL	# for no traits
trait_fn = 2			# for a randomly-generated traits file with 2 states
#trait_fn = fn			# for a given traits file


param_iter = 1
simnum = 1
#for (param_iter in 1:nrow(dej_params))
 for (param_iter in 1:1)
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

		simdir = slashslash(paste(master_wd, "/sims/", "ps", ps_count_txt, "_sim", simcount_txt, "/", sep=""))
		simdir


		txt = paste("\nSimulation #", simcount, " in ", simdir, "/", "\n", sep="")
		cat(txt)


		dir.create(path=simdir, showWarnings=FALSE)
		setwd(simdir)
		list.files(simdir)


		#######################################################
		# Put parameters into simulation
		#######################################################
		# Set up the BD / SSE parameters
		SSEmodel = NULL
		SSEmodel$brate = dej_params$brate[param_iter]
		SSEmodel$drate = dej_params$drate[param_iter]
		SSEmodel$rangesize_b_exponent = dej_params$b_exp[param_iter]
		SSEmodel$rangesize_d_exponent = dej_params$d_exp[param_iter]
		SSEmodel$dej_params = dej_params

		print(param_iter)
		print(dej_params[param_iter,])

		sim_function <- function(dej_params, param_iter, number_of_trait_states, trait_fn, numsims, simcount)
			{
			# Set up the BioGeoBEARS_run_object -- DEFAULTS, OR
			# SIMILAR TO AN ACTUAL ANALYSIS!! 
			# e.g. DEC+J parameters
			BioGeoBEARS_run_object = define_BioGeoBEARS_run()
			BioGeoBEARS_run_object$max_range_size = 4
		
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dej_params$d[param_iter]
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dej_params$d[param_iter]

			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = dej_params$e[param_iter]
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = dej_params$e[param_iter]

			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = dej_params$j[param_iter]
			BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = dej_params$j[param_iter]

			# You gotta re-calculate with the new params;
			# And, add trait
			# Input the simulation values of the trait parameters, if needed
			number_of_trait_states = 2
			trait_fn=number_of_trait_states
			SSEsim_inputs = SSEsim_setup_inputs(SSEmodel=SSEmodel, BioGeoBEARS_run_object=BioGeoBEARS_run_object, trait_fn=trait_fn)
		
			SSEsim_inputs$COO_weights_columnar

			# Run the inputs now that the BioGeoBEARS_run_object is complete
			SSEmodel=SSEsim_inputs$SSEmodel
			BioGeoBEARS_run_object=SSEsim_inputs$BioGeoBEARS_run_object
			BioGeoBEARS_run_object$trait_Pmat_txt

			names(SSEsim_inputs)
			names(SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object)

			SSEsim_inputs$BioGeoBEARS_run_object$trait_Pmat_txt
			SSEsim_inputs$Qmat
			tail(SSEsim_inputs$Qmat)
			SSEsim_inputs$Qmat[14:20,]
			SSEsim_inputs$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
			SSEsim_inputs$COO_weights_columnar
		
		
		
			# Do a simulation
			#testwd = "/drives/SkyDrive/_________thesis/_doc2/ch2_submission/2014-02-11_reviews/testsim/"
			SSEsim_inputs=SSEsim_inputs;
		
			# 4 areas, 16 states; states 1 and 17 are null geographic range
# 			if (simnum < (numsims/2))
# 				{
# 				rootstate = 18; # Start in area 1, trait state 2
# 				} else {
# 				rootstate = 2;
# 				}
			rootstate = 18; # Start trait state 2 (lower dispersal rate), in area 1 ("K")
			                # (but, have higher rate of exit from this state)

			time_stop=1000;
			taxa_stop=50; # 50
			seed = seedval + simcount;
			printlevel=0;
			testwd=simdir
			SSEsim_results = SSEsim_run(SSEsim_inputs=SSEsim_inputs, rootstate=rootstate, time_stop=time_stop, taxa_stop=taxa_stop, seed=seed, printlevel=printlevel, testwd=testwd)
			names(SSEsim_results)

			#names(SSEsim_results)
			#SSEsim_results$trynum
			#SSEsim_results$simtr
			#dej_params[param_iter,]
			#plot(SSEsim_results$simtr)

			# Process the results
			fossils_older_than=0.001
			printlevel=2
			write_files=TRUE
			SSEsim_results_processed = SSEsim_to_files(SSEsim_results, simdir=simdir, fossils_older_than=fossils_older_than, printlevel=printlevel, write_files=write_files)
			SSEsim_results_processed
			names(SSEsim_results_processed)
			SSEsim_results_processed$tipranges_observed
			SSEsim_results_processed$traits_observed
			names(SSEsim_results_processed)
		
			save(SSEsim_results_processed, file="SSEsim_results_processed.Rdata")


			dim(SSEsim_results_processed$table_of_cladogenetic_events_translated)
			dim(SSEsim_results_processed$table_of_cladogenetic_events_observed)
			dim(SSEsim_results_processed$table_of_range_change_events_translated)


			#####################################################################################################
			# Tabulate the percentage time each part of the tree was in traitstate #1 or traitstate #2
			#####################################################################################################
			# Raw simulation results -- one edge per trait-state
			num_geog_ranges = length(SSEsim_results_processed$SSEsim_results_raw$states_list)
			edge_ranges = SSEsim_results_processed$SSEsim_results_raw$edge_ranges
			edge_length = SSEsim_results_processed$SSEsim_results_raw$edge_length
			
			trait1_TF = edge_ranges <= num_geog_ranges
			trait2_TF = edge_ranges > num_geog_ranges
			time_spent_in_trait1 = sum(edge_length[trait1_TF])
			time_spent_in_trait2 = sum(edge_length[trait2_TF])
			total_time = time_spent_in_trait1 + time_spent_in_trait2
			pct_spent_in_trait1 = round(100*(time_spent_in_trait1 / total_time), digits=1)
			pct_spent_in_trait2 = round(100*(time_spent_in_trait2 / total_time), digits=1)
			time_spent_by_traitstate = c(time_spent_in_trait1, time_spent_in_trait2, total_time, pct_spent_in_trait1, pct_spent_in_trait2)
			time_spent_by_traitstate_names = c("time_spent_in_trait1", "time_spent_in_trait2", "total_time", "pct_spent_in_trait1", "pct_spent_in_trait2")
			names(time_spent_by_traitstate) = time_spent_by_traitstate_names
			time_spent_by_traitstate
			
			#####################################################################################################
			# Tabulate the numbers of anagenetic events, both dispersal and traits, in the FULL SIMULATED tree
			#####################################################################################################
			simtr_ana_event_counts = paste0(SSEsim_results_processed$table_of_range_change_events_translated$event_type, ";", SSEsim_results_processed$table_of_range_change_events_translated$trait_events)
			simtr_ana_event_counts_table = table(simtr_ana_event_counts)
			names(simtr_ana_event_counts_table)
		
			# 3 anagenetic event types in traitstate #1, 3 anagenetic event types in traitstate #2; 2 trait changes
			simtr_ana_event_counts_vector = NULL
		
			simtr_ana_event_counts_vector$sim_t12 = unname(simtr_ana_event_counts_table["trait change (t);1->2"])
			simtr_ana_event_counts_vector$sim_t21 = unname(simtr_ana_event_counts_table["trait change (t);2->1"])
			simtr_ana_event_counts_vector$sim_d_t1 = unname(simtr_ana_event_counts_table["expansion (d);1->1"])
			simtr_ana_event_counts_vector$sim_d_t2 = unname(simtr_ana_event_counts_table["expansion (d);2->2"])
			simtr_ana_event_counts_vector$sim_e_t1 = unname(simtr_ana_event_counts_table["contraction (e);1->1"])
			simtr_ana_event_counts_vector$sim_e_t2 = unname(simtr_ana_event_counts_table["contraction (e);2->2"])
			simtr_ana_event_counts_vector$sim_a_t1 = unname(simtr_ana_event_counts_table["range-switching (a);1->1"])
			simtr_ana_event_counts_vector$sim_a_t2 = unname(simtr_ana_event_counts_table["range-switching (a);2->2"])
			simtr_ana_event_counts_vector = unlist(simtr_ana_event_counts_vector)
		
			#####################################################################################################
			# Tabulate the numbers of anagenetic events, both dispersal and traits, in the OBSERVED tree
			#####################################################################################################
			obstr_ana_event_counts = paste0(SSEsim_results_processed$table_of_range_change_events_observed$event_type, ";", SSEsim_results_processed$table_of_range_change_events_observed$trait_events)
			obstr_ana_event_counts_table = table(obstr_ana_event_counts)
			names(obstr_ana_event_counts_table)
		
			# 3 anagenetic event types in traitstate #1, 3 anagenetic event types in traitstate #2; 2 trait changes
			obstr_ana_event_counts_vector = NULL
		
			obstr_ana_event_counts_vector$obs_t12 = unname(obstr_ana_event_counts_table["trait change (t);1->2"])
			obstr_ana_event_counts_vector$obs_t21 = unname(obstr_ana_event_counts_table["trait change (t);2->1"])
			obstr_ana_event_counts_vector$obs_d_t1 = unname(obstr_ana_event_counts_table["expansion (d);1->1"])
			obstr_ana_event_counts_vector$obs_d_t2 = unname(obstr_ana_event_counts_table["expansion (d);2->2"])
			obstr_ana_event_counts_vector$obs_e_t1 = unname(obstr_ana_event_counts_table["contraction (e);1->1"])
			obstr_ana_event_counts_vector$obs_e_t2 = unname(obstr_ana_event_counts_table["contraction (e);2->2"])
			obstr_ana_event_counts_vector$obs_a_t1 = unname(obstr_ana_event_counts_table["range-switching (a);1->1"])
			obstr_ana_event_counts_vector$obs_a_t2 = unname(obstr_ana_event_counts_table["range-switching (a);2->2"])
			obstr_ana_event_counts_vector = unlist(obstr_ana_event_counts_vector)

		
		
			#####################################################################################################
			# Tabulate the numbers of cladogenetic events, within each traitstate, in the FULL SIMULATED tree
			#####################################################################################################
			simtr_clado_event_counts = paste0(SSEsim_results_processed$table_of_cladogenetic_events_translated$event_type, ";", SSEsim_results_processed$table_of_cladogenetic_events_translated$trait_events)
			simtr_clado_event_counts_table = table(simtr_clado_event_counts)
			names(simtr_clado_event_counts_table)
		
			# 4 clado event types in traitstate #1, 4 clado event types in traitstate #2
			simtr_clado_event_counts_vector = NULL
			# jump dispersal
			simtr_clado_event_counts_vector$sim_j_t1 = unname(simtr_clado_event_counts_table["founder (j);1->1,1"])
			simtr_clado_event_counts_vector$sim_j_t2 = unname(simtr_clado_event_counts_table["founder (j);2->2,2"])
			# narrow sympatry
			simtr_clado_event_counts_vector$sim_y_t1 = unname(simtr_clado_event_counts_table["sympatry (y);1->1,1"])
			simtr_clado_event_counts_vector$sim_y_t2 = unname(simtr_clado_event_counts_table["sympatry (y);2->2,2"])
			# subset sympatry
			simtr_clado_event_counts_vector$sim_s_t1 = unname(simtr_clado_event_counts_table["subset (s);1->1,1"])
			simtr_clado_event_counts_vector$sim_s_t2 = unname(simtr_clado_event_counts_table["subset (s);2->2,2"])
			# vicariance
			simtr_clado_event_counts_vector$sim_v_t1 = unname(simtr_clado_event_counts_table["vicariance (v);1->1,1"])
			simtr_clado_event_counts_vector$sim_v_t2 = unname(simtr_clado_event_counts_table["vicariance (v);2->2,2"])
			simtr_clado_event_counts_vector = unlist(simtr_clado_event_counts_vector)
		
		
			#####################################################################################################
			# Tabulate the numbers of cladogenetic events, within each traitstate, in the OBSERVED tree
			#####################################################################################################
			obstr_clado_event_counts = paste0(SSEsim_results_processed$table_of_cladogenetic_events_observed$event_type, ";", SSEsim_results_processed$table_of_cladogenetic_events_observed$trait_events)
			obstr_clado_event_counts_table = table(obstr_clado_event_counts)
			names(obstr_clado_event_counts_table)

			# 4 clado event types in traitstate #1, 4 clado event types in traitstate #2
			obstr_clado_event_counts_vector = NULL
			# jump dispersal
			obstr_clado_event_counts_vector$obs_j_t1 = unname(obstr_clado_event_counts_table["founder (j);1->1,1"])
			obstr_clado_event_counts_vector$obs_j_t2 = unname(obstr_clado_event_counts_table["founder (j);2->2,2"])
			# narrow sympatry
			obstr_clado_event_counts_vector$obs_y_t1 = unname(obstr_clado_event_counts_table["sympatry (y);1->1,1"])
			obstr_clado_event_counts_vector$obs_y_t2 = unname(obstr_clado_event_counts_table["sympatry (y);2->2,2"])
			# subset sympatry
			obstr_clado_event_counts_vector$obs_s_t1 = unname(obstr_clado_event_counts_table["subset (s);1->1,1"])
			obstr_clado_event_counts_vector$obs_s_t2 = unname(obstr_clado_event_counts_table["subset (s);2->2,2"])
			# vicariance
			obstr_clado_event_counts_vector$obs_v_t1 = unname(obstr_clado_event_counts_table["vicariance (v);1->1,1"])
			obstr_clado_event_counts_vector$obs_v_t2 = unname(obstr_clado_event_counts_table["vicariance (v);2->2,2"])	
			obstr_clado_event_counts_vector = unlist(obstr_clado_event_counts_vector)
	
			# Simulation tree stats
			nsp_sim = length(SSEsim_results_processed$simtr$tip.label)
			nsp_obs = length(SSEsim_results_processed$simtr_observed$tip.label)
			
			simtrtable = prt(SSEsim_results_processed$simtr, printflag=FALSE)
			tree_ht_sim = max(simtrtable$time_bp, na.rm=TRUE)
			tree_len_sim = sum(SSEsim_results_processed$simtr$edge.length)

			obstrtable = prt(SSEsim_results_processed$simtr_observed, printflag=FALSE)
			tree_ht_obs = max(obstrtable$time_bp, na.rm=TRUE)
			tree_len_obs = sum(SSEsim_results_processed$simtr_observed$edge.length)
			
			treestats = c(nsp_sim, tree_ht_sim, tree_len_sim, nsp_obs, tree_ht_obs, tree_len_obs)
			treestats_names = c("nsp_sim", "tree_ht_sim", "tree_len_sim", "nsp_obs", "tree_ht_obs", "tree_len_obs")
			names(treestats) = treestats_names
			
			sim_summary = c(treestats, time_spent_by_traitstate, simtr_ana_event_counts_vector, obstr_ana_event_counts_vector, simtr_clado_event_counts_vector, obstr_clado_event_counts_vector)

			# Set NAs to 0
			sim_summary[is.na(sim_summary)] = 0
			
			return(sim_summary)
			} # END sim_function
		
		# Try the above-specified simulation with a try() function
		sim_summary = try(expr=sim_function(dej_params, param_iter, number_of_trait_states, trait_fn, numsims,simcount))
		
		
		
		# Error catch
		if (class(sim_summary) == "try-error")
			{
			time_spent_by_traitstate_names = c("time_spent_in_trait1", "time_spent_in_trait2", "total_time", "pct_spent_in_trait1", "pct_spent_in_trait2")

			sim_summary_names = c("nsp_sim", "tree_ht_sim", "tree_len_sim", "nsp_obs", "tree_ht_obs", "tree_len_obs", time_spent_by_traitstate_names, "t12", "t21", "d_t1", "d_t2", "e_t1", "e_t2", "a_t1", "a_t2", "t12", "t21", "d_t1", "d_t2", "e_t1", "e_t2", "a_t1", "a_t2", "j_t1", "j_t2", "y_t1", "y_t2", "s_t1", "s_t2", "v_t1", "v_t2", "j_t1", "j_t2", "y_t1", "y_t2", "s_t1", "s_t2", "v_t1", "v_t2")
			sim_summary = rep(NA, times=length(sim_summary_names))
			names(sim_summary) = sim_summary_names
			}
		
		
		# The input simulation parameters
		seed = seedval + simcount
		preface = c(param_iter, simnum, simcount, seed, dej_params[param_iter,])
		names(preface) = c("param_iter", "simnum", "simcount", "seed", names(dej_params))
		
		# Plot the results
		# plot(SSEsim_results$simtr, label.offset=0.2)
		# axisPhylo()
		# title("SSE sim, with species labels")
		# nodelabels()
		# tiplabels()


		tmp_results = c(preface, sim_summary)
		tmp_results_mat = as.matrix(tmp_results, nrow=1)
		tmp_results_mat = as.data.frame(tmp_results_mat, stringsAsFactors=FALSE)
		
		outfn = slashslash(paste0(master_wd, "/", "sim_summary.txt"))
		if ( (param_iter == 1) && (simnum == 1) )
			{
			write.table(x=tmp_results, file=outfn, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
			} else {
			write.table(x=tmp_results, file=outfn, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
			}
		moref(outfn)
		

		# Print the results
		txt = paste("\nSimulation finished. Printing sim_summary.", "\n", sep="")
		cat(txt)
		print(sim_summary)

		} # end 100 simulations
	} # end dej_params


# stop("Manual ending after simulations...")


 for (param_iter in 1:1)
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

		simdir = slashslash(paste(master_wd, "/sims/", "ps", ps_count_txt, "_sim", simcount_txt, sep=""))
		simdir


		txt = paste("\nINFERENCE on simulation #", simcount, " in ", simdir, "/", "\n", sep="")
		cat(txt)

		setwd(simdir)


		#######################################################
		# Inference
		#######################################################
		max_range_size = 4


		#######################################################
		# Traits-only model -- 1 rate
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = 1
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "tipranges_observed_1area.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		# Set up DEC model, but set all rates to 0 (data are 1 invariant area)
		# (nothing to do; defaults)

		# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
		BioGeoBEARS_run_object

		# This contains the model object
		BioGeoBEARS_run_object$BioGeoBEARS_model_object

		# This table contains the parameters of the model 
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0.0

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.0

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.0

		# Set up BAYAREALIKE model
		# No subset sympatry
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

		# No vicariance
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

		# No jump dispersal/founder-event speciation
		# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
		# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

		# Adjust linkage between parameters
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

		# Only sympatric/range-copying (y) events allowed, and with 
		# exact copying (both descendants always the same size as the ancestor)
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		geog_values = getranges_from_LagrangePHYLIP("tipranges_observed_1area.data")

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


		# Look at the params table
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "t12"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 1

		# No multipliers on geog (set m1 and m2 to 1)
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "desc"] = "trait-based dispersal rate multipliers m1"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "desc"] = "trait-based dispersal rate multipliers m2"

		# Run this to check inputs. Read the error messages if you get them!
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		# For a slow analysis, run once, then set runslow=FALSE to just 
		# load the saved result.
		runslow = TRUE
		resfn = "sim_traitsOnly_1rate_v1.Rdata"
		if (runslow)
				{
				res = bears_optim_run(BioGeoBEARS_run_object)
				res    

				save(res, file=resfn)
				resTrait_1rate = res
				} else {
				# Loads to "res"
				load(resfn)
				resTrait_1rate = res
				} # END if (runslow)



		#######################################################
		# Traits-only model -- 2 rates
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = 1
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "tipranges_observed_1area.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		# Set up DEC model, but set all rates to 0 (data are 1 invariant area)
		# (nothing to do; defaults)

		# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
		BioGeoBEARS_run_object

		# This contains the model object
		BioGeoBEARS_run_object$BioGeoBEARS_model_object

		# This table contains the parameters of the model 
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0.0

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.0

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.0

		# Set up BAYAREALIKE model
		# No subset sympatry
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

		# No vicariance
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

		# No jump dispersal/founder-event speciation
		# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
		# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

		# Adjust linkage between parameters
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

		# Only sympatric/range-copying (y) events allowed, and with 
		# exact copying (both descendants always the same size as the ancestor)
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


		# Look at the params table
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 1

		# No multipliers on geog (set m1 and m2 to 1)
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1.0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "desc"] = "trait-based dispersal rate multipliers m1"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "desc"] = "trait-based dispersal rate multipliers m2"

		# Run this to check inputs. Read the error messages if you get them!
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		# For a slow analysis, run once, then set runslow=FALSE to just 
		# load the saved result.
		runslow = TRUE
		resfn = "sim_traitsOnly_2rates_v1.Rdata"
		if (runslow)
				{
				res = bears_optim_run(BioGeoBEARS_run_object)
				res    

				save(res, file=resfn)
				resTrait_2rates = res
				} else {
				# Loads to "res"
				load(resfn)
				resTrait_2rates = res
				} # END if (runslow)



		#######################################################
		# Run DEC (on geography only)
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		#tr = read.tree(BioGeoBEARS_run_object$trfn)

		
		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
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
		# Run DEC+J (on geography only)
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		#tr = read.tree(BioGeoBEARS_run_object$trfn)
		
		dstart = resDEC$outputs@params_table["d","est"]
		estart = resDEC$outputs@params_table["e","est"]
		jstart = 0.0001
		
		# Add j as a free parameter
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart


		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)
		
		print("Printing warnings: 'warnings()':")
		print(warnings())
		
		runslow = TRUE
		resfn = "DECj_inf.Rdata"
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






		#######################################################
		# Run DEC + t12 + t21 + m2, starting from DEC-geog and 2rates
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


		# Look at the params table
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
		
		
		# Starting values from ML results of simpler run
		t12_start = resTrait_2rates$outputs@params_table["t12","est"]
		t21_start = resTrait_2rates$outputs@params_table["t21","est"]
		m2_start = 1
		dstart = resDEC$outputs@params_table["d","est"]
		estart = max(c(resDEC$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
		jstart = 0.0001

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart


		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0.1


		# Set 0/1 multipliers on dispersal rate
		# For flightlessness (m2), max multiplier is 1, and
		# fix to a small value, or estimate
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10




		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		runslow = TRUE
		resfn = "DEC+t12+t21+m2_inf.Rdata"
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)
			resDEC_t12_t21_m2 = res
			} else {
			# Loads to "res"
			load(resfn)
			resDEC_t12_t21_m2 = res
			}




		#######################################################
		# Run DECj + t12 + t21 + m2, starting from DECj-geog and 2rates
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



		# Starting values from ML results of simpler run
		t12_start = resTrait_2rates$outputs@params_table["t12","est"]
		t21_start = resTrait_2rates$outputs@params_table["t21","est"]
		m2_start = 1
		dstart = resDECj$outputs@params_table["d","est"]
		estart = max(c(resDECj$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
		jstart = resDECj$outputs@params_table["j","est"]


		# Set up DEC+J model
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

		# Add j as a free parameter
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

		# Crash fix
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13

		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0.1


		# Set 0/1 multipliers on dispersal rate
		# For flightlessness (m2), max multiplier is 1, and
		# fix to a small value, or estimate
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10




		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		resfn = "DECJ+t12+t21+m2_inf.Rdata"
		runslow = TRUE
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)

			resDECj_t12_t21_m2 = res
			} else {
			# Loads to "res"
			load(resfn)
			resDECj_t12_t21_m2 = res
			}





		#######################################################
		# Run DECj + t12 + t21 + m2, starting from DEC + t12 + t21 + m2
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



		# Starting values from ML results of simpler run
		t12_start = resDEC_t12_t21_m2$outputs@params_table["t12","est"]
		t21_start = resDEC_t12_t21_m2$outputs@params_table["t21","est"]
		m2_start = resDEC_t12_t21_m2$outputs@params_table["m2","est"]
		dstart = resDEC_t12_t21_m2$outputs@params_table["d","est"]
		estart = max(c(resDEC_t12_t21_m2$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
		jstart = 0.0001


		# Set up DEC+J model
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

		# Add j as a free parameter
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

		# Crash fix
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13

		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0.1


		# Set 0/1 multipliers on dispersal rate
		# For flightlessness (m2), max multiplier is 1, and
		# fix to a small value, or estimate
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10



		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		resfn = "DECJ+t12+t21+m2_rep2_inf.Rdata"
		runslow = TRUE
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)

			resDECj_t12_t21_m2_rep2 = res
			} else {
			# Loads to "res"
			load(resfn)
			resDECj_t12_t21_m2_rep2 = res
			}







		#######################################################
		# Run DEC - no M!(independent check)
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


		# Look at the params table
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0.007
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0.1


		# Set 0/1 multipliers on dispersal rate
		# For flightlessness (m2), max multiplier is 1, and
		# fix to a small value, or estimate
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 10

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10




		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		runslow = TRUE
		resfn = "DEC+t12+t21+noM2_inf.Rdata"
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)
			resDEC_noM = res
			} else {
			# Loads to "res"
			load(resfn)
			resDEC_noM = res
			}



		#######################################################
		# Run DECj - no M! (independent check)
		#######################################################
		BioGeoBEARS_run_object = define_BioGeoBEARS_run()
		BioGeoBEARS_run_object$print_optim = TRUE
		BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
		BioGeoBEARS_run_object$max_range_size = max_range_size
		BioGeoBEARS_run_object$num_cores_to_use=1
		BioGeoBEARS_run_object$cluster_already_open = cluster_already_open
		BioGeoBEARS_run_object$use_optimx="GenSA"
		BioGeoBEARS_run_object$speedup=TRUE
		BioGeoBEARS_run_object$geogfn = "geog.data"
		BioGeoBEARS_run_object$trfn = "tree.newick"
		BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
		BioGeoBEARS_run_object$return_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
		BioGeoBEARS_run_object$calc_ancprobs = TRUE
		BioGeoBEARS_run_object$on_NaN_error = -1000000

		tr = read.tree(BioGeoBEARS_run_object$trfn)
		#plot(tr); axisPhylo()

		trait_fn = "traits.data"
		trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
		trait_values

		# Add the traits data and model
		BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)

		# Set up DEC+J model

		# Starting values from ML results of simpler run
		t12_start = resDEC_noM$outputs@params_table["t12","est"]
		t21_start = resDEC_noM$outputs@params_table["t21","est"]
		m2_start = resDEC_noM$outputs@params_table["m2","est"]
		dstart = resDEC_noM$outputs@params_table["d","est"]
		estart = resDEC_noM$outputs@params_table["e","est"]
		jstart = 0.0001


		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart


		# Add j as a free parameter
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart



		#######################################################
		# Manual modifications of trait-based model
		#######################################################
		# Edit t12 and t21 rates

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "free"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = t21_start
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0.00001
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0.1


		# Set 0/1 multipliers on dispersal rate
		# For flightlessness (m2), max multiplier is 1, and
		# fix to a small value, or estimate
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 10

		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
		BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 10


		BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
		check_BioGeoBEARS_run(BioGeoBEARS_run_object)

		resfn = "DECJ+t12+t21+noM2_inf.Rdata"
		runslow = TRUE
		if (runslow)
			{
			res = bears_optim_run(BioGeoBEARS_run_object)
			res    

			save(res, file=resfn)

			resDECj_noM = res
			} else {
			# Loads to "res"
			load(resfn)
			resDECj_noM = res
			}
		
		
		param_names = c("lnL", "d", "e", "j", "t12", "t21", "m1", "m2")

		Trait_1rate_results = c(
		resTrait_1rate$total_loglikelihood,
		resTrait_1rate$output@params_table["d", "est"], 
		resTrait_1rate$output@params_table["e", "est"], 
		resTrait_1rate$output@params_table["j", "est"], 
		resTrait_1rate$output@params_table["t12", "est"], 
		resTrait_1rate$output@params_table["t21", "est"], 
		resTrait_1rate$output@params_table["m1", "est"], 
		resTrait_1rate$output@params_table["m2", "est"]
		)
		names(Trait_1rate_results) = paste("Trait_1rate_", param_names, sep="")

		Trait_2rates_results = c(
		resTrait_2rates$total_loglikelihood,
		resTrait_2rates$output@params_table["d", "est"], 
		resTrait_2rates$output@params_table["e", "est"], 
		resTrait_2rates$output@params_table["j", "est"], 
		resTrait_2rates$output@params_table["t12", "est"], 
		resTrait_2rates$output@params_table["t21", "est"], 
		resTrait_2rates$output@params_table["m1", "est"], 
		resTrait_2rates$output@params_table["m2", "est"]
		)
		names(Trait_2rates_results) = paste("Trait_2rates_", param_names, sep="")

		DEC_results = c(
		resDEC$total_loglikelihood,
		resDEC$output@params_table["d", "est"], 
		resDEC$output@params_table["e", "est"], 
		resDEC$output@params_table["j", "est"], 
		resDEC$output@params_table["t12", "est"], 
		resDEC$output@params_table["t21", "est"], 
		resDEC$output@params_table["m1", "est"], 
		resDEC$output@params_table["m2", "est"]
		)
		names(DEC_results) = paste("DEC_", param_names, sep="")


		DECj_results = c(
		resDECj$total_loglikelihood,
		resDECj$output@params_table["d", "est"], 
		resDECj$output@params_table["e", "est"], 
		resDECj$output@params_table["j", "est"], 
		resDECj$output@params_table["t12", "est"], 
		resDECj$output@params_table["t21", "est"], 
		resDECj$output@params_table["m1", "est"], 
		resDECj$output@params_table["m2", "est"]
		)
		names(DECj_results) = paste("DECj_", param_names, sep="")

		DEC_t12_t21_m2_results = c(
		resDEC_t12_t21_m2$total_loglikelihood,
		resDEC_t12_t21_m2$output@params_table["d", "est"], 
		resDEC_t12_t21_m2$output@params_table["e", "est"], 
		resDEC_t12_t21_m2$output@params_table["j", "est"], 
		resDEC_t12_t21_m2$output@params_table["t12", "est"], 
		resDEC_t12_t21_m2$output@params_table["t21", "est"], 
		resDEC_t12_t21_m2$output@params_table["m1", "est"], 
		resDEC_t12_t21_m2$output@params_table["m2", "est"]
		)
		names(DEC_t12_t21_m2_results) = paste("DEC_t12_t21_m2_", param_names, sep="")

		DECj_t12_t21_m2_results = c(
		resDECj_t12_t21_m2$total_loglikelihood,
		resDECj_t12_t21_m2$output@params_table["d", "est"], 
		resDECj_t12_t21_m2$output@params_table["e", "est"], 
		resDECj_t12_t21_m2$output@params_table["j", "est"], 
		resDECj_t12_t21_m2$output@params_table["t12", "est"], 
		resDECj_t12_t21_m2$output@params_table["t21", "est"], 
		resDECj_t12_t21_m2$output@params_table["m1", "est"], 
		resDECj_t12_t21_m2$output@params_table["m2", "est"]
		)
		names(DECj_t12_t21_m2_results) = paste("DECj_t12_t21_m2_", param_names, sep="")

		DECj_t12_t21_m2_rep2_results = c(
		resDECj_t12_t21_m2_rep2$total_loglikelihood,
		resDECj_t12_t21_m2_rep2$output@params_table["d", "est"], 
		resDECj_t12_t21_m2_rep2$output@params_table["e", "est"], 
		resDECj_t12_t21_m2_rep2$output@params_table["j", "est"], 
		resDECj_t12_t21_m2_rep2$output@params_table["t12", "est"], 
		resDECj_t12_t21_m2_rep2$output@params_table["t21", "est"], 
		resDECj_t12_t21_m2_rep2$output@params_table["m1", "est"], 
		resDECj_t12_t21_m2_rep2$output@params_table["m2", "est"]
		)
		names(DECj_t12_t21_m2_rep2_results) = paste("DECj_t12_t21_m2_rep2_", param_names, sep="")

		
		DEC_noM_results = c(
		resDEC_noM$total_loglikelihood,
		resDEC_noM$output@params_table["d", "est"], 
		resDEC_noM$output@params_table["e", "est"], 
		resDEC_noM$output@params_table["j", "est"], 
		resDEC_noM$output@params_table["t12", "est"], 
		resDEC_noM$output@params_table["t21", "est"], 
		resDEC_noM$output@params_table["m1", "est"], 
		resDEC_noM$output@params_table["m2", "est"]
		)
		names(DEC_noM_results) = paste("DEC_noM_", param_names, sep="")

		DECj_noM_results = c(
		resDECj_noM$total_loglikelihood,
		resDECj_noM$output@params_table["d", "est"], 
		resDECj_noM$output@params_table["e", "est"], 
		resDECj_noM$output@params_table["j", "est"], 
		resDECj_noM$output@params_table["t12", "est"], 
		resDECj_noM$output@params_table["t21", "est"], 
		resDECj_noM$output@params_table["m1", "est"], 
		resDECj_noM$output@params_table["m2", "est"]
		)
		names(DECj_noM_results) = paste("DECj_noM_", param_names, sep="")
	
		tmp_results = c(dej_params[param_iter,], Trait_1rate_results, Trait_2rates_results, DEC_results, DECj_results, DEC_t12_t21_m2_results, DECj_t12_t21_m2_results, DECj_t12_t21_m2_rep2_results, DEC_noM_results, DECj_noM_results)
		tmp_results_mat = as.matrix(tmp_results, nrow=1)
		tmp_results_mat = as.data.frame(tmp_results_mat, stringsAsFactors=FALSE)
		
		outfn = slashslash(paste0(master_wd, "/", "params_sim_inferred.txt"))
		if ( (param_iter == 1) && (simnum == 1) )
			{
			write.table(x=tmp_results, file=outfn, append=FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
			} else {
			write.table(x=tmp_results, file=outfn, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
			}
		moref(outfn)
		

		# Print the results
		txt = paste("\nSimulation finished. DEC LnL=", resDEC$total_loglikelihood, ", DEC+J LnL=", resDECj$total_loglikelihood, "\n", sep="")
		cat(txt)
		cat("Model params:\n")
		print(dej_params[param_iter,])



		} # end 100 simulations
	} # end dej_params


if (.Platform$GUI != "AQUA" && num_cores_to_use > 1)
	{
	cat("\n\nStopping cluster with ", num_cores_to_use, " cores.\n\n", sep="")
	stopCluster(cluster_already_open)
	}

stop("End of inference runs")


