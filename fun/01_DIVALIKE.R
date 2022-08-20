####### Writing began in August 2022 by Thales de Lima
####### Written during the course Introduction to Scientific Computation
####### Written as part of the final evaluation of this same course
####### Written based on the template script of the package BioGeoBEARS (Matzke, Nicholas J. (2018))
####### Serrapilheira ICTP/SAIRF QBIO program

####### This script set the functions for the ancestral range estimation
####### The model used is an adapted likelyhood version of
####### the DIVA (dispersal-vicariance) model (Ronquist (1997).

####### I have tried to modify it as less as possible from the original template
####### Therefore, only the main inputs are called as arguments of the function,
####### while all other parameters are left as fixed, in the default option
####### (i.e., as in template script)

####### The template script and many useful explanations and details about the package
######## and its models can be found here:
####### http://phylo.wikidot.com/biogeobears#script

####### The package GitHub page with some informations can be found here:
####### https://github.com/nmatzke/BioGeoBEARS



### Library -------------------------------------------------------------
library(optimx)
library(GenSA)
library(FD)
library(snow)
library(parallel)
library(ape)

library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)



### Function for DIVA model without founder-event speciation -------------------------

DIVALIKE <- function(trfn, geogfn, max_range_size) {
  #######################################################
  # Run DIVALIKE
  #######################################################
  BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu,
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change

  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
  # if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
  # problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table

  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

  # Set up DIVALIKE model
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  runslow = TRUE
  if (runslow) {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    resDIVALIKE = res
  } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKE = res
  }
  return(resDIVALIKE)
}




### Function for DIVA model with founder-event speciation ----------------------------

DIVALIKE_J <- function(trfn, geogfn, max_range_size,
                       resDIVALIKE) {

  #######################################################
  # Run DIVALIKE+J
  #######################################################
  resDIVALIKE <- resDIVALIKE
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu,
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change

  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
  # if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
  # problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table

  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

  # Set up DIVALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = 0.0001

  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

  check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  runslow = TRUE
  if (runslow)
  {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res <-bears_optim_run(BioGeoBEARS_run_object)
    res


    resDIVALIKEj <- res
  } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKEj = res
  }
  return(resDIVALIKEj)
}
