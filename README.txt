# This document describes how to run the code which produced the simulation results in
#  the 2020 Survey Methodology paper Clark (2020): Model-assisted sample design is minimax for model-based prediction
#
# The code was run on a unix high performance system with a large number of cores.
# It could be adapted to run on a PC but would take 1-2 weeks or longer to run.
#
# The files AV_bound_allk.sh, Rjob_AVbound_k.txt, AV_single_NCI230119.R, AV_preamble_220119.R and AV_analyse210119.R should all be in the current directory.
# You should change the directory (paperdir) that results are saved to in AV_analyse210119.R (they currently save to my google drive)
# From the directory containing the above scripts, run
     AV_bound_allk.sh
# This will (eventually) produce the files simres1.rdata, …, simres184.rdata
# These can then be analysed by opening R and running the R code in AV_analyse210119.R. 
# Text files containing the bodies of LaTex tables will be produced, also some plot files.



