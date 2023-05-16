#!/usr/bin/env bash

# ngsF-HMM.sh - a wrapper to help ensure convergence of ngsF-HMM replicates.
# Execute with whatever arguments you would give to ngsF-HMM and it will run
# replicates until they converge or a maximum number is reached

# Depends on: ngsF-HMM (case sensitive) in $PATH

# Author: Zachary J. Nolen
# Last updated: 2023-05-16

# This version is modified to be compatible with a snakemake pipeline, see 
# previous commits for a more general version.

#####################
### Housekeeping ####
#####################

(# Define function to make printing to stderr easy
echoerr() { echo "$@" 1>&2; }

# Make it so stuff fails more often
set -eo pipefail

################################
## Import Snakemake Arguments ##
################################

geno="${snakemake_input[beagle]}"
pos="${snakemake_output[pos]}"
zcat $geno | awk '{print $1}' | sed 's/\(.*\)_/\1\t/' | tail -n +2 > $pos
nind="${snakemake_params[nind]}"
nsites=$(cat $pos | wc -l)
outpre="${snakemake_params[out]}"
threads="${snakemake[threads]}"
TMPDIR="${snakemake_resources[tmpdir]}"

passedargs="--geno $geno --pos $pos --n_ind $nind --n_sites $nsites --lkl --n_threads $threads"

# log what happened
echoerr "The following options will be passed to ngsF-HMM, please make sure"
echoerr "they are correct:"
echoerr "$passedargs --out $outpre"
echoerr

#####################
# Setup directories #
#####################

outfile=$(basename $outpre)
outdir=$(dirname $outpre)

# temp directory
temp=$TMPDIR/ngsF-HMM_opt_$outfile
mkdir -p $temp

# out directory
mkdir -p $outdir

# Make sure output and temp directories aren't the same. There is a very low 
# chance they would be, but why risk it...
if [ $temp == $outdir ]; then
	echoerr "ERROR: Temporary directory is the same as output directory."
	echoerr 'Please change variable $TMPDIR or -o/-outfiles path so that'
	echoerr "they do not match."
	exit -1
fi

# log what happened
echoerr "Temporary files will be written to:"
echoerr $temp
echoerr
echoerr "Highest likelihood replicate will be written with this prefix:"
echoerr $outpre
echoerr

#####################
##### Set params ####
#####################

# set maximum number of replicates to perform before stopping if not 
# converged
if [ -z $reps ] || [ "$reps" = "None" ]; then
	reps=100
else
	echoerr 'WARNING: $reps has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

# set maximum range between top reps in likelihood units to qualify reps
# as reaching convergance
if [ -z $thresh ] || [ "$thresh" = "None" ]; then
	thresh=2
else
	echoerr 'WARNING: $thresh has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

# set the minimum number of top replicates that must be within the value 
# $thresh to qualify as converged
if [ -z $conv ] || [ "$conv" = "None" ]; then
	conv=3
else
	echoerr 'WARNING: $conv has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

# log what got set
echoerr "Beginning ngsF-HMM replicates, will perform replicate runs until"
echoerr "the top $conv runs are within $thresh log-likelihood units of each"
echoerr "other or until $reps replicates have finished. If you'd like to"
echoerr 'change these convergence criteria, edit the $reps, $thresh, or'
echoerr '$conv variables in the wrapper or export them before starting.'
echoerr

# Maybe one day I will add resuming a terminated run, but I don't have a good
# way to save the best likelihood from the previous run without it being a 
# false positive of convergence. Maybe temp names?

# if [ $resume == "1" ]; then
# 	echoerr "Picking up from previous interrupted run at the last completed"
# 	echoerr "replicate. WARNING: I won't check to make sure your inputs or"
# 	echoerr "options are the same since the last run, so please be sure of"
# 	echoerr "this or start over by removing the -resume option."
# 	echoerr
# fi

#####################
## Run the wrapper ##
#####################

# set log file where replicate likelihoods will be stored and assessed from
templog=$temp/${outfile}_optimization_wrapper.log
permlog=$outdir/${outfile}_optimization_wrapper.log

# initialize the log and make the temporary directory for kept results
> $permlog
> $templog
mkdir -p $temp/bestrep

# run replicates of ngsF-HMM until we hit either convergence or $reps
for i in $(seq 1 $reps); do
	# check if loop should break when convergence is reached
	if (( "$i" > "20" )); then
		# sort likelihoods, see if top three are within $thresh of each other
		# and break if so
		diff=$(sort -k2gr $templog | head -n 3 | awk '{print $2}' | \
				tr '\n' ' ' | awk '{print $1-$3}')
		if (( $(echo $diff $thresh | awk '{print $1 <= $2}') )); then
			break
		fi
	fi

	# make a folder for the replicate results
	reptemp=$temp/iter$i
	mkdir -p $reptemp

	# run replicate
	echoerr "Beginning replicate $i..."
	ngsF-HMM $passedargs --out $reptemp/$outfile --seed $RANDOM

	# add likelihood to log
	like=$(head -n 1 $reptemp/$outfile.indF)
	echo "$i	$like" >> $templog
	cp $templog $permlog

	# if better than previous reps, store results
	bestrep=$(sort -k2gr $templog | head -n 1 | awk '{print $1}')
	if [ "$bestrep" == "$i" ]; then
		mv $reptemp/* $temp/bestrep/
	fi

	# delete replicate folder
	rm -r $reptemp

done

if (( $(echo $diff $thresh | awk '{print $1 <= $2}') )); then
	echoerr
	echoerr "Optimization complete! Top $conv replicates have converged with"
	echoerr "$diff log-likelihood unit difference."
	echoerr
else
	echoerr
	echoerr "WARNING: Optimization reached max reps ($reps) without converging"
	echoerr "with a difference of $diff log-likelihood units between the top"
	echoerr "$conv results. Consider increasing the number of replicates or"
	echoerr "performing the analysis only for lower values of K."
	echoerr
fi

echoerr "Moving best likelihood results to output folder..."
mv $temp/bestrep/* $outdir/

echoerr "Done!"
echoerr

echoerr "Final checks to make sure you have the right files..."
bestlike=$(sort -k2gr $templog | head -n 1 | awk '{print $2}')
yourlike=$(head -n 1 $outpre.indF)

if (( $(awk 'BEGIN{ print "'$bestlike'"=="'$yourlike'" }') )); then
	echoerr "Best likelihood results successfully written to output directory!"
else
	echoerr "ERROR: Something went wrong, try comparing the values of the"
	echoerr "optimization_wrapper.log file and the .log file of the ngsF-HMM"
	echoerr "run."
	echoerr "best like = $bestlike; output like = $yourlike"
	exit -1
fi) 2> "${snakemake_log[0]}"