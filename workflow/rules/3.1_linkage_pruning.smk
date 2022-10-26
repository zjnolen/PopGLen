localrules:

rule ngsLD_estLD:
	input:
		beagle=rules.angsd_doGlf2.output.beagle,
		bamlist=rules.angsd_makeBamlist.output
	output:
		ld=results+"/genotyping/pruned_beagle/ngsLD/"+dataset+
			"_{population}{dp}_chunk{chunk}.ld.gz",
		pos=results+"/genotyping/pruned_beagle/ngsLD/"+dataset+
			"_{population}{dp}_chunk{chunk}.pos"
	log:
		logs + "/ngsLD/estLD/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		"library://james-s-santangelo/ngsld/ngsld:1.1.1"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*720
	shell:
		r"""
		zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos} 2> {log}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		if [[ $nsites == 0 ]]; then
			touch {output.ld}
		else
			ngsLD --geno {input.beagle} --n_ind $nind --n_sites $nsites \
				--pos {output.pos} --probs --n_threads {threads} \
				| gzip > {output.ld} 2>> {log}
		fi
		"""


rule ngsLD_prune_sites:
	input:
		ld=rules.ngsLD_estLD.output.ld,
		pos=rules.ngsLD_estLD.output.pos
	output:
		sites=results+"/genotyping/pruned_beagle/ngsLD/"+dataset+
			"{population}{dp}_chunk{chunk}_pruned.sites"
	log:
		logs + "/ngsLD/prune_sites/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	conda:
		"../envs/pruning.yaml"
	threads: lambda wildcards, attempt: attempt*5
	resources:
		time=lambda wildcards, attempt: attempt*1440
	shell:
		"""
		nsites=$(cat {input.pos} | wc -l)

		if [[ $nsites == 0 ]]; then
			touch {output.sites}
		else
			workflow/scripts/prune_ngsLD.py --input {input.ld} --max_dist 50000 \
				--min_weight 0.1 --output {output.sites} 2> {log}
		fi
		"""

rule prune_chunk_beagle:
	input:
		beagle=rules.angsd_doGlf2.output.beagle,
		sites=expand(results+"/genotyping/pruned_beagle/ngsLD/"+dataset+
			"{population}{{dp}}_chunk{{chunk}}_pruned.sites",
			population="all")
	output:
		pruned=temp(results+"/genotyping/pruned_beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}_pruned.beagle"),
		prunedgz=results+"/genotyping/pruned_beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}_pruned.beagle.gz"
	log:
		logs+"/ngsLD/prune_beagle/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	threads: lambda wildcards, attempt: attempt*10
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		r"""
		set +o pipefail;
		zcat {input.beagle} | head -n 1 > {output.pruned} 2> {log}
		
		join -t $'\t' <(sort -k1,1 {input.sites}) <(zcat {input.beagle} | sort -k1,1) | \
			sed 's/_/\t/' | sort -k1,1 -k2,2n | sed 's/\t/_/' \
			>> {output.pruned} 2>> {log}

		gzip -c {output.pruned} > {output.prunedgz} 2>> {log}

		Nsites=$(cat {input.sites} | wc -l | awk '{{print $1+1}}') &>> {log}
		NsitesB=$(zcat {output.prunedgz} | wc -l) &>> {log}

		echo "Sites searched for: $Nsites" &>> {log}
		echo "Sites in pruned beagle: $NsitesB" &>> {log}

		# No longer works since pruned site list isn't custom to population.
		# Would fail for any population that happens to be missing a site in
		# the pruned dataset in all individuals.
		# if [ $Nsites = $NsitesB ]; then
		# 	echo "Pruning successful!" &>> {log}
		# else
		# 	echo "Number of sites in pruned beagle differ from sites file." \
		# 		&>> {log}
		# 	exit 1
		# fi
		"""

rule merge_pruned_beagles:
	input:
		pruned=lambda w: expand(results+"/genotyping/pruned_beagle/chunk/"+
			dataset+"_{{population}}{{dp}}_chunk{chunk}_pruned.beagle.gz", 
			chunk=chunklist),
		full=results+"/genotyping/beagle/"+dataset+
			"_{population}{dp}.beagle.gz"
	output:
		beagle=results+"/genotyping/pruned_beagle/"+dataset+
			"_{population}{dp}_pruned.beagle.gz"
	log:
		logs + "/ngsLD/{population}{dp}_merge_pruned.log"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		echo "cat file order:"
		echo {input.pruned} | tr ' ' '\n'

		echo "Printing header to beagle..."
		set +o pipefail;
		zcat {input.full} | head -n 1 | gzip > {output.beagle}

		echo "Adding requested chunks to beagle..."
		for f in {input.pruned}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.beagle} 2>> {log}
		done
		"""