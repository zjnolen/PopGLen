ruleorder: merge_pruned_beagles > prune_chrom_beagle
localrules: merge_pruned_beagles

rule ngsLD_estLD:
	input:
		beagle=rules.angsd_chrom_beagle.output.beagle,
		bamlist=rules.angsd_makeBamlist.output
	output:
		ld=results+"/ngsLD/{population}_chr{chrom}.ld.gz",
		pos=results+"/ngsLD/{population}_chr{chrom}.pos"
	log:
		logs + "/ngsLD/estLD/{population}_chr{chrom}.log"
	container:
		"library://james-s-santangelo/ngsld/ngsld:1.1.1"
	threads: lambda wildcards, attempt: attempt
	resources:
		time="12:00:00"
	shell:
		r"""
		zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos} 2> {log}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1 + 1}}')

		ngsLD --geno {input.beagle} --n_ind $nind --n_sites $nsites \
			--pos {output.pos} --probs --n_threads {threads} \
			| gzip > {output.ld} 2>> {log}
		"""

rule ngsLD_prune_sites:
	input:
		ld=rules.ngsLD_estLD.output.ld
	output:
		sites=results+"/ngsLD/{population}_chr{chrom}_pruned.sites"
	log:
		logs + "/ngsLD/prune_sites/{population}_chr{chrom}.log"
	conda:
		"../envs/pruning.yaml"
	threads: lambda wildcards, attempt: attempt*4
	resources:
		time="24:00:00"
	shell:
		"""
		workflow/scripts/prune_ngsLD.py --input {input.ld} \
			--output {output.sites} 2> {log}
		"""

rule prune_chrom_beagle:
	input:
		beagle=rules.angsd_chrom_beagle.output.beagle,
		sites=rules.ngsLD_prune_sites.output.sites
	output:
		pruned=temp(results+"/angsd/beagle/pruned/chrom/{population}_chr{chrom}_pruned.beagle"),
		prunedgz=results+"/angsd/beagle/pruned/chrom/{population}_chr{chrom}_pruned.beagle.gz"
	log:
		logs + "/ngsLD/prune_beagle/{population}_chr{chrom}.log"
	threads: lambda wildcards, attempt: attempt
	resources:
		time="04:00:00"
	shell:
		r"""
		set +o pipefail;
		zcat {input.beagle} | head -n 1 > {output.pruned} 2> {log}

		while read pos; do
			zcat {input.beagle} | grep "$pos	" >> {output.pruned} 2>> {log}
		done < {input.sites}

		gzip -c {output.pruned} > {output.prunedgz} 2> {log}

		Nsites=$(cat {input.sites} | wc -l | awk '{{print $1+1}}') &>> {log}
		NsitesB=$(zcat {output.prunedgz} | wc -l) &>> {log}

		if [ $Nsites = $NsitesB ]; then
			echo "Pruning successful!" &>> {log}
		else
			echo "Number of sites in pruned beagle differ from sites file." \
				&>> {log}
			exit 1
		fi
		"""

rule merge_pruned_beagles:
	input:
		pruned=lambda w: expand(results+"/angsd/beagle/pruned/chrom/" \
			"{{population}}_chr{chrom}_pruned.beagle.gz", chrom=get_contigs()),
		full=results+"/angsd/beagle/{population}_genome.beagle.gz"
	output:
		beagle=results+"/angsd/beagle/pruned/{population}_genome_pruned.beagle.gz"
	log:
		logs + "/ngsLD/{population}_genome_merge_pruned_beagle.log"
	shell:
		r"""
		echo "cat file order:"
		echo {input.pruned} | tr ' ' '\n'

		echo "Printing header to beagle..."
		set +o pipefail;
		zcat {input.full} | head -n 1 | gzip > {output.beagle}

		echo "Adding requested chromosomes to beagle..."
		for f in {input.pruned}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.beagle} 2>> {log}
		done
		"""