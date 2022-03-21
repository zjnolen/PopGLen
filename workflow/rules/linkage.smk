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
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time="24:00:00"
	shell:
		"""
		workflow/scripts/prune_ngsLD.py --input {input.ld} > {output.sites} \
			2> {log}
		"""

rule prune_chrom_beagle:
	input:
		beagle=rules.angsd_chrom_beagle.output.beagle,
		sites=rules.ngsLD_prune_sites.output.sites
	output:
		pruned=results+"/angsd/beagle/pruned/chrom/{population}_chr{chrom}_pruned.beagle.gz",
		sites=temp(results+"/ngsLD/prune_sites/{population}_chr{chrom}_pruned.sites.tab")
	log:
		logs + "/ngsLD/prune_beagle/{population}_chr{chrom}.log"
	threads: lambda wildcards, attempt: attempt*2
	shell:
		r"""
		echo "The beagle file produced by this will have no header." > {log}
		echo "It will be added in when the files are cat together." >> {log}

		# Add trailing tab to sites file ensuring only exact matches are
		# kept.
		sed "s/$/\t/g" {input.sites} > {output.sites} 2>> {log}

		zgrep -f {output.sites} {input.beagle} | gzip > {output.pruned} \
			2>> {log}
		
		Nsites=$(cat {input.sites} | wc -l) &>> {log}
		NsitesB=$(zcat {output.pruned} | wc -l) &>> {log}

		if [ $Nsites = $NsitesB ]; then
			echo "Pruning successful!" &>> {log}
		else
			echo "Number of sites in pruned beagle differ from sites file." \
				&>> {log}
			echo "Please check to make sure sites file has format" \
				"chrom_position\t for each site, one per line." &>> {log}
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

		cat {input.pruned} >> {output.beagle}
		"""