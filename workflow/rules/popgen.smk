rule ngsLD_estLD:
	input:
		beagle=results + "/angsd/beagle/{population}_{chrom}_md{miss}.beagle.gz",
		bamlist=results + "/angsd/bamlists/{population}.bamlist"
	output:
		ld=results + "/ngsLD/{population}_{chrom}_md{miss}.ld",
		pos=results + "/ngsLD/{population}_{chrom}_md{miss}.pos"
	log:
		"logs/ngsLD/{population}_{chrom}_md{miss}.log"
	container:
		"library://james-s-santangelo/ngsld/ngsld:1.1.1"
	threads: 2
	resources:
		time="12:00:00"
	shell:
		r"""
		zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos} 2> {log}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1 + 1}}')

		ngsLD --geno {input.beagle} --n_ind $nind --n_sites $nsites \
			--pos {output.pos} --probs --out {output.ld} &>> {log}
		"""

rule ngsLD_prune_sites:
	input:
		ld=results + "/ngsLD/{population}_{chrom}_md{miss}.ld"
	output:
		sites=results + "/ngsLD/{population}_{chrom}_md{miss}_pruned.sites"
	log:
		"logs/ngsLD/{population}_{chrom}_md{miss}_prune.log"
	conda:
		"../envs/pruning.yaml"
	threads: lambda wildcards, input: ((input.size_mb * 0.8) /6800)+1
	resources:
		time="24:00:00"
	shell:
		"workflow/scripts/prune_ngsLD.py --input {input.ld} --output {output.sites} &> {log}"

rule prune_chrom_beagle:
	input:

	output:
		beagle=results + "/beagle/{population}_genome_md{miss}_pruned.beagle.gz"
	log:
		"logs/ngsLD/{population}_genome_md{miss}_prune.log"
	shell:
		"zcat {input.beagle} | tail -n +2"