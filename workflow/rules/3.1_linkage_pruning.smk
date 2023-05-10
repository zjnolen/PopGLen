rule ngsLD_estLD:
	input:
		beagle="results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.beagle.gz",
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist"
	output:
		ld=temp("results/datasets/{dataset}/beagles/pruned/ngsLD/{dataset}{ref}_{population}{dp}_chunk{chunk}.ld.gz"),
		pos="results/datasets/{dataset}/beagles/pruned/ngsLD/{dataset}{ref}_{population}{dp}_chunk{chunk}.pos"
	log:
		"logs/{dataset}/ngsLD/estLD/{dataset}{ref}_{population}{dp}_chunk{chunk}.log"
	container:
		"library://james-s-santangelo/ngsld/ngsld:1.1.1"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*720
	shell:
		r"""
		(zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		if [[ $nsites == 0 ]]; then
			touch {output.ld}
		else
			ngsLD --geno {input.beagle} --n_ind $nind --n_sites $nsites \
				--pos {output.pos} --probs --n_threads {threads} \
				| gzip > {output.ld}
		fi) 2> {log}
		"""


rule ngsLD_prune_sites:
	input:
		ld="results/datasets/{dataset}/beagles/pruned/ngsLD/{dataset}{ref}_{population}{dp}_chunk{chunk}.ld.gz",
		pos="results/datasets/{dataset}/beagles/pruned/ngsLD/{dataset}{ref}_{population}{dp}_chunk{chunk}.pos"
	output:
		sites="results/datasets/{dataset}/beagles/pruned/ngsLD/{dataset}{ref}_{population}{dp}_chunk{chunk}_pruned.sites"
	log:
		"logs/{dataset}/ngsLD/prune_sites/{dataset}{ref}_{population}{dp}_chunk{chunk}.log"
	conda:
		"../envs/pruning.yaml"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*1440
	shell:
		"""
		(nsites=$(cat {input.pos} | wc -l)

		if [[ $nsites == 0 ]]; then
			touch {output.sites}
		else
			workflow/scripts/prune_ngsLD.py --input {input.ld} --max_dist 50000 \
				--min_weight 0.1 --output {output.sites}
		fi) 2> {log}
		"""

rule prune_chunk_beagle:
	input:
		beagle="results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.beagle.gz",
		sites="results/datasets/{dataset}/beagles/pruned/ngsLD/{dataset}{ref}_all{dp}_chunk{chunk}_pruned.sites"
	output:
		prunedgz="results/datasets/{dataset}/beagles/pruned/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}_pruned.beagle.gz"
	log:
		"logs/{dataset}/ngsLD/prune_beagle/{dataset}{ref}_{population}{dp}_chunk{chunk}.log"
	conda:
		"../envs/shell.yaml"
	shadow: "minimal"
	threads: lambda wildcards, attempt: attempt*10
	params:
		pruned="results/datasets/{dataset}/beagles/pruned/chunks/{dataset}{ref}_{population}{dp}_chunk{chunk}_pruned.beagle"
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		r"""
		(set +o pipefail;
		zcat {input.beagle} | head -n 1 > {params.pruned}
		
		join -t $'\t' <(sort -k1,1 {input.sites}) <(zcat {input.beagle} | sort -k1,1) | \
			sed 's/_/\t/' | sort -k1,1 -k2,2n | sed 's/\t/_/' \
			>> {params.pruned}

		gzip -c {params.pruned} > {output.prunedgz}

		Nsites=$(cat {input.sites} | wc -l | awk '{{print $1+1}}')
		NsitesB=$(zcat {output.prunedgz} | wc -l)

		echo "Sites searched for: $Nsites"
		echo "Sites in pruned beagle: $NsitesB") &> {log}
		"""

rule merge_pruned_beagles:
	input:
		pruned=lambda w: expand("results/datasets/{{dataset}}/beagles/pruned/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_pruned.beagle.gz", chunk=chunklist)
	output:
		beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_{population}{dp}_pruned.beagle.gz"
	log:
		"logs/{dataset}/ngsLD/merge_pruned/{dataset}.{ref}_{population}{dp}_merge_pruned.log"
	wildcard_constraints:
		population="|".join(
			["all"] +
			[i for i in samples.index.tolist()] +
			[i for i in samples.population.values.tolist()] +
			[i for i in samples.depth.values.tolist()]
			)
	conda:
		"../envs/shell.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		(echo "cat file order:"
		echo {input.pruned} | tr ' ' '\n'

		echo "Printing header to beagle..."
		set +o pipefail;
		zcat {input.pruned} | head -n 1 | gzip > {output.beagle}

		echo "Adding requested chunks to beagle..."
		for f in {input.pruned}; do
			zcat $f | tail -n +2
		done | gzip >> {output.beagle}) &> {log}
		"""

rule remove_excl_pca_admix:
	input:
		"results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_pruned.beagle.gz"
	output:
		"results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all_excl_pca-admix{dp}_pruned.beagle.gz"
	log:
		"logs/{dataset}/ngsLD/excl_pca_admix_beagle/{dataset}.{ref}_all_excl_pca-admix{dp}.log"
	conda:
		"../envs/shell.yaml"
	params:
		remcols=get_excl_ind_cols
	shell:
		"""
		zcat {input} | cut -f{params.remcols} --complement | gzip > {output} 2> {log}
		"""