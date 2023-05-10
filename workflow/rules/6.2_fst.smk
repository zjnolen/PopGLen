rule realSFS_fst_index:
	input:
		saf1="results/datasets/{dataset}/safs/{dataset}.{ref}_{population1}{dp}.saf.idx",
		saf1_others=multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population1}{dp}.saf",".pos.gz",".gz"),
		saf2="results/datasets/{dataset}/safs/{dataset}.{ref}_{population2}{dp}.saf.idx",
		saf2_others=multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population2}{dp}.saf",".pos.gz",".gz"),
		sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}.sfs"
	output:
		fstidx="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}.fst.idx"
	container:
		angsd_container
	log:
		"logs/{dataset}/realSFS/fst/index/{dataset}.{ref}_{population1}-{population2}{dp}.log"
	params:
		out=lambda w, output: output.fstidx.removesuffix('.fst.idx'),
		fst=config["params"]["fst"]["whichFst"]
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS fst index -whichFst {params.fst} \
			{input.saf1} {input.saf2} -sfs {input.sfs} \
			-fstout {params.out} &> {log}
		"""

rule realSFS_fst_stats:
	input:
		fstidx="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}.fst.idx"
	output:
		fstglob="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}.fst.global"
	container:
		angsd_container
	log:
		"logs/{dataset}/realSFS/fst/stats/{dataset}.{ref}_{population1}-{population2}{dp}.log"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		realSFS fst stats {input.fstidx} | \
			awk '{{print "{wildcards.population1}\t{wildcards.population2}\t"\
			$1"\t"$2}}' > {output.fstglob} 2> {log}
		"""

# rule realSFS_fst_stats2:
# 	input:
# 		fstidx="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}.fst.idx"
# 	output:
# 		fstwin="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{population1}-{population2}{dp}.fst.window"
# 	container:
# 		angsd_container
# 	log:
# 		"logs/{dataset}/realSFS/fst/stats2/{dataset}.{ref}_{population1}-{population2}{dp}.log"
# 	shell:
# 		"""
# 		realSFS fst stats2 {input.fstidx} -win 50000 -step 10000 \
# 			> {output.fstwin} 2> {log}
# 		"""

rule aggregate_global_fst:
	input:
		get_fst
	output:
		"results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{unit}pairs{dp}.fst.sum"
	log:
		"logs/{dataset}/realSFS/fst/aggregate/{dataset}.{ref}_{unit}pairs{dp}.log"
	conda:
		"../envs/shell.yaml"
	wildcard_constraints:
		unit="ind|pop"
	shell:
		"""
		(echo "pop1\tpop2\tunweight.fst\tweight.fst" > {output}
		cat {input} >> {output}) 2> {log}
		"""

rule plot_fst:
	input:
		"results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_{unit}pairs{dp}.fst.sum"
	output:
		report("results/datasets/{dataset}/plots/fst/{dataset}.{ref}_{unit}pairs{dp}_fst.pdf",
			category="Fst",
			labels={
				"Topic":"Pairwise Fst",
				"Unit":"{unit}",
				"Subsampling":"{dp}",
				"Type":"heatmap"
			})
	log:
		"logs/{dataset}/realSFS/fst/plot/{dataset}.{ref}_{unit}pairs{dp}.log"
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_fst.R"