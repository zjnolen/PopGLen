rule realSFS_1dSFS:
	input:
		saf="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}.saf.idx",
		others=multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}.saf",".pos.gz",".gz")
	output:
		sfs=ensure("results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}.sfs",
			non_empty=True)
	container:
		angsd_container
	log:
		"logs/{dataset}/realSFS/1dSFS/{dataset}.{ref}_{population}{dp}.log"
	params:
		fold=config["params"]["realsfs"]["fold"]
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS {input.saf} -fold {params.fold} -P {threads} \
			> {output.sfs} 2> {log}
		"""

rule realSFS_2dSFS:
	input:
		saf1="results/datasets/{dataset}/safs/{dataset}.{ref}_{population1}{dp}.saf.idx",
		saf1_others=multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population1}{dp}.saf",".pos.gz",".gz"),
		saf2="results/datasets/{dataset}/safs/{dataset}.{ref}_{population2}{dp}.saf.idx",
		saf2_others=multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population2}{dp}.saf",".pos.gz",".gz")
	output:
		sfs=ensure("results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population1}-{population2}{dp}.sfs", non_empty=True)
	container:
		angsd_container
	log:
		"logs/{dataset}/realSFS/2dSFS/{dataset}.{ref}_{population1}-{population2}{dp}.log"
	wildcard_constraints:
		population1="|".join(
			[i for i in samples.index.tolist()] +
			[i for i in samples.population.values.tolist()]
			),
		population2="|".join(
			[i for i in samples.index.tolist()] +
			[i for i in samples.population.values.tolist()]
			)
	params:
		fold=config["params"]["realsfs"]["fold"]
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*180
	shell:
		"""
		realSFS {input.saf1} {input.saf2} -fold {params.fold} \
			-P {threads} > {output.sfs} 2> {log}
		"""

rule plot_heterozygosity:
	input:
		sfs=expand("results/datasets/{{dataset}}/analyses/sfs/{{dataset}}.{{ref}}_{sample}{{dp}}.sfs",
			sample=samples.index),
		popfile="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list"
	output:
		"results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all{dp}_heterozygosity.tsv",
		report("results/datasets/{dataset}/plots/heterozygosity/{dataset}.{ref}_all{dp}_heterozygosity.pdf",
			category="Heterozygosity",
			labels={
				"Topic":"Heterozygosity",
				"Subsampling":"{dp}",
				"Type":"boxplot"
			})
	log:
		"logs/{dataset}/heterozygosity/{dataset}.{ref}_all{dp}_calc-plot.log"
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_heterozygosity.R"