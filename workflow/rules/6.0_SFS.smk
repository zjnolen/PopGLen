rule realSFS_1dSFS:
	input:
		saf=rules.realSFS_catsaf.output
	output:
		sfs=ensure(results + "/analyses/sfs/"+dataset+"_{population}{dp}.sfs",
			non_empty=True)
	container:
		angsd_container
	log:
		logs + "/realSFS/1dSFS/"+dataset+"_{population}{dp}.log"
	params:
		fold=config["params"]["realsfs"]["fold"]
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS {input.saf[0]} -fold {params.fold} -P {threads} \
			> {output.sfs} 2> {log}
		"""

rule realSFS_2dSFS:
	input:
		saf1=multiext(results+"/genotyping/saf/"+dataset+
			"_{population1}{dp}.saf",".idx",".pos.gz",".gz"),
		saf2=multiext(results+"/genotyping/saf/"+dataset+
			"_{population2}{dp}.saf",".idx",".pos.gz",".gz")
	output:
		sfs=ensure(results + "/analyses/sfs/"+dataset+
			"_{population1}-{population2}{dp}.sfs", non_empty=True)
	container:
		angsd_container
	log:
		logs + "/realSFS/2dSFS/"+dataset+"_{population1}-{population2}{dp}.log"
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
		realSFS {input.saf1[0]} {input.saf2[0]} -fold {params.fold} \
			-P {threads} > {output.sfs} 2> {log}
		"""

rule plot_heterozygosity:
	input:
		sfs=expand(results+
			"/analyses/sfs/"+dataset+"_{sample}{{dp}}.sfs",
			sample=samples.index),
		popfile=results+"/genotyping/pop_lists/"+dataset+"_all.indiv.list"
	output:
		results+"/analyses/heterozygosity/"+dataset+
			"_all{dp}_heterozygosity.tsv",
		report(results+"/plots/heterozygosity/"+dataset+
			"_all{dp}_heterozygosity.pdf",
			category="Heterozygosity",
			labels={
				"Topic":"Heterozygosity",
				"Subsampling":"{dp}",
				"Type":"boxplot"
			})
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_heterozygosity.R"