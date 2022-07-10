rule realSFS_1dSFS:
	input:
		saf=rules.realSFS_catsaf.output
	output:
		sfs=results + "/analyses/sfs/"+dataset+"_{population}{dp}.sfs"
	container:
		"docker://zjnolen/angsd:0.937"
	log:
		logs + "/realSFS/1dSFS/"+dataset+"_{population}{dp}.log"
	params:
		fold=config["params"]["angsd"]["fold"]
	threads: lambda wildcards, attempt: attempt*5
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS {input.saf} -fold {params.fold} -P {threads} \
			> {output.sfs} 2> {log}
		"""

rule realSFS_2dSFS:
	input:
		saf1=results+"/genotyping/saf/"+dataset+
			"_{population1}{dp}.saf.idx",
		saf2=results+"/genotyping/saf/"+dataset+
			"_{population2}{dp}.saf.idx"
	output:
		sfs=results + "/analyses/sfs/"+dataset+
			"_{population1}-{population2}{dp}.sfs"
	container:
		"docker://zjnolen/angsd:0.937"
	log:
		logs + "/realSFS/2dSFS/"+dataset+"_{population1}-{population2}{dp}.log"
	params:
		fold=config["params"]["angsd"]["fold"]
	threads: lambda wildcards, attempt: attempt*5
	resources:
		time=lambda wildcards, attempt: attempt*180
	shell:
		"""
		realSFS {input.saf1} {input.saf2} -fold {params.fold} -P {threads} \
			> {output.sfs} 2> {log}
		"""

localrules: plot_heterozygosity

rule plot_heterozygosity:
	input:
		sfs=expand(results+
			"/analyses/sfs/"+dataset+"_{sample}{{dp}}.sfs",
			sample=samples.index),
		popfile=results+"/genotyping/pop_lists/"+dataset+"_all.indiv.list"
	output:
		report(results+"/plots/heterozygosity/"+dataset+
			"_all{dp}_heterozygosity.pdf",
			category="Heterozygosity")
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_heterozygosity.R"