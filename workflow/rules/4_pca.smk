rule pca_pcangsd:
	input:
		beagle=results+"/genotyping/pruned_beagle/"+dataset+
			"_{population}{dp}_pruned.beagle.gz"
	output:
		cov=results+"/analyses/pcangsd/"+dataset+"_{population}{dp}.cov"
	log:
		logs + "/pcangsd/"+dataset+"_{population}{dp}.log"
	container:
		pcangsd_container
	params:
		prefix=results + "/analyses/pcangsd/"+dataset+"_{population}{dp}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		pcangsd -b {input.beagle} -o {params.prefix} &> {log}
		"""

localrules: plot_pca

rule plot_pca:
	input:
		rules.pca_pcangsd.output.cov,
		rules.popfile.output.inds
	output:
		report(
			results+"/plots/pca/"+dataset+
				"_{population}{dp}_pc{xpc}-{ypc}.pdf",
			category="PCA",
			labels={
				"Topic":"PCA",
				"PCs":"PC{xpc}-PC{ypc}",
				"Subsampling":"{dp}",
				"Type":"scatterplot"
			})
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_PCA.R"