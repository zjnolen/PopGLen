rule pca_pcangsd:
	input:
		beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_pruned.beagle.gz"
	output:
		cov="results/datasets/{dataset}/analyses/pcangsd/{dataset}.{ref}_all{dp}.cov"
	log:
		"logs/{dataset}/pcangsd/{dataset}.{ref}_all{dp}.log"
	container:
		pcangsd_container
	params:
		prefix=lambda w, output: os.path.splitext(output.cov)[0]
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		pcangsd -b {input.beagle} -o {params.prefix} &> {log}
		"""

rule plot_pca:
	input:
		"results/datasets/{dataset}/analyses/pcangsd/{dataset}.{ref}_all{dp}.cov",
		"results/datasets/{dataset}/poplists/{dataset}_all.indiv.list"
	output:
		report(
			"results/datasets/{dataset}/plots/pca/{dataset}.{ref}_all{dp}_pc{xpc}-{ypc}.svg",
			category="PCA",
			labels={
				"Topic":"PCA",
				"PCs":"PC{xpc}-PC{ypc}",
				"Subsampling":"{dp}",
				"Type":"scatterplot"
			})
	log:
		"logs/{dataset}/pcangsd/{dataset}.{ref}_all{dp}_pc{xpc}-{ypc}_plot.log"
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_PCA.R"