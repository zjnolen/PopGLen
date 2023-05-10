rule realSFS_saf2theta:
	input:
		safidx="results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}.saf.idx",
		saf_others=multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}.saf",".pos.gz",".gz"),
		sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{population}{dp}.sfs"
	output:
		thetasidx="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}.thetas.idx",
		thetas="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}.thetas.gz"
	container:
		angsd_container
	log:
		"logs/{dataset}/realSFS/saf2theta/{dataset}.{ref}_{population}{dp}.log"
	params:
		out=lambda w, output: output.thetas.removesuffix('.thetas.gz'),
		fold=config["params"]["realsfs"]["fold"]
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS saf2theta {input.safidx} -sfs {input.sfs} -fold {params.fold} \
			-outname {params.out} &> {log}
		"""

rule thetaStat:
	input:
		thetasidx="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}.thetas.idx",
		thetas="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}.thetas.gz"
	output:
		thetas="results/datasets/{dataset}/analyses/thetas/{dataset}.{ref}_{population}{dp}.thetaWindows.pestPG"
	container:
		angsd_container
	log:
		"logs/{dataset}/thetaStat/{dataset}.{ref}_{population}{dp}.log"
	params:
		out=lambda w, output: os.path.splitext(output.thetas)[0],
		win_size=config["params"]["thetas"]["win_size"],
		win_step=config["params"]["thetas"]["win_step"],
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		thetaStat do_stat {input.thetasidx} -win {params.win_size} -step {params.win_step} \
			-outnames {params.out} &> {log}
		"""

rule plot_thetas:
	input:
		expand("results/datasets/{{dataset}}/analyses/thetas/{{dataset}}.{{ref}}_{population}{{dp}}.thetaWindows.pestPG",
			population=pop_list)
	output:
		report(expand("results/datasets/{{dataset}}/plots/thetas/{{dataset}}.{{ref}}_all{{dp}}.{stat}.pdf",
			stat=["watterson","pi","tajima"]),
			category="Thetas")
	log:
		"logs/{dataset}/thetaStat/{dataset}.{ref}_all{dp}_plot.log"
	conda:
		"../envs/r.yaml"
	params:
		popnames=pop_list,
		outpre=lambda w, output: output[0].removesuffix('.watterson.pdf')
	script:
		"../scripts/plot_thetas.R"