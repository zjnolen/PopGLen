rule realSFS_saf2theta:
	input:
		safidx=rules.realSFS_catsaf.output,
		sfs=rules.realSFS_1dSFS.output.sfs
	output:
		thetasidx=results+"/analyses/thetas/"+dataset+
			"_{population}{dp}.thetas.idx",
		thetas=results+"/analyses/thetas/"+dataset+
			"_{population}{dp}.thetas.gz"
	container:
		angsd_container
	log:
		logs+"/realSFS/saf2theta/"+dataset+"_{population}{dp}.log"
	params:
		out=results + "/analyses/thetas/"+dataset+
			"_{population}{dp}",
		fold=config["params"]["angsd"]["fold"]
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS saf2theta {input.safidx} -sfs {input.sfs} -fold {params.fold} \
			-outname {params.out}
		"""

rule thetaStat:
	input:
		thetas=rules.realSFS_saf2theta.output.thetasidx
	output:
		thetas=results+"/analyses/thetas/"+dataset+
			"_{population}{dp}.thetaWindows.pestPG"
	container:
		angsd_container
	log:
		logs + "/thetaStat/"+dataset+"_{population}{dp}.log"
	params:
		out=results+"/analyses/thetas/"+dataset+
			"_{population}{dp}.thetaWindows"
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		thetaStat do_stat {input.thetas} -win 50000 -step 10000 \
			-outnames {params.out}
		"""

localrules: plot_thetas

rule plot_thetas:
	input:
		expand(results+"/analyses/thetas/"+dataset+
			"_{population}{{dp}}.thetaWindows.pestPG",
			population=pop_list)
	output:
		report(expand(results+"/plots/thetas/"+dataset+"_all{{dp}}.{stat}.pdf",
			stat=["watterson","pi","tajima"]),
			category="Thetas")
	conda:
		"../envs/r.yaml"
	params:
		popnames=pop_list,
		outpre=results+"/plots/thetas/"+dataset+"_all{dp}"
	script:
		"../scripts/plot_thetas.R"