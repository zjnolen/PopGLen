rule angsd_saf2theta:
	input:
		safidx=rules.angsd_doSaf.output.safidx,
		sfs=rules.angsd_realSFS.output.sfs
	output:
		thetasidx=results + "/angsd/thetas/{population}.thetas.idx",
		thetas=results + "/angsd/thetas/{population}.thetas.gz"
	conda:
		"../envs/angsd.yaml"
	log:
		logs + "/angsd/saf2theta/{population}.log"
	params:
		out_prefix=results + "/angsd/thetas/{population}",
		fold=config["params"]["angsd"]["fold"]
	resources:
		time="12:00:00"
	shell:
		"""
		realSFS saf2theta {input.safidx} -sfs {input.sfs} -fold {params.fold} \
			-outname {params.out_prefix}
		"""

rule angsd_thetaStat:
	input:
		thetas=rules.angsd_saf2theta.output.thetasidx
	output:
		thetas=results + "/angsd/thetas/{population}.thetaWindows.pestPG"
	conda:
		"../envs/angsd.yaml"
	log:
		logs + "/angsd/thetaStat/{population}.log"
	params:
		out_prefix=results + "/angsd/thetas/{population}.thetaWindows"
	resources:
		time="12:00:00"
	shell:
		"""
		thetaStat do_stat {input.thetas} -win 50000 -step 10000 \
			-outnames {params.out_prefix}
		"""