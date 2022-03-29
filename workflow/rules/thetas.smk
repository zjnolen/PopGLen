rule angsd_saf2theta:
	input:
		safidx=rules.angsd_doSaf.output.safidx,
		sfs=rules.angsd_realSFS.output.sfs,
		sites=get_sites_file
	output:
		thetasidx=results + "/angsd/thetas/{population}{sites}.thetas.idx",
		thetas=results + "/angsd/thetas/{population}{sites}.thetas.gz"
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	log:
		logs + "/angsd/saf2theta/{population}{sites}.log"
	params:
		out_prefix=results + "/angsd/thetas/{population}{sites}",
		fold=config["params"]["angsd"]["fold"],
		sites=get_sites_option
	resources:
		time="12:00:00"
	shell:
		"""
		realSFS saf2theta {input.safidx} -sfs {input.sfs} -fold {params.fold} \
			{params.sites} -outname {params.out_prefix}
		"""

rule angsd_thetaStat:
	input:
		thetas=rules.angsd_saf2theta.output.thetasidx
	output:
		thetas=results + "/angsd/thetas/{population}{sites}.thetaWindows.pestPG"
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	log:
		logs + "/angsd/thetaStat/{population}{sites}.log"
	params:
		out_prefix=results + "/angsd/thetas/{population}{sites}.thetaWindows"
	resources:
		time="12:00:00"
	shell:
		"""
		thetaStat do_stat {input.thetas} -win 50000 -step 10000 \
			-outnames {params.out_prefix}
		"""