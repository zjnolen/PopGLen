rule pca_pcangsd:
	input:
		beagle=results+"/genotyping/pruned_beagle/"+dataset+
			"_{population}{dp}_pruned.beagle.gz"
	output:
		cov=results+"/analyses/pcangsd/"+dataset+"_{population}{dp}.cov",
		kin=results+"/analyses/pcangsd/"+dataset+"_{population}{dp}.kinship.npy"
	log:
		logs + "/pcangsd/"+dataset+"_{population}{dp}.log"
	params:
		prefix=results + "/analyses/pcangsd/"+dataset+"_{population}{dp}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		module load bioinfo-tools
		module load PCAngsd/0.982

		pcangsd.py -b {input.beagle} -o {params.prefix} &> {log}
		"""