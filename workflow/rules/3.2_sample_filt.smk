rule ngsrelate:
	input:
		beagle=rules.merge_beagle.output.beagle,
		bamlist=rules.angsd_makeBamlist.output,
		inds=rules.popfile.output.inds
	output:
		relate=results+"/analyses/ngsrelate/"+dataset+
			"_{population}{dp}_relate.tsv",
		samples=results+"/analyses/ngsrelate/"+dataset+
			"_{population}{dp}_samples.list"
	log:
		logs + "/ngsrelate/"+dataset+"_{population}{dp}.log"
	threads: lambda wildcards, attempt: attempt*4
	resources:
		time=lambda wildcards, attempt: attempt*360
	shell:
		r"""
		module load bioinfo-tools NgsRelate

		nsites=$(zcat {input.beagle} | tail -n +2 | wc -l) 2>> {log}

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}') 2>> {log}

		echo "nsites nind" >> {log}

		echo $nsites $nind >> {log}

		cut -f1 {input.inds} | tail -n +2 > {output.samples} 2>> {log}

		ngsrelate -G {input.beagle} -n $nind -L $nsites -O {output.relate} \
			-z {output.samples} 2>> {log}
		"""
	