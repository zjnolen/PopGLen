rule ngsrelate:
	input:
		beagle=rules.merge_beagle.output.beagle,
		bamlist=rules.angsd_makeBamlist.output
	output:
		results+"/analyses/ngsrelate/"+dataset+"_{population}{dp}_relate.tsv"
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

		ngsrelate -G {input.beagle} -n $nind -L $nsites -O {output} 2>> {log}
		"""
	