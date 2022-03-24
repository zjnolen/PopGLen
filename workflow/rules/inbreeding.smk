rule ngsf_hmm:
	input:
		beagle=rules.merge_pruned_beagles.output.beagle,
		bamlist=rules.angsd_makeBamlist.output
	output:
		#indf=results + "/inbreeding/{population}.indF",
		idb=results + "/inbreeding/{population}.ibd",
		indF=results + "/inbreeding/{population}.indF",
		pos=results + "/inbreeding/{population}.pos"
	log:
		logs + "/ngsF-HMM/{population}.log"
	params:
		out_prefix=results + "/inbreeding/{population}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time="06:00:00"
	shell:
		r"""
		module load bioinfo-tools
		module load ngsF-HMM/20200720-2df9690

		zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos} 2> {log}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1 + 1}}')

		ngsF-HMM.sh --geno {input.beagle} --n_ind $nind \
			--n_sites $nsites --pos {output.pos} --lkl \
			--out {params.out_prefix} &>> {log}
		"""

