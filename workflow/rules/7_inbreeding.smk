localrules: convert_ibd

rule ngsf_hmm:
	input:
		beagle=rules.merge_pruned_beagles.output.beagle,
		bamlist=rules.angsd_makeBamlist.output
	output:
		#indf=results + "/inbreeding/{population}.indF",
		ibd=results+"/analyses/ngsF-HMM/"+dataset+
			"_{population}{dp}.ibd",
		indF=results+"/analyses/ngsF-HMM/"+dataset+
			"_{population}{dp}.indF",
		pos=results+"/analyses/ngsF-HMM/"+dataset+
			"_{population}{dp}.pos"
	log:
		logs + "/ngsF-HMM/"+dataset+"_{population}{dp}.log"
	params:
		out=results + "/analyses/ngsF-HMM/"+dataset+
			"_{population}{dp}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*360
	shell:
		r"""
		module load bioinfo-tools
		module load ngsF-HMM/20200720-2df9690

		zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos} 2> {log}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')

		ngsF-HMM.sh --geno {input.beagle} --n_ind $nind \
			--n_sites $nsites --pos {output.pos} --lkl \
			--out {params.out} &>> {log}
		"""

rule convert_ibd:
	input:
		ibd=rules.ngsf_hmm.output.ibd,
		pos=rules.ngsf_hmm.output.pos,
		inds=rules.popfile.output.inds
	output:
		roh=results + "/analyses/ngsF-HMM/"+dataset+"_{population}{dp}.roh"
	log:
		logs + "/ngsF-HMM/"+dataset+"_{population}{dp}_convert_ibd.log"
	shell:
		"""
		perl workflow/scripts/convert_ibd_mod.pl --ind_file {input.inds} \
			--pos_file {input.pos} --ibd_pos_file {input.ibd} > {output.roh} \
			2> {log}
		"""