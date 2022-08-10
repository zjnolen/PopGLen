localrules: convert_ibd, plot_froh

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
	threads: lambda wildcards, attempt: attempt*10
	resources:
		time=lambda wildcards, attempt: attempt*2880
	shell:
		r"""
		module load bioinfo-tools
		module load ngsF-HMM/20200720-2df9690

		zcat {input.beagle} | awk '{{print $1}}' | sed 's/\(.*\)_/\1\t/' \
			| tail -n +2 > {output.pos} 2> {log}
		
		nsites=$(cat {output.pos} | wc -l)

		nind=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')

		export TMPDIR={resources.tmpdir}

		workflow/scripts/ngsF-HMM.sh --geno {input.beagle} --n_ind $nind \
			--n_sites $nsites --pos {output.pos} --lkl \
			--out {params.out} --n_threads {threads} &>> {log}
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
		inds=$(basename {input.inds})
		tail -n +2 {input.inds} > {resources.tmpdir}/$inds.tmp

		perl workflow/scripts/convert_ibd_mod.pl --pos_file {input.pos} \
			--ind_file {resources.tmpdir}/$inds.tmp \
			--ibd_pos_file {input.ibd} > {output.roh} 2> {log}
		"""

rule plot_froh:
	input:
		roh=expand(results+"/analyses/ngsF-HMM/"+dataset+
			"_{population}{{dp}}.roh", population=pop_list),
		inds=results+"/genotyping/pop_lists/"+dataset+"_all.indiv.list",
		autos=rules.sexlink_sum.output.sum
	output:
		report(expand(results+"/plots/inbreeding/"+dataset+
			"_all{{dp}}.{stat}.pdf",
			stat=["froh","rohreg"]),
			category="Inbreeding")
	conda:
		"../envs/r.yaml"
	params:
		popnames=pop_list,
		outpre=results+"/plots/inbreeding/"+dataset+"_all{dp}"
	script:
		"../scripts/plot_Froh.R"