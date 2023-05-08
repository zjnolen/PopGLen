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
	container:
		ngsf_hmm_container
	params:
		out=results + "/analyses/ngsF-HMM/"+dataset+
			"_{population}{dp}"
	threads: lambda wildcards, attempt: attempt*10
	resources:
		time=lambda wildcards, attempt: attempt*2880
	shell:
		r"""
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
		roh=results+"/analyses/ngsF-HMM/"+dataset+"_{population}{dp}.roh"
	log:
		logs + "/ngsF-HMM/"+dataset+"_{population}{dp}_convert_ibd.log"
	container:
		ngsf_hmm_container
	shadow:
		"copy-minimal"
	shell:
		"""
		tail -n +2 {input.inds} > {input.inds}.tmp 2> {log}

		# convert_ibd.pl drops warnings when chromosomes have non-numeric 
		# names, and may not handle them well. As a work around, a numeric 
		# index of chromosome names is created, then the names replaced with 
		# these numbers. Then returns their names after the conversion is 
		# complete.

		# first, create the index file
		n_contigs=$(awk '{{print $1}}' {input.pos} | uniq | wc -l) 2>> {log}
		awk '{{print $1}}' {input.pos} | uniq > {input.pos}.contigs 2>> {log}
		seq $n_contigs > {input.pos}.contigs.idx 2>> {log}
		
		# create index based pos file, adapted from:
		# Glenn Jackman - https://stackoverflow.com/a/7198895
		awk '
    		FILENAME == ARGV[1] {{ listA[$1] = FNR; next }}
    		FILENAME == ARGV[2] {{ listB[FNR] = $1; next }}
    		{{
        		for (i = 1; i <= NF; i++) {{
            		if ($1 in listA) {{
               			$1 = listB[listA[$i]]
            		}}
        		}}
        		print
    		}}' \
		{input.pos}.contigs {input.pos}.contigs.idx {input.pos} \
			> {input.pos}.tmp 2>> {log}

		convert_ibd.pl --pos_file {input.pos}.tmp \
			--ind_file {input.inds}.tmp	--ibd_pos_file {input.ibd} \
			> {output.roh}.tmp 2>> {log}
		
		# Revert back to original contig name:
		awk '
    		FILENAME == ARGV[1] {{ listA[$1] = FNR; next }}
    		FILENAME == ARGV[2] {{ listB[FNR] = $1; next }}
    		{{
        		for (i = 1; i <= NF; i++) {{
            		if ($1 in listA) {{
               			$1 = listB[listA[$i]]
            		}}
        		}}
        		print
    		}}' \
		{input.pos}.contigs.idx {input.pos}.contigs {output.roh}.tmp \
			> {output.roh} 2>> {log}
		"""

def get_auto_sum(wildcards):
	if config["reference"]["sex-linked"] or config["reference"]["exclude"] \
		or config["reference"]["mito"]:
		return REF_DIR+"/beds/"+REF_NAME+"_excl.bed.sum"
	else:
		return REF_DIR+"/beds/"+REF_NAME+"_genome.bed"

rule plot_froh:
	input:
		roh=expand(results+"/analyses/ngsF-HMM/"+dataset+
			"_{population}{{dp}}.roh", population=pop_list),
		inds=results+"/genotyping/pop_lists/"+dataset+"_all.indiv.list",
		autos=get_auto_sum
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