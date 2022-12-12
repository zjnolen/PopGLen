# Per individual estimates of endogenous content - primarily useful for 
# historical samples or when using references that are not specific to 
# the species of interest

def get_endo_cont_stat(wildcards):
    # Gets bam file for endogenous content calculation
    s = wildcards.sample
    if s in samples.index[samples.time == "modern"]:
        return "results/mapping/mapped/"+s+".paired.flagstat"
    elif s in samples.index[samples.time == "historical"]:
        return "results/mapping/mapped/"+s+".merged.flagstat"

rule endo_cont:
    input:
        get_endo_cont_stat
    output:
        "results/mapping/qc/endogenous_content/{sample}.endo"
    shell:
        r"""
        total=$(grep -E "^[0-9]+ \+ [0-9]+ in total" {input} \
            | awk '{{print $1}}')
        mapped=$(grep -E "^[0-9]+ \+ [0-9]+ mapped" {input} \
            | awk '{{print $1}}')
        primary=$(grep -E "^[0-9]+ \+ [0-9]+ primary mapped" {input} \
            | awk '{{print $1}}')
        
        echo $total $mapped $primary {wildcards.sample} | \
            awk '{{printf "%s\t%.3f\t%.3f\n",$4,$2/$1*100,$3/$1*100}}' \
            > {output}
        """

rule compile_endo_cont:
	input:
		lambda w: expand("results/mapping/qc/endogenous_content/{sample}.endo",
			sample=get_samples_from_pop("all"))
	output:
		results+"/qc/endogenous_content/"+dataset+"_all.endo.tsv"
	resources:
		time=lambda wildcards, attempt: attempt*15
	shell:
		"""
		echo "sample	perc.endo	perc.prim.endo" > {output}
		cat {input} >> {output}
		"""

# Pairwise individual relatedness with R0, R1, KING-robust kinship 
# method from Waples et al. 2019, MolEcol

rule est_kinship_stats:
	input:
		sfs=results + "/analyses/sfs/"+dataset+
			"_{ind1}-{ind2}{dp}.sfs"
	output:
		results+"/analyses/kinship/waples2019/"+dataset+
			"_{ind1}-{ind2}{dp}.kinship"
	log:
		logs + "/kinship/waples2019/"+dataset+"_{ind1}-{ind2}{dp}_kinship.log"
	wildcard_constraints:
		ind1="|".join(
			[i for i in samples.index.tolist()]
			),
		ind2="|".join(
			[i for i in samples.index.tolist()]
			)
	conda:
		"../envs/r.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*15
	script:
		"../scripts/kinship.R"

def get_kinship(wildcards):
	combos = list(itertools.combinations(samples.index, 2))
	# sort inds alphebetically, this ensures that should new inds be added
	# after generating some SFS, the reordering of the combinations won't
	# lead to generating identical SFS with the individuals swapped
	combos = [sorted(pair) for pair in combos]
	ind1 = [pair[0] for pair in combos]
	ind2 = [pair[1] for pair in combos]
	return expand(results+"/analyses/kinship/waples2019/"+dataset+
				"_{ind1}-{ind2}"+wildcards.dp+".kinship",
				zip, ind1=ind1, ind2=ind2)

rule compile_kinship_stats:
	input:
		get_kinship
	output:
		results+"/analyses/kinship/waples2019/"+dataset+"_all{dp}.kinship"
	resources:
		time=lambda wildcards, attempt: attempt*15
	shell:
		"""
		echo "ind1	ind2	R0	R1	KING" > {output}
		cat {input} >> {output}
		"""

rule kinship_table_html:
	input:
		results+"/analyses/kinship/waples2019/"+dataset+"_all{dp}.kinship"
	output:
		report(results+"/analyses/kinship/waples2019/"+dataset+
					"_all{dp}.kinship.html",
				category="Kinship",
				labels={
					"Topic":"Waples et al. 2019 R0,R1,KING",
					"Type":"Table"
				})
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"

# Individual depth histograms

rule ind_unfiltered_depth:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		anc=REF,
		ref=REF
	output:
		sample_hist="results/mapping/qc/ind_depth/unfiltered/"+dataset+ \
			"_{population}{dp}.depthSample",
		global_hist=temp("results/mapping/qc/ind_depth/unfiltered/"+dataset+ \
			"_{population}{dp}.depthGlobal"),
		arg="results/mapping/qc/ind_depth/unfiltered/"+dataset+ \
			"_{population}{dp}.arg"
	log:
		"logs/mapping/ind_depth/unfiltered/"+dataset+"_{population}{dp}.log"
	container:
		angsd_container
	params:
		out="results/mapping/qc/ind_depth/unfiltered/"+dataset+ \
			"_{population}{dp}",
		maxdepth=config["params"]["angsd"]["maxdepth"]
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
			-bam {input.bamlist} -nThreads {threads} -out {params.out} &> {log}
		"""

rule ind_filtered_depth:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		anc=REF,
		ref=REF,
		sites=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites.idx"
	output:
		sample_hist=results+"/qc/ind_depth/filtered/"+dataset+ \
			"_{population}{dp}.depthSample",
		global_hist=temp(results+"/qc/ind_depth/filtered/"+dataset+ \
			"_{population}{dp}.depthGlobal"),
		arg=results+"/qc/ind_depth/filtered/"+dataset+ \
			"_{population}{dp}.arg"
	log:
		logs+"/ind_depth/filtered/"+dataset+"_{population}{dp}.log"
	container:
		angsd_container
	params:
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/qc/ind_depth/filtered/"+dataset+ \
			"_{population}{dp}",
		maxdepth=config["params"]["angsd"]["maxdepth"]
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
			-bam {input.bamlist} -ref {input.ref} -nThreads {threads} \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-sites {input.sites} -out {params.out} &> {log}
		"""

def get_total_bed(wildcards):
	if wildcards.prefix == "results/mapping/qc/ind_depth/unfiltered/":
		return REF_DIR+"/beds/"+REF_NAME+"_genome.bed"
	elif wildcards.prefix == results+"/qc/ind_depth/filtered/":
		return results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.bed"

rule summarize_ind_depth:
	input:
		sample_hist="{prefix}"+dataset+"_{sample}{dp}.depthSample",
		bed=get_total_bed
	output:
		sample_summ="{prefix}"+dataset+"_{sample}{dp}.depth.sum"
	conda:
		"../envs/r.yaml"
	wildcard_constraints:
		prefix="results/mapping/qc/ind_depth/unfiltered/|"+
			results+"/qc/ind_depth/filtered/"
	threads: lambda wildcards, attempt: attempt
	script:
		"../scripts/calc_depth.R"

def get_depth_header(wildcards):
	if wildcards.prefix == "results/mapping/qc/ind_depth/unfiltered/":
		return "genome"
	elif wildcards.prefix == results+"/qc/ind_depth/filtered/":
		return "filt"

rule merge_ind_depth:
	input:
		depth=lambda w: expand("{{prefix}}"+dataset+
			"_{sample}{{dp}}.depthSample",
			sample=get_samples_from_pop("all")),
		summary=lambda w: expand("{{prefix}}"+dataset+ \
			"_{sample}{{dp}}.depth.sum",
			sample=get_samples_from_pop("all"))
	output:
		"{prefix}"+dataset+"_all{dp}.depth",
		"{prefix}"+dataset+"_all{dp}.depth.sum"
	params:
		header=get_depth_header
	shell:
		"""
		cat {input.depth} > {output[0]}
		echo "sample	{params.header}.depth.mean	{params.header}.depth.stdev" > {output[1]}
		cat {input.summary} >> {output[1]}
		"""

def get_sample_qcs(wildcards):
	inputs = [results+"/genotyping/pop_lists/"+dataset+
					"_all.indiv.list",
				"results/mapping/qc/ind_depth/unfiltered/"+dataset+
					"_all{dp}.depth.sum",
				results+"/qc/ind_depth/filtered/"+dataset+
					"_all{dp}.depth.sum"]
	if config["analyses"]["endogenous_content"]:
		inputs.append(
			results+"/qc/endogenous_content/"+dataset+"_all.endo.tsv")
	return inputs

localrules: combine_sample_qc

rule combine_sample_qc:
	input:
		get_sample_qcs
	output:
		results+"/qc/"+dataset+"_all{dp}.sampleqc.tsv"
	shadow: "copy-minimal"
	shell:
		r"""
		for i in {input}; do
			head -n 1 $i > ${{i}}.tmp
			tail -n +2 $i | sort -k1 >> ${{i}}.tmp
		done

		cut -d '	' -f 1 {input[0]}.tmp > {output}
		
		for i in {input}; do
			cut -d '	' -f 2- ${{i}}.tmp | paste {output} - > {output}.tmp
			mv {output}.tmp {output}
		done
		"""

rule sample_qc_summary:
	input:
		results+"/qc/"+dataset+"_all{dp}.sampleqc.tsv"
	output:
		report(results+"/qc/"+dataset+"_all{dp}.sampleqc.html",
				category="Quality Control",
				labels={
					"Topic":"Sample QC",
					"Type":"Table"
				})
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"

rule ngsrelate:
	input:
		beagle=rules.merge_beagle.output.beagle,
		bamlist=rules.angsd_makeBamlist.output,
		inds=rules.popfile.output.inds
	output:
		relate=results+"/analyses/kinship/ngsrelate/"+dataset+
			"_{population}{dp}_relate.tsv",
		samples=results+"/analyses/kinship/ngsrelate/"+dataset+
			"_{population}{dp}_samples.list"
	log:
		logs + "/kinship/ngsrelate/"+dataset+"_{population}{dp}.log"
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

rule ngsrelate_summary:
	input:
		results+"/analyses/kinship/ngsrelate/"+dataset+
			"_{population}{dp}_relate.tsv"
	output:
		report(results+"/analyses/kinship/ngsrelate/"+dataset+
				"_{population}{dp}_relate.html",
				category="Kinship",
				labels={
					"Topic":"NgsRelate",
					"Type":"Table"
				})
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"