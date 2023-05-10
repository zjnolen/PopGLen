# Pairwise individual relatedness with R0, R1, KING-robust kinship 
# method from Waples et al. 2019, MolEcol

rule est_kinship_stats:
	input:
		sfs="results/datasets/{dataset}/analyses/sfs/{dataset}.{ref}_{ind1}-{ind2}{dp}.sfs"
	output:
		"results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_{ind1}-{ind2}{dp}.kinship"
	log:
		"logs/{dataset}/kinship/waples2019/{dataset}.{ref}_{ind1}-{ind2}{dp}_kinship.log"
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

rule compile_kinship_stats:
	input:
		get_kinship
	output:
		"results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_all{dp}.kinship"
	log:
		"logs/{dataset}/kinship/waples2019/{dataset}.{ref}_all{dp}_compile-stats.log"
	conda:
		"../envs/shell.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*15
	shell:
		"""
		(echo "ind1	ind2	R0	R1	KING" > {output}
		cat {input} >> {output}) 2> {log}
		"""

rule kinship_table_html:
	input:
		"results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_all{dp}.kinship"
	output:
		report("results/datasets/{dataset}/analyses/kinship/waples2019/{dataset}.{ref}_all{dp}.kinship.html",
				category="Kinship",
				labels={
					"Topic":"Waples et al. 2019 R0,R1,KING",
					"Type":"Table"
				})
	log:
		"logs/{dataset}/kinship/waples2019/{dataset}.{ref}_all{dp}_tsv2html.log"
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"

# Individual depth histograms

rule ind_unfiltered_depth:
	input:
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		ref="results/ref/{ref}/{ref}.fa"
	output:
		sample_hist="results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}.depthSample",
		global_hist=temp("results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}.depthGlobal"),
		arg="results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}.arg"
	log:
		"logs/mapping/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}.log"
	container:
		angsd_container
	params:
		out=lambda w, output: os.path.splitext(output.arg)[0],
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
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		ref="results/ref/{ref}/{ref}.fa",
		sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites",
		idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites.idx"
	output:
		sample_hist="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}.depthSample",
		global_hist=temp("results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}.depthGlobal"),
		arg="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}.arg"
	log:
		"logs/{dataset}/ind_depth/filtered/{dataset}.{ref}_{population}{dp}.log"
	container:
		angsd_container
	params:
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=lambda w, output: os.path.splitext(output.arg)[0],
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

rule summarize_ind_depth:
	input:
		sample_hist="{prefix}{dataset}.{ref}_{sample}{dp}.depthSample",
		bed=get_total_bed
	output:
		sample_summ="{prefix}{dataset}.{ref}_{sample}{dp}.depth.sum"
	log:
		"logs/summarize_ind_depth/{prefix}{dataset}.{ref}_{sample}{dp}.log"
	conda:
		"../envs/r.yaml"
	wildcard_constraints:
		prefix="results/mapping/qc/ind_depth/unfiltered/|"+
			"results/datasets/{dataset}/qc/ind_depth/filtered/"
	threads: lambda wildcards, attempt: attempt
	script:
		"../scripts/calc_depth.R"

rule merge_ind_depth:
	input:
		depth=lambda w: expand("{{prefix}}{{dataset}}.{{ref}}_{sample}{{dp}}.depthSample",
			sample=get_samples_from_pop("all")),
		summary=lambda w: expand("{{prefix}}{{dataset}}.{{ref}}_{sample}{{dp}}.depth.sum",
			sample=get_samples_from_pop("all"))
	output:
		dep="{prefix}{dataset}.{ref}_all{dp}.depth",
		sum="{prefix}{dataset}.{ref}_all{dp}.depth.sum"
	log:
		"logs/merge_depth/{prefix}{dataset}.{ref}_all{dp}.log"
	conda:
		"../envs/shell.yaml"
	params:
		header=get_depth_header
	shell:
		"""
		(cat {input.depth} > {output.dep}
		echo "sample	{params.header}.depth.mean	{params.header}.depth.stdev" > {output.sum}
		cat {input.summary} >> {output.sum}) 2> {log}
		"""

localrules: combine_sample_qc

rule combine_sample_qc:
	input:
		unpack(get_sample_qcs)
	output:
		"results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv"
	log:
		"log/datasets/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}.log"
	conda:
		"../envs/shell.yaml"
	shadow: "minimal"
	shell:
		r"""
		(for i in {input}; do
			head -n 1 $i > ${{i}}.tmp
			tail -n +2 $i | sort -k1 >> ${{i}}.tmp
		done

		cut -d '	' -f 1 {input.inds}.tmp > {output}
		
		for i in {input}; do
			cut -d '	' -f 2- ${{i}}.tmp | paste {output} - > {output}.tmp
			mv {output}.tmp {output}
		done) 2> {log}
		"""

rule sample_qc_summary:
	input:
		"results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv"
	output:
		report("results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.html",
				category="Quality Control",
				labels={
					"Topic":"Sample QC",
					"Type":"Table"
				})
	log:
		"logs/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}_tsv2html.log"
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"

rule ngsrelate:
	input:
		beagle="results/datasets/{dataset}/beagles/pruned/{dataset}.{ref}_all{dp}_pruned.beagle.gz",
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_all{dp}.bamlist",
		inds="results/datasets/{dataset}/poplists/{dataset}_all.indiv.list"
	output:
		relate="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_relate.tsv",
		samples="results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_samples.list"
	log:
		"logs/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}.log"
	container:
		ngsrelate_container
	threads: lambda wildcards, attempt: attempt*4
	params:
		nind=lambda w: len(get_samples_from_pop("all"))
	resources:
		time=lambda wildcards, attempt: attempt*360
	shell:
		r"""
		(nsites=$(zcat {input.beagle} | tail -n +2 | wc -l)
		echo "nsites nind"
		echo $nsites {params.nind}
		cut -f1 {input.inds} | tail -n +2 > {output.samples}
		ngsRelate -G {input.beagle} -n {params.nind} -L $nsites -O {output.relate} \
			-z {output.samples}) &> {log}
		"""

rule ngsrelate_summary:
	input:
		"results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_relate.tsv"
	output:
		report("results/datasets/{dataset}/analyses/kinship/ngsrelate/{dataset}.{ref}_all{dp}_relate.html",
				category="Kinship",
				labels={
					"Topic":"NgsRelate",
					"Type":"Table"
				})
	log:
		"results/{dataset}/kinship/ngsrelate/{dataset}.{ref}_all{dp}_tsv2html.log"
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"