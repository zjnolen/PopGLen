localrules: angsd_makeBamlist
ruleorder: angsd_merge_beagle > angsd_chrom_beagle

rule angsd_makeBamlist:
	input:
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		results + "/angsd/bamlists/{population}.bamlist"
	shell:
		"""
        (readlink -f {input.bams} | perl -pe 'chomp if eof') > {output}
        """

rule angsd_doSaf:
	input:
		bamlist=results + "/angsd/bamlists/{population}.bamlist",
		bai=get_bamlist_bais,
		ref=genome_file(),
		anc=genome_file(),
		fai=genome_file() + ".fai"
	output:
		saf=results + "/angsd/saf/{population}.saf.gz",
		safidx=results + "/angsd/saf/{population}.saf.idx",
		arg=results + "/angsd/saf/{population}.arg"
	log:
		"logs/angsd/{population}_saf.log"
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	params:
		extra=config["params"]["angsd"]["extra"],
		gl_model=config["params"]["angsd"]["gl_model"],
		out_prefix=results + "/angsd/saf/{population}"
	resources:
		time="12:00:00"
	threads: 2
	shell:
		"""
		angsd -doSaf 1 -bam {input.bamlist} -GL {params.gl_model} \
			-anc {input.anc} -ref {input.ref} -nThreads {threads} \
			{params.extra} -out {params.out_prefix} 2> {log}
		"""

rule angsd_realSFS:
	input:
		safidx=rules.angsd_doSaf.output.safidx
	output:
		sfs=(results + "/angsd/sfs/{population}_fold" 
			+ str(config["params"]["angsd"]["fold"]) + ".sfs")
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	log:
		"logs/angsd/{population}_sfs.log"
	params:
		fold=config["params"]["angsd"]["fold"]
	shell:
		"""
		realSFS {input.safidx} -fold {params.fold} -P {threads} 2> {log} \
			> {output.sfs}
		"""

rule angsd_realSFS_boot:
	input:
		safidx=rules.angsd_doSaf.output.safidx
	output:
		sfs=results + "/angsd/sfs/{population}_fold" \
			+ str(config["params"]["angsd"]["fold"]) + ".bootsfs"
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	log:
		"logs/angsd/{population}_bootsfs.log"
	params:
		fold=config["params"]["angsd"]["fold"],
		boots=config["params"]["angsd"]["sfsboot"]
	resources:
		time="24:00:00"
	shell:
		"""
		realSFS {input.safidx} -fold {params.fold} -P {threads} \
			-bootstrap {params.boots} 2> {log} > {output.sfs}
		"""

rule angsd_chrom_beagle:
	input:
		bamlist=results + "/angsd/bamlists/{population}.bamlist",
		ref=genome_file(),
		fai=genome_file() + ".fai"
	output:
		beagle=results + "/angsd/beagle/{population}_{chrom}_md{miss}.beagle.gz",
		arg=results + "/angsd/beagle/{population}_{chrom}_md{miss}.arg"
	log:
		"logs/angsd/beagle/{population}_{chrom}_md{miss}.log"
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	params:
		extra=config["params"]["angsd"]["extra"],
		gl_model=config["params"]["angsd"]["gl_model"],
		out_prefix=results + "/angsd/beagle/{population}_{chrom}_md{miss}",
		pval=config["params"]["angsd"]["snp_pval"]
	resources:
		time="48:00:00"
	shell:
		r"""
		minInd=$(cat {input.bamlist} \
			| wc -l \
			| awk '{{print ($1+1)*((100-{wildcards.miss})/100)}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
		
		minDP=$(cat {input.bamlist} \
			| wc -l \
			| awk '{{print ($1+1)*2}}')
		
		angsd -GL {params.gl_model} -doGlf 2 -doMaf 1 -SNP_pval {params.pval} \
			-bam {input.bamlist} -nThreads {threads} -out {params.out_prefix} \
			-doMajorMinor 1 -r {wildcards.chrom} -ref {input.ref} \
			-minInd $minInd {params.extra} -P {threads} -doCounts 1 \
			-minMaf 0.05 -setMinDepth 130 -setMaxDepth 2080 &> {log}
		"""

rule angsd_merge_beagle:
	input:
		aggregate_beagles
	output:
		beagle=results + "/angsd/beagle/{population}_genome_md{miss}.beagle.gz"
	log:
		"logs/angsd/beagle/aggregate_{population}_md{miss}.log"
	shell:
		r"""
		echo {input} | tr ' ' '\n' > {log}
		cat {input} > {output.beagle} 2>> {log}
		"""