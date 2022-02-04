localrules: angsd_makeBamlist

rule angsd_makeBamlist:
	input:
		bams=get_bamlist_bams
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
	threads: 4
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
	threads: 4
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
	threads: 4
	resources:
		time="24:00:00"
	shell:
		"""
		realSFS {input.safidx} -fold {params.fold} -P {threads} \
			-bootstrap {params.boots} 2> {log} > {output.sfs}
		"""