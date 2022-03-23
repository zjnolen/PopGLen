localrules: angsd_makeBamlist

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
		bamlist=rules.angsd_makeBamlist.output,
		ref=genome_file(),
		anc=genome_file(),
		fai=genome_file() + ".fai",
		popDP=results+"/depth/{population}.depthMean"
	output:
		saf=results + "/angsd/saf/{population}.saf.gz",
		safidx=results + "/angsd/saf/{population}.saf.idx",
		arg=results + "/angsd/saf/{population}.arg"
	log:
		logs + "/angsd/doSaf/{population}.log"
	conda:
		"../envs/angsd.yaml"
	params:
		extra=config["params"]["angsd"]["extra"],
		gl_model=config["params"]["angsd"]["gl_model"],
		miss=get_miss_data_prop,
		out_prefix=results + "/angsd/saf/{population}"
	resources:
		time="12:00:00"
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		minInd=$(echo $nInd \
			| awk '{{print $1*(1-{params.miss})}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
		minDP=$(echo $nInd | awk '{{print $1*2}}')
		maxDP=$(cat {input.popDP} | awk '{{print 3*$1}}')

		angsd -doSaf 1 -bam {input.bamlist} -GL {params.gl_model} \
			-anc {input.anc} -ref {input.ref} -nThreads {threads} \
			{params.extra} -out {params.out_prefix} -doCounts 1 \
			-setMinDepth $minDP -setMaxDepth $maxDP -minInd $minInd \
			-setMinDepthInd 3 &> {log}
		"""

rule angsd_realSFS:
	input:
		safidx=rules.angsd_doSaf.output.safidx
	output:
		sfs=(results + "/angsd/sfs/{population}.sfs")
	conda:
		"../envs/angsd.yaml"
	log:
		logs + "/angsd/sfs/{population}.log"
	params:
		fold=config["params"]["angsd"]["fold"]
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		realSFS {input.safidx} -fold {params.fold} -P {threads} 2> {log} \
			> {output.sfs}
		"""

rule angsd_realSFS_boot:
	input:
		safidx=rules.angsd_doSaf.output.safidx
	output:
		sfs=results + "/angsd/sfs/{population}.bootsfs"
	conda:
		"../envs/angsd.yaml"
	log:
		logs + "/angsd/sfs/{population}_boot.log"
	params:
		fold=config["params"]["angsd"]["fold"],
		boots=config["params"]["angsd"]["sfsboot"]
	resources:
		time="24:00:00"
	threads: lambda wildcards, attempt: attempt
	shell:
		"""
		realSFS {input.safidx} -fold {params.fold} -P {threads} \
			-bootstrap {params.boots} 2> {log} > {output.sfs}
		"""

rule angsd_chrom_beagle:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		ref=genome_file(),
		fai=rules.samtools_faidx.output,
		popDP=results+"/depth/{population}.depthMean"
	output:
		beagle=results + "/angsd/beagle/chrom/{population}_chr{chrom}.beagle.gz",
		arg=results + "/angsd/beagle/chrom/{population}_chr{chrom}.arg"
	log:
		logs + "/angsd/beagle/chrom/{population}_chr{chrom}.log"
	conda:
		"../envs/angsd.yaml"
	params:
		extra=config["params"]["angsd"]["extra"],
		gl_model=config["params"]["angsd"]["gl_model"],
		miss=get_miss_data_prop,
		pval=config["params"]["angsd"]["snp_pval"],
		out_prefix=results + "/angsd/beagle/chrom/{population}_chr{chrom}"
	resources:
		time="48:00:00"
	shell:
		r"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		minInd=$(echo $nInd \
			| awk '{{print $1*(1-{params.miss})}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
		minDP=$(echo $nInd | awk '{{print $1*2}}')
		maxDP=$(cat {input.popDP} | awk '{{print 3*$1}}')

		angsd -GL {params.gl_model} -doGlf 2 -doMaf 1 -SNP_pval {params.pval} \
			-bam {input.bamlist} -nThreads {threads} -out {params.out_prefix} \
			-doMajorMinor 1 -r {wildcards.chrom} -ref {input.ref} \
			-minInd $minInd {params.extra} -doCounts 1 -minMaf 0.05 \
			-setMinDepth $minDP -setMinDepthInd 3 -setMaxDepth $maxDP &> {log}
		"""

rule angsd_merge_beagle:
	input:
		lambda w: expand(results+"/angsd/beagle/chrom/{{population}}_chr{chrom}.beagle.gz", chrom=get_contigs())
	output:
		beagle=results + "/angsd/beagle/{population}_genome.beagle.gz"
	log:
		logs + "/angsd/beagle/aggregate_{population}.log"
	shell:
		r"""
		set +o pipefail;
		zcat {input[0]} | head -n 1 | gzip > {output} 2> {log}

		for f in {input}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.beagle} 2>> {log}
		done
		"""