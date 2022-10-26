localrules: angsd_makeBamlist, popfile

rule angsd_makeBamlist:
	input:
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		results + "/genotyping/bamlists/"+dataset+"_{population}{dp}.bamlist"
	shell:
		"""
        (readlink -f {input.bams} | perl -pe 'chomp if eof') > {output}
        """

rule popfile:
	input:
		bamlist=results+"/genotyping/bamlists/"+dataset+
			"_{population}.bamlist"
	output:
		inds=results+"/genotyping/pop_lists/"+dataset+
			"_{population}.indiv.list"
	run:
		inds = get_samples_from_pop(wildcards.population)
		samples.loc[inds].to_csv(output.inds, sep="\t", quoting=csv.QUOTE_NONE,
			header = True, index = False)

rule angsd_sites_index:
	input:
		"{prefix}.sites"
	output:
		multiext("{prefix}.sites",".idx",".bin")
	log:
		logs + "/angsd_sites_index/{prefix}.log"
	container:
		angsd_container
	shell:
		"""
		angsd sites index {input} 2> {log}
		"""

def get_snpset(wildcards):
	pop = wildcards.population
	if pop == "all":
		return [results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
				results+"/genotyping/filters/beds/"+dataset+
					"{dp}_filts.sites.idx"]
	else:
		return [results+"/genotyping/filters/beds/"+dataset+"{dp}_snps.sites",
				results+"/genotyping/filters/beds/"+dataset+
					"{dp}_snps.sites.idx"]

def get_popopts(wildcards):
	pop = wildcards.population
	if pop == "all":
		return "-doMajorMinor 1 -SNP_pval "+ \
				str(config["params"]["angsd"]["snp_pval"])+" -minMaf "+ \
				str(config["params"]["angsd"]["min_maf"])
	else:
		return "-doMajorMinor 3"

rule angsd_doGlf2:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		anc=REF,
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=get_snpset
	output:
		beagle=results+"/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.beagle.gz",
		maf=results+"/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.mafs.gz"
	log:
		logs + "/angsd/doGlf2/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		gl_model=config["params"]["angsd"]["gl_model"],
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		# miss=get_miss_data_prop,
		popopts=get_popopts,
		out=results + "/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*720
	shell:
		"""
		# nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')

		angsd -doGlf 2 -bam {input.bamlist} -GL {params.gl_model} \
			{params.popopts} -doMaf 1 -ref {input.ref} -nThreads {threads}  \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-sites {input.sites[0]} -rf {input.regions} \
			-out {params.out} &> {log}
		"""

rule merge_beagle:
	input:
		lambda w: expand(results+"/genotyping/beagle/chunk/"+dataset+
			"_{{population}}{{dp}}_chunk{chunk}.beagle.gz", chunk=chunklist)
	output:
		beagle=results+"/genotyping/beagle/"+dataset+
			"_{population}{dp}.beagle.gz"
	log:
		logs + "/angsd/doGlf2/"+dataset+
			"_{population}{dp}_mergebeag.log"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		set +o pipefail;
		zcat {input[0]} | head -n 1 | gzip > {output} 2> {log}

		for f in {input}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.beagle} 2>> {log}
		done
		"""

rule merge_maf:
	input:
		lambda w: expand(results+"/genotyping/beagle/chunk/"+dataset+
			"_{{population}}{{dp}}_chunk{chunk}.mafs.gz", chunk=chunklist)
	output:
		maf=results+"/genotyping/mafs/"+dataset+
			"_{population}{dp}.mafs.gz"
	log:
		logs + "/angsd/doGlf2/"+dataset+
			"_{population}{dp}_mergemafs.log"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		set +o pipefail;
		zcat {input[0]} | head -n 1 | gzip > {output.maf} 2> {log}

		for f in {input}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.maf} 2>> {log}
		done
		"""

rule snpset:
	input:
		results+"/genotyping/mafs/"+dataset+"_all{dp}.mafs.gz"
	output:
		results+"/genotyping/filters/beds/"+dataset+"{dp}_snps.sites"
	shell:
		"""
		zcat {input} | tail -n +2 | cut -f1-4 > {output}
		"""

rule angsd_doSaf:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		anc=REF,
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites.idx"
	output:
		saf=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.saf.gz",
		safidx=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.saf.idx",
		arg=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.arg",
	log:
		logs + "/angsd/doSaf/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		gl_model=config["params"]["angsd"]["gl_model"],
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	resources:
		time=lambda wildcards, attempt: attempt*180
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')

		angsd -doSaf 1 -bam {input.bamlist} -GL {params.gl_model} \
			-ref {input.ref} -anc {input.anc} -nThreads {threads} \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-sites {input.sites} -rf {input.regions} -out {params.out} &> {log}
		"""

rule realSFS_catsaf:
	input:
		safs=lambda w: expand(results+"/genotyping/saf/chunk/"+
			dataset+"_{{population}}{{dp}}_chunk{chunk}.saf.idx",
			chunk=chunklist)
	output:
		results+"/genotyping/saf/"+dataset+
			"_{population}{dp}.saf.idx"
	log:
		logs+"/realSFS/cat/"+dataset+"_{population}{dp}.log"
	container:
		angsd_container
	params:
		out=results+"/genotyping/saf/"+dataset+"_{population}{dp}"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		realSFS cat {input.safs} -P {threads} -outnames {params.out} 2> {log}
		"""

rule angsd_haplocall:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites.idx"
	output:
		results+"/genotyping/haploid_calls/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.haplo.gz"
	log:
		logs+"/angsd/doHaploCall/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/genotyping/haploid_calls/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	resources:
		time=lambda wildcards, attempt: attempt*240
	shell:
		"""
		angsd -doHaploCall 1 -bam {input.bamlist} -ref {input.ref} \
			-nThreads {threads} -doCounts 1 -minMinor 1 {params.extra} \
			-minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
			-rf {input.regions} -out {params.out} &> {log}
		"""

rule merge_haplocall:
	input:
		haplos=lambda w: expand(results+"/genotyping/haploid_calls/chunk/"+
			dataset+"_{{population}}{{dp}}_chunk{chunk}.haplo.gz",
			chunk=chunklist)
	output:
		results+"/genotyping/haploid_calls/"+dataset+
			"_{population}{dp}.haplo.gz"
	log:
		logs+"/angsd/doHaploCall/"+dataset+"_{population}{dp}_merge.log"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		set +o pipefail;
		zcat {input.haplos[0]} | head -n 1 | gzip > {output} 2> {log}

		for i in {input.haplos}; do
			zcat $i | tail -n +2 | gzip >> {output}
		done 2>> {log}
		"""

rule angsd_doIBS:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		sites=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites.idx"
	output:
		results+"/genotyping/single-read-sampling/"+dataset+
			"_{population}{dp}.ibs.gz",
		results+"/genotyping/single-read-sampling/"+dataset+
			"_{population}{dp}.ibsMat"
	log:
		logs+"/angsd/doIBS/"+dataset+"_{population}{dp}.log"
	container:
		angsd_container
	params:
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/genotyping/single-read-sampling/"+dataset+
			"_{population}{dp}"
	resources:
		time=lambda wildcards, attempt: attempt*1440
	shell:
		"""
		angsd -doIBS 1 -bam {input.bamlist} -nThreads {threads} -doCounts 1 \
			-minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
			-makeMatrix 1 -out {params.out} &> {log}
		"""