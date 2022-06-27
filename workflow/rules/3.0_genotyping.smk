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
		"docker://zjnolen/angsd:0.937"
	shell:
		"""
		angsd sites index {input} 2> {log}
		"""

rule angsd_doGlf2:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		anc=REF,
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=results+"/genotyping/filters/beds/"+dataset+"_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"_filts.sites.idx"
	output:
		beagle=results+"/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.beagle.gz",
		maf=results+"/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.mafs.gz"
	log:
		logs + "/angsd/doGlf2/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		"docker://zjnolen/angsd:0.937"
	params:
		gl_model=config["params"]["angsd"]["gl_model"],
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		miss=get_miss_data_prop,
		pval=config["params"]["angsd"]["snp_pval"],
		maf=config["params"]["angsd"]["min_maf"],
		out=results + "/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*180
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		minInd=$(echo $nInd \
			| awk '{{print $1*(1-{params.miss})}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')

		angsd -doGlf 2 -bam {input.bamlist} -GL {params.gl_model} \
			-doMajorMinor 1 -doMaf 1 -SNP_pval {params.pval} \
			-minMaf {params.maf} -ref {input.ref} -nThreads {threads}  \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-minInd $minInd -sites {input.sites} -rf {input.regions} \
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


rule angsd_doSaf:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		anc=REF,
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=results+"/genotyping/filters/beds/"+dataset+"_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"_filts.sites.idx"
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
		"docker://zjnolen/angsd:0.937"
	params:
		gl_model=config["params"]["angsd"]["gl_model"],
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		miss=get_miss_data_prop,
		out=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	resources:
		time=lambda wildcards, attempt: attempt*180
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		minInd=$(echo $nInd \
			| awk '{{print $1*(1-{params.miss})}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')

		angsd -doSaf 1 -bam {input.bamlist} -GL {params.gl_model} \
			-ref {input.ref} -anc {input.anc} -nThreads {threads} \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-minInd $minInd -sites {input.sites} -rf {input.regions} \
			-setMinDepthInd 5 -doCounts 1 -out {params.out} &> {log}
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
		"docker://zjnolen/angsd:0.937"
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
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=results+"/genotyping/filters/beds/"+dataset+"_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"_filts.sites.idx"
	output:
		results+"/genotyping/haploid_calls/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.haplo.gz"
	log:
		logs+"/angsd/doHaploCall/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		"docker://zjnolen/angsd:0.937"
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
		sites=results+"/genotyping/filters/beds/"+dataset+"_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"_filts.sites.idx"
	output:
		results+"/genotyping/single-read-sampling/"+dataset+
			"_{population}{dp}.ibs.gz",
		results+"/genotyping/single-read-sampling/"+dataset+
			"_{population}{dp}.ibsMat"
	log:
		logs+"/angsd/doIBS/"+dataset+"_{population}{dp}.log"
	container:
		"docker://zjnolen/angsd:0.937"
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

# rule angsd_doGlf4:
# 	input:
# 		ref=REF,
# 		bamlist=rules.angsd_makeBamlist.output,
# 		sites=results+"/genotyping/filters/beds/"+dataset+"_filts.sites",
# 		idx=results+"/genotyping/filters/beds/"+dataset+"_filts.sites.idx"
# 	output:
# 		glf=results+"/genotyping/glf4/chrom/"+dataset+
# 			"_{population}{dp}_chr{chr}.glf.gz",
# 		nind=results+"/genotyping/glf4/chrom/"+dataset+
# 			"_{population}{dp}_chr{chr}.glf.nInd"
# 	log:
# 		logs + "/genotyping/doGlf4/"+dataset+"_{population}{dp}_chr{chr}.log"
# 	container:
# 		"docker://zjnolen/angsd:0.937"
# 	params:
# 		gl_model=config["params"]["angsd"]["gl_model"],
# 		extra=config["params"]["angsd"]["extra"],
# 		mapQ=config["mapQ"],
# 		baseQ=config["baseQ"],
# 		miss=get_miss_data_prop,
# 		out=results + "/genotyping/glf4/chrom/"+dataset+
# 			"_{population}{dp}_chr{chr}"
# 	threads: lambda wildcards, attempt: attempt*2
# 	resources:
# 		time=lambda wildcards, attempt: attempt*300
# 	shell:
# 		r"""
# 		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
# 		minInd=$(echo $nInd \
# 			| awk '{{print $1*(1-{params.miss})}}' \
# 			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
# 		echo $nInd > {output.nind}

# 		angsd -doGlf 1 -bam {input.bamlist} -GL {params.gl_model} \
# 			-ref {input.ref} -nThreads {threads} {params.extra} \
# 			-minMapQ {params.mapQ} -minQ {params.baseQ} -minInd $minInd \
# 			-sites {input.sites} -r {wildcards.chr} -out {params.out} &> {log}
# 		"""

# rule angsd_doGlf2:
# 	input:
# 		glf=rules.angsd_doGlf4.output.glf,
# 		nind=rules.angsd_doGlf4.output.nind,
# 		fai=REF+".fai"
# 	output:
# 		beagle=results+"/genotyping/beagle/chrom/"+dataset+
# 			"_{population}{dp}_chr{chr}.beagle.gz"
# 	log:
# 		logs + "/angsd/doGlf2/"+dataset+"_{population}{dp}_chr{chr}.log"
# 	container:
# 		"docker://zjnolen/angsd:0.937"
# 	params:
# 		pval=config["params"]["angsd"]["snp_pval"],
# 		maf=config["params"]["angsd"]["min_maf"],
# 		out=results + "/genotyping/beagle/chrom/"+dataset+
# 			"_{population}{dp}_chr{chr}"
# 	threads: lambda wildcards, attempt: attempt
# 	resources:
# 		time="02:00:00"
# 	shell:
# 		"""
# 		nInd=$(cat {input.nind})

# 		angsd -glf10_text {input.glf} -doGlf 2 -nInd $nInd -fai {input.fai} \
# 			-nThreads {threads} -doMaf 1 -doMajorMinor 1 \
# 			-SNP_pval {params.pval} -minMaf {params.maf}  \
# 			-out {params.out} &> {log}
# 		"""

# rule merge_glf:
# 	input:
# 		glf=lambda w: expand(results+"/genotyping/glf4/chrom/"+dataset+
# 			"_{{population}}{{dp}}_chr{chr}.glf.gz", chr=filt_chroms()),
# 		nind=lambda w: expand(results+"/genotyping/glf4/chrom/"+dataset+
# 			"_{{population}}{{dp}}_chr{chr}.glf.nInd", chr=filt_chroms())
# 	output:
# 		glf=results+"/genotyping/glf4/"+dataset+"_{population}{dp}.glf.gz",
# 		nind=results+"/genotyping/glf4/"+dataset+"_{population}{dp}.glf.nInd"
# 	log:
# 		logs + "/angsd/doGlf4/"+dataset+
# 			"_{population}{dp}_mergeglf.log"
# 	shell:
# 		"""
# 		zcat {input.glf} | gzip > {output.glf} 2> {log}
# 		cat {input.nind[0]} > {output.nind} 2> {log}
# 		"""

# rule merge_beagle:
# 	input:
# 		lambda w: expand(results+"/genotyping/beagle/chrom/"+dataset+
# 			"_{{population}}{{dp}}_chr{chr}.beagle.gz", chr=filt_chroms())
# 	output:
# 		beagle=results+"/genotyping/beagle/"+dataset+
# 			"_{population}{dp}.beagle.gz"
# 	log:
# 		logs + "/angsd/doGlf2/"+dataset+
# 			"_{population}{dp}_mergebeag.log"
# 	shell:
# 		r"""
# 		set +o pipefail;
# 		zcat {input[0]} | head -n 1 | gzip > {output} 2> {log}

# 		for f in {input}; do
# 			zcat $f | tail -n +2 | gzip | cat >> {output.beagle} 2>> {log}
# 		done
# 		"""

# rule angsd_doIBS:
# 	input:
# 		bamlist=rules.angsd_makeBamlist.output,
# 		sites=rules.combine_beds.output.sit,
# 		sitesidx=rules.combine_beds.output.sit+".idx",
# 		ref=REF
# 	output:
# 		ibs=results+"/genotyping/single-read-sampling/"+dataset+
# 			"_{population}{dp}.ibs.gz",
# 		mat=results+"/genotyping/single-read-sampling/"+dataset+
# 			"_{population}{dp}.ibsMat"
# 	log:
# 	 	logs + "/angsd/doIBS/"+dataset+"_{population}{dp}.log"
# 	container:
# 		"docker://zjnolen/angsd:0.937"
# 	params:
# 		extra=config["params"]["angsd"]["extra"],
# 		mapQ=config["mapQ"],
# 		baseQ=config["baseQ"],
# 		miss=get_miss_data_prop,
# 		out=results+"/genotyping/single-read-sampling/"+dataset+
# 			"_{population}{dp}"
# 	resources:
# 		time=lambda wildcards, attempt: 2880*attempt
# 	shell:
# 		r"""
# 		angsd -doIBS 1 -makeMatrix 1 -bam {input.bamlist} -doCounts 1 \
# 			-ref {input.ref} -nThreads {threads} {params.extra} \
# 			-minMapQ {params.mapQ} -minQ {params.baseQ} \
# 			-sites {input.sites} -out {params.out} &> {log}
# 		"""