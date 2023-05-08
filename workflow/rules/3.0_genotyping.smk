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
		glf=results+"/genotyping/glf/chunk/"+dataset+ \
			"_{population}{dp}_chunk{chunk}.glf.gz",
		fai=REF+".fai",
		sites=get_snpset
	output:
		beagle=temp(results+"/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.beagle.gz"),
		maf=temp(results+"/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.mafs.gz")
	log:
		logs + "/angsd/doGlf2/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		popopts=get_popopts,
		nind=lambda w: len(get_samples_from_pop(w.population)),
		out=results + "/genotyping/beagle/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*720
	shell:
		"""
		angsd -doGlf 2 -glf10_text {input.glf} {params.popopts} -doMaf 1 \
			-nThreads {threads} -sites {input.sites[0]} -fai {input.fai} \
			-nInd {params.nind} -out {params.out} &>> {log}
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
		glf=results+"/genotyping/glf/chunk/"+dataset+ \
			"_{population}{dp}_chunk{chunk}.glf.gz",
		fai=REF+".fai",
		anc=REF
	output:
		saf=temp(results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.saf.gz"),
		safidx=temp(results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.saf.idx"),
		safpos=temp(results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.saf.pos.gz"),
		arg=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}.arg",
	log:
		logs + "/angsd/doSaf/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		nind=lambda w: len(get_samples_from_pop(w.population)),
		out=results+"/genotyping/saf/chunk/"+dataset+
			"_{population}{dp}_chunk{chunk}"
	resources:
		time=lambda wildcards, attempt: attempt*180
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		angsd -doSaf 1 -glf10_text {input.glf} -anc {input.anc} \
			-nThreads 1 -fai {input.fai} -nInd {params.nind} \
			-out {params.out} &> {log}
		"""

rule realSFS_catsaf:
	input:
		safs=lambda w: expand(results+"/genotyping/saf/chunk/"+
			dataset+"_{{population}}{{dp}}_chunk{chunk}.saf.idx",
			chunk=chunklist),
		safgz=lambda w: expand(results+"/genotyping/saf/chunk/"+
			dataset+"_{{population}}{{dp}}_chunk{chunk}.saf.gz",
			chunk=chunklist),
		safposgz=lambda w: expand(results+"/genotyping/saf/chunk/"+
			dataset+"_{{population}}{{dp}}_chunk{chunk}.saf.pos.gz",
			chunk=chunklist)
	output:
		multiext(results+"/genotyping/saf/"+dataset+
			"_{population}{dp}.saf",".idx",".pos.gz",".gz")
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
		realSFS cat {input.safs} -P 1 -outnames {params.out} 2> {log}
		"""

rule angsd_doIBS:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		sites=results+"/genotyping/filters/beds/"+dataset+"{dp}_snps.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"{dp}_snps.sites.idx"
	output:
		results+"/genotyping/IBSmatrix/"+dataset+
			"_{population}{dp}.ibs.gz",
		results+"/genotyping/IBSmatrix/"+dataset+
			"_{population}{dp}.ibsMat",
		results+"/genotyping/IBSmatrix/"+dataset+
			"_{population}{dp}.arg"
	log:
		logs+"/angsd/doIBS/"+dataset+"_{population}{dp}.log"
	container:
		angsd_container
	params:
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/genotyping/IBSmatrix/"+dataset+
			"_{population}{dp}"
	threads: 8
	resources:
		time=lambda wildcards, attempt: attempt*2880
	shell:
		"""
		angsd -doIBS 1 -bam {input.bamlist} -nThreads {threads} -doCounts 1 \
			-minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
			-doMajorMinor 3 -makeMatrix 1 -out {params.out} &> {log}
		"""

rule angsd_doGlf4:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		ref=REF,
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		sites=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
		idx=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites.idx"
	output:
		glf=results+"/genotyping/glf/chunk/"+dataset+ \
			"_{population}{dp}_chunk{chunk}.glf.gz",
		arg=results+"/genotyping/glf/chunk/"+dataset+ \
			"_{population}{dp}_chunk{chunk}.arg"
	log:
		logs+"/angsd/doGlf4/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	wildcard_constraints:
		population="all"
	container:
		angsd_container
	params:
		gl_model=config["params"]["angsd"]["gl_model"],
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/genotyping/glf/chunk/"+dataset+ \
			"_{population}{dp}_chunk{chunk}"
	resources:
		time=lambda wildcards, attempt: attempt*360
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		angsd -doGlf 4 -bam {input.bamlist} -GL {params.gl_model} \
			-ref {input.ref} -nThreads {threads} {params.extra} \
			-minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
			-rf {input.regions} -out {params.out} &> {log}
		"""

# rule catglf:
# 	input:
# 		glfs=lambda w: expand(results+"/genotyping/glf/chunk/"+
# 			dataset+"_all{{dp}}_chunk{chunk}.glf.gz",
# 			chunk=chunklist)
# 	output:
# 		results+"/genotyping/glf/"+dataset+
# 			"_all{dp}.glf.gz"
# 	log:
# 		logs+"/angsd/doGlf4/cat/"+dataset+"_all{dp}.log"
# 	resources:
# 		time=lambda wildcards, attempt: attempt*60
# 	shell:
# 		"""
# 		cat {input.glfs} > {output} 2> {log}
# 		"""

rule sampleglf:
	input:
		glf=results+"/genotyping/glf/chunk/"+dataset+
			"_all{dp}_chunk{chunk}.glf.gz"
	output:
		glf=results+"/genotyping/glf/chunk/"+dataset+
			"_{sample}{dp}_chunk{chunk}.glf.gz"
	log:
		logs+"/angsd/sampleglf/"+dataset+"_{sample}{dp}_chunk{chunk}.log"
	params:
		start=lambda w: samples.index.values.tolist().index(w.sample),
		end=lambda w: str(samples.index.values.tolist().index(w.sample)+1)
	shell:
		"""
		zcat {input.glf} | cut -f1-2,{params.start}3-{params.end}2 | \
			gzip > {output} 2> {log}
		"""

rule popglf:
	input:
		sample_glfs=lambda w: expand(results+"/genotyping/glf/chunk/"+dataset+
			"_{sample}{{dp}}_chunk{{chunk}}.glf.gz", 
			sample=get_samples_from_pop(w.population))
	output:
		glf=results+"/genotyping/glf/chunk/"+dataset+ \
			"_{population}{dp}_chunk{chunk}.glf.gz"
	log:
		logs+"/angsd/sampleglf/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	wildcard_constraints:
		population="|".join(
			[i for i in samples.index.tolist()] +
			[i for i in samples.population.values.tolist()] +
			[i for i in samples.depth.values.tolist()]
			)
	params:
		tmpfile=dataset+"_{population}{dp}_chunk{chunk}.glf"
	resources:
		time=360
	shell:
		"""
		zcat {input.sample_glfs[0]} | cut -f1-2 > {resources.tmpdir}/{params.tmpfile} 2> {log}
		for i in {input.sample_glfs}; do
			echo "Adding $i to glf..." >> {log}
			zcat $i | cut -f3-12 | paste -d '	' {resources.tmpdir}/{params.tmpfile} - > {resources.tmpdir}/{params.tmpfile}.tmp 2>> {log}
			mv {resources.tmpdir}/{params.tmpfile}.tmp {resources.tmpdir}/{params.tmpfile} 2>> {log}
		done
		echo "Gzipping final glf..." >> {log}
		gzip -c {resources.tmpdir}/{params.tmpfile} > {output.glf} 2>> {log}
		"""