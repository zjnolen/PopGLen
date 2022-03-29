localrules: angsd_makeBamlist, popfile, angsd_make_auto_rf, merge_pop_site_lists

# input and filtering file creation

rule popfile:
	output:
		inds=results + "/pop_lists/{population}.indiv.list"
	run:
		inds = get_samples_from_pop(wildcards.population)
		with open(output.inds, "w") as out:
			for item in inds:
				out.write("%s\n" % item)

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

rule angsd_make_auto_rf:
	input:
		genome_file() + ".fai"
	output:
		rf=results+"/angsd/regions/autosomes.regions"
	run:
		autos = get_autos()
		with open(output.rf, "w") as out:
			for item in autos:
				out.write("%s\n" % item)

rule get_pop_site_lists:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		rf=rules.angsd_make_auto_rf.output.rf,
		popDP=results+"/depth/{population}.depthMean",
		ref=genome_file()
	output:
		posgz=results+"/sites/{population}_chr{chrom}_filtautos.pos.gz",
		sites=results+"/sites/{population}_chr{chrom}_filtautos.sites",
		bed=results+"/sites/{population}_chr{chrom}_filtautos.bed"
	params:
		out_prefix=results+"/sites/{population}_chr{chrom}_filtautos",
		extra=config["params"]["angsd"]["extra"],
		miss=get_miss_data_prop
	log:
		logs + "/angsd/sites/{population}_chr{chrom}_filtautos.log"
	conda:
		"../envs/angsd.yaml"
	resources:
		time="12:00:00"
	shell:
		r"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		minInd=$(echo $nInd \
			| awk '{{print $1*(1-{params.miss})}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
		minDP=$(awk '{{print $2}}' {input.popDP})
		maxDP=$(awk '{{print $3}}' {input.popDP})

		angsd -bam {input.bamlist} -ref {input.ref} -nThreads {threads} \
			{params.extra} -out {params.out_prefix} -doCounts 1 \
			-setMinDepth $minDP -setMaxDepth $maxDP -minInd $minInd \
			-setMinDepthInd 3 -r {wildcards.chrom} -dumpCounts 1 &> {log}
        
        zcat {output.posgz} | tail -n +2 | awk '{{print $1"\t"$2}}' \
            > {output.sites} 2>> {log}
		
		awk '{{print $1"\t"$2"\t"$2+1}}' {output.sites} > {output.bed} 2>> {log}
        """

rule get_site_intersects:
	input:
		bed=get_intersect_inputs
	output:
		sites=results+"/sites/{type}_chr{chrom}_autos_intersect.sites",
		bed=results+"/sites/{type}_chr{chrom}_autos_intersect.bed"
	log:
		logs + "/bedops/{type}_chr{chrom}_autos_intersect.log"
	conda:
		"../envs/bedtools.yaml"
	resources:
		time="06:00:00"
	shell:
		r"""
		cat {input.bed[0]} > {output.bed} 2> {log}

		for b in {input.bed}; do
			a={output.bed}
			wc -l $a >> {log}
			echo "intersecting $b" >> {log}
			bedtools intersect -a $a -b $b -wa -sorted > {output.bed}.tmp 2>> {log}
			mv {output.bed}.tmp {output.bed} 2>> {log}
		done
		
		awk '{{print $1"\t"$2}}' {output.bed} > {output.sites} 2>> {log}
		"""

rule merge_pop_site_lists:
	input:
		sites=lambda w: expand(results+"/sites/{{type}}_chr{chrom}_autos_intersect.sites", chrom=get_autos()),
		bed=lambda w: expand(results+"/sites/{{type}}_chr{chrom}_autos_intersect.bed", chrom=get_autos())
	output:
		sites=results+"/sites/{type}_autos_intersect.sites",
		bed=results+"/sites/{type}_autos_intersect.bed"
	wildcard_constraints:
		type="population|sample"
	log:
		logs + "/angsd/sites/{type}_merge_sites.log"
	shell:
		"""
		cat {input.sites} > {output.sites} 2> {log}
		cat {input.bed} > {output.bed} 2>> {log}
		"""

rule angsd_sites_index:
	input:
		"{prefix}.sites"
	output:
		"{prefix}.sites.idx"
	conda:
		"../envs/angsd.yaml"
	shell:
		"""
		angsd sites index {input}
		"""

rule angsd_doSaf:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		ref=genome_file(),
		anc=genome_file(),
		fai=genome_file() + ".fai",
		popDP=results+"/depth/{population}.depthMean",
		rf=rules.angsd_make_auto_rf.output.rf
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
		minDP=$(awk '{{print $2}}' {input.popDP})
		maxDP=$(awk '{{print $3}}' {input.popDP})

		angsd -doSaf 1 -bam {input.bamlist} -GL {params.gl_model} \
			-anc {input.anc} -ref {input.ref} -nThreads {threads} \
			{params.extra} -out {params.out_prefix} -doCounts 1 \
			-setMinDepth $minDP -setMaxDepth $maxDP -minInd $minInd \
			-setMinDepthInd 3 -rf {input.rf} &> {log}
		"""

rule angsd_realSFS:
	input:
		safidx=rules.angsd_doSaf.output.safidx,
		sites=get_sites_file
	output:
		sfs=results + "/angsd/sfs/{population}{sites}.sfs"
	container:
		"library://james-s-santangelo/angsd/angsd:0.933"
	log:
		logs + "/angsd/sfs/{population}{sites}.log"
	params:
		fold=config["params"]["angsd"]["fold"],
		sites=get_sites_option
	threads: lambda wildcards, attempt: attempt*2
	shell:
		"""
		realSFS {input.safidx} -fold {params.fold} -P {threads} \
			{params.sites} > {output.sfs} 2> {log} 
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
		minDP=$(awk '{{print $2}}' {input.popDP})
		maxDP=$(awk '{{print $3}}' {input.popDP})

		angsd -GL {params.gl_model} -doGlf 2 -doMaf 1 -SNP_pval {params.pval} \
			-bam {input.bamlist} -nThreads {threads} -out {params.out_prefix} \
			-doMajorMinor 1 -r {wildcards.chrom} -ref {input.ref} \
			-minInd $minInd {params.extra} -doCounts 1 -minMaf 0.05 \
			-setMinDepth $minDP -setMinDepthInd 3 -setMaxDepth $maxDP &> {log}
		"""

rule angsd_merge_beagle:
	input:
		lambda w: expand(results+"/angsd/beagle/chrom/{{population}}_chr{chrom}.beagle.gz", chrom=get_autos())
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