localrules: angsd_makeBamlist, popfile

rule angsd_makeBamlist:
	input:
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		"results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist"
	log:
		"logs/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.log"
	conda:
		"../envs/shell.yaml"
	shell:
		"""
        (readlink -f {input.bams} | perl -pe 'chomp if eof') > {output} 2> {log}
        """

rule popfile:
	output:
		inds="results/datasets/{dataset}/poplists/{dataset}_{population}.indiv.list"
	log:
		"logs/{dataset}/poplists/{dataset}_{population}_makelist.log"
	conda:
		"../envs/python.yaml"
	params:
		samplelist= samples,
		inds = lambda w: get_samples_from_pop(w.population)
	script:
		"../scripts/make_popfile.py"

rule angsd_sites_index:
	input:
		"results/datasets/{dataset}/filters/{prefix}.sites"
	output:
		multiext("results/datasets/{dataset}/filters/{prefix}.sites",".idx",".bin")
	log:
		"logs/{dataset}/angsd/sites_index/{prefix}.log"
	container:
		angsd_container
	shell:
		"""
		angsd sites index {input} 2> {log}
		"""

rule angsd_doGlf4:
	input:
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		ref="results/ref/{ref}/{ref}.fa",
		regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
		sites="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites",
		idx="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites.idx"
	output:
		glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.glf.gz",
		arg="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.arg"
	log:
		"logs/{dataset}/angsd/doGlf4/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
	wildcard_constraints:
		population="all"
	container:
		angsd_container
	params:
		gl_model=config["params"]["angsd"]["gl_model"],
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=lambda w, output: os.path.splitext(output.arg)[0]
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

rule sampleglf:
	input:
		glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_all{dp}_chunk{chunk}.glf.gz"
	output:
		glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{sample}{dp}_chunk{chunk}.glf.gz"
	log:
		"logs/{dataset}/angsd/sampleglf/{dataset}.{ref}_{sample}{dp}_chunk{chunk}.log"
	conda:
		"../envs/shell.yaml"
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
		sample_glfs=lambda w: expand("results/datasets/{{dataset}}/glfs/chunks/{{dataset}}.{{ref}}_{sample}{{dp}}_chunk{{chunk}}.glf.gz", 
			sample=get_samples_from_pop(w.population))
	output:
		glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.glf.gz"
	log:
		"logs/{dataset}/angsd/popglf/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
	wildcard_constraints:
		population="|".join(
			[i for i in samples.index.tolist()] +
			[i for i in samples.population.values.tolist()] +
			[i for i in samples.depth.values.tolist()]
			)
	conda:
		"../envs/shell.yaml"
	params:
		tmpfile="{dataset}.{ref}_{population}{dp}_chunk{chunk}.glf"
	resources:
		time=360
	shell:
		"""
		(zcat {input.sample_glfs[0]} | cut -f1-2 > {resources.tmpdir}/{params.tmpfile}
		for i in {input.sample_glfs}; do
			echo "Adding $i to glf..."
			zcat $i | cut -f3-12 | paste -d '	' {resources.tmpdir}/{params.tmpfile} - > {resources.tmpdir}/{params.tmpfile}.tmp
			mv {resources.tmpdir}/{params.tmpfile}.tmp {resources.tmpdir}/{params.tmpfile}
		done
		echo "Gzipping final glf..."
		gzip -c {resources.tmpdir}/{params.tmpfile} > {output.glf}) &> {log}
		"""

rule angsd_doGlf2:
	input:
		glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.glf.gz",
		fai="results/ref/{ref}/{ref}.fa.fai",
		sites=get_snpset
	output:
		beagle=temp("results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.beagle.gz"),
		maf=temp("results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.mafs.gz"),
		arg="results/datasets/{dataset}/beagles/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.arg"
	log:
		"logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		popopts=get_popopts,
		nind=lambda w: len(get_samples_from_pop(w.population)),
		out=lambda w, output: os.path.splitext(output.arg)[0]
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*720
	shell:
		"""
		angsd -doGlf 2 -glf10_text {input.glf} {params.popopts} -doMaf 1 \
			-nThreads {threads} -sites {input.sites[0]} -fai {input.fai} \
			-nInd {params.nind} -out {params.out} &> {log}
		"""

rule merge_beagle:
	input:
		lambda w: expand("results/datasets/{{dataset}}/beagles/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.beagle.gz", chunk=chunklist)
	output:
		beagle="results/datasets/{dataset}/beagles/{dataset}.{ref}_{population}{dp}.beagle.gz"
	log:
		"logs/{dataset}/angsd/doGlf2/{dataset}.{ref}_{population}{dp}_merge-beagle.log"
	conda:
		"../envs/shell.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		(set +o pipefail;
		zcat {input} | head -n 1 | gzip > {output}

		for f in {input}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.beagle}
		done) 2> {log}
		"""

rule merge_maf:
	input:
		lambda w: expand("results/datasets/{{dataset}}/beagles/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.mafs.gz", chunk=chunklist)
	output:
		maf="results/datasets/{dataset}/mafs/{dataset}.{ref}_{population}{dp}.mafs.gz"
	log:
		"logs/{dataset}/angsd/doGlf2/{dataset}{ref}_{population}{dp}_merge-mafs.log"
	conda:
		"../envs/shell.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		(set +o pipefail;
		zcat {input} | head -n 1 | gzip > {output.maf}

		for f in {input}; do
			zcat $f | tail -n +2 | gzip | cat >> {output.maf}
		done) 2> {log}
		"""

rule snpset:
	input:
		"results/datasets/{dataset}/mafs/{dataset}.{ref}_all{dp}.mafs.gz"
	output:
		"results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_snps.sites"
	log:
		"logs/{dataset}/filters/snps/{dataset}.{ref}{dp}_snps.log"
	conda:
		"../envs/shell.yaml"
	shell:
		"""
		zcat {input} | tail -n +2 | cut -f1-4 > {output} 2> {log}
		"""

rule angsd_doSaf:
	input:
		glf="results/datasets/{dataset}/glfs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.glf.gz",
		fai="results/ref/{ref}/{ref}.fa.fai",
		anc="results/ref/{ref}/{ref}.fa"
	output:
		saf=temp("results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.saf.gz"),
		safidx=temp("results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.saf.idx"),
		safpos=temp("results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.saf.pos.gz"),
		arg="results/datasets/{dataset}/safs/chunks/{dataset}.{ref}_{population}{dp}_chunk{chunk}.arg",
	log:
		"logs/{dataset}/angsd/doSaf/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		nind=lambda w: len(get_samples_from_pop(w.population)),
		out=lambda w, output: os.path.splitext(output.arg)[0]
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
		safs=lambda w: expand("results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.saf.idx",
			chunk=chunklist),
		safgz=lambda w: expand("results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.saf.gz",
			chunk=chunklist),
		safposgz=lambda w: expand("results/datasets/{{dataset}}/safs/chunks/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.saf.pos.gz",
			chunk=chunklist)
	output:
		multiext("results/datasets/{dataset}/safs/{dataset}.{ref}_{population}{dp}.saf",".idx",".pos.gz",".gz")
	log:
		"logs/{dataset}/realSFS/cat/{dataset}.{ref}_{population}{dp}.log"
	container:
		angsd_container
	params:
		out=lambda w, output: output[0].removesuffix('.saf.idx')
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		realSFS cat {input.safs} -P 1 -outnames {params.out} 2> {log}
		"""

rule angsd_doIBS:
	input:
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
		bams=get_bamlist_bams,
		bais=get_bamlist_bais,
		sites="results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_snps.sites",
		idx="results/datasets/{dataset}/filters/snps/{dataset}.{ref}{dp}_snps.sites.idx"
	output:
		ibs="results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_{population}{dp}.ibs.gz",
		ibsmat="results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_{population}{dp}.ibsMat",
		arg="results/datasets/{dataset}/analyses/IBS/{dataset}.{ref}_{population}{dp}.arg"
	log:
		"logs/{dataset}/angsd/doIBS/{dataset}{ref}_{population}{dp}.log"
	container:
		angsd_container
	params:
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=lambda w, output: os.path.splitext(output.arg)[0]
	threads: 8
	resources:
		time=lambda wildcards, attempt: attempt*2880
	shell:
		"""
		angsd -doIBS 1 -bam {input.bamlist} -nThreads {threads} -doCounts 1 \
			-minMapQ {params.mapQ} -minQ {params.baseQ} -sites {input.sites} \
			-doMajorMinor 3 -makeMatrix 1 -out {params.out} &> {log}
		"""