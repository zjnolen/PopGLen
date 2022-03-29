localrules: merge_plink, tped2ped

rule angsd_doPlink:
	input:
		bamlist=rules.angsd_makeBamlist.output,
		ref=genome_file(),
		fai=rules.samtools_faidx.output,
		popDP=results+"/depth/{population}.depthMean"
	output:
		tfam=temp(results + "/angsd/plink/chrom/{population}_chr{chrom}.tfam"),
		tped=temp(results + "/angsd/plink/chrom/{population}_chr{chrom}.tped")
	log:
		logs + "/angsd/plink/chrom/{population}_chr{chrom}.log"
	conda:
		"../envs/angsd.yaml"
	params:
		extra=config["params"]["angsd"]["extra"],
		gl_model=config["params"]["angsd"]["gl_model"],
		miss=get_miss_data_prop,
		pval=config["params"]["angsd"]["snp_pval"],
		out_prefix=results + "/angsd/plink/chrom/{population}_chr{chrom}"
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

		angsd -GL {params.gl_model} -doMaf 1 -SNP_pval {params.pval} \
			-doPlink 2 -doPost 1 -doGeno -4 -postCutoff 0.99 \
			-bam {input.bamlist} -nThreads {threads} -out {params.out_prefix} \
			-doMajorMinor 1 -r {wildcards.chrom} -ref {input.ref} \
			-minInd $minInd {params.extra} -doCounts 1 -minMaf 0.05 \
			-setMinDepth $minDP -setMinDepthInd 3 -setMaxDepth $maxDP &> {log}
		"""

rule merge_plink:
	input:
		tped=lambda w: expand(results+"/angsd/plink/chrom/{{population}}_chr{chrom}.tped", chrom=get_autos()),
		tfam=lambda w: expand(results+"/angsd/plink/chrom/{{population}}_chr{chrom}.tfam", chrom=get_autos())
	output:
		tped=temp(results + "/angsd/plink/{population}_genome.tped"),
		tfam=temp(results + "/angsd/plink/{population}_genome.tfam")
	run:
		from snakemake.shell import shell
		open(output.tped, 'w').close()
		with open(output.tped, 'a') as out:
			for i in range(len(input.tped)):
				if i < 200:
					chr = i + 1
					with open(input.tped[i]) as f:
						print(re.sub('^(.*?) ', str(chr) + ' ', f.read(), 
							flags=re.MULTILINE), file=out)
				else:
					print("Skipping chromosomes after first 200...")
		shell(
			"sed -i '/^$/d' {output.tped}"
		)
		shell(
			"cp {input.tfam[0]} {output.tfam}"
		)

rule tped2ped:
	input:
		rules.merge_plink.output.tped,
		rules.merge_plink.output.tfam
	output:
		ped=results + "/angsd/plink/{population}_genome.ped",
		map=results + "/angsd/plink/{population}_genome.map",
		nosex=temp(results + "/angsd/plink/{population}_genome.nosex")
	log:
		logs + "/plink/{population}_tped2ped.log"
	conda:
		"../envs/plink.yaml"
	params:
		inprefix=results + "/angsd/plink/{population}_genome",
		outprefix=results + "/angsd/plink/{population}_genome"
	shell:
		"""
		plink --tfile {params.inprefix} --recode \
			--out {params.outprefix} &> {log}
		"""