# Creates a filtered list of autosomal sites for analyses to filter the 
# analyses to. This filtering limits analyses to autosomal scaffolds and 
# removes sites with low mappability, low complexity, excessively high 
# or low depth, or excess heterozygosity. This filtering regime was 
# adapted from Pečnerová et al. 2021 (Current Biology).

localrules: genome_bed, smallscaff_bed, sexlink_bed, mito_bed, \
genmap_filt_bed, excess_hetero_bed, genome_sum, smallscaff_sum, \
depth_sum, genmap_filt_sum, repeat_sum, sexlink_sum, excess_hetero_sum, \
excess_hetero_bed, combine_beds

#################################################
# Reference-based filters (dataset independent) #
#################################################

# Whole genome - Use index to get bed file of all sites

rule genome_bed:
	input:
		fai=REF+".fai"
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_genome.bed"
	shell:
		r"""
		awk -v OFS='\t' '{{print $1, "0", $2}}' {input.fai} > {output.bed}
		"""

rule genome_sum:
	input:
		bed=REF_DIR+"/beds/"+REF_NAME+"_genome.bed"
	output:
		sum=REF_DIR+"/beds/"+REF_NAME+"_genome.bed.sum"
	shell:
		r"""
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.bed})
		echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > \
			{output.sum}
		"""

# Scaffold size - remove scaffolds smaller than a threshold

rule smallscaff_bed:
	input:
		genbed=rules.genome_bed.output.bed
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_scaff"+
			str(config["reference"]["min_size"])+".bed"
	params:
		minsize=config["reference"]["min_size"]
	shell:
		r"""
		awk '$3 < {params.minsize}' {input.genbed} > {output.bed}
		"""

rule smallscaff_sum:
	input:
		sum=rules.genome_sum.output.sum,
		bed=rules.smallscaff_bed.output.bed
	output:
		sum=rules.smallscaff_bed.output.bed+".sum"
	params:
		minsize=config["reference"]["min_size"]
	shell:
		r"""
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Scaffolds<{params.minsize}bp\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

# Autosomal scaffolds - User supplied

rule sexlink_bed:
	input:
		rules.genome_bed.output.bed
	output:
		exclbed=REF_DIR+"/beds/"+REF_NAME+"_excl.bed",
		xzbed=REF_DIR+"/beds/"+REF_NAME+"_XZ.bed"
	params:
		xz=config["reference"]["XZ"],
		excl=config["reference"]["exclude"]+config["reference"]["mito"]
	run:
		shell("printf '%s\\n' {params.xz} | "
			"grep -f - {input[0]} > {output.xzbed}")
		shell("printf '%s\\n' {params.xz} {params.excl} | "
			"grep -f - {input[0]} > {output.exclbed}")

rule sexlink_sum:
	input:
		sum=rules.genome_sum.output.sum,
		bed=rules.sexlink_bed.output.exclbed
	output:
		sum=rules.sexlink_bed.output.exclbed+".sum"
	shell:
		r"""
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Autosomes\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

rule mito_bed:
	input:
		genbed=rules.genome_bed.output.bed
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_mito.bed"
	params:
		mito=config["reference"]["mito"]
	shell:
		"""
		grep {params.mito} {input.genbed} > {output.bed}
		"""

# Mappability - Genmap

rule genmap_index:
	input:
		ref=REF
	output:
		directory(REF_DIR+"/genmap/index"),
		multiext(REF_DIR+"/genmap/index/index",
			".ids.concat",".ids.limits",".info.concat",
			".info.limits",".lf.drp",".lf.drp.sbl",".lf.drs",".lf.drv",
			".lf.drv.sbl",".lf.pst",".rev.lf.drp",".rev.lf.drp.sbl",
			".rev.lf.drs",".rev.lf.drv",".rev.lf.drv.sbl",".rev.lf.pst",
			".sa.ind",".sa.len",".sa.val",".txt.concat",".txt.limits"
			)
	log:
		"logs/reffilt/genmap/index/"+REF_NAME+".log"
	conda:
		"../envs/genmap.yaml"
	params:
		out=REF_DIR+"/genmap/index"
	shell:
		"""
		# genmap index annoyingly fails if directory already exists,
		# delete it to keep it happy
		rm -r {output[0]}

		genmap index -F {input.ref} -I {params.out} &> {log}
		"""

rule genmap_map:
	input:
		index=rules.genmap_index.output
	output:
		bed=REF_DIR+"/genmap/map/"+REF_NAME+".bedgraph"
	log:
		"logs/reffilt/genmap/map/"+REF_NAME+".log"
	conda:
		"../envs/genmap.yaml"
	params:
		index=rules.genmap_index.params.out,
		K=config["params"]["genmap"]["K"],
		E=config["params"]["genmap"]["E"],
		out=REF_DIR+"/genmap/map/"+REF_NAME
	threads: lambda wildcards, attempt: attempt
	resources: 
		time=lambda wildcards, attempt: attempt*360
	shell:
		"""
		genmap map -K {params.K} -E {params.E} -I {params.index} \
			-O {params.out} -bg &> {log}
		"""

rule genmap_filt_bed:
	input:
		bed=rules.genmap_map.output.bed
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed"
	log:
		"logs/reffilt/genmap/map/"+REF_NAME+"_lowmap.log"
	params:
		thresh=config["params"]["genmap"]["map_thresh"]
	shell:
		"""
		awk '$4 < {params.thresh}' {input.bed} > {output.bed} 2> {log}
		"""

rule genmap_filt_sum:
	input:
		sum=rules.genome_sum.output.sum,
		bed=rules.genmap_filt_bed.output.bed
	output:
		sum=rules.genmap_filt_bed.output.bed+".sum"
	params:
		map=config["params"]["genmap"]["map_thresh"]
	shell:
		r"""
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Mappability<{params.map}\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

# Repeat filtering - RepeatModeler & RepeatMasker

rule repeat_builddatabase:
	input:
		ref=REF
	output:
		multiext(REF_DIR+"/repeatmodeler/"+REF_NAME+".",
				"nhr","nin","nnd","nni","nog","nsq","translation")
	conda:
		"../envs/repeatmasker.yaml"
	log:
		"logs/reffilt/repeatmodeler/builddatabase/"+REF_NAME+".log"
	params:
		db=REF_DIR+"/repeatmodeler/"+REF_NAME
	shell:
		"""
		BuildDatabase -name {params.db} {input.ref} &> {log}
		"""

rule repeatmodeler:
	input:
		database=rules.repeat_builddatabase.output
	output:
		fa=REF_DIR+"/repeatmodeler/"+REF_NAME+"-families.fa",
		stk=REF_DIR+"/repeatmodeler/"+REF_NAME+"-families.stk",
		log=REF_DIR+"/repeatmodeler/"+REF_NAME+"-rmod.log"
	log:
		"logs/reffilt/repeatmodeler/repeatmodeler/"+REF_NAME+".log"
	conda:
		"../envs/repeatmasker.yaml"
	params:
		db=rules.repeat_builddatabase.params.db,
		ref=REF_NAME
	threads: 10
	resources:
		time="7-00:00:00"
	shadow: "copy-minimal"
	shell:
		"""
		RepeatModeler -database {params.db} -pa {threads} &> {log}
		"""

## Get repeatmasker inputs
repmaskin = []
repmaskin.append(REF)
if config["analyses"]["repeatmasker"]["local_lib"]:
	repmaskin.append(config["analyses"]["repeatmasker"]["local_lib"])
elif config["analyses"]["repeatmasker"]["build_lib"]:
	repmaskin.append(rules.repeatmodeler.output.fa)

rule repeatmasker:
	input:
		repmaskin
	output:
		gff=REF_DIR+"/repeatmasker/"+os.path.basename(REF)+".out.gff"
	log:
		"logs/reffilt/repeatmasker/"+REF_NAME+".log"
	conda:
		"../envs/repeatmasker.yaml"
	params:
		out=REF_DIR+"/repeatmasker/",
		lib="-species '"+config["analyses"]["repeatmasker"]["dfam_lib"]+"'" \
				if config["analyses"]["repeatmasker"]["dfam_lib"] \
				else "-lib "+repmaskin[1],
		gff=os.path.basename(REF)+".out.gff"
	threads: 5
	resources:
		time="12:00:00"
	shadow: "copy-minimal"
	shell:
		"""
		RepeatMasker -pa {threads} {params.lib} -gff -x -no_is \
			-dir {params.out} {input[0]} &> {log}
		"""

rule repeat_sum:
	input:
		sum=rules.genome_sum.output.sum,
		gff=rules.repeatmasker.output.gff
	output:
		sum=rules.repeatmasker.output.gff+".sum"
	shell:
		r"""
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$5-$4-1}}END{{print SUM}}' {input.gff})
		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Repeats\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

#################################################
#   Dataset-based filters (dataset dependent)   #
#################################################

# Global sequencing depth - ANGSD

rule angsd_depth:
	input:
		bamlist=results+"/genotyping/bamlists/"+dataset+"_{set}.bamlist",
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		ref=REF
	output:
		posgz=results+"/genotyping/filters/depthfilt/"+dataset+ \
				"_{set}_chunk{chunk}.pos.gz",
		hist=results+"/genotyping/filters/depthfilt/"+dataset+ \
				"_{set}_chunk{chunk}.depthGlobal"
	log:
		logs+ "/depthfilt/"+dataset+"_{set}_chunk{chunk}.log"
	container:
		"docker://zjnolen/angsd:0.937"
	params:
		out=results+
			"/genotyping/filters/depthfilt/"+dataset+"_{set}_chunk{chunk}"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time="04:00:00"
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
        maxDP=$(echo 100 $nInd | awk '{{print $1 * $2}}')

		angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
			-doCounts 1 -dumpCounts 1 -doDepth 1 \
			-maxDepth $maxDP -out {params.out} &> {log}
		"""

rule combine_depths:
	input:
		lambda w: expand(results+"/genotyping/filters/depthfilt/"+dataset+
						"_{{set}}_chunk{chunk}.depthGlobal",
						chunk=chunklist)
	output:
		results+"/genotyping/filters/depthfilt/"+dataset+
			"_{set}.depthGlobal"
	shell:
		"""
		cat {input} > {output}
		"""

rule summarize_depths:
	input:
		rules.combine_depths.output
	output:
		results+"/genotyping/filters/depthfilt/"+dataset+
			"_{set}_depth.summary"
	conda:
		"../envs/r.yaml"
	params:
		lower=0.025,
		upper=0.975
	threads: lambda wildcards, attempt: attempt*2
	script:
		"../scripts/depth.R"

rule depth_bed:
	input:
		genbed=rules.genome_bed.output.bed,
		quants=rules.summarize_depths.output,
		pos=lambda w: expand(results+"/genotyping/filters/depthfilt/"+dataset+\
						"_{{set}}_chunk{chunk}.pos.gz",
						chunk=chunklist)
	output:
		results+"/genotyping/filters/depthfilt/"+dataset+ \
			"_{set}_depthextremes.bed"
	conda:
		"../envs/bedtools.yaml"
	shell:
		r"""
		lower=$(awk '{{print $2}}' {input.quants})
		upper=$(awk '{{print $3}}' {input.quants})
		for i in {input.pos}; do
			zcat $i | tail -n +2 | \
			awk -v lower=$lower -v upper=$upper '$3 > lower && $3 < upper'
		done | \
		awk '{{print $1"\t"$2-1"\t"$2}}' > {output}
		bedtools merge -i {output} > {output}.tmp
		bedtools subtract -a {input.genbed} -b {output}.tmp > {output}
		rm {output}.tmp
		"""

rule combine_depth_bed:
	input:
		expand(results+"/genotyping/filters/depthfilt/"+dataset+"_{set}" \
			"_depthextremes.bed",set=["all"]+list(set(samples.depth.values)))
	output:
		results+"/genotyping/filters/beds/"+dataset+"_depthfilt.bed"
	log:
		logs+"/depthfilt/"+dataset+"_combinebed.log"
	conda:
		"../envs/bedtools.yaml"
	resources:
		time="04:00:00"
	shell:
		"""
		cat {input} > {output}.tmp
		sort -k1,1 -k2,2n {output}.tmp > {output}.tmp.sort
		rm {output}.tmp

		bedtools merge -i {output}.tmp.sort > {output}
		rm {output}.tmp.sort
		"""

rule depth_sum:
	input:
		sum=rules.genome_sum.output.sum,
		bed=rules.combine_depth_bed.output
	output:
		sum=rules.combine_depth_bed.output[0]+".sum"
	shell:
		r"""
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Depth\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

# # Excess heterozygosity - PCAngsd

# rule excess_hetero_beagle:
# 	input:
# 		bamlist=results+"/genotyping/bamlists/"+dataset+"_all.bamlist",
# 		regions=REF_DIR+"/beds/chunk{chunk}_"+config["chunk_size"]+"bp.rf",
# 		ref=REF
# 	output:
# 		beagle=results+"/genotyping/filters/heterofilt/beagle/chunk/"+dataset+
# 				"_chunk{chunk}_heterofilt.beagle.gz"
# 	log:
# 		logs + "/heterofilt/beagle/"+dataset+"_chunk{chunk}_heterofilt.log"
# 	container:
# 		"docker://zjnolen/angsd:0.937"
# 	params:
# 		gl_model=config["params"]["angsd"]["gl_model"],
# 		extra=config["params"]["angsd"]["extra"],
# 		mapQ=config["mapQ"],
# 		baseQ=config["baseQ"],
# 		pval=config["params"]["angsd"]["snp_pval"],
# 		maf=config["params"]["angsd"]["min_maf"],
# 		out=results+"/genotyping/filters/heterofilt/beagle/chunk/"+dataset+
# 			"_chunk{chunk}_heterofilt"
# 	threads: lambda wildcards, attempt: attempt*2
# 	resources:
# 		time="06:00:00"
# 	shell:
# 		"""
# 		angsd -doGlf 2 -bam {input.bamlist} -GL {params.gl_model} \
# 			-ref {input.ref} -nThreads {threads} {params.extra} \
# 			-doMaf 1 -doMajorMinor 1 -minMapQ {params.mapQ} \
# 			-minQ {params.baseQ} -SNP_pval {params.pval} -minMaf {params.maf} \
# 			-rf {input.regions} -out {params.out} &> {log}
# 		"""

# rule excess_hetero_beagle_merge:
# 	input:
# 		lambda w: expand(results+"/genotyping/filters/heterofilt/beagle/" \
# 						"chunk/"+dataset+"_chunk{chunk}_heterofilt.beagle.gz",
# 						chunk=chunklist)
# 	output:
# 		beagle=results+"/genotyping/filters/heterofilt/beagle/"+dataset+ \
# 				"_heterofilt.beagle.gz"
# 	log:
# 		logs + "/heterofilt/beagle/"+dataset+"_heterofilt.log"
# 	shell:
# 		r"""
# 		set +o pipefail;
# 		zcat {input[0]} | head -n 1 | gzip > {output} 2> {log}

# 		for f in {input}; do
# 			zcat $f | tail -n +2 | gzip >> {output.beagle} 2>> {log}
# 		done
# 		"""

# rule excess_hetero_pcangsd:
# 	input:
# 		beagle=rules.excess_hetero_beagle_merge.output.beagle
# 	output:
# 		cov=results+"/genotyping/filters/heterofilt/"+dataset+
# 			"_heterofilt.cov",
# 		inb=results+"/genotyping/filters/heterofilt/"+dataset+
# 			"_heterofilt.inbreed.sites.npy",
# 		lrt=results+"/genotyping/filters/heterofilt/"+dataset+
# 			"_heterofilt.lrt.sites.npy",
# 		sites=results+"/genotyping/filters/heterofilt/"+dataset+
# 			"_heterofilt.sites"
# 	log:
# 		logs + "/heterofilt/pcangsd/"+dataset+"_heterofilt.log"
# 	params:
# 		out=results+"/genotyping/filters/heterofilt/"+dataset+"_heterofilt"
# 	resources:
# 		time="04:00:00"
# 	threads: 4
# 	shell:
# 		r"""
# 		module load bioinfo-tools
# 		module load PCAngsd/0.982

# 		pcangsd.py -threads {threads} -b {input.beagle} -inbreedSites \
# 			-sites_save -o {params.out} &> {log}
		
# 		cp {output.sites} {output.sites}.bak
# 		sed -i -r 's/(.*)_/\1\t/' {output.sites}
# 		"""

# rule excess_hetero_bed:
# 	input:
# 		inb=rules.excess_hetero_pcangsd.output.inb,
# 		lrt=rules.excess_hetero_pcangsd.output.lrt,
# 		sites=rules.excess_hetero_pcangsd.output.sites,
# 		genbed=rules.genome_bed.output.bed
# 	output:
# 		bed=results+"/genotyping/filters/beds/"+dataset+"_heterofilt.bed",
# 		sites=results+"/genotyping/filters/beds/"+dataset+
# 			"_heterofilt.filt.sites"
# 	params:
# 		hwe=1e-6,
# 		F=-0.95,
# 		win=50000
# 	run:
# 		import pandas as pd
# 		import numpy as np
# 		genbed = pd.read_csv(input.genbed, sep = "\t", header = None,
# 			names = ["chr","start","stop"])
# 		df = pd.read_csv(input.sites, sep = "\t", header = None, 
# 			names = ["chr","pos"])
# 		inb = np.load(input.inb)
# 		lrt = np.load(input.lrt)
# 		df['inb'] = inb.tolist()
# 		df['lrt'] = lrt.tolist()
# 		df = pd.merge(genbed,df)
# 		df['max'] = df['stop']
# 		df = df[(df['inb'] < params.F) & (df['lrt'] < params.hwe)]
# 		df['start'] = df['pos'] - params.win
# 		df['stop'] = df['pos'] + params.win
# 		df['start'][df['start'] < 0] = 0
# 		df['stop'][df['stop'] > df['max']] = df['max']
# 		df.to_csv(output.sites, sep = "\t", columns=['chr','pos'], 
# 			header=False, index=False)
# 		df.to_csv(output.bed, sep = "\t", columns=['chr','start','stop'], 
# 			header=False, index=False)

# rule excess_hetero_sum:
# 	input:
# 		sum=rules.genome_sum.output.sum,
# 		bed=rules.excess_hetero_bed.output.bed
# 	output:
# 		sum=rules.excess_hetero_bed.output.bed+".sum"
# 	shell:
# 		r"""
# 		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {input.bed})
# 		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
# 			awk '{{print "Excess heterozygosity\t"$2-$1"\t"($2-$1)/$2*100}}' \
# 			> {output.sum}
# 		"""

#################################################
#              Combine all filters              #
#################################################

# Joining reference filters - BEDTools

def get_bed_filts(wildcards):
	bedin = []
	bedsum = []
	if config["reference"]["min_size"] and config["reference"]["min_size"] > 0:
		bedin.append(REF_DIR+"/beds/"+REF_NAME+"_scaff"+
					 str(config["reference"]["min_size"])+".bed")
		bedsum.append(REF_DIR+"/beds/"+REF_NAME+"_scaff"+
					 str(config["reference"]["min_size"])+".bed.sum")
	if config["reference"]["XZ"] or config["reference"]["exclude"]:
		bedin.append(REF_DIR+"/beds/"+REF_NAME+"_excl.bed")
		bedsum.append(REF_DIR+"/beds/"+REF_NAME+"_excl.bed.sum")
	if config["analyses"]["genmap"]:
		bedin.append(REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed")
		bedsum.append(REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed.sum")
	if config["analyses"]["repeatmasker"]:
		bedin.append(REF_DIR+"/repeatmasker/"+os.path.basename(REF)+".out.gff")
		bedsum.append(REF_DIR+"/repeatmasker/"+os.path.basename(REF)+
			".out.gff.sum")
	if config["analyses"]["extreme_depth"]:
		bedin.append(results+"/genotyping/filters/beds/"+dataset+
					 "_depthfilt.bed")
		bedsum.append(results+"/genotyping/filters/beds/"+dataset+
					 "_depthfilt.bed.sum")
	# if config["analyses"]["excess_hetero"]:
	# 	bedin.append(results+"/genotyping/filters/beds/"+dataset+
	# 				 "_heterofilt.bed")
	# 	bedsum.append(results+"/genotyping/filters/beds/"+dataset+
	# 				 "_heterofilt.bed.sum")
	return {'gen': REF_DIR+"/beds/"+REF_NAME+"_genome.bed", 
			'sum': REF_DIR+"/beds/"+REF_NAME+"_genome.bed.sum", 
			'filt': bedin, 'sums': bedsum}

checkpoint combine_beds:
	input:
		unpack(get_bed_filts)
	output:
		bed=results+"/genotyping/filters/beds/"+dataset+"_filts.bed",
		lis=results+"/genotyping/filters/beds/"+dataset+"_filts.list",
		sit=results+"/genotyping/filters/beds/"+dataset+"_filts.sites",
		sum=results+"/genotyping/filters/beds/"+dataset+"_filts.sum"
	log:
		logs + "/reffilt/combine.log"
	conda:
		"../envs/bedtools.yaml"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time="04:00:00"
	shell:
		r"""
		printf '%s\n' {input.filt} > {output.lis} 2> {log}
		cat {input.gen} > {output.bed} 2>> {log}

		echo "Name	Length(bp)	Percent" > {output.sum}
		cat {input.sum} >> {output.sum}

		for i in {input.filt}; do
			bedtools subtract -a {output.bed} -b $i > {output.bed}.tmp
			mv {output.bed}.tmp {output.bed}
		done 2>> {log}

		for i in {input.sums}; do
			cat $i >> {output.sum}
		done 2>> {log}

		filtlen=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2-1}}END{{print SUM}}' \
			{output.bed})
		echo $filtlen $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Combined	"$1"	"$1/$2*100}}' >> {output.sum}

		awk '{{print $1"\t"$2+1"\t"$3}}' {output.bed} > {output.sit}.tmp \
			2>> {log}
		sort -k1 {output.sit}.tmp > {output.sit}
		rm {output.sit}.tmp
		"""