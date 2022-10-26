# Creates a filtered list of autosomal sites for analyses to filter the 
# analyses to. This filtering limits analyses to autosomal scaffolds and 
# removes sites with low mappability, low complexity, excessively high 
# or low depth, or excess heterozygosity. This filtering regime was 
# adapted from Pečnerová et al. 2021 (Current Biology).

localrules: genome_bed, smallscaff_bed, sexlink_bed, mito_bed, \
genmap_filt_bed, combine_beds, repeat_sum

#################################################
# Reference-based filters (dataset independent) #
#################################################

# Whole genome - Use index to get bed file of all sites

rule genome_bed:
	input:
		fai=REF+".fai"
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_genome.bed",
		sum=REF_DIR+"/beds/"+REF_NAME+"_genome.bed.sum"
	shell:
		r"""
		# generate bed
		awk -v OFS='\t' '{{print $1, "0", $2}}' {input.fai} > {output.bed}

		# summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > \
			{output.sum}
		"""

# Scaffold size - remove scaffolds smaller than a threshold

rule smallscaff_bed:
	input:
		genbed=rules.genome_bed.output.bed,
		gensum=rules.genome_bed.output.sum
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_scaff"+
			str(config["reference"]["min_size"])+".bed",
		sum=REF_DIR+"/beds/"+REF_NAME+"_scaff"+
			str(config["reference"]["min_size"])+".bed.sum"
	params:
		minsize=config["reference"]["min_size"]
	shell:
		r"""
		# generate bed
		awk '$3 < {params.minsize}' {input.genbed} > {output.bed}

		# summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
			awk '{{print "Scaffolds<{params.minsize}bp\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

# Autosomal scaffolds - User supplied

rule sexlink_bed:
	input:
		genbed=rules.genome_bed.output.bed,
		gensum=rules.genome_bed.output.sum
	output:
		exclbed=REF_DIR+"/beds/"+REF_NAME+"_excl.bed",
		sum=REF_DIR+"/beds/"+REF_NAME+"_excl.bed.sum",
		xzbed=REF_DIR+"/beds/"+REF_NAME+"_XZ.bed"
	params:
		xz=config["reference"]["XZ"],
		excl=config["reference"]["exclude"]+config["reference"]["mito"]
	shell:
		r"""
		# generate beds
		printf '%s\n' {params.xz} | grep -f - {input.genbed} > {output.xzbed}
		printf '%s\n' {params.xz} {params.excl} | grep -f - {input.genbed} > \
			{output.exclbed}
		
		# summarize exclbed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' \
			{output.exclbed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
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
		genbed=rules.genmap_map.output.bed,
		gensum=rules.genome_bed.output.sum
	output:
		bed=REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed",
		sum=REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed.sum"
	log:
		"logs/reffilt/genmap/map/"+REF_NAME+"_lowmap.log"
	params:
		thresh=config["params"]["genmap"]["map_thresh"]
	shell:
		r"""
		# generate bed
		awk '$4 < {params.thresh}' {input.genbed} > {output.bed} 2> {log}

		#summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' \
			{output.bed}) 2>> {log}
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | awk \
			'{{print "Mappability<{params.thresh}\t"$2-$1"\t"($2-$1)/$2*100}}' \
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
		time=10080
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
		time=720
	shadow: "copy-minimal"
	shell:
		"""
		RepeatMasker -pa {threads} {params.lib} -gff -x -no_is \
			-dir {params.out} {input[0]} &> {log}
		"""

rule repeat_sum:
	input:
		sum=rules.genome_bed.output.sum,
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
		bamlist=results+"/genotyping/bamlists/"+dataset+
			"_{population}{dp}.bamlist",
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+"bp.rf",
		ref=REF,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		posgz=temp(results+"/genotyping/filters/depthfilt/"+dataset+ \
				"_{population}{dp}_chunk{chunk}.pos.gz"),
		hist=temp(results+"/genotyping/filters/depthfilt/"+dataset+ \
				"_{population}{dp}_chunk{chunk}.depthGlobal"),
		samphist=temp(results+"/genotyping/filters/depthfilt/"+dataset+ \
				"_{population}{dp}_chunk{chunk}.depthSample"),
		arg=temp(results+"/genotyping/filters/depthfilt/"+dataset+ \
				"_{population}{dp}_chunk{chunk}.arg")
	log:
		logs+ "/depthfilt/"+dataset+"_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		out=results+"/genotyping/filters/depthfilt/"+dataset+ \
			"_{population}{dp}_chunk{chunk}"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*360
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
        maxDP=$(echo 1000 $nInd | awk '{{print $1 * $2}}')

		angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
			-doCounts 1 -dumpCounts 1 -doDepth 1 -maxDepth $maxDP \
			-out {params.out} &> {log}
		"""

rule combine_depths:
	input:
		lambda w: expand(results+"/genotyping/filters/depthfilt/"+dataset+
						"_{{population}}{{dp}}_chunk{chunk}.depthGlobal",
						chunk=chunklist)
	output:
		results+"/genotyping/filters/depthfilt/"+dataset+ \
			"_{population}{dp}.depthGlobal"
	shell:
		"""
		cat {input} > {output}
		"""

rule summarize_depths:
	input:
		rules.combine_depths.output
	output:
		results+"/genotyping/filters/depthfilt/"+dataset+
			"_{population}{dp}_depth.summary"
	conda:
		"../envs/r.yaml"
	params:
		lower=0.025,
		upper=0.975
	threads: lambda wildcards, attempt: attempt*2
	script:
		"../scripts/depth_extremes.R"

rule depth_bed:
	input:
		genbed=rules.genome_bed.output.bed,
		quants=rules.summarize_depths.output,
		pos=lambda w: expand(results+"/genotyping/filters/depthfilt/"+dataset+\
						"_{{population}}{{dp}}_chunk{chunk}.pos.gz",
						chunk=chunklist)
	output:
		results+"/genotyping/filters/depthfilt/"+dataset+ \
			"_{population}{dp}_depthextremes.bed"
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
		beds=expand(results+"/genotyping/filters/depthfilt/"+dataset+
			"_{population}{{dp}}_depthextremes.bed",population=["all"]+
			list(set(samples.depth.values))),
		gensum=rules.genome_bed.output.sum
	output:
		bed=results+"/genotyping/filters/beds/"+dataset+"{dp}_depthfilt.bed",
		sum=results+"/genotyping/filters/beds/"+dataset+"{dp}_depthfilt.bed.sum"
	log:
		logs+"/depthfilt/"+dataset+"{dp}_combinebed.log"
	conda:
		"../envs/bedtools.yaml"
	shadow: "copy-minimal"
	resources:
		time=240
	shell:
		"""
		# combine beds
		cat {input.beds} > {output.bed}.tmp
		sort -k1,1 -k2,2n {output.bed}.tmp > {output.bed}.tmp.sort
		rm {output.bed}.tmp

		bedtools merge -i {output.bed}.tmp.sort > {output.bed}
		rm {output.bed}.tmp.sort

		# summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
			awk '{{print "Depth\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}
		"""

# Missing data across dataset - ANGSD

rule angsd_missdata:
	input:
		bamlist=results+"/genotyping/bamlists/"+dataset+"_{population}{dp}.bamlist",
		regions=REF_DIR+"/beds/chunk{chunk}_"+str(config["chunk_size"])+
			"bp.rf",
		ref=REF,
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		posgz=temp(results+"/genotyping/filters/missdata/"+dataset+ \
				"_{population}{dp}_chunk{chunk}_over{miss}.pos.gz"),
		arg=temp(results+"/genotyping/filters/missdata/"+dataset+ \
				"_{population}{dp}_chunk{chunk}_over{miss}.arg")
	log:
		logs+ "/missdata/"+dataset+"_{population}{dp}_chunk{chunk}_over{miss}.log"
	container:
		angsd_container
	params:
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=results+"/genotyping/filters/missdata/"+dataset+
			"_{population}{dp}_chunk{chunk}_over{miss}"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*360
	shell:
		"""
		nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
		minInd=$(echo $nInd \
			| awk '{{print $1*{wildcards.miss}}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
		
		angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
			-ref {input.ref} -doCounts 1 -dumpCounts 1 -minInd $minInd \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-out {params.out} &> {log}
		"""

rule missdata_bed:
	input:
		pos=lambda w: expand(results+"/genotyping/filters/missdata/"+dataset+
						"_{{population}}{{dp}}_chunk{chunk}_over{{miss}}.pos.gz",
						chunk=chunklist),
		genbed=rules.genome_bed.output.bed,
		gensum=rules.genome_bed.output.sum
	output:
		bed=results+"/genotyping/filters/missdata/"+dataset+
			"_{population}{dp}_under{miss}.bed",
		sum=results+"/genotyping/filters/missdata/"+dataset+
			"_{population}{dp}_under{miss}.bed.sum",
		tmp=temp(results+"/genotyping/filters/missdata/"+dataset+
			"_{population}{dp}_under{miss}.bed.tmp")
	conda:
		"../envs/bedtools.yaml"
	shell:
		r"""
		# generate bed
		> {output.tmp}
		for i in {input.pos}; do
			zcat $i | tail -n +2 >> {output.tmp}
		done
		
		awk '{{print $1"\t"$2-1"\t"$2}}' {output.tmp} > {output.bed}
		bedtools merge -i {output.bed} > {output.tmp}
		bedtools subtract -a {input.genbed} -b {output.tmp} > {output.bed}
		
		# summarize bed
		perc=$(echo {wildcards.miss} | awk '{{print $1*100}}')
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) $perc | \
			awk '{{print "Min "$3"% data ({wildcards.population})\t"$2-$1"\t" \
			($2-$1)/$2*100}}' > {output.sum}
		"""

#################################################
#              Combine all filters              #
#################################################

# Joining reference filters - BEDTools

def get_bed_filts(wildcards):
	bedin = []
	bedsum = []
	# add minimum size filter if set
	if config["reference"]["min_size"] and config["reference"]["min_size"] > 0:
		bedin.append(REF_DIR+"/beds/"+REF_NAME+"_scaff"+
					 str(config["reference"]["min_size"])+".bed")
		bedsum.append(REF_DIR+"/beds/"+REF_NAME+"_scaff"+
					 str(config["reference"]["min_size"])+".bed.sum")
	# add sex chromosome filter if set
	if config["reference"]["XZ"] or config["reference"]["exclude"]:
		bedin.append(REF_DIR+"/beds/"+REF_NAME+"_excl.bed")
		bedsum.append(REF_DIR+"/beds/"+REF_NAME+"_excl.bed.sum")
	# add mappability filter if set
	if config["analyses"]["genmap"]:
		bedin.append(REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed")
		bedsum.append(REF_DIR+"/beds/"+REF_NAME+"_lowmap.bed.sum")
	# add repeat filter if set
	if config["analyses"]["repeatmasker"]:
		bedin.append(REF_DIR+"/repeatmasker/"+os.path.basename(REF)+".out.gff")
		bedsum.append(REF_DIR+"/repeatmasker/"+os.path.basename(REF)+
			".out.gff.sum")
	# add global depth extremes filter if set
	if config["analyses"]["extreme_depth"]:
		bedin.append(results+"/genotyping/filters/beds/"+dataset+
					 "{dp}_depthfilt.bed")
		bedsum.append(results+"/genotyping/filters/beds/"+dataset+
					 "{dp}_depthfilt.bed.sum")
	# add dataset level missing data filter if set
	if config["analyses"]["dataset_missing_data"]:
		bedin.append(results+"/genotyping/filters/missdata/"+dataset+
			"_all{dp}_under"+
			str(config["analyses"]["dataset_missing_data"])+".bed")
		bedsum.append(results+"/genotyping/filters/missdata/"+dataset+
			"_all{dp}_under"+
			str(config["analyses"]["dataset_missing_data"])+".sum")
	# add population level missing data filter if set
	if config["analyses"]["population_missing_data"]:
		[bedin.append(i) for i in \
			expand(results+"/genotyping/filters/missdata/"+dataset+
					"_{population}{{dp}}_under"+
				str(config["analyses"]["population_missing_data"])+".bed",
				population=pop_list)]
		[bedsum.append(i) for i in \
			expand(results+"/genotyping/filters/missdata/"+dataset+
					"_{population}{{dp}}_under"+
				str(config["analyses"]["population_missing_data"])+".sum",
				population=pop_list)]	
	# if config["analyses"]["excess_hetero"]:
	# 	bedin.append(results+"/genotyping/filters/beds/"+dataset+
	# 				 "_heterofilt.bed")
	# 	bedsum.append(results+"/genotyping/filters/beds/"+dataset+
	# 				 "_heterofilt.bed.sum")
	return {'gen': REF_DIR+"/beds/"+REF_NAME+"_genome.bed", 
			'sum': REF_DIR+"/beds/"+REF_NAME+"_genome.bed.sum", 
			'filt': bedin, 'sums': bedsum}

rule combine_beds:
	input:
		unpack(get_bed_filts)
	output:
		bed=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.bed",
		lis=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.list",
		sit=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sites",
		sum=results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sum"
	log:
		logs + "/reffilt/combine{dp}.log"
	conda:
		"../envs/bedtools.yaml"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=240
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

rule filter_summary:
	input:
		results+"/genotyping/filters/beds/"+dataset+"{dp}_filts.sum"
	output:
		report(results+"/genotyping/filters/"+dataset+"{dp}_filts.html",
				category="Quality Control",
				labels={
					"Topic":"Site filters",
					"Type":"Table"
				})
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"