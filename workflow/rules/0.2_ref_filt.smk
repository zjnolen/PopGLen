# Creates a filtered list of autosomal sites for analyses to filter the 
# analyses to. This filtering limits analyses to autosomal scaffolds and 
# removes sites with low mappability, low complexity, excessively high 
# or low depth, or excess heterozygosity. This filtering regime was 
# adapted from Pečnerová et al. 2021 (Current Biology).

localrules: genome_bed, smallscaff_bed, sexlink_bed, genmap_filt_bed, combine_beds, repeat_sum

#################################################
# Reference-based filters (sample independent) #
#################################################

# Whole genome - Use index to get bed file of all sites

rule genome_bed:
	input:
		fai="results/ref/{ref}/{ref}.fa.fai"
	output:
		bed="results/ref/{ref}/beds/genome.bed",
		sum="results/ref/{ref}/beds/genome.bed.sum"
	log:
		"logs/ref/genome_bed/{ref}.log"
	conda:
		"../envs/shell.yaml"
	shell:
		r"""
		(# generate bed
		awk -v OFS='\t' '{{print $1, "0", $2}}' {input.fai} > {output.bed}

		# summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > \
			{output.sum}) 2> {log}
		"""

# Scaffold size - remove scaffolds smaller than a threshold

rule smallscaff_bed:
	input:
		genbed="results/ref/{ref}/beds/genome.bed",
		gensum="results/ref/{ref}/beds/genome.bed.sum"
	output:
		bed="results/datasets/{dataset}/filters/small_scaffs/{ref}_scaff{size}bp.bed",
		sum="results/datasets/{dataset}/filters/small_scaffs/{ref}_scaff{size}bp.bed.sum"
	log:
		"logs/{dataset}/filters/smallscaff/{ref}_scaff{size}bp.log"
	conda:
		"../envs/shell.yaml"
	params:
		minsize="{size}"
	shell:
		r"""
		(# generate bed
		awk '$3 < {params.minsize}' {input.genbed} > {output.bed}

		# summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
			awk '{{print "Scaffolds<{params.minsize}bp\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}) 2> {log}
		"""

# Autosomal scaffolds - User supplied

rule sexlink_bed:
	input:
		genbed="results/ref/{ref}/beds/genome.bed",
		gensum="results/ref/{ref}/beds/genome.bed.sum"
	output:
		exclbed="results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed",
		sum="results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_excl.bed.sum",
		sexbed="results/datasets/{dataset}/filters/sex-link_mito_excl/{ref}_sex-linked.bed"
	log:
		"logs/{dataset}/filters/sex-link_mito_excl/{ref}.log"
	conda:
		"../envs/shell.yaml"
	params:
		sex=config["reference"]["sex-linked"],
		excl=config["reference"]["exclude"],
		mito=config["reference"]["mito"]
	shell:
		r"""
		(# generate beds
		printf '%s\n' {params.sex} | grep -f - {input.genbed} > {output.sexbed}
		printf '%s\n' {params.sex} {params.excl} {params.mito} | grep -f - {input.genbed} > \
			{output.exclbed}
		
		# summarize exclbed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' \
			{output.exclbed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
			awk '{{print "Autosomes\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}) 2> {log}
		"""

# Mappability - Genmap

rule genmap_index:
	input:
		ref="results/ref/{ref}/{ref}.fa"
	output:
		fold=directory("results/ref/{ref}/genmap/index"),
		files=multiext("results/ref/{ref}/genmap/index/index",
			".ids.concat",".ids.limits",".info.concat",
			".info.limits",".lf.drp",".lf.drp.sbl",".lf.drs",".lf.drv",
			".lf.drv.sbl",".lf.pst",".rev.lf.drp",".rev.lf.drp.sbl",
			".rev.lf.drs",".rev.lf.drv",".rev.lf.drv.sbl",".rev.lf.pst",
			".sa.ind",".sa.len",".sa.val",".txt.concat",".txt.limits"
			)
	log:
		"logs/ref/genmap/index/{ref}.log"
	conda:
		"../envs/genmap.yaml"
	threads: lambda wildcards, attempt: attempt
	shell:
		"""
		# genmap index annoyingly fails if directory already exists,
		# delete it to keep it happy
		rm -r {output.fold} 2> {log}

		genmap index -F {input.ref} -I {output.fold} &>> {log}
		"""

rule genmap_map:
	input:
		fold=rules.genmap_index.output.fold,
		files=rules.genmap_index.output.files
	output:
		bed="results/ref/{ref}/genmap/map/{ref}.bedgraph"
	log:
		"logs/ref/genmap/map/{ref}.log"
	conda:
		"../envs/genmap.yaml"
	params:
		K=config["params"]["genmap"]["K"],
		E=config["params"]["genmap"]["E"],
		out=lambda w, output: os.path.splitext(output.bed)[0]
	threads: lambda wildcards, attempt: attempt
	resources: 
		time=lambda wildcards, attempt: attempt*360
	shell:
		"""
		genmap map -K {params.K} -E {params.E} -I {input.fold} \
			-O {params.out} -bg &> {log}
		"""

rule genmap_filt_bed:
	input:
		genbed="results/ref/{ref}/genmap/map/{ref}.bedgraph",
		gensum="results/ref/{ref}/beds/genome.bed.sum"
	output:
		bed="results/datasets/{dataset}/filters/lowmap/{ref}_lowmap.bed",
		sum="results/datasets/{dataset}/filters/lowmap/{ref}_lowmap.bed.sum"
	log:
		"logs/{dataset}/filters/lowmap/{ref}_lowmap.log"
	conda:
		"../envs/shell.yaml"
	params:
		thresh=config["params"]["genmap"]["map_thresh"]
	shell:
		r"""
		# generate bed
		(awk '$4 < {params.thresh}' {input.genbed} > {output.bed}

		#summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' \
			{output.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | awk \
			'{{print "Mappability<{params.thresh}\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}) 2> {log}
		"""

# Repeat filtering - RepeatModeler & RepeatMasker

rule repeat_builddatabase:
	input:
		ref="results/ref/{ref}/{ref}.fa"
	output:
		multiext("results/ref/{ref}/repeatmodeler/{ref}.",
				"nhr","nin","nnd","nni","nog","nsq","translation")
	conda:
		"../envs/repeatmasker.yaml"
	log:
		"logs/ref/repeatmodeler/builddatabase/{ref}.log"
	params:
		db=lambda w, output: os.path.splitext(output[0])[0]
	shell:
		"""
		BuildDatabase -name {params.db} {input.ref} &> {log}
		"""

rule repeatmodeler:
	input:
		database=rules.repeat_builddatabase.output
	output:
		fa="results/ref/{ref}/repeatmodeler/{ref}-families.fa",
		stk="results/ref/{ref}/repeatmodeler/{ref}-families.stk",
		log="results/ref/{ref}/repeatmodeler/{ref}-rmod.log"
	log:
		"logs/ref/repeatmodeler/repeatmodeler/{ref}.log"
	conda:
		"../envs/repeatmasker.yaml"
	params:
		db=lambda w, input: os.path.splitext(input[0])[0],
		ref="{ref}"
	threads: 10
	resources:
		time=10080
	shadow: "minimal"
	shell:
		"""
		RepeatModeler -database {params.db} -pa {threads} &> {log}
		"""

rule repeatmasker:
	input:
		unpack(get_repmaskin)
	output:
		gff="results/ref/{ref}/repeatmasker/{ref}.fa.out.gff"
	log:
		"logs/ref/repeatmasker/{ref}.log"
	conda:
		"../envs/repeatmasker.yaml"
	params:
		out=lambda w, output: os.path.dirname(output.gff),
		libpre="-species" if config["analyses"]["repeatmasker"]["dfam_lib"] else "-lib",
		lib=lambda w, input: f"'{config['analyses']['repeatmasker']['dfam_lib']}'" if config["analyses"]["repeatmasker"]["dfam_lib"] else input.lib
	threads: 5
	resources:
		time=720
	shadow: "shallow"
	shell:
		"""
		RepeatMasker -pa {threads} {params.libpre} {params.lib} -gff -x -no_is \
			-dir {params.out} {input.ref} &> {log}
		"""

rule repeat_sum:
	input:
		sum="results/ref/{ref}/beds/genome.bed.sum",
		gff="results/ref/{ref}/repeatmasker/{ref}.fa.out.gff"
	output:
		sum="results/ref/{ref}/repeatmasker/{ref}.fa.out.gff.sum"
	log:
		"logs/ref/repeatmasker/summarize_gff/{ref}.log"
	conda:
		"../envs/shell.yaml"
	shell:
		r"""
		(len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$5-$4-1}}END{{print SUM}}' {input.gff})
		echo $len $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Repeats\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}) &> {log}
		"""

#################################################
#   Sample-based filters (sample dependent)   #
#################################################

# Global sequencing depth - ANGSD

rule angsd_depth:
	input:
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
		regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
		ref="results/ref/{ref}/{ref}.fa",
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		posgz=temp("results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.pos.gz"),
		hist=temp("results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.depthGlobal"),
		samphist=temp("results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.depthSample"),
		arg=temp("results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.arg")
	log:
		"logs/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_chunk{chunk}.log"
	container:
		angsd_container
	params:
		out=lambda w, output: os.path.splitext(output.arg)[0]
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*720
	shell:
		"""
		(nInd=$(cat {input.bamlist} | wc -l | awk '{{print $1+1}}')
        maxDP=$(echo 1000 $nInd | awk '{{print $1 * $2}}')

		angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
			-doCounts 1 -dumpCounts 1 -doDepth 1 -maxDepth $maxDP \
			-out {params.out}) 2> {log}
		"""

rule combine_depths:
	input:
		lambda w: expand("results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.depthGlobal",
						chunk=chunklist)
	output:
		"results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}.depthGlobal"
	log:
		"logs/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_combined.log"
	conda:
		"../envs/shell.yaml"
	shell:
		"""
		cat {input} > {output} 2> {log}
		"""

rule summarize_depths:
	input:
		"results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}.depthGlobal"
	output:
		"results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth.summary"
	log:
		"logs/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth_extremes.log"
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
		genbed="results/ref/{ref}/beds/genome.bed",
		quants="results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_depth.summary",
		pos=lambda w: expand("results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}.pos.gz",
						chunk=chunklist)
	output:
		"results/datasets/{dataset}/filters/depth/{dataset}.{ref}_{population}{dp}_extreme-depth.bed"
	log:
		"logs/{dataset}/filters/depth/bed/{dataset}.{ref}_{population}{dp}.log"
	conda:
		"../envs/bedtools.yaml"
	shell:
		r"""
		(lower=$(awk '{{print $2}}' {input.quants})
		upper=$(awk '{{print $3}}' {input.quants})
		for i in {input.pos}; do
			zcat $i | tail -n +2 | \
			awk -v lower=$lower -v upper=$upper '$3 > lower && $3 < upper'
		done | \
		awk '{{print $1"\t"$2-1"\t"$2}}' > {output}
		bedtools merge -i {output} > {output}.tmp
		bedtools subtract -a {input.genbed} -b {output}.tmp > {output}
		rm {output}.tmp) 2> {log}
		"""

rule combine_depth_bed:
	input:
		beds=expand("results/datasets/{{dataset}}/filters/depth/{{dataset}}.{{ref}}_{population}{{dp}}_extreme-depth.bed",population=["all"]+list(set(samples.depth.values))),
		gensum="results/ref/{ref}/beds/genome.bed.sum"
	output:
		bed="results/datasets/{dataset}/filters/depth/{dataset}.{ref}{dp}_extreme-depth.bed",
		sum="results/datasets/{dataset}/filters/depth/{dataset}.{ref}{dp}_extreme-depth.bed.sum"
	log:
		"logs/{dataset}/filters/depth/bed/{dataset}.{ref}{dp}_combine-bed.log"
	conda:
		"../envs/bedtools.yaml"
	shadow: "minimal"
	resources:
		time=240
	shell:
		"""
		# combine beds
		(cat {input.beds} > {output.bed}.tmp
		sort -k1,1 -k2,2n {output.bed}.tmp > {output.bed}.tmp.sort
		rm {output.bed}.tmp

		bedtools merge -i {output.bed}.tmp.sort > {output.bed}
		rm {output.bed}.tmp.sort

		# summarize bed
		len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' {output.bed})
		echo $len $(awk -F "\t" '{{print $2}}' {input.gensum}) | \
			awk '{{print "Depth\t"$2-$1"\t"($2-$1)/$2*100}}' \
			> {output.sum}) 2> {log}
		"""

# Missing data across dataset - ANGSD

rule angsd_missdata:
	input:
		bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
		regions="results/datasets/{dataset}/filters/chunks/{ref}_chunk{chunk}.rf",
		ref="results/ref/{ref}/{ref}.fa",
		bams=get_bamlist_bams,
		bais=get_bamlist_bais
	output:
		posgz=temp("results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.pos.gz"),
		arg=temp("results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.arg")
	log:
		"logs/{dataset}/filters/angsd_missdata/{dataset}.{ref}_{population}{dp}_chunk{chunk}_over{miss}.log"
	container:
		angsd_container
	params:
		nind=lambda w: len(get_samples_from_pop(w.population)),
		extra=config["params"]["angsd"]["extra"],
		mapQ=config["mapQ"],
		baseQ=config["baseQ"],
		out=lambda w, output: os.path.splitext(output.arg)[0]
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=lambda wildcards, attempt: attempt*360
	shell:
		"""
		(minInd=$(echo {params.nind} \
			| awk '{{print $1*{wildcards.miss}}}' \
			| awk '{{print int($1) + ( $1!=int($1) && $1>=0 )}}')
		
		angsd -bam {input.bamlist} -nThreads {threads} -rf {input.regions} \
			-ref {input.ref} -doCounts 1 -dumpCounts 1 -minInd $minInd \
			{params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
			-out {params.out}) 2> {log}
		"""

rule missdata_bed:
	input:
		pos=lambda w: expand("results/datasets/{{dataset}}/filters/missdata/{{dataset}}.{{ref}}_{{population}}{{dp}}_chunk{chunk}_over{{miss}}.pos.gz",
						chunk=chunklist),
		genbed="results/ref/{ref}/beds/genome.bed",
		gensum="results/ref/{ref}/beds/genome.bed.sum"
	output:
		bed="results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_under{miss}.bed",
		sum="results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_under{miss}.bed.sum",
		tmp=temp("results/datasets/{dataset}/filters/missdata/{dataset}.{ref}_{population}{dp}_under{miss}.bed.tmp")
	log:
		"logs/{dataset}/filters/missdata_bed/{dataset}.{ref}_{population}{dp}_under{miss}.log"
	conda:
		"../envs/bedtools.yaml"
	shell:
		r"""
		# generate bed
		(> {output.tmp}
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
			($2-$1)/$2*100}}' > {output.sum}) 2> {log}
		"""

#################################################
#              Combine all filters              #
#################################################

rule combine_beds:
	input:
		unpack(get_bed_filts)
	output:
		bed="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.bed",
		lis="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.list",
		sit="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sites",
		sum="results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sum"
	log:
		"logs/{dataset}/filters/combine/{dataset}.{ref}{dp}_combine_beds.log"
	conda:
		"../envs/bedtools.yaml"
	threads: lambda wildcards, attempt: attempt*2
	resources:
		time=240
	shell:
		r"""
		(printf '%s\n' {input.filt} > {output.lis}
		cat {input.gen} > {output.bed}

		echo "Name	Length(bp)	Percent" > {output.sum}
		cat {input.sum} >> {output.sum}

		for i in {input.filt}; do
			bedtools subtract -a {output.bed} -b $i > {output.bed}.tmp
			mv {output.bed}.tmp {output.bed}
		done

		for i in {input.sums}; do
			cat $i >> {output.sum}
		done

		filtlen=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2-1}}END{{print SUM}}' \
			{output.bed})
		echo $filtlen $(awk -F "\t" '{{print $2}}' {input.sum}) | \
			awk '{{print "Combined	"$1"	"$1/$2*100}}' >> {output.sum}

		awk '{{print $1"\t"$2+1"\t"$3}}' {output.bed} > {output.sit}.tmp
		sort -V {output.sit}.tmp > {output.sit}
		rm {output.sit}.tmp) 2> {log}
		"""

rule filter_summary:
	input:
		"results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.sum"
	output:
		report("results/datasets/{dataset}/filters/combined/{dataset}.{ref}{dp}_filts.html",
				category="Quality Control",
				labels={
					"Topic":"Site filters",
					"Type":"Table"
				})
	log:
		"logs/{dataset}/filters/combine/{dataset}.{ref}{dp}_tsv2html.log"
	conda:
		"../envs/r-rectable.yaml"
	script:
		"../scripts/tsv2html.Rmd"