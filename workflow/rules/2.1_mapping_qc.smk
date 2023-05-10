rule samtools_flagstat:
	input:
		"results/mapping/mapped/{prefix}.bam"
	output:
		"results/mapping/mapped/{prefix}.flagstat"
	log:
		"logs/mapping/samtools/flagstat/{prefix}.log"
	conda:
		"../envs/samtools.yaml"
	shell:
		"""
		samtools flagstat {input} > {output} 2> {log}
		"""

rule qualimap:
	input:
		unpack(get_final_bam)
	output:
		html="results/mapping/qc/qualimap/{sample}.{ref}/qualimapReport.html",
		txt="results/mapping/qc/qualimap/{sample}.{ref}/genome_results.txt"
	params:
		out=lambda w, output: os.path.dirname(output.html)
	conda:
		"../envs/qualimap.yaml"
	log:
		"logs/mapping/qualimap/{sample}.{ref}.log"
	resources:
		time=360
	shell:
		"""
		(unset DISPLAY

		qualimap bamqc --java-mem-size={resources.mem_mb}M -bam {input.bam} \
			-outdir {params.out}) &> {log}
		"""

# Per individual estimates of endogenous content - primarily useful for 
# historical samples or when using references that are not specific to 
# the species of interest

rule endo_cont:
	input:
		get_endo_cont_stat
	output:
		"results/mapping/qc/endogenous_content/{sample}.{ref}.endo"
	conda:
		"../envs/shell.yaml"
	log:
		"logs/mapping/endogenous_content/{sample}.{ref}.log"
	shell:
		r"""
		(total=$(grep -E "^[0-9]+ \+ [0-9]+ in total" {input} \
			| awk '{{print $1}}')
		mapped=$(grep -E "^[0-9]+ \+ [0-9]+ mapped" {input} \
			| awk '{{print $1}}')
		primary=$(grep -E "^[0-9]+ \+ [0-9]+ primary mapped" {input} \
			| awk '{{print $1}}')
		
		echo $total $mapped $primary {wildcards.sample} | \
			awk '{{printf "%s\t%.3f\t%.3f\n",$4,$2/$1*100,$3/$1*100}}' \
			> {output}) 2> {log}
		"""

rule compile_endo_cont:
	input:
		lambda w: expand("results/mapping/qc/endogenous_content/{sample}.{{ref}}.endo", sample=get_samples_from_pop("all"))
	output:
		"results/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_all.endo.tsv"
	log:
		"logs/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_compile-endocont.log"
	conda:
		"../envs/shell.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*15
	shell:
		"""
		(echo "sample	perc.endo	perc.prim.endo" > {output}
		cat {input} >> {output}) 2> {log}
		"""