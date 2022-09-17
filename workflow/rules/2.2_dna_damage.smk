rule damageprofiler:
	input:
		bam="results/mapping/{sample}.rmdup.realn.bam",
		ref=REF
	output:
		multiext("results/mapping/qc/damageprofiler/{sample}/",
				"5pCtoT_freq.txt","3pGtoA_freq.txt","Length_plot.pdf",
				"DamagePlot_five_prime.svg","DamagePlot.pdf",
				"DamagePlot_three_prime.svg","DamageProfiler.log",
				"lgdistribution.txt","edit_distance.svg","edit_distance.pdf",
				"editDistance.txt","Length_plot_combined_data.svg",
				"Length_plot_forward_reverse_separated.svg",
				"misincorporation.txt","5p_freq_misincorporations.txt",
				"3p_freq_misincorporations.txt","DNA_comp_genome.txt",
				"DNA_composition_sample.txt","dmgprof.json")
	log:
		"logs/mapping/damageprofiler/{sample}.log"
	conda:
		"../envs/damageprofiler.yaml"
	params:
		out="results/mapping/qc/damageprofiler/{sample}"
	threads: lambda wildcards, attempt: attempt
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		"""
		damageprofiler -Xmx{resources.mem_mb}m -i {input.bam} -r {input.ref} \
			-o {params.out} &> {log}
		"""