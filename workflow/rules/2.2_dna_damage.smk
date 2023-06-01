# Rules for visualizing and correcting signatures of post-mortem DNA damage


rule damageprofiler:
    """
    Estimates varous metrics related to post-mortem DNA damage. Informative 
    rather than corrective.
    """
    input:
        bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        multiext(
            "results/mapping/qc/damageprofiler/{sample}.{ref}/",
            "5pCtoT_freq.txt",
            "3pGtoA_freq.txt",
            "Length_plot.pdf",
            "DamagePlot_five_prime.svg",
            "DamagePlot.pdf",
            "DamagePlot_three_prime.svg",
            "DamageProfiler.log",
            "lgdistribution.txt",
            "edit_distance.svg",
            "edit_distance.pdf",
            "editDistance.txt",
            "Length_plot_combined_data.svg",
            "Length_plot_forward_reverse_separated.svg",
            "misincorporation.txt",
            "5p_freq_misincorporations.txt",
            "3p_freq_misincorporations.txt",
            "DNA_comp_genome.txt",
            "DNA_composition_sample.txt",
            "dmgprof.json",
        ),
    log:
        "logs/mapping/damageprofiler/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/damageprofiler/{sample}.{ref}.log"
    conda:
        "../envs/damageprofiler.yaml"
    params:
        out=lambda w, output: os.path.dirname(output[0]),
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        damageprofiler -Xmx{resources.mem_mb}m -i {input.bam} -r {input.ref} \
            -o {params.out} &> {log}
        """


rule mapDamage2_rescaling:
    """
    Estimates various metrics related to post-mortem DNA damage and rescales base 
    quality scores to correct for damage.
    """
    input:
        bam="results/mapping/bams/{sample}.{ref}.rmdup.realn.bam",
        ref="results/ref/{ref}/{ref}.fa",
    output:
        outdir=directory("results/mapping/qc/mapdamage/{sample}.{ref}/"),
        rescaled="results/mapping/bams/{sample}.{ref}.rmdup.realn.rescaled.bam",
    log:
        "logs/mapping/mapdamage/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/mapdamage/{sample}.{ref}.log"
    container:
        mapdamage_container
    threads: 4
    resources:
        runtime=1440,
    params:
        tmp="results/mapping/qc/mapdamage/{sample}.{ref}/{sample}.{ref}.rmdup.realn.rescaled.bam",
    shell:
        """
        (mapDamage -i {input.bam} -r {input.ref} -d {output.outdir} --rescale
        mv {params.tmp} {output.rescaled}) &> {log}
        """
