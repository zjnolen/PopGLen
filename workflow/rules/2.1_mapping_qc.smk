rule samtools_flagstat:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.flagstat"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule samtools_idxstats:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai"
    output:
        "{prefix}.idxstats"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

rule qualimap:
    input:
        "results/mapping/bams/{sample}{dp}.rmdup.realn.bam"
    output:
        "results/mapping/qc/qualimap/{sample}{dp}/qualimapReport.html",
        "results/mapping/qc/qualimap/{sample}{dp}/genome_results.txt"
    params:
        out="results/mapping/qc/qualimap/{sample}{dp}"
    conda:
        "../envs/qualimap.yaml"
    log:
        "logs/mapping/qualimap/{sample}{dp}.log"
    resources:
        time=360
    shell:
        """
        unset DISPLAY

        qualimap bamqc --java-mem-size={resources.mem_mb}M -bam {input} \
            -outdir {params.out} 2> {log}
        """