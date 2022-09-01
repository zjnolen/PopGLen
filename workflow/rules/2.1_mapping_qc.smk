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
        results+"/mapping/{sample}.rmdup.realn.bam"
    output:
        results+"/qc/qualimap/{sample}/qualimapReport.html",
        results+"/qc/qualimap/{sample}/genome_results.txt"
    params:
        out=results+"/qc/qualimap/{sample}"
    conda:
        "../envs/qualimap.yaml"
    log:
        logs + "/qualimap/{sample}.log"
    resources:
        time=360
    shell:
        """
        unset DISPLAY

        qualimap bamqc --java-mem-size={resources.mem_mb}M -bam {input} \
            -outdir {params.out} 2> {log}
        """