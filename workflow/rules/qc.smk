localrules: pop_meanDepth

rule samtools_flagstat:
    input:
        results+"/mapped/{sample}.bam"
    output:
        results+"/samtools/stats/{sample}.flagstat"
    wrapper:
        "v1.0.0/bio/samtools/flagstat"

rule endorspy:
    input:
        results+"/samtools/stats/{sample}.flagstat"
    output:
        results+"/endorspy/{sample}_endogenous_dna_mqc.json"
    params:
        outprefix=results+"/endorspy/{sample}"
    conda:
        "../envs/endorspy.yaml"
    shell:
        "endorspy -o json -n {params.outprefix} {input}"

rule qualimap:
    input:
        results+"/dedup/{sample}.bam"
    output:
        results+"/qualimap/{sample}/qualimapReport.html",
        results+"/qualimap/{sample}/genome_results.txt"
    params:
        outdir=results+"/qualimap/{sample}"
    conda:
        "../envs/qualimap.yaml"
    log:
        logs + "/qualimap/{sample}_dedup.log"
    resources:
        time="06:00:00"
    shell:
        """
        qualimap bamqc --java-mem-size={resources.mem_mb}M -bam {input} \
            -outdir {params.outdir} 2> {log}
        """

rule angsd_doDepth:
    input:
        results+"/angsd/bamlists/{population}.bamlist"
    output:
        depthSample=results+"/depth/{population}_chr{chrom}.depthSample",
        depthGlobal=results+"/depth/{population}_chr{chrom}.depthGlobal"
    params:
        out_prefix=results+"/depth/{population}_chr{chrom}"
    log:
        logs + "angsd/depth/{population}_chr{chrom}.log"
    conda:
        "../envs/angsd.yaml"
    resources:
        time=lambda wildcards, attempt: attempt*360
    shell:
        """
        angsd -bam {input} -doDepth 1 -doCounts 1 -r {wildcards.chrom} \
            -maxDepth 200 -out {params.out_prefix} &> {log}
        """

rule pop_meanDepth:
    input:
        lambda w: expand(results+"/depth/{{population}}_chr{chrom}.depthGlobal", chrom=get_contigs())
    output:
        results+"/depth/{population}.depthMean"
    shell:
        """
        cat {input} | workflow/scripts/pop_depth.py > {output}
        """