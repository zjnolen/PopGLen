rule mapdamage:
    input:
        ref=genome_file(),
        bam=results+"/bam_dedup/{sample}.mem.bam"
    output:
        log=results+"/mapping/mapdamage/{sample}/Runtime_log.txt",  # output folder is infered from this file, so it needs to be the same folder for all output files
        GtoA3p=results+"/mapping/mapdamage/{sample}/3pGtoA_freq.txt",
        CtoT5p=results+"/mapping/mapdamage/{sample}/5pCtoT_freq.txt",
        dnacomp=results+"/mapping/mapdamage/{sample}/dnacomp.txt",
        frag_misincorp=results+"/mapping/mapdamage/{sample}/Fragmisincorporation_plot.pdf",
        len=results+"/mapping/mapdamage/{sample}/Length_plot.pdf",
        lg_dist=results+"/mapping/mapdamage/{sample}/lgdistribution.txt",
        misincorp=results+"/mapping/mapdamage/{sample}/misincorporation.txt"
    log:
        "logs/mapdamage/{sample}.log"
    resources:
        time="04:00:00"
    wrapper:
        "0.84.0/bio/mapdamage2"

rule samtools_flagstat:
    input:
        intermediate+"/mapping/{sample}.mem.bam"
    output:
        results+"/samtools/stats/{sample}.mem.flagstat"
    wrapper:
        "v1.0.0/bio/samtools/flagstat"

rule endorspy:
    input:
        results+"/samtools/stats/{sample}.mem.flagstat"
    output:
        results+"/endorspy/{sample}_mem_endogenous_dna_mqc.json"
    params:
        outprefix=results+"/endorspy/{sample}_mem"
    conda:
        "../envs/endorspy.yaml"
    shell:
        "endorspy -o json -n {params.outprefix} {input}"

rule qualimap:
    input:
        results+"/dedup/{sample}.mem.bam"
    output:
        results+"/qualimap/{sample}_mem_dedup/qualimapReport.html",
        results+"/qualimap/{sample}_mem_dedup/genome_results.txt"
    params:
        outdir=results+"/qualimap/{sample}_mem_dedup"
    conda:
        "../envs/qualimap.yaml"
    log:
        "logs/qualimap/{sample}_mem_dedup.log"
    resources:
        time="06:00:00",
        mem_mb=5120
    shell:
        """
        qualimap bamqc --java-mem-size={resources.mem_mb}M -bam {input} \
            -outdir {params.outdir} 2> {log}
        """