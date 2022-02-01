rule qualimap:
    input:
        results+"/bam_dedup/{sample}.mem.bam"
    output:
        results+"/mapping/qualimap/{sample}_mem/qualimapReport.html"
    params:
        outdir = results+"/mapping/qualimap/{sample}_mem"
    shell:
        "qualimap bamqc "
        "-bam {input} "
        "-outdir {params.outdir}"
        
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