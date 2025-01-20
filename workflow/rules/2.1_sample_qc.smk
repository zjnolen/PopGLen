# Rules for various assessments of samples, mainly regarding alignment quality


localrules:
    combine_sample_qc,
    endo_cont,
    compile_endo_cont,
    merge_ibs_ref_bias,
    merge_ind_depth,


rule samtools_flagstat:
    """
    Estimate statistics for samtools flags. Needed to assess mapping percentage
    """
    input:
        "results/{prefix}.bam",
    output:
        "results/{prefix}.flagstat",
    container:
        samtools_container
    log:
        "logs/mapping/samtools/flagstat/{prefix}.log",
    benchmark:
        "benchmarks/mapping/samtools/flagstat/{prefix}.log"
    resources:
        runtime="1h",
    group:
        "endocont"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """


rule qualimap:
    """
    Estimate general mapping statistics for each final sample bam file.
    """
    input:
        unpack(get_final_bam),
    output:
        fold=directory("results/mapping/qc/qualimap/{sample}.{ref}"),
        rep="results/mapping/qc/qualimap/{sample}.{ref}/qualimapReport.html",
        txt="results/mapping/qc/qualimap/{sample}.{ref}/genome_results.txt",
    container:
        qualimap_container
    log:
        "logs/mapping/qualimap/{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/qualimap/{sample}.{ref}.log"
    resources:
        runtime="6h",
    shell:
        """
        qualimap bamqc -bam {input.bam} --java-mem-size={resources.mem_mb}M \
            -outdir {output.fold} 2> {log}
        """


rule qualimap_userprovided:
    """
    Estimate general mapping statistics for each user-provided sample bam file.
    """
    input:
        unpack(get_final_bam),
    output:
        fold=directory(
            "results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}"
        ),
        rep="results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}/qualimapReport.html",
        txt="results/datasets/{dataset}/qc/user-provided-bams/qualimap/{sample}.{ref}/genome_results.txt",
    container:
        qualimap_container
    log:
        "logs/mapping/qualimap/{dataset}.{sample}.{ref}.log",
    benchmark:
        "benchmarks/mapping/qualimap/{dataset}.{sample}.{ref}.log"
    resources:
        runtime="6h",
    shell:
        """
        qualimap bamqc -bam {input.bam} --java-mem-size={resources.mem_mb}M \
            -outdir {output.fold} 2> {log}
        """


rule qualimap_multiqc:
    """
    Merge Qualimap outputs into a MultiQC report.
    """
    input:
        multiqc_input_qualimap,
    output:
        report(
            "results/datasets/{dataset}/qc/qualimap/qualimap_all.{ref}_mqc.html",
            category="00 Quality Control",
            subcategory="5 Qualimap",
            labels={"Type": "MultiQC Report"},
        ),
    log:
        "logs/mapping/qualimap/{dataset}.{ref}_mqc.log",
    container:
        multiqc_container
    params:
        extra="--cl-config \"extra_fn_clean_exts: ['.rmdup', '.clip']\" "
        '--cl-config "qualimap_config: { '
        'general_stats_coverage: [1,2,3,5,10,15] }"',
    resources:
        runtime="1h",
    shell:
        """
        multiqc {params.extra} --no-data-dir \
            --filename {output} {input} 2> {log}
        """


rule endo_cont:
    """
    Estimate the proportion of reads mapping to the reference as a proxy for
    endogenous content.
    """
    input:
        unpack(get_endo_cont_stat),
    output:
        endo="results/datasets/{dataset}/qc/endogenous_content/{dataset}.{sample}.{ref}.endo",
    container:
        shell_container
    log:
        "logs/datasets/{dataset}/qc/endogenous_content/{dataset}.{sample}.{ref}.log",
    benchmark:
        "benchmarks/datasets/{dataset}/qc/endogenous_content/{dataset}.{sample}.{ref}.log"
    resources:
        runtime="15m",
    group:
        "endocont"
    script:
        "../scripts/calc_endocont.sh"


rule compile_endo_cont:
    """
    Merge per sample endogenous content estimates into a single table.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/qc/endogenous_content/{{dataset}}.{sample}.{{ref}}.endo",
            sample=get_popfile_inds(w),
        ),
    output:
        "results/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_{population}{dp}.endo.tsv",
    wildcard_constraints:
        population="all",
    log:
        "logs/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_{population}{dp}_compile-endocont.log",
    benchmark:
        "benchmarks/datasets/{dataset}/qc/endogenous_content/{dataset}.{ref}_{population}{dp}_compile-endocont.log"
    container:
        shell_container
    resources:
        runtime="15m",
    shell:
        """
        (printf "sample\tperc.collapsed.map\tperc.uncollapsed.map\tperc.total.map\n" > {output}
        cat {input} >> {output}) 2> {log}
        """


rule ind_unfiltered_depth:
    """
    Estimate unfiltered sample depth, only removing reads using default ANGSD filters
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
    output:
        sample_hist="results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.depthSample",
        global_hist=temp(
            "results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.depthGlobal"
        ),
        arg="results/mapping/qc/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.arg",
    log:
        "logs/mapping/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.log",
    benchmark:
        "benchmarks/mapping/ind_depth/unfiltered/{dataset}.{ref}_{population}{dp}_allsites-unfilt.log"
    container:
        angsd_container
    params:
        out=lambda w, output: os.path.splitext(output.arg)[0],
        maxdepth=config["params"]["angsd"]["maxdepth"],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
    shell:
        """
        angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
            -bam {input.bamlist} -nThreads {threads} -out {params.out} &> {log}
        """


rule ind_mapQ_baseQ_depth:
    """
    Estimate sample depth with only mapping and base quality filters. Also
    includes universal 'extra' params user provides for ANGSD. This makes it
    compatible with using -trim
    """
    input:
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
    output:
        sample_hist="results/mapping/qc/ind_depth/mapq-baseq-filtered/{dataset}.{ref}_{population}{dp}_allsites-mapq{mapq}-baseq{baseq}-filt.depthSample",
        global_hist=temp(
            "results/mapping/qc/ind_depth/mapq-baseq-filtered/{dataset}.{ref}_{population}{dp}_allsites-mapq{mapq}-baseq{baseq}-filt.depthGlobal"
        ),
        arg="results/mapping/qc/ind_depth/mapq-baseq-filtered/{dataset}.{ref}_{population}{dp}_allsites-mapq{mapq}-baseq{baseq}-filt.arg",
    log:
        "logs/mapping/ind_depth/mapq-baseq-filtered/{dataset}.{ref}_{population}{dp}_allsites-mapq{mapq}-baseq{baseq}-filt.log",
    benchmark:
        "benchmarks/mapping/ind_depth/mapq-baseq-filtered/{dataset}.{ref}_{population}{dp}_allsites-mapq{mapq}-baseq{baseq}-filt.log"
    container:
        angsd_container
    params:
        out=lambda w, output: os.path.splitext(output.arg)[0],
        extra=config["params"]["angsd"]["extra"],
        maxdepth=config["params"]["angsd"]["maxdepth"],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
    shell:
        """
        angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
            -minMapQ {wildcards.mapq} -minQ {wildcards.baseq} {params.extra} \
            -bam {input.bamlist} -nThreads {threads} -out {params.out} &> {log}
        """


rule ind_filtered_depth:
    """
    Estimate depth at positions using filters for main workflow. This describes the
    depth of positions that will go into likelihood calculations
    """
    input:
        unpack(filt_depth),
        bamlist="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
    output:
        sample_hist="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.depthSample",
        global_hist=temp(
            "results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.depthGlobal"
        ),
        arg="results/datasets/{dataset}/qc/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.arg",
    log:
        "logs/{dataset}/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/ind_depth/filtered/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    params:
        extra=config["params"]["angsd"]["extra"],
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
        maxdepth=config["params"]["angsd"]["maxdepth"],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 60,
    shell:
        """
        angsd -doDepth 1 -doCounts 1 -maxDepth {params.maxdepth} \
            -bam {input.bamlist} -ref {input.ref} -nThreads {threads} \
            {params.extra} -minMapQ {params.mapQ} -minQ {params.baseQ} \
            -sites {input.sites} -out {params.out} &> {log}
        """


rule summarize_ind_depth:
    """
    Get mean and standard deviation of depth from individual depth distributions
    """
    input:
        sample_hist="{prefix}{dataset}.{ref}_{sample}{dp}_{group}.depthSample",
        filtersum=get_filter_sum,
    output:
        sample_summ="{prefix}{dataset}.{ref}_{sample}{dp}_{group}.depth.sum",
    log:
        "logs/summarize_ind_depth/{prefix}{dataset}.{ref}_{sample}{dp}_{group}.log",
    benchmark:
        "benchmarks/summarize_ind_depth/{prefix}{dataset}.{ref}_{sample}{dp}_{group}.log"
    container:
        r_container
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime="1h",
    script:
        "../scripts/calc_depth.R"


rule merge_ind_depth:
    """
    Combine depth summaries for all individuals
    """
    input:
        depth=lambda w: expand(
            "{{prefix}}{{dataset}}.{{ref}}_{sample}{{dp}}_{{group}}.depthSample",
            sample=get_popfile_inds(w),
        ),
        summary=lambda w: expand(
            "{{prefix}}{{dataset}}.{{ref}}_{sample}{{dp}}_{{group}}.depth.sum",
            sample=get_popfile_inds(w),
        ),
    output:
        dep="{prefix}{dataset}.{ref}_{population}{dp}_{group}.depth",
        sum="{prefix}{dataset}.{ref}_{population}{dp}_{group}.depth.sum",
    wildcard_constraints:
        population="all",
    log:
        "logs/merge_depth/{prefix}{dataset}.{ref}_{population}{dp}_{group}.log",
    benchmark:
        "benchmarks/merge_depth/{prefix}{dataset}.{ref}_{population}{dp}_{group}.log"
    container:
        shell_container
    resources:
        runtime="15m",
    shell:
        """
        (cat {input.depth} > {output.dep}
        printf "sample\t{wildcards.group}.depth.mean\t{wildcards.group}.depth.stdev\n" \
            > {output.sum}
        cat {input.summary} >> {output.sum}) 2> {log}
        """


rule combine_sample_qc:
    """
    Compile summary table of all sample QC measures
    """
    input:
        unpack(get_sample_qcs),
    output:
        "results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv",
    log:
        "logs/datasets/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}.log",
    benchmark:
        "benchmarks/datasets/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}.log"
    container:
        shell_container
    shadow:
        "minimal"
    resources:
        runtime="15m",
    shell:
        """
        (for i in {input}; do
            head -n 1 $i > ${{i}}.tmp
            tail -n +2 $i | sort -k1 >> ${{i}}.tmp
        done

        cut -d '\t' -f 1 {input.inds}.tmp > {output}

        for i in {input}; do
            cut -d '\t' -f 2- ${{i}}.tmp | paste {output} - > {output}.tmp
            mv {output}.tmp {output}
        done) 2> {log}
        """


rule sample_qc_summary:
    """
    Convert sample QC summary table to html for report
    """
    input:
        "results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.tsv",
    output:
        report(
            "results/datasets/{dataset}/qc/{dataset}.{ref}_all{dp}.sampleqc.html",
            category="00 Quality Control",
            subcategory="6 Sample depth and endogenous content",
            labels=lambda w: {**dp_report(w), "Type": "Table"},
        ),
    log:
        "logs/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/combine_sample_qc/{dataset}.{ref}{dp}_tsv2html.log"
    container:
        r_container
    shadow:
        "minimal"
    resources:
        runtime="15m",
    script:
        "../scripts/tsv2html.R"


rule ibs_ref_bias_nofilts:
    """
    Calculates the average identity by state to the reference for all sites,
    even those filtered by the sites filter. Useful to determine if batches
    of samples have reference biases. Compare to filtered version to see if
    site filters reduce this bias.
    """
    input:
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
    output:
        ibs=temp(
            "results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_allsites-unfilt.ibs.gz"
        ),
        arg="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_allsites-unfilt.arg",
        stats="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_allsites-unfilt.refibs.tsv",
    wildcard_constraints:
        population="|".join(samples.index),
    log:
        "logs/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_allsites-unfilt.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_allsites-unfilt.log"
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        mindepthind=config["params"]["angsd"]["mindepthind"],
        extra=config["params"]["angsd"]["extra"],
        trans=get_trans,
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        r"""
        (angsd -doIBS 1 -bam {input.bam} -ref {input.ref} -nThreads {threads} \
            -doMajorMinor 4 {params.extra} -GL {params.gl_model} \
            -minMapQ {params.mapQ} -setMinDepthInd {params.mindepthind} \
            -minQ {params.baseQ} -doCounts 1 -output01 1 \
            -rmTrans {params.trans} -out {params.out}
        
        ibs=$(zcat {output.ibs} | tail -n+2 | \
            awk '{{ sum += $5 }} END {{ if (NR > 0) print 1 - (sum / NR) }}')
        printf "{wildcards.population}\t$ibs\n" > {output.stats}) 2> {log}
        """


rule ibs_ref_bias_filts:
    """
    Calculates the average identity by state to the reference for all sites,
    even those filtered by the sites filter. Useful to determine if batches
    of samples have reference biases. Compare to filtered version to see if
    site filters reduce this bias.
    """
    input:
        unpack(filt_depth),
        bam="results/datasets/{dataset}/bamlists/{dataset}.{ref}_{population}{dp}.bamlist",
        bams=get_bamlist_bams,
        bais=get_bamlist_bais,
        ref="results/ref/{ref}/{ref}.fa",
        reffai="results/ref/{ref}/{ref}.fa.fai",
    output:
        ibs=temp(
            "results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{sites}-filts.ibs.gz"
        ),
        arg="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{sites}-filts.arg",
        stats="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{sites}-filts.refibs.tsv",
    wildcard_constraints:
        population="|".join(samples.index),
    log:
        "logs/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_{sites}-filts.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_{sites}-filts.log"
    container:
        angsd_container
    params:
        gl_model=config["params"]["angsd"]["gl_model"],
        mindepthind=config["params"]["angsd"]["mindepthind"],
        extra=config["params"]["angsd"]["extra"],
        trans=get_trans,
        mapQ=config["mapQ"],
        baseQ=config["baseQ"],
        out=lambda w, output: os.path.splitext(output.arg)[0],
    threads: lambda wildcards, attempt: attempt
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
    shell:
        r"""
        (angsd -doIBS 1 -bam {input.bam} -ref {input.ref} -nThreads {threads} \
            -doMajorMinor 4 {params.extra} -GL {params.gl_model} \
            -minMapQ {params.mapQ} -setMinDepthInd {params.mindepthind} \
            -minQ {params.baseQ} -doCounts 1 -output01 1 \
            -rmTrans {params.trans} -sites {input.sites} -out {params.out}
        
        ibs=$(zcat {output.ibs} | tail -n+2 | \
            awk '{{ sum += $5 }} END {{ if (NR > 0) print 1 - (sum / NR) }}')
        printf "{wildcards.population}\t$ibs\n" > {output.stats}) 2> {log}
        """


rule merge_ibs_ref_bias:
    """
    Combines all sample IBS ref bias results into a single table.
    """
    input:
        lambda w: expand(
            "results/datasets/{{dataset}}/qc/ibs_refbias/{{dataset}}.{{ref}}_{population}{{dp}}_{{filts}}.refibs.tsv",
            population=get_popfile_inds(w),
        ),
    output:
        ibs=temp(
            "results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{filts}.refibs.merged.tsv"
        ),
    log:
        "logs/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_{filts}_merge.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_{filts}_merge.log"
    container:
        shell_container
    resources:
        runtime="15m",
    group:
        "combine_ibsrefbias"
    shell:
        r"""
        printf "sample\tibs.to.ref\n" > {output}
        cat {input} >> {output}
        """


rule plot_ibs_ref_bias:
    """
    Create a boxplot comparing groups for reference bias. Makes a plot for all
    available groupings in the sample list: population, time, and depth. Does
    not delve into user defined columns.
    """
    input:
        ibs="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_all{dp}_{filts}.refibs.merged.tsv",
        pops="results/datasets/{dataset}/poplists/{dataset}_all{dp}.indiv.list",
    output:
        pop_plot=report(
            "results/datasets/{dataset}/plots/ibs_refbias/{dataset}.{ref}_all{dp}_{filts}.population.pdf",
            category="00 Quality Control",
            subcategory="7 Sample IBS to reference (ref bias)",
            labels=lambda w: {
                **dp_report(w),
                "Filter": "{filts}",
                "Type": "Boxplot",
                "Grouping": "Population",
            },
        ),
        tim_plot=report(
            "results/datasets/{dataset}/plots/ibs_refbias/{dataset}.{ref}_all{dp}_{filts}.time.pdf",
            category="00 Quality Control",
            subcategory="7 Sample IBS to reference (ref bias)",
            labels=lambda w: {
                **dp_report(w),
                "Filter": "{filts}",
                "Type": "Boxplot",
                "Grouping": "Time",
            },
        ),
        dep_plot=report(
            "results/datasets/{dataset}/plots/ibs_refbias/{dataset}.{ref}_all{dp}_{filts}.depth.pdf",
            category="00 Quality Control",
            subcategory="7 Sample IBS to reference (ref bias)",
            labels=lambda w: {
                **dp_report(w),
                "Filter": "{filts}",
                "Type": "Boxplot",
                "Grouping": "Depth",
            },
        ),
        table="results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_all{dp}_{filts}.refibs.tsv",
    log:
        "logs/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_all{dp}_{filts}_plot.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_all{dp}_{filts}_plot.log"
    params:
        plotpre=lambda w, output: output["pop_plot"].removesuffix(".population.pdf"),
    container:
        r_container
    resources:
        runtime="30m",
    group:
        "combine_ibsrefbias"
    script:
        "../scripts/plot_ref_bias.R"


rule ibs_ref_bias_table_html:
    """
    Convert sample QC summary table to html for report
    """
    input:
        "results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{filts}.refibs.tsv",
    output:
        report(
            "results/datasets/{dataset}/qc/ibs_refbias/{dataset}.{ref}_{population}{dp}_{filts}.refibs.html",
            category="00 Quality Control",
            subcategory="7 Sample IBS to reference (ref bias)",
            labels=lambda w: {**dp_report(w), "Filter": "{filts}", "Type": "Table"},
        ),
    log:
        "logs/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_{filts}_tsv2html.log",
    benchmark:
        "benchmarks/{dataset}/angsd/ibs_ref_bias/{dataset}.{ref}_{population}{dp}_{filts}_tsv2html.log"
    container:
        r_container
    shadow:
        "minimal"
    resources:
        runtime="15m",
    group:
        "combine_ibsrefbias"
    script:
        "../scripts/tsv2html.R"
