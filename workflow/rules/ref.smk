localrules: get_genome, chromosome_list

rule get_genome:
    output:
        "resources/reference/" + os.path.basename(config['reference']['fasta'])
    log:
        "logs/get_genome/get_genome.log"
    params:
        fasta_url=config['reference']
    shell:
        """
        wget -P resources/reference -o {log} {params.fasta_url}
        """

rule gunzip_genome:
    input:
        genome_file() + ".gz"
    output:
        protected(genome_file())
    shell:
        "gunzip {input}"

rule bwa_index:
    input:
        genome_file()
    output:
        protected(multiext(genome_file(),".amb",".ann",".bwt",".pac",".sa"))
    log:
        "logs/bwa_index/genome.log"
    resources:
        time="02:00:00"
    wrapper:
        "0.84.0/bio/bwa/index"

rule samtools_faidx:
    input:
        genome_file()
    output:
        protected(genome_file() + ".fai")
    wrapper:
        "v1.0.0/bio/samtools/faidx"

checkpoint chromosome_list:
    input:
        genome_file()
    output:
        protected(genome_file() + ".chromlist")
    shell:
        """
        grep ">" {input} | awk '{{print $1}}' | cut -c2- > {output}
        """