localrules: get_genome

rule get_genome:
    output:
        "resources/reference/"+os.path.basename(config['reference'])
    log:
        "logs/reference/get_genome.log"
    params:
        fasta_url=config['reference']
    shell:
        """
        mkdir -p resources/reference
        cd resources/reference

        wget {params.fasta_url}
        """

rule gunzip_genome:
    input:
        genome_file()+".gz"
    output:
        genome_file()
    shell:
        "gunzip {input}"

rule bwa_index:
    input:
        genome_file()
    output:
        genome_file()+".amb",
        genome_file()+".ann",
        genome_file()+".bwt",
        genome_file()+".pac",
        genome_file()+".sa"
    log:
        "logs/bwa_index/genome.log"
    resources:
        time="02:00:00"
    wrapper:
        "0.84.0/bio/bwa/index"