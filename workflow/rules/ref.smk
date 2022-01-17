localrules: get_genome

rule get_genome:
    output:
        genome_file()
    log:
        "logs/reference/get_genome.log"
    run:
        if pd.isna(config['reference']['fasta_path']):
            # Download reference to resources folder
            import urllib.request
            fasta_url = config['reference']['fasta_url']
            urllib.request.urlretrieve(fasta_url, genome_file())
        else:
            # Link specified genome file to expected location
            src = os.path.abspath(config['reference']['fasta_path'])
            dst = os.path.abspath(genome_file())
            os.symlink(src, dst)

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