# Per individual estimates of endogenous content - primarily useful for 
# historical samples or when using references that are not specific to 
# the species of interest

def get_endo_cont_stat(wildcards):
    # Gets bam file for endogenous content calculation
    s = wildcards.sample
    if s in samples.index[samples.time == "modern"]:
        return "results/mapping/mapped/"+s+".paired.flagstat"
    elif s in samples.index[samples.time == "historical"]:
        return "results/mapping/mapped/"+s+".merged.flagstat"

rule endo_cont:
    input:
        get_endo_cont_stat
    output:
        "results/mapping/qc/endogenous_content/{sample}.endo"
    shell:
        r"""
        total=$(grep -E "^[0-9]+ \+ [0-9]+ in total" {input} \
            | awk '{{print $1}}')
        mapped=$(grep -E "^[0-9]+ \+ [0-9]+ mapped" {input} \
            | awk '{{print $1}}')
        primary=$(grep -E "^[0-9]+ \+ [0-9]+ primary mapped" {input} \
            | awk '{{print $1}}')
        
        echo $total $mapped $primary {wildcards.sample} | \
            awk '{{printf "%s\t%.3f\t%.3f\n",$4,$2/$1*100,$3/$1*100}}' \
            > {output}
        """

rule compile_endo_cont:
	input:
		lambda w: expand("results/mapping/qc/endogenous_content/{sample}.endo",
			sample=get_samples_from_pop("all"))
	output:
		results+"/qc/endogenous_content/"+dataset+"_all.endo.tsv"
	resources:
		time=lambda wildcards, attempt: attempt*15
	shell:
		"""
		echo "sample	perc.endo	perc.prim.endo" > {output}
		cat {input} >> {output}
		"""

# Pairwise individual relatedness with R0, R1, KING-robust kinship 
# method from Waples et al. 2019, MolEcol

rule est_kinship_stats:
	input:
		sfs=results + "/analyses/sfs/"+dataset+
			"_{ind1}-{ind2}{dp}.sfs"
	output:
		results+"/analyses/kinship/"+dataset+
			"_{ind1}-{ind2}{dp}.kinship"
	log:
		logs + "/kinship/"+dataset+"_{ind1}-{ind2}{dp}_kinship.log"
	conda:
		"../envs/r.yaml"
	resources:
		time=lambda wildcards, attempt: attempt*15
	script:
		"../scripts/kinship.R"

def get_kinship(wildcards):
	combos = list(itertools.combinations(samples.index, 2))
	# sort inds alphebetically, this ensures that should new inds be added
	# after generating some SFS, the reordering of the combinations won't
	# lead to generating identical SFS with the individuals swapped
	combos = [sorted(pair) for pair in combos]
	ind1 = [pair[0] for pair in combos]
	ind2 = [pair[1] for pair in combos]
	return expand(results+"/analyses/kinship/"+dataset+
				"_{ind1}-{ind2}"+wildcards.dp+".kinship",
				zip, ind1=ind1, ind2=ind2)

rule compile_kinship_stats:
	input:
		get_kinship
	output:
		results+"/analyses/kinship/"+dataset+"_all{dp}.kinship"
	resources:
		time=lambda wildcards, attempt: attempt*15
	shell:
		"""
		echo "ind1	ind2	R0	R1	KING" > {output}
		cat {input} >> {output}
		"""

		echo "nsites nind" >> {log}

		echo $nsites $nind >> {log}

		cut -f1 {input.inds} | tail -n +2 > {output.samples} 2>> {log}

		ngsrelate -G {input.beagle} -n $nind -L $nsites -O {output.relate} \
			-z {output.samples} 2>> {log}
		"""
	