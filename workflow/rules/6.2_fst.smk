rule realSFS_fst_index:
	input:
		saf1=rules.realSFS_2dSFS.input.saf1,
		saf2=rules.realSFS_2dSFS.input.saf2,
		sfs=rules.realSFS_2dSFS.output.sfs
	output:
		fstidx=results+"/analyses/fst/"+dataset+
			"_{population1}-{population2}{dp}.fst.idx"
	container:
		angsd_container
	log:
		logs+"/realSFS/fst/index/"+dataset+
			"_{population1}-{population2}{dp}.log"
	params:
		out=results+"/analyses/fst/"+dataset+"_{population1}-{population2}{dp}"
	resources:
		time=lambda wildcards, attempt: attempt*120
	shell:
		"""
		realSFS fst index {input.saf1} {input.saf2} -sfs {input.sfs} \
			-fstout {params.out} 2> {log}
		"""

rule realSFS_fst_stats:
	input:
		fstidx=rules.realSFS_fst_index.output.fstidx
	output:
		fstglob=results+"/analyses/fst/"+dataset+
			"_{population1}-{population2}{dp}.fst.global"
	container:
		angsd_container
	log:
		logs+"/realSFS/fst/stats/"+dataset+
			"_{population1}-{population2}{dp}.log"
	resources:
		time=lambda wildcards, attempt: attempt*60
	shell:
		r"""
		realSFS fst stats {input.fstidx} | \
			awk '{{print "{wildcards.population1}\t{wildcards.population2}\t"\
			$1"\t"$2}}' > {output.fstglob} 2> {log}
		"""

rule realSFS_fst_stats2:
	input:
		fstidx=rules.realSFS_fst_index.output.fstidx
	output:
		fstwin=results+"/analyses/fst/"+dataset+
			"_{population1}-{population2}{dp}.fst.window"
	container:
		angsd_container
	log:
		logs+"/realSFS/fst/stats2/"+dataset+
			"_{population1}-{population2}{dp}.log"
	shell:
		"""
		realSFS fst stats2 {input.fstidx} -win 50000 -step 10000 \
			> {output.fstwin} 2> {log}
		"""

def get_fst(wildcards):
	if wildcards.unit == "ind":
		unit = samples.index
	elif wildcards.unit == "pop":
		unit = pop_list
	combos = list(itertools.combinations(unit, 2))
	# sort pops alphebetically, this ensures that should new pops be added
	# after generating some SFS, the reordering of the combinations won't
	# lead to generating identical SFS with the populations swapped
	combos = [sorted(pair) for pair in combos]
	pop1 = [pair[0] for pair in combos]
	pop2 = [pair[1] for pair in combos]
	return expand(results+"/analyses/fst/"+dataset+
				"_{population1}-{population2}"+wildcards.dp+".fst.global",
				zip, population1=pop1, population2=pop2)

rule aggregate_global_fst:
	input:
		get_fst
	output:
		results+"/analyses/fst/"+dataset+"_{unit}pairs{dp}.fst.sum"
	log:
		logs+"/realSFS/fst/aggregate/"+dataset+"_{unit}pairs{dp}.log"
	wildcard_constraints:
		unit="ind|pop"
	shell:
		"""
		echo "pop1\tpop2\tunweight.fst\tweight.fst" > {output} 2> {log}
		cat {input} >> {output} 2>> {log}
		"""

localrules: plot_fst

rule plot_fst:
	input:
		results+"/analyses/fst/"+dataset+"_{unit}pairs{dp}.fst.sum"
	output:
		report(results+"/plots/fst/"+dataset+"_{unit}pairs{dp}_fst.pdf",
			category="Fst")
	conda:
		"../envs/r.yaml"
	script:
		"../scripts/plot_fst.R"