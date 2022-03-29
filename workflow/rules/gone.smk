rule gone:
	input:
		ped=rules.tped2ped.output.ped,
		map=rules.tped2ped.output.map
	output:
		output=results + "/GONE/OUTPUT_{population}_cMMb{rr}_hc{hc}",
		Ne=results + "/GONE/Output_Ne_{population}_cMMb{rr}_hc{hc}",
		d2=results + "/GONE/Output_d2_{population}_cMMb{rr}_hc{hc}"
	log:
		logs + "/GONE/{population}_cMMb{rr}_hc{hc}.log"
	params:
		gonedir=results + "/GONE/"
	threads: 20
	resources:
		time="12:00:00"
	shell:
		"""
		date
		mkdir -p {params.gonedir}

		mkdir -p {resources.tmpdir}/{rule}.{wildcards}/
		cp -r bin/GONE/* {resources.tmpdir}/{rule}.{wildcards}/
		cp {input.ped} {resources.tmpdir}/{rule}.{wildcards}/
		cp {input.map} {resources.tmpdir}/{rule}.{wildcards}/
		
		workdir=$(pwd -P)
		cd {resources.tmpdir}/{rule}.{wildcards}/
		ls
		sed -i -e 's/cMMb=1/cMMb={wildcards.rr}/g' \
			-e 's/hc=0.5/hc={wildcards.hc}/g' INPUT_PARAMETERS_FILE

		bash script_GONE.sh {wildcards.population}_genome

		cp OUTPUT_{wildcards.population}_genome ${{workdir}}/{output.output}
		cp Output_Ne_{wildcards.population}_genome ${{workdir}}/{output.Ne}
		cp Output_d2_{wildcards.population}_genome ${{workdir}}/{output.d2}
		date
		"""