import os

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "%s/config.yaml" % SNAKEMAKE_DIR

if not os.path.exists("log"):
    os.makedirs("log")

localrules: all, clean

rule filter_sun_repeats:
    input: "suns.bed", config["repeats"]
    output: "suns.no_repeats_36bp_flanking.bed"
    params: sge_opts="-l h_rt=1:0:0 -l mfree=4G"
    shell:
        "bedtools window -a {input[0]} -b {input[1]} -w 36 -v > {output}"

rule get_suns:
    input: "merged/suns.sort.bed", "sunks.bed"
    output: "suns.bed"
    params: sge_opts="-l h_rt=1:0:0 -l mfree=4G"
    shell:
        "intersectBed -a {input[0]} -b {input[1]} -sorted > {output}"

rule sort_mismatches:
    input: dynamic("suns/suns.{part}.bed")
    output: "merged/suns.sort.bed"
    params: sge_opts="-l h_rt=1:0:0"
    shell:
        "sort -m -k 1,1 -k 2,2n {input} > {output}"

rule get_mismatches:
    input: wgac="split/unique_wgac{part}",
           fasta="regions.fasta"
    output: "suns/suns.{part}.bed"
    shadow: True
    params: sge_opts="-l h_rt=5:0:0 -soft -l ssd=True"
    shell:
        """python find_suns.py {input.wgac} {input.fasta} . --local
        sort -k 1,1 -k 2,2n suns.bed > {output}"""


rule split_wgac:
    input: "unique_wgac.bed"
    output: dynamic("split/unique_wgac{part}")
    params: sge_opts="-l h_rt=1:0:0", ofprefix="split/unique_wgac", size=config["size"]
    shell:
        "split -l {params.size} -d -a 4 {input} {params.ofprefix}"

rule regions_fasta:
    input: bed="regions.bed", ref=config["reference"]
    output: "regions.fasta"
    params: sge_opts="-l h_rt=1:0:0"
    shell:
        "bedtools getfasta -fi {input.ref} -bed {input.bed} -fo {output} -s"

rule regions_bed:
    input: "unique_wgac.bed"
    output: "regions.bed"
    params: sge_opts="-l h_rt=1:0:0"
    shell:
        """awk 'OFS="\\t" {{ print $1,$2,$3,$1"_"$2"_"$3,0,$4; print $5,$6,$7,$5"_"$6"_"$7,0,$8 }}' {input} | \
        sort -k 1,1 -k 2,2n | uniq > {output}"""

rule unique_wgac:
    input: "wgac.bed"
    output: "unique_wgac.bed"
    params: sge_opts="-l h_rt=1:0:0"
    shell:
        "python get_unique_wgac.py {input} > {output}"

rule wgac_bed:
    input: config["wgac"]
    output: "wgac.bed"
    params: sge_opts="-l h_rt=1:0:0"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            for line in infile:
                chr1, s1, e1, name, score, strand1, chr2, s2, e2 = line.rstrip().split()[0:9]
                if strand1 in ["_", "-"]:
                    strand2 = "-"
                else:
                    strand2 = "+"
                strand1 = "+"
                print(chr1, s1, e1, strand1, chr2, s2, e2, strand2, sep="\t", file=outfile)

rule sunks_bed:
    input: config["sunks"]
    output: "sunks.bed"
    params: sge_opts="-l h_rt=1:0:0"
    shell:
        "zcat {input} | mergeBed -i stdin -d 0 > {output}"

rule clean:
    shell:
        "rm -f *.bed *.fasta *.tab"
