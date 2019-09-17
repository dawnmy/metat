from os import path

# Load paths for input and output from config
configfile: "config/config.yaml"

dataset = config["dataset"]
fq_dir = config["fq_dir"]
proj_dir = path.join(config["out_dir"], dataset)
ref = config["ref"]
paired = config["paired"]
suffix = config["suffix"]
threads = config["threads"]
cd = path.dirname(path.abspath(__name__))

out_data_dir = path.join(proj_dir, "data")
qc_dir = path.join(out_data_dir, "qc_fq")
reports_dir = path.join(proj_dir, "reports")
results_dir = path.join(proj_dir, "results")

wildcard_constraints:
    sample = "[^\.\/]+"

if paired:
    samples, = glob_wildcards(fq_dir + "/{sample, [^\.\/]+}"+suffix[0])
else:
    raise RuntimeError("The single end data is not supported so far!")

def get_pe_fq(wc):
    return [path.join(fq_dir, wc.sample + end) for end in suffix]


# def get_se_fq(wc):
#     return path.join(fq_dir, wc.sample + ".fq.gz")

def sortme_db(cd):
    rrna, = glob_wildcards(cd + "/data/rRNA/fasta/{rrna}.fasta")
    dbs = []
    for i in rrna:
        fa = path.join(cd, "data/rRNA/fasta/{}.fasta".format(i))
        idx = path.join(cd, "data/rRNA/idx/{}".format(i))
        dbs.append(",".join([fa, idx]))
    return ":".join(dbs)

sortme_dbs = sortme_db(cd)

rule all:
    input:
        count = expand(
            results_dir + "/count/{dataset}.count", dataset=dataset)

rule build_idx:
    input:
        ref = ref
    output:
        fa_idx = ref + ".fai",
        bwa_idx = ref + ".bwt"
    conda:
        "config/conda.metat.yaml"
    shell:
        """
        samtools faidx {input.ref}
        bwa index {input.ref}
        """

rule sortme_idx:
    input:
        ref = ref
    output:
        touch(cd + "/data/rRNA/idx.done")
    params:
        bin_dir = path.join(cd, "libs/sortmerna"),
        dbs = sortme_dbs
    shell:
        """
        {params.bin_dir}/indexdb --ref {params.dbs}
        """

rule fastp:
        input:
            get_pe_fq
        output:
            or1 = qc_dir + "/{sample}.qc_1.fq",
            or2 = qc_dir + "/{sample}.qc_2.fq",
            html = reports_dir + "/fastp/{sample}.qc.report.html"
        threads: 8
        conda:
            "config/conda.metat.yaml"
        shell:
            """
            fastp --in1 {input[0]} --in2 {input[1]} -o {output.or1} -O {output.or2} \
                -5 20 -3 20 -l 30 -h {output.html} -w {threads}
            """


rule sortmerna:
    input:
        r1 = rules.fastp.output.or1,
        r2 = rules.fastp.output.or2,
        idx_done = cd + "/data/rRNA/idx.done"
    output:
        merged_in = qc_dir + "/{sample}.merged.fq",
        merged_out = qc_dir + "/mrna/{sample}.merged.fq",
        rrna = qc_dir + "/rrna/{sample}.merged.fq",
        mr1 = qc_dir + "/mrna/{sample}.mrna_1.fq",
        mr2 = qc_dir + "/mrna/{sample}.mrna_2.fq"
    params:
        bin_dir = path.join(cd, "libs/sortmerna"),
        dbs = sortme_dbs
    log:
        reports_dir + "/sortmerna/{sample}.derrna.log"
    benchmark:
        reports_dir + "/benchmarks/{sample}.sortmerna.txt"
    threads: threads
    shell:
        """
        {params.bin_dir}/merge-paired-reads.sh {input.r1} {input.r2} {output.merged_in}
        {params.bin_dir}/sortmerna --ref {params.dbs} \
            --reads {output.merged_in} --fastx --aligned {output.rrna} \
            --other {output.merged_out} --paired_out -a {threads} --log -v >> {log}
        {params.bin_dir}/unmerge-paired-reads.sh {output.merged_out} {output.mr1} {output.mr2}
        """

# The host contamination sequences removal if needed
# rule host_removal:


rule bwa:
    input:
        r1 = rules.sortmerna.output.mr1,
        r2 = rules.sortmerna.output.mr2,
        ref = ref,
        ref_idx = rules.build_idx.output.bwa_idx
    output:
        out_data_dir + "/bam/{sample}.bam"
    conda:
        "config/conda.metat.yaml"
    benchmark:
        reports_dir + "/benchmarks/{sample}.bwa.txt"
    log:
        reports_dir + "/bwa/{sample}.log"
    threads: threads
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} |\
            samtools view -@ {threads} -Shb - 2>> {log} |\
            samtools sort -@ {threads} -m 20G - -o {output} >> {log} 2>&1
        """

rule bam_filter:
    input:
        rules.bwa.output
    output:
        out_data_dir + "/bam/{sample}.filtered.bam"
    conda:
        "config/conda.metat.yaml"
    log:
        reports_dir + "/samtools/{sample}.log"
    threads: threads
    shell:
        """
        samtools view -@ {threads} -Shb -F 4 -q 10 {input} 2>> {log} |\
            samtools sort -@ {threads} -m 20G - -o {output} >> {log} 2>&1
        samtools index {output}
        """

rule samcount:
    input:
        bam = rules.bam_filter.output
    output:
        results_dir + "/count/{sample}.samcount"
    conda:
        "config/conda.metat.yaml"
    threads: threads
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

rule cat_count:
    input:
        expand(
            results_dir + "/count/{sample}.samcount", sample=samples)
    output:
        results_dir + "/count/{dataset}.count"
    conda:
        "config/conda.metat.yaml"
    params:
        cols = ",".join(["1"] + [str(3 * i)
                                 for i in range(1, len(samples) + 1)])
    threads: threads
    shell:
        """
        cat <(printf "gene\t";sed 's/[^ ]\+\/\|\.samcount//g;s/ \+/\t/g' <(echo {input})) \
            <(csvtk join -TtHf 1 {input}|csvtk cut -TtHf {params.cols}) > {output}
        """
