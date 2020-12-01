import os

from shutil import copyfile
from snakemake.utils import min_version

min_version("5.4.3")

# Get parent directory of snakefile
SNAKEDIR = os.path.dirname(workflow.snakefile)
# Get parent directory of config file
CONFDIR = SNAKEDIR
if workflow.configfiles:
    conf_dir = workflow.configfiles[0]
    if os.path.isabs(conf_dir):
        CONFDIR = path.dirname(conf_dir)

WORKDIR = "."
if config.get("workdir_top"):
    if os.path.isabs(config["workdir_top"]):
        WORKDIR = os.path.join(config.get("workdir_top", ""))
    else:
        WORKDIR = os.path.join(CONFDIR, config["workdir_top"])

workdir: WORKDIR

print("Working directory: {}".format(WORKDIR))

pipeline_version="2.0.2"

#########################
###### PARAMETERS #######
#########################


# INPUT BAM folder
bam = None
if config.get("bam"):
    bam = os.path.join(CONFDIR, config["bam"])

# INPUT FASTQ folder
FQ_INPUT_DIRECTORY = []
if not bam:
    if not "input_fastq" in config:
        print("\"input_fastq\" not specified in config file. Exiting...")

    FQ_INPUT_DIRECTORY = os.path.join(CONFDIR, config["input_fastq"])
    if not os.path.exists(FQ_INPUT_DIRECTORY):
        print("Could not find {}".format(FQ_INPUT_DIRECTORY))

    MAPPED_BAM = "{sample}/alignment/{sample}_lra.bam"
else:
    MAPPED_BAM = bam
    if not os.path.exists(bam):
        print("Could not find {}".format(bam))


# Input reference FASTA
FA_REF = os.path.join(CONFDIR, config["reference_fasta"])
if not os.path.exists(FA_REF):
    print("Could not find {}".format(FA_REF))

# Reference index name
FA_REF_INDEX = FA_REF + ".gli"

# Parameter: sample_name
sample = config.get('sample_name', "sv_sample01")

# Parameter: target_bed
target_bed = ""
if config.get("target_bed"):
    target = config["target_bed"]
    if not os.path.isabs(target):
        target = os.path.join(CONFDIR, target)
    if os.path.exists(target):
        target_bed = target
        print("Using {} as target file".format(target_bed))
    else:
        print("Target BED {} not found. Continuing without target".format(target))

thread_n = config.get("threads", 30)

#########################
######## RULES ##########
#########################

rule call:
    input:
        expand("{name}/sv_calls/{name}_cutesv_filtered.vcf.gz", name=sample),
        expand("{name}/version.txt", name=sample)

rule qc:
    input:
        expand("{name}/qc", name=sample),
        expand("{name}/version.txt", name=sample)

rule eval:
    input:
        expand("{name}/evaluation_summary.json", name=sample),
        expand("{name}/version.txt", name=sample)


rule print_version:
    output:
        "{name}/version.txt"
    params:
        version = pipeline_version
    shell:
        "echo {params.version} > {output}"


rule index_lra:
   input:
       REF = FA_REF
   output:
       INDEX = FA_REF_INDEX
   conda: "env.yml"
   threads: thread_n
   shell:
       "lra index -ONT {input}"

rule map_lra:
   input:
       FQ = FQ_INPUT_DIRECTORY,
       REF = FA_REF,
       INDEX = FA_REF_INDEX,
   output:
       BAM = "{sample}/alignment/{sample}_lra.bam",
       BAI = "{sample}/alignment/{sample}_lra.bam.bai"
   conda: "env.yml"
   threads: thread_n
   benchmark: "{sample}/benchmarks/map_lra_{sample}.time"
   shell:
       "catfishq -r {input.FQ} | seqtk seq -A - | lra align -ONT -t {threads} {input.REF} - -p s | samtools addreplacerg -r \"@RG\tID:{sample}\tSM:{sample}\" - | samtools sort -@ {threads} -T {sample} -O BAM -o {output.BAM} - && samtools index -@ {threads} {output.BAM}"


rule call_cutesv:
    input:
        BAM = MAPPED_BAM,
        REF = FA_REF,
    output:
        VCF = "{sample}/sv_calls/{sample}_cutesv_tmp.vcf"
    params:
        min_size = config.get("min_sv_length", 30),
        max_size = config.get("max_sv_length", 100000),
        min_read_support = 2,
        min_read_length = config.get("min_read_length", 1000),
        min_mq = config.get("min_read_mapping_quality", 20),
    conda: "env.yml"
    threads: thread_n
    benchmark: "{sample}/benchmarks/call_cutesv_{sample}.time"
    shell:
        "cuteSV -t {threads} --min_size {params.min_size} --max_size  {params.max_size} -S {sample} --retain_work_dir --report_readid --min_support {params.min_read_support} --genotype {input.BAM} {input.REF} {output.VCF} {sample}/sv_calls/ "


rule filter_vcf:
    input:
         MOS = "{sample}/depth",
         VCF = rules.call_cutesv.output.VCF,
    output:
         VCF = temp("{sample}/sv_calls/{sample}_cutesv_filtered_tmp.vcf")
    params:
        min_sv_length = config.get("min_sv_length", 30),
        max_sv_length = config.get("max_sv_length", 100000),
        target_bed = config.get("target_bed", None),
        sv_types = config.get("sv_type", "DEL INS"),
    conda: "env.yml"
    wrapper:
         f"file:{CONFDIR}/wrappers/filter"

rule sort_vcf:
    input:
        VCF = rules.filter_vcf.output.VCF
    output:
        VCF = temp("{sample}/sv_calls/{sample}_cutesv_filtered.vcf")
    conda: "env.yml"
    shell:
         "vcfsort {input.VCF} > {output.VCF}"


rule index_vcf:
    input:
         VCF = rules.sort_vcf.output.VCF
    output:
         VCF = "{sample}/sv_calls/{sample}_cutesv_filtered.vcf.gz"
    conda: "env.yml"
    shell:
         "cat {input.VCF} | bgziptabix {output.VCF}"


rule nanoplot_qc:
    input:
        BAM = MAPPED_BAM
    output:
        DIR = directory("{sample}/qc")
    params:
        sample = sample
    conda: "env.yml"
    threads: thread_n
    shell:
        "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.sample}_ --N50 --title {params.sample} --downsample 100000"


rule calc_depth:
    input:
         BAM = MAPPED_BAM,
    output:
        DIR = directory("{sample}/depth"),
    params:
        BED = config.get("target", "1000000")
    conda: "env.yml"
    threads: thread_n
    shell:
        "mkdir -p {output.DIR}; mosdepth -x -t {threads} -n -b {params.BED} {output.DIR}/{sample} {input.BAM}"


rule download_hg002_truthset:
    output:
        VCF = "HG002_SVs_Tier1_v0.6.vcf.gz",
        TBI = "HG002_SVs_Tier1_v0.6.vcf.gz.tbi",
        BED = "HG002_SVs_Tier1_v0.6.bed"
    conda: "env.yml"
    shell:
        "wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz && wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi && wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed"


rule intersect_target_highconf:
    input:
        TRUTH_BED = rules.download_hg002_truthset.output.BED,
        VCF = rules.index_vcf.output.VCF
    output:
        BED = "{sample}/eval_high_conf.bed"
    conda: "env.yml"
    shell:
        "bedtools intersect -a {input.TRUTH_BED} -b {input.VCF} -u > {output.BED}"

rule eval_reformat:
    input:
        VCF = rules.index_vcf.output.VCF,
    output:
        VCF = "{sample}/sv_calls/{sample}_cutesv_filtered_eval.vcf.gz",
    conda: "env.yml"
    shell:
        "zcat {input.VCF} | sed 's/SVTYPE=DUP/SVTYPE=INS/g' |  bcftools view -i '(SVTYPE = \"INS\" || SVTYPE = \"DEL\")' | bgziptabix {output.VCF}"


rule eval_vcf:
    input:
        VCF = rules.eval_reformat.output.VCF,
        REF = FA_REF,
        TRUTH_VCF = rules.download_hg002_truthset.output.VCF,
        TRUTH_BED = rules.intersect_target_highconf.output.BED
    output:
        dir = directory("{sample}/evaluation/"),
        summary = "{sample}/evaluation_summary.json"
    conda: "env.yml"
    shell:
        """
        truvari bench --passonly -b {input.TRUTH_VCF} --includebed {input.TRUTH_BED} --pctsim 0 -c {input.VCF} -f {input.REF} -o {sample}/eval/ -o {output.dir}
        cp {output.dir}/summary.txt {output.summary}
        """
