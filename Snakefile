import os
import glob
import gzip

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

configfile: os.path.join(CONFDIR, "config.yml")

config['threads'] = 120

if os.path.isabs(config["workdir_top"]):
    WORKDIR = os.path.join(config["workdir_top"], config["pipeline"])
else:
    WORKDIR = os.path.join(CONFDIR, config["workdir_top"], config["pipeline"])

if os.path.isabs(config["resdir"]):
    RESDIR =  config["resdir"]
else:
    RESDIR = os.path.join(CONFDIR, config["resdir"])


workdir: os.path.join(config["workdir_top"], config["pipeline"])
include: "snakelib/utils.snake"

print("Working directory: {}".format(WORKDIR))

#########################
### HELPER FUNCTIONS ####
#########################

def find_file_in_folder(folder, pattern="*.fastq", single=False):
    if os.path.isfile(folder):
        return folder
    files = []
    for file in glob.glob(os.path.join(folder, pattern)):
       files.append(file)

    if len(files) == 0:
        print("Could not find {} files in {}".format(pattern, folder))

    if single:
        if len(files) > 1:
            print("Warning: Multiple {} files found in {}".format(pattern, folder))
        return files[0]
    else:
        return files

def get_param_from_file(filename, default=0):
    try:
        with open(filename, "r") as fh:
            return int(float(fh.readline().strip()))
    except Exception as e:
        print(e)
        print("Could not read parameter from {}. Falling back to default: {}".format(filename, default))
        return default


#########################
###### PARAMETERS #######
#########################


# INPUT BAM folder
bam_folder = None
if "bam_folder" in config:
    bam_folder = os.path.join(CONFDIR, config["bam_folder"])

# INPUT FASTQ folder
FQ_INPUT_DIRECTORY = []
if not bam_folder:
    if not "input_fastq" in config:
        print("\"input_fastq\" not specified in config file. Exiting...")

    FQ_INPUT_DIRECTORY = os.path.join(CONFDIR, config["input_fastq"])
    if not os.path.exists(FQ_INPUT_DIRECTORY):
        print("Could not find {}".format(FQ_INPUT_DIRECTORY))

    MAPPED_BAM = "{sample}/alignment/{sample}_minimap2.bam"
else:
    MAPPED_BAM = find_file_in_folder(bam_folder, "*.bam", single=True)


# Input reference FASTA
FA_REF = os.path.join(CONFDIR, config["reference_fastq"])


# Parameter: sample_name
sample = "sv_sample01"
if "sample_name" in config:
    sample = config['sample_name']


# Parameter: target
target_bed = "{sample}/target.bed"
if "target" in config:
    target = os.path.join(CONFDIR, config["target"])
    if os.path.exists(target):
        # copyfile(target, "{}/target.bed".format(sample))
        target_bed = target
    else:
        print("Target BED {} not found. Continuing without target".format(target))


#########################
######## RULES ##########
#########################

rule all:
    input:
        expand("{name}/sv_calls/{name}_sniffles_filtered.vcf.gz", name=sample),
        expand("{name}/qc", name=sample),

rule call:
    input:
        expand("{name}/sv_calls/{name}_sniffles_filtered.vcf.gz", name=sample),

rule qc:
    input:
        expand("{name}/qc", name=sample),

rule json:
    input:
        expand("{name}/sv_calls/{name}_sniffles_filtered.vcf.gz", name=sample),
        expand("{name}/json/{name}_sniffles_filtered.json", name=sample),

rule eval:
    input:
        expand("{name}/eval/summary.txt", name=sample),

rule init:
    output: "init"
    conda: "env.yml"
    shell:
        "pip install {SNAKEDIR}/scripts &> {output}"

rule index_minimap2:
   input:
       REF = FA_REF
   output:
       "{sample}/index/minimap2.idx"
   threads: config['threads']
   conda: "env.yml"
   shell:
       "minimap2 -t {threads} -ax map-ont --MD -Y {input.REF} -d {output}"

rule map_minimap2:
   input:
       FQ = FQ_INPUT_DIRECTORY,
       IDX = rules.index_minimap2.output
   output:
       BAM = "{sample}/alignment/{sample}_minimap2.bam",
       BAI = "{sample}/alignment/{sample}_minimap2.bam.bai"
   params:
       min_qscore = config["min_qscore"] if "min_qscore" in config else 6,
       min_read_length = config["min_read_length"] if "min_read_length" in config else 1000,
       sort_threads = max(1, (max(1, config["threads"]) * 0.1))
   conda: "env.yml"
   threads: config["threads"]
   shell:
       "cat_fastq {input.FQ} | minimap2 -t {threads} -ax map-ont --MD -Y {input.IDX} - | samtools sort -@ {params.sort_threads} -o {output.BAM} - && samtools index -@ {threads} {output.BAM}"


rule bed_from_bam:
    input:
        BAM = MAPPED_BAM,
        SETUP = "init"
    output:
        "{sample}/target.bed"
    params:
        filter = "_Un _random"
    conda: "env.yml"
    shell:
        "bamref2bed -b {input.BAM} -f {params.filter} > {output}"


rule call_sniffles:
    input:
        BAM = MAPPED_BAM,
    output:
        VCF = "{sample}/sv_calls/{sample}_sniffles_tmp.vcf"
    params:
        read_support = 3,
        min_read_length = config['min_read_length'] if 'min_read_length' in config else 1000,
        min_mq = config['min_read_mapping_quality'] if 'min_read_mapping_quality' in config else 20
    conda: "env.yml"
    threads: config["threads"]
    shell:
        "sniffles -m {input.BAM} -v {output.VCF} -s {params.read_support} -r {params.min_read_length} -q {params.min_mq} --genotype --report_read_strands"


#rule filter_region:
#    input:
#        VCF = rules.call_sniffles.output.VCF,
#        BED = rules.bed_from_bam.output,
#        SETUP = "init"
#    output:
#        VCF = temp("{sample}/sv_calls/{sample}_sniffles_region_filtered.vcf")
#    conda: "env.yml"
#    shell:
#        "bedtools intersect -header -u -a {input.VCF} -b {input.BED} > {output.VCF}"


rule reformat_vcf:
    input:
         VCF = rules.call_sniffles.output.VCF,
         SETUP = "init"
    output:
         VCF = "{sample}/sv_calls/{sample}_sniffles.vcf"
    conda: "env.yml"
    shell:
         "sniffles-edit --ins-length --check --vcf-version -v {input.VCF} -o {output.VCF}"


rule filter_vcf:
    input:
         VCF = rules.reformat_vcf.output.VCF,
         RS = "{sample}/parameter_min_read_support.tsv",
         SETUP = "init"
    output:
         VCF = temp("{sample}/sv_calls/{sample}_sniffles_filtered_tmp.vcf")
    params:
        min_sv_length = config['min_sv_length'] if "min_sv_length" in config else 50,
        max_sv_length = config['max_sv_length'] if "max_sv_length" in config else 400000,
        strand_support = config['advanced_strand_support'] if "advanced_strand_support" in config else 0.05,
        sv_types = "DEL INS DUP",
        min_af = config['advanced_min_af'] if "advanced_min_af" in config else 0.15
    conda: "env.yml"
    shell:
         "sniffles-filter -v {input.VCF} -m `cat {input.RS}` -t {params.sv_types}  --strand-support {params.strand_support} -l {params.min_sv_length} --min-af {params.min_af} --max-length {params.max_sv_length} -o {output.VCF}"

rule sort_vcf:
    input:
        VCF = rules.filter_vcf.output.VCF
    output:
        VCF = temp("{sample}/sv_calls/{sample}_sniffles_filtered.vcf")
    conda: "env.yml"
    shell:
         "vcfsort {input.VCF} > {output.VCF}"


rule index_vcf:
    input:
         VCF = rules.sort_vcf.output.VCF
    output:
         VCF = "{sample}/sv_calls/{sample}_sniffles_filtered.vcf.gz"
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
    threads: config["threads"]
    conda: "env.yml"
    shell:
        "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.sample}_ --N50 --title {params.sample} --downsample 100000"


rule calc_depth:
    input:
         BAM = MAPPED_BAM,
         BED = target_bed
    output:
        DIR = directory("{sample}/depth"),
    threads: config["threads"]
    conda: "env.yml"
    shell:
         "mkdir -p {output.DIR}; mosdepth -x -t {threads} -n -b {input.BED} {output.DIR}/{sample} {input.BAM}"


rule auto_read_support:
    input:
         "{sample}/depth"
    output:
         "{sample}/parameter_min_read_support.tsv"
    #conda: "env.yml"
    run:
        with open(output[0], "w") as out:
            if "min_read_support" not in config or config["min_read_support"] == "auto":
                mosdepth_file = os.path.join(input[0], "{}.regions.bed.gz".format(wildcards.sample))
                with gzip.open(mosdepth_file, "r") as fh:
                    sum_depth = 0
                    count_depth = 0
                    total_size = 0
                    for line in fh:
                        if not line:
                            continue
                        cols = line.strip().split(b"\t")
                        sum_depth += float(cols[3]) * int(cols[2])
			total_size += int(cols[2])
                        count_depth += 1
                print(sum_depth, total_size, count_depth)
                min_rs = round((sum_depth / total_size) * 0.3)
            else:
                min_rs = config["min_read_support"]

            if min_rs < 5:
                print("Min read support < 5 not allowed. Falling back to minimum of 5.")
                min_rs = 5
            print(min_rs, file=out)


rule telemetry:
    input:
         BAM = MAPPED_BAM,
         VCF = rules.sort_vcf.output.VCF,
         DEPTH = rules.calc_depth.output.DIR,
         RS = rules.auto_read_support.output,
         SETUP = "init"
    output:
          "{sample}/json/{sample}_sniffles_filtered.json"
    conda: "env.yml"
    shell:
         "sniffles-telemetry -v {input.VCF} -b {input.BAM} -d {input.DEPTH}/{sample}.regions.bed.gz -s {input.RS} > {output}"


rule download_hg002_truthset:
    output:
        VCF = "HG002_SVs_Tier1_v0.6.vcf.gz",
        TBI = "HG002_SVs_Tier1_v0.6.vcf.gz.tbi",
        BED = "HG002_SVs_Tier1_v0.6.bed"
    conda: "env.yml"
    shell:
        "wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz && wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi && wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed"


rule download_truvari:
    output:
        "truvari.py"
    conda: "env.yml"
    shell:
         "wget https://raw.githubusercontent.com/spiralgenetics/truvari/develop/truvari.py"


rule intersect_target_highconf:
    input:
        TRUTH_BED = rules.download_hg002_truthset.output.BED,
        TARGET = target_bed
    output:
        BED = "{sample}/eval_high_conf.bed"
    conda: "env.yml"
    shell:
        "bedtools intersect -a {input.TRUTH_BED} -b {input.TARGET} > {output.BED}"


rule eval_vcf:
    input:
        SCRIPT = "truvari.py",
        VCF = rules.index_vcf.output.VCF,
        REF = FA_REF,
        TRUTH_VCF = rules.download_hg002_truthset.output.VCF,
        TRUTH_BED = rules.intersect_target_highconf.output.BED
    output:
        "{sample}/eval/summary.txt",
    conda: "env.yml"
    shell:
        "rmdir {sample}/eval/ && python {input.SCRIPT} --passonly -b {input.TRUTH_VCF} --includebed {input.TRUTH_BED} --pctsize 0 --pctsim 0 -t -c {input.VCF} -f {input.REF} -o {sample}/eval/"
