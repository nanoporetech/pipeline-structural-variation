import os
import glob
import gzip

from shutil import copyfile
from snakemake.utils import min_version

min_version("5.4.3")

configfile: "config.yml"

# Get parent directory of snakefile
SNAKEDIR = os.path.dirname(workflow.snakefile)

# Get parent directory of config file
CONFDIR = SNAKEDIR
if workflow.configfiles:
    conf_dir = workflow.configfiles[0]
    if os.path.isabs(conf_dir):
        CONFDIR = path.dirname(conf_dir)

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
FQ_INPUT_DIRECTORY = None
if not bam_folder:
    if not "input_fastq" in config:
        print("\"input_fastq\" not specified in config file. Exiting...")

    FQ_INPUT_DIRECTORY = os.path.join(CONFDIR, config["input_fastq"])
    if not os.path.exists(FQ_INPUT_DIRECTORY):
        print("Could not find {}".format(FQ_INPUT_DIRECTORY))

    MAPPED_BAM = "alignment/{sample}_minimap2.bam"
else:
    MAPPED_BAM = find_file_in_folder(bam_folder, "*.bam", single=True)


# Input reference FASTA
FA_REF = os.path.join(CONFDIR, config["reference_fastq"])


# Parameter: sample_name
sample = "sv_sample01"
if "sample_name" in config:
    sample = config['sample_name']


# Parameter: target
target_bed = "target_{sample}.bed"
if "target" in config:
    target = os.path.join(CONFDIR, config["target"])
    if os.path.exists(target):
        copyfile(target, "target_{}.bed".format(sample))
    else:
        print("Target BED {} not found. Continuing without target".format(target))



#########################
###### PREPROCESS #######
#########################

# Get FASTQ files from input directory
FASTQ_FILES = []
if FQ_INPUT_DIRECTORY:
    FASTQ_FILES = find_file_in_folder(FQ_INPUT_DIRECTORY)


#########################
######## RULES ##########
#########################

rule all:
    input:
        expand("sv_calls/{name}_sniffles_filtered.vcf.gz", name=sample),
        expand("qc/{name}/", name=sample)

rule call:
    input:
        expand("sv_calls/{name}_sniffles_filtered.vcf.gz", name=sample)

rule qc:
    input:
        expand("qc/{name}/", name=sample)

rule json:
    input:
        expand("json/{name}_sniffles_filtered.json", name=sample)


rule map_minimap2:
   input:
       FQ = FASTQ_FILES,
       REF = FA_REF
   output:
       BAM = MAPPED_BAM,
       BAI = MAPPED_BAM + ".bai"
   params:
       min_qscore = config["min_qscore"] if "min_qscore" in config else 6,
       min_read_length = config["min_read_length"] if "min_read_length" in config else 1000
   conda: "env.yml"
   threads: config["threads"]
   shell:
       "cat {input.FQ} | NanoFilt -q {params.min_qscore} -l {params.min_read_length} | minimap2 -t {threads} -ax map-ont --MD -Y {input.REF} - | samtools sort -o {output.BAM} - && samtools index -@ {threads} {output.BAM}"


rule bed_from_bam:
    input:
        BAM = rules.map_minimap2.output.BAM
    output:
        "target_{sample}.bed"
    params:
        filter = "_Un _random"
    conda: "env.yml"
    shell:
        "bamref2bed -b {input.BAM} -f {params.filter} > {output}"


rule call_sniffles:
    input:
        BAM = rules.map_minimap2.output.BAM,
    output:
        VCF = temp("sv_calls/{sample}_sniffles_tmp.vcf")
    params:
        read_support = 3,
        min_read_length = config['min_read_length'] if 'min_read_length' in config else 1000,
        min_mq = config['min_read_mapping_quality'] if 'min_read_mapping_quality' in config else 20
    conda: "env.yml"
    threads: config["threads"]
    shell:
        "sniffles -m {input.BAM} -v {output.VCF} -s {params.read_support} -r {params.min_read_length} -q {params.min_mq} --genotype --report_read_strands"


rule filter_region:
    input:
        VCF = rules.call_sniffles.output.VCF,
        BED = rules.bed_from_bam.output
    output:
        VCF = temp("sv_calls/{sample}_sniffles_region_filtered.vcf")
    conda: "env.yml"
    shell:
        "bedtools intersect -header -u -a {input.VCF} -b {input.BED} > {output.VCF}"


rule reformat_vcf:
    input:
         VCF = rules.filter_region.output.VCF
    output:
         VCF = "sv_calls/{sample}_sniffles.vcf"
    conda: "env.yml"
    shell:
         "sniffles-edit --ins-length --check --vcf-version -v {input.VCF} -o {output.VCF}"


rule filter_vcf:
    input:
         VCF = rules.reformat_vcf.output.VCF,
         RS = "{sample}_parameter_min_read_support.tsv"
    output:
         VCF = temp("sv_calls/{sample}_sniffles_filtered_tmp.vcf")
    params:
        min_read_support = lambda wildcards: get_param_from_file("{}_parameter_min_read_support.tsv".format(wildcards.sample), default=10),
        min_sv_length = config['min_sv_length'] if "min_sv_length" in config else 50,
        max_sv_length = config['max_sv_length'] if "max_sv_length" in config else 400000,
        strand_support = config['advanced_strand_support'] if "advanced_strand_support" in config else 0.05,
        sv_types = "DEL INS DUP",
        min_af = config['advanced_min_af'] if "advanced_min_af" in config else 0.15
    conda: "env.yml"
    shell:
         "sniffles-filter -v {input.VCF} -m {params.min_read_support} -t {params.sv_types}  --strand-support {params.strand_support} -l {params.min_sv_length} --min-af {params.min_af} --max-length {params.max_sv_length} -o {output.VCF}"

rule sort_vcf:
    input:
        VCF = rules.filter_vcf.output.VCF
    output:
        VCF = temp("sv_calls/{sample}_sniffles_filtered.vcf")
    conda: "env.yml"
    shell:
         "vcfsort {input.VCF} > {output.VCF}"


rule index_vcf:
    input:
         VCF = rules.sort_vcf.output.VCF
    output:
         VCF = "sv_calls/{sample}_sniffles_filtered.vcf.gz"
    conda: "env.yml"
    shell:
         "cat {input.VCF} | bgziptabix {output.VCF}"


rule nanoplot_qc:
    input:
        BAM = rules.map_minimap2.output.BAM
    output:
        DIR = directory("qc/{sample}")
    params:
        sample = sample
    threads: config["threads"]
    conda: "env.yml"
    shell:
        "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.sample}_ --N50 --title {params.sample}" # --downsample 1000


rule calc_depth:
    input:
         BAM = rules.map_minimap2.output.BAM,
         BED = target_bed
    output:
        DIR = directory("depth/{sample}/"),
    threads: config["threads"]
    conda: "env.yml"
    shell:
         "mosdepth -x -t {threads} -n -b {input.BED} {output.DIR}/{sample} {input.BAM}"


rule auto_read_support:
    input:
         "depth/{sample}/"
    output:
        "{sample}_parameter_min_read_support.tsv"
    conda: "env.yml"
    run:
        with open(output[0], "w") as out:
            if "min_read_support" not in config or config["min_read_support"] == "auto":
                mosdepth_file = os.path.join(input[0], "{}.regions.bed.gz".format(wildcards.sample))
                with gzip.open(mosdepth_file, "r") as fh:
                    sum_depth = 0
                    count_depth = 0
                    for line in fh:
                        if not line:
                            continue
                        cols = line.strip().split(b"\t")
                        sum_depth += float(cols[3])
                        count_depth += 1
                min_rs = round((sum_depth / count_depth) * 0.33)
            else:
                min_rs = config["min_read_support"]

            if min_rs < 5:
                print("Min read support < 5 not allowed. Falling back to minimum of 5.")
                min_rs = 5
            print(min_rs, file=out)


rule telemetry:
    input:
         BAM = rules.map_minimap2.output.BAM,
         VCF = rules.sort_vcf.output.VCF,
         DEPTH = rules.calc_depth.output.DIR,
         RS = rules.auto_read_support.output
    output:
          "json/{sample}_sniffles_filtered.json"
    conda: "env.yml"
    shell:
         "sniffles-telemetry -v {input.VCF} -b {input.BAM} -d {input.DEPTH}/{sample}.regions.bed.gz -s {input.RS} > {output}"
