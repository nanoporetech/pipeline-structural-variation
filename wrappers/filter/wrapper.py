import os
import gzip
from snakemake.shell import shell

threshold_lookup = ['0'] + ['2'] * 10 + ['3'] * 9 + ['5'] * 20 + ['8'] * 100

min_rs = 10
min_rs_limit = 2

try:
    mosdepth_file = os.path.join(snakemake.input.MOS, "{}.regions.bed.gz".format(snakemake.wildcards.sample))
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
    avg_depth = sum_depth / total_size
    print("Avg. depth found: {}".format(avg_depth))

    avg_depth = min(avg_depth, len(threshold_lookup) - 1)

    if "min_read_support" not in snakemake.config or snakemake.config["min_read_support"] == "auto":
        min_rs = int(threshold_lookup[round(avg_depth)])

    else:
        min_rs = snakemake.config["min_read_support"]

    if min_rs < min_rs_limit:
        print("Min read support {} not allowed. "
              "Falling back to minimum of {}.".format(min_rs, min_rs_limit))
        min_rs = min_rs_limit
    print("Setting min. read support to {} for "
          "overall read depth of {}".format(min_rs,
                                            round((sum_depth / total_size))))
except Exception as e:
    print("Could not auto detect reads support. "
          "Falling back to {}".format(min_rs))


sv_types = ["DEL", "INS", "DUP", "INV"]
bcf_string=""

types = list(set(snakemake.params.sv_types.split(" ")).intersection(sv_types))

print("Keeping SV when SVTYPE = {} and SVLEN > {} and < {} ".format(' or '.join(map(str, types)),snakemake.params.min_sv_length, snakemake.params.max_sv_length ))

for type in types:
    if type != types[-1]:
        bcf_string = str(bcf_string) + "SVTYPE = \"" + type + "\" || "
    else:
        bcf_string = str(bcf_string) + "SVTYPE = \"" + type + "\" "


if snakemake.params.target_bed is None:
    shell('bcftools view -i \'({bcf_string}) && ABS(SVLEN) > {snakemake.params.min_sv_length} && ABS(SVLEN) < {snakemake.params.max_sv_length} && INFO/RE >= {min_rs}\' {snakemake.input.VCF} > {snakemake.output.VCF}')
else:
    shell('bcftools view -T {snakemake.params.target_bed} -i \'({bcf_string}) && ABS(SVLEN) > {snakemake.params.min_sv_length} && ABS(SVLEN) < {snakemake.params.max_sv_length} && INFO/RE >= {min_rs}\' {snakemake.input.VCF} > {snakemake.output.VCF}')

