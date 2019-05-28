import os
import gzip

min_rs = 10
min_rs_limit = 3

try:
    if "min_read_support" not in snakemake.config or snakemake.config["min_read_support"] == "auto":
        mosdepth_file = os.path.join(snakemake.input[0], "{}.regions.bed.gz".format(snakemake.wildcards.sample))
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
        min_rs = round((sum_depth / total_size) * 0.20)
    else:
        min_rs = snakemake.config["min_read_support"]

    if min_rs < min_rs_limit:
        print("Min read support < 5 not allowed. "
              "Falling back to minimum of {}.".format(min_rs_limit))
        min_rs = min_rs_limit
    print("Setting min. read support to {} for "
          "overall read depth of {}".format(min_rs,
                                            round((sum_depth / total_size))))
except Exception as e:
    print("Could not auto detect reads support. "
          "Falling back to {}".format(min_rs))

with open(snakemake.output[0], "w") as out:
    print(min_rs, file=out)
