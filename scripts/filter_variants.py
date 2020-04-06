import os

for file in os.listdir("/Users/johnnavarro/AS_variants_list"):
    if file.endswith("_variants.txt"):
        with open(file) as fin:
            outfile = file.split(".")[0] + "filtered.txt"
            with open(outfile, "w+") as fout:
                for line in fin:
                    if float(line.split("\t")[5]) > 0.01:
                        fout.write(line)
