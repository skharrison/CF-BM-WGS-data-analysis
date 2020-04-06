import os
import pandas as pd

column_names = ["Sample", "Chromosome", "Position", "Ref Allele", "Alt Allele", "Allele Frequency"]
filenames = []
samplenames = []
df_list = []
for file in os.listdir("/Users/johnnavarro/AS_variants_list"):
    if file.endswith("_variantsfiltered.txt"):
        filenames.append(file)
for f in sorted(filenames):
    specific_names = []
    for val in column_names[0:4]:
        specific_names.append(val)
    for val in column_names[4:]:
        val = val + "_" + f.split("_")[0]
        specific_names.append(val)
    df_list.append(pd.read_csv(f,sep="\t",header=None,names=specific_names))
dfs2 = df_list[0].drop("Sample",axis=1)
sortednames = sorted(filenames)
i=1
for df in df_list[1:]:
    outfile = "AS_list" + str(i) + ".txt"
    dfs2 = pd.merge(dfs2,df.drop(["Sample"],axis=1),on=["Chromosome","Position","Ref Allele"], how="outer")
    dfs2.to_csv(outfile,sep="\t")
    i+=1