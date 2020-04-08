from functools import reduce
import os
import pandas as pd

column_names = ["Sample", "Chromosome", "Position", "Ref Allele", "Alt Allele", "Allele Frequency"]
df_list = []
for file in os.listdir("/Users/johnnavarro/AS_variants_list"):
    if file.endswith("_variants.txt"):
        df_list.append(pd.read_csv(file,sep="\t",header=None,names=column_names))
# combinethem = lambda d1,d2: pd.merge(d1,d2,on="Position", how="outer")
# dfs1 = reduce(combinethem,df_list)
dfs2 = df_list[0]
i=1
for df in df_list[1:]:
    print(i)
    dfs2 = pd.merge(dfs2,df,on="Position", how="outer")
    i+=1
# dfs1.to_csv("AS_list1.txt",sep="\t")
dfs2.to_csv("AS_list2.txt",sep="\t")
