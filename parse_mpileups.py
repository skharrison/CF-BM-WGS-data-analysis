#!/usr/bin/python

import re
import os
for file in os.listdir("/home/jnavarro/AS2Calls_test/ASpileups"):
    if file.endswith(".mpileup") == True:
        output_text = ""
        with open(file) as fin:   
            for line in fin:
                snps = []
                indels = []
                ifhas = False
                acount = 0
                ccount = 0
                tcount = 0
                gcount = 0
                sample = os.path.splitext(file)[0]
                contig = line.split("\t")[0]
                pos = line.split("\t")[1]
                refallele = line.split("\t")[2]
                length = int(line.split("\t")[3])
                for part in re.split("[.,^]",line.split("\t")[4]):
                    if "-" in part or "+" in part:
                        indels.append(part)
                        continue
                    elif "A" in part or "T" in part or "C" in part or "G" in part or "a" in part or "t" in part or "c" in part or "g" in part:
                        snps.append(part)
                        for base in part:
                            if base == "A" or base == "a":
                                acount+=1
                            elif base == "T" or base == "t":
                                tcount+=1
                            elif base == "C" or base == "c":
                                ccount+=1
                            elif base == "G" or base == "g":
                                gcount+=1
                allcounts = [acount,tcount,ccount,gcount]
                try:
                    allfreqs = [acount/length,tcount/length,ccount/length,gcount/length]
                except ZeroDivisionError:
                    continue
                totalreads = max(allcounts)
                maxfreq = max(allfreqs)
                altposition = allfreqs.index(maxfreq)
                altallele = ""
                if altposition == 0:
                    altallele = "A"
                elif altposition == 1:
                    altallele = "T"
                elif altposition == 2:
                    altallele = "C"
                elif altposition == 3:
                    altallele = "G"
                if maxfreq == 0.0:
                    continue
                output_text += sample + "\t" + contig + "\t" + str(pos) + "\t" + str(refallele) + "\t" + str(altallele) + "\t" + str(maxfreq) + "\t" + str(totalreads) + "\t" + str(length) + "\n"
        fileout = sample + "_variants_withreads.txt"
        print(file)
        with open(fileout, "w+") as fout:
            fout.write(output_text)
    
