#!/usr/bin/python

import os
for file in os.listdir("/home/jnavarro/AS2Calls_test/ASpileups"):
    if file.endswith(".mpileup") == True:
        print(file)
        output_text = ""
        with open(file) as fin:   
            for line in fin:
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
                bases = line.split("\t")[4]
                for i in range(length):
                    try:
                        if bases[i] == "A" or bases[i] == "a":
                            if bases[i-2] != "-" or bases[i-1] != "+":
                                acount +=1
                                ifhas = True
                        elif bases[i] == "C" or bases[i] == "c":
                            if bases[i-2] != "-" or bases[i-1] != "+":
                                ccount +=1
                                ifhas = True
                        elif bases[i] == "T" or bases[i] == "t":
                            if bases[i-2] != "-" or bases[i-1] != "+":
                                tcount += 1
                                ifhas = True
                        elif bases[i] == "G" or bases[i] == "g":
                            if bases[i-2] != "-" or bases[i-1] != "+":
                                gcount += 1
                                ifhas = True
                    except IndexError:
                        continue
                if ifhas:
                    lstcounts = [acount,ccount,tcount,gcount]
                    lststuff = [acount/length,ccount/length,tcount/length,gcount/length]
                    totalreads = max(lstcounts)
                    allelefreq = max(lststuff)
                    altposition = lststuff.index(allelefreq)
                    altallele = ""
                    if altposition == 0:
                        altallele = "A"
                    elif altposition == 1:
                        altallele = "C"
                    elif altposition == 2:
                        altallele = "T"
                    elif altposition == 3:
                        altallele = "G"
                    output_text += sample + "\t" + contig + "\t" + str(pos) + "\t" + refallele + "\t" + altallele + "\t" + str(allelefreq) + "\t" + str(totalreads) + "\t" +  str(length) + "\n"
        fileout = sample + "_variants_withreads.txt"
        with open(fileout, "w+") as fout:
            fout.write(output_text)



                


                
    
            
