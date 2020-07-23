import json                                                                                                                                                                                                
JSON=open("mitokons_annot.json")
dataannot = json.load(JSON)
JSON.close()

JSON=open("mitokons_inter.json")
datainter = json.load(JSON)
JSON.close()


BED = open("mitokons.bed","w")
BED.write("#chrom\tstart\tend\tcons_eukaryota\tcons_metazoa\tcons_bilateria\tcons_vertebrata\tcons_tetrapoda\tcons_amniota\tcons_mammalia\tcons_eutheria\tcons_euarchontoglires\tcons_primates\tjs_eukaryota\tjs_metazoa\tjs_bilateria\tjs_vertebrata\tjs_tetrapoda\tjs_amniota\tjs_mammalia\tjs_eutheria\tjs_euarchontoglires\tjs_primates\n")

for pos_nucl in dataannot:
    if not pos_nucl=="lst_truncated_codon_pos":
        cons_eukaryota = []
        cons_metazoa = []
        cons_bilateria = []
        cons_vertebrata = []
        cons_tetrapoda = []
        cons_amniota = []
        cons_mammalia = []
        cons_eutheria = []
        cons_euarchontoglires = []
        cons_primates = []
        js_eukaryota = []
        js_metazoa = []
        js_bilateria = []
        js_vertebrata = []
        js_tetrapoda = []
        js_amniota = []
        js_mammalia = []
        js_eutheria = []
        js_euarchontoglires = []
        js_primates = []
        for i in range(len(dataannot[pos_nucl]["mitocog_id"])):
            pos_aa = str(dataannot[pos_nucl]["pos_aa"][i])
            gene = dataannot[pos_nucl]["gene"][i]
            mitocog_id = dataannot[pos_nucl]["mitocog_id"][i]
            if dataannot[pos_nucl]["aa"][i]!="*" and mitocog_id in datainter:
                cons_eukaryota.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Eukaryota"])))
                cons_metazoa.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Metazoa"])))
                cons_bilateria.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Bilateria"])))
                cons_vertebrata.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Vertebrata"])))
                cons_tetrapoda.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Tetrapoda"])))
                cons_amniota.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Amniota"])))
                cons_mammalia.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Mammalia"])))
                cons_eutheria.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Eutheria"])))
                cons_euarchontoglires.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Euarchontoglires"])))
                cons_primates.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["%cons"]["Primates"])))
                js_eukaryota.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Eukaryota"],2)))
                js_metazoa.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Metazoa"],2)))
                js_bilateria.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Bilateria"],2)))
                js_vertebrata.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Vertebrata"],2)))
                js_tetrapoda.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Tetrapoda"],2)))
                js_amniota.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Amniota"],2)))
                js_mammalia.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Mammalia"],2)))
                js_eutheria.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Eutheria"],2)))
                js_euarchontoglires.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Euarchontoglires"],2)))
                js_primates.append(gene+":"+str(round(datainter[mitocog_id][pos_aa]["JS"]["Primates"],2)))
        if len(cons_eukaryota)>0:
            line = "MT\t"+str(int(pos_nucl)-1)+"\t"+pos_nucl
            line+="\t"+str("/").join(cons_eukaryota)
            line+="\t"+str("/").join(cons_metazoa)
            line+="\t"+str("/").join(cons_bilateria)
            line+="\t"+str("/").join(cons_vertebrata)
            line+="\t"+str("/").join(cons_tetrapoda)
            line+="\t"+str("/").join(cons_amniota)
            line+="\t"+str("/").join(cons_mammalia)
            line+="\t"+str("/").join(cons_eutheria)
            line+="\t"+str("/").join(cons_euarchontoglires)
            line+="\t"+str("/").join(cons_primates)
            line+="\t"+str("/").join(js_eukaryota)
            line+="\t"+str("/").join(js_metazoa)
            line+="\t"+str("/").join(js_bilateria)
            line+="\t"+str("/").join(js_vertebrata)
            line+="\t"+str("/").join(js_tetrapoda)
            line+="\t"+str("/").join(js_amniota)
            line+="\t"+str("/").join(js_mammalia)
            line+="\t"+str("/").join(js_eutheria)
            line+="\t"+str("/").join(js_euarchontoglires)
            line+="\t"+str("/").join(js_primates)+"\n"
            if len(dataannot[pos_nucl]["mitocog_id"])==1: BED.write(line.replace(gene+":","").replace("-1000.0","na"))
            else: BED.write(line.replace("-1000.0","na"))

BED.close()