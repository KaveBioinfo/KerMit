import os
import subprocess   
import datetime                                                                                                                                               
import time
import operator
import shutil
from KerMit_display import *
try: from cyvcf2 import VCF
except ImportError: sys.exit("""You need cyvcf2! install it from http://pypi.python.org/pypi/cyvcf2 or run pip install cyvcf2""")
try: import xlsxwriter
except ImportError: sys.exit("""You need xlsxwriter! install it from http://pypi.python.org/pypi/xlsxwriter or run pip install xlsxwriter""")


#***** VALIDATE VCF *****#
def validateVCF(dicoInit,pathVcf):
    # Report folder
    path_dir_validate = dicoInit["pathTmp"]+"/kermit_inputvalidateVCF"
    if os.path.isdir(path_dir_validate): shutil.rmtree(path_dir_validate)
    os.mkdir(path_dir_validate)
    # Create process
    cmd_validate = dicoInit["pathVcfvalidator"]+" -o "+path_dir_validate+" --require-evidence -l error < "+pathVcf
    process = subprocess.Popen([cmd_validate], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    # retrieve statuts & summary file path
    statut = ""
    boolvalid = False
    reportPath = ""
    lst_errors = []
    for line in err.decode('utf-8').split("\n"):
        if line.__contains__("Summary report written to"): reportPath = line.split(": ")[1]
        elif line.__contains__("the input file is not valid"): statut = "not valid"
        elif line.__contains__("the input file is valid"): statut = "valid"
    if statut=="" or reportPath=="": lst_errors.append(["program error"])
    elif statut=="not valid":
        IN = open(reportPath,'r') ; lst_lines = IN.read().split("\n") ; IN.close()
        for line in lst_lines:
            if "Error:" in line: lst_errors.append(line.split("Error: ")[1])
    else: boolvalid = True
    shutil.rmtree(path_dir_validate)
    return (boolvalid,lst_errors)



#***** VEP REFORMAT *****#
def vep_reformat(dicoInit):
    VCFvep = open(dicoInit["pathVEPvcf"],'r')
    lstLines = VCFvep.read().split("\n")
    VCFvep.close()
    VCFformated = open(dicoInit["pathVEPformatedvcf"],'w')
    for line in lstLines:
        lineNew = ""
        if line=="": continue
        # Modify vep header
        elif "INFO=<ID=CSQ" in line: lineNew = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Protein_position|Amino_acids\">"
        # Modify variant features
        elif line[0]!="#": 
            lstInfo = line.split("\t")[7].split(";")
            for info in lstInfo:
                fieldName = info.split("=")[0]
                fieldValue = info.split("=")[1]
                if fieldName=="CSQ":
                    lstNewField = []
                    for alleleField in fieldValue.split(","):
                        lstVEPfield = alleleField.split("|")
                        lstNewField.append(lstVEPfield[0]+"|"+lstVEPfield[1]+"|"+lstVEPfield[2]+"|"+lstVEPfield[3]+"|"+lstVEPfield[14]+"|"+lstVEPfield[15])
                    newfieldValue = ",".join(lstNewField)
                    lineNew = line.replace(fieldValue,newfieldValue).replace("chrM","MT")
        else: lineNew = line.replace("chrM","MT")
        lineNew = lineNew.replace("3_prime_UTR_variant","OTHERS")
        lineNew = lineNew.replace("coding_sequence_variant","CDS")
        lineNew = lineNew.replace("frameshift_variant","FS")
        lineNew = lineNew.replace("inframe_deletion","INFDel")
        lineNew = lineNew.replace("inframe_insertion","INFIns")
        lineNew = lineNew.replace("intergenic_variant","INT")
        lineNew = lineNew.replace("missense_variant","MIS")
        lineNew = lineNew.replace("non_coding_transcript_exon_variant","RNA")
        lineNew = lineNew.replace("protein_altering_variant","OTHERS")
        lineNew = lineNew.replace("start_lost","STARTlost")
        lineNew = lineNew.replace("stop_gained","STOPgain")
        lineNew = lineNew.replace("stop_lost","STOPlost")
        lineNew = lineNew.replace("stop_retained_variant","STOPret")
        lineNew = lineNew.replace("synonymous_variant","SYN")
        lineNew = lineNew.replace("transcript_ablation","OTHERS")
        VCFformated.write(lineNew+"\n")
    VCFformated.close()



#***** CHECK data and download if required (or --dicoInit["update"]) *****#
def check_and_download(dicoInit):
    # Required data
    pathClinvar = dicoInit["pathDirData"]+"/clinvar_MT.vcf.gz"
    pathHmtdbPatho = dicoInit["pathDirData"]+"/HMTDB_var_pa_tot.tsv.gz"
    pathHmtdbHealthy = dicoInit["pathDirData"]+"/HMTDB_var_tot.tsv.gz"
    pathMitomapPatho = dicoInit["pathDirData"]+"/MITOMAP_disease.vcf.gz"
    pathMitomapPoly = dicoInit["pathDirData"]+"/MITOMAP_polymorphisms.vcf.gz"
    pathMitotip = dicoInit["pathDirData"]+"/MITOTIP_scores.vcf.gz"
    pathMitimpact = dicoInit["pathDirData"]+"/MitImpact_db.vcf.gz"
    #***** ClinVar *****#
    if dicoInit["update"]==True or not os.path.isfile(pathClinvar) or os.path.getsize(pathClinvar)==0:
        printcolor(" Download ClinVar\n","0","greenlighter",False,dicoInit["colorBool"])
        try : os.remove(pathClinvar) ; os.remove(pathClinvar+".tbi")
        except: pass
        # Download
        cmd_download = "wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -P "+dicoInit["pathTmp"]
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","ClinVar VCF download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Download tbi
        cmd_download = "wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi -P "+dicoInit["pathTmp"]
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","ClinVar VCF index download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Filter MT chromosom
        cmd_filter = "bcftools filter -r MT -o "+pathClinvar+" -Oz "+dicoInit["pathTmp"]+"/clinvar.vcf.gz"
        process = subprocess.Popen([cmd_filter], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","ClinVar bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Tabix filtered file
        cmd_tabix = "tabix -f -p vcf "+pathClinvar
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","ClinVar tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])
    #***** Hmtdb *****#
    if dicoInit["update"]==True or not os.path.isfile(pathHmtdbPatho) or not os.path.isfile(pathHmtdbHealthy) or os.path.getsize(pathHmtdbPatho)==0 or os.path.getsize(pathHmtdbHealthy)==0:
        printcolor(" Download Hmtdb\n","0","greenlighter",False,dicoInit["colorBool"])
        try : os.remove(pathHmtdbPatho) ; os.remove(pathHmtdbHealthy) ; os.remove(pathHmtdbPatho+".tbi") ; os.remove(pathHmtdbHealthy+".tbi")
        except: pass
        # Download
        cmd_download = "wget -q https://www.hmtdb.uniba.it/download_file?dataset=var_tot.zip -P "+dicoInit["pathTmp"]
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Healthy download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_download = "wget -q https://www.hmtdb.uniba.it/download_file?dataset=var_pa_tot.zip -P "+dicoInit["pathTmp"]
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Patho download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Unzip healthy
        cmd_unzip = "unzip -p "+dicoInit["pathTmp"]+"/download_file?dataset=var_tot.zip > "+dicoInit["pathTmp"]+"/var_tot.csv"
        process = subprocess.Popen([cmd_unzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Healthy unzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_unzip = "unzip -p "+dicoInit["pathTmp"]+"/download_file?dataset=var_pa_tot.zip > "+dicoInit["pathTmp"]+"/var_pa_tot.csv"
        process = subprocess.Popen([cmd_unzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Patho unzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Remove header + replace comma by tab + Remove gap and convert position to int + add MT first column
        IN = open(dicoInit["pathTmp"]+"/var_tot.csv",'r') ; lst_lines = IN.read().split("\n") ; IN.close()
        OUT = open(dicoInit["pathTmp"]+"/var_tot_modif.csv",'w')
        for line in lst_lines:
            if line!="" and not "Position" in line: # no header
                split_line = line.split(",")
                if ".0" in split_line[0]: OUT.write("MT\t"+split_line[0].replace(".0","")+"\t"+split_line[1]+"\n")
        OUT.close()
        IN = open(dicoInit["pathTmp"]+"/var_pa_tot.csv",'r') ; lst_lines = IN.read().split("\n") ; IN.close()
        OUT = open(dicoInit["pathTmp"]+"/var_pa_tot_modif.csv",'w')
        for line in lst_lines:
            if line!="" and not "Position" in line: # no header
                split_line = line.split(",")
                if ".0" in split_line[0]: OUT.write("MT\t"+split_line[0].replace(".0","")+"\t"+split_line[1]+"\n")
        OUT.close()
        # Bgzip
        cmd_bgzip = "bgzip -c "+dicoInit["pathTmp"]+"/var_tot_modif.csv > "+pathHmtdbHealthy
        process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Healthy bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_bgzip = "bgzip -c "+dicoInit["pathTmp"]+"/var_pa_tot_modif.csv > "+pathHmtdbPatho
        process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Patho bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Touch
        cmd_touch = "touch -cr "+dicoInit["pathTmp"]+"/download_file?dataset=var_tot.zip "+pathHmtdbHealthy
        process = subprocess.Popen([cmd_touch], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Healthy touch\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_touch = "touch -cr "+dicoInit["pathTmp"]+"/download_file?dataset=var_tot.zip "+pathHmtdbPatho
        process = subprocess.Popen([cmd_touch], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Patho touch\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Tabix
        cmd_tabix = "tabix -b 2 -e 2 "+pathHmtdbHealthy
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Healthy tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_tabix = "tabix -b 2 -e 2 "+pathHmtdbPatho
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","Hmtdb Patho tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])
    #***** MITOMAP *****#
    if dicoInit["update"]==True or not os.path.isfile(pathMitomapPatho) or not os.path.isfile(pathMitomapPoly) or os.path.getsize(pathMitomapPatho)==0 or os.path.getsize(pathMitomapPoly)==0:
        printcolor(" Download MITOMAP\n","0","greenlighter",False,dicoInit["colorBool"])
        try : os.remove(pathMitomapPatho) ; os.remove(pathMitomapPoly) ; os.remove(pathMitomapPatho+".tbi") ; os.remove(pathMitomapPoly+".tbi")
        except: pass
        # Download
        cmd_download = "wget -q https://mitomap.org/cgi-bin/disease.cgi?format=vcf -O "+pathMitomapPatho.replace(".gz","")
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MITOMAP diseases download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_download = "wget -q https://mitomap.org/cgi-bin/polymorphisms.cgi?format=vcf -O "+pathMitomapPoly.replace(".gz","")
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MITOMAP polymorphisms download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Correct malformed & redundant field in MITOMAP diseases
        VCF = open(pathMitomapPatho.replace(".gz",""),'r')
        lst_lines = VCF.read().split("\n")
        VCF.close()
        VCF = open(pathMitomapPatho.replace(".gz",""),'w')
        for line in lst_lines:
            if line!="":
                newLine = line
                if line[0]!="#":
                    for field in line.split("\t")[7].split(";"):
                        key = field.split("=")[0]
                        value = field.split("=")[1]
                        if key=="Disease":
                            newValue = ", ".join(set(value.replace("|",",").replace("/",",").replace("-,-",",").replace(",-",",").replace("-"," ").split(",")))
                            newLine = newLine.replace(key+"="+value,key+"="+newValue)
                        elif key=="DiseaseStatus":
                            set_status = set()
                            if "cfrm" in value.lower(): set_status.add("Cfrm")
                            if "reported" in value.lower(): set_status.add("Reported")
                            if "conflicting" in value.lower(): set_status.add("Conflict")
                            if "unclear" in value.lower(): set_status.add("Unclear")
                            if "marker" in value.lower() or "population" in value.lower(): set_status.add("HG-marker")
                            if "synergistic" in value.lower(): set_status.add("Synergistic")
                            if "secondary" in value.lower(): set_status.add("Secondary")
                            newValue = ",".join(set_status)
                            newLine = newLine.replace(key+"="+value,key+"="+newValue)
                VCF.write(newLine+"\n")
        VCF.close()
        # Bgzip
        cmd_bgzip = "bgzip "+pathMitomapPatho.replace(".gz","")
        process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MITOMAP diseases bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_bgzip = "bgzip "+pathMitomapPoly.replace(".gz","")
        process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MITOMAP polymorphisms bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Tabix
        cmd_tabix = "tabix -p vcf "+pathMitomapPatho
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MITOMAP diseases tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        cmd_tabix = "tabix -p vcf "+pathMitomapPoly
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MITOMAP polymorphisms tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])
    #***** MitoTip *****#
    if dicoInit["update"]==True or not os.path.isfile(pathMitotip) or os.path.getsize(pathMitotip)==0:
        printcolor(" Download MitoTip\n","0","greenlighter",False,dicoInit["colorBool"])
        try : os.remove(pathMitotip) ; os.remove(pathMitotip+".tbi")
        except: pass
        # Download
        cmd_download = "wget -q https://mitomap.org/downloads/mitotip_scores.txt -O "+dicoInit["pathTmp"]+"/mitotip.tsv"
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitoTip download\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Retrieve mtDNA sequence
        FASTA = open(dicoInit["pathFasta"],'r')
        seq = "".join(FASTA.read().split("\n")[1:])
        FASTA.close()
        # Convert to VCF
        filedate = datetime.datetime.strptime(time.ctime(os.path.getmtime(dicoInit["pathTmp"]+"/mitotip.tsv")), "%a %b %d %H:%M:%S %Y").strftime("%Y%m%d")
        headerVCF = "##fileformat=VCFv4.2\n"\
                    "##fileDate="+filedate+"\n"\
                    "##source=https://mitomap.org/downloads/mitotip_scores.txt\n"\
                    "##reference=http://www.ncbi.nlm.nih.gov/nuccore/251831106\n"\
                    "##contig=<ID=MT,length=16569,assembly=B37>\n"\
                    "##INFO=<ID=MTs,Number=.,Type=String,Description=\"MitoTip score\">\n"\
                    "##INFO=<ID=MTq,Number=.,Type=String,Description=\"MitoTip quartile\">\n"\
                    "#CHROM  POS ID  REF ALT QUAL    FILTER  INFO\n"
        TSV = open(dicoInit["pathTmp"]+"/mitotip.tsv",'r')
        lst_lines = TSV.read().split("\n")
        TSV.close()
        VCF = open(pathMitotip.replace(".gz",""),"w")
        VCF.write(headerVCF)
        dico_line = {}
        for line in lst_lines:
            if line!="" and not "Position" in line:
                split_line = line.split("\t")
                pos = split_line[0]
                ref = split_line[1]
                alt = split_line[2]
                # for deletion search after nucl
                if alt==":":
                    alt = seq[int(pos)]
                    ref = seq[int(pos)-1]+alt
                score = split_line[3]
                quartile = split_line[4]
                vcfLine =  "MT\t"+pos+"\t.\t"+ref+"\t"+alt+"\t.\tPASS\tMTs="+score+";MTq="+quartile+"\n"
                try: dico_line[int(pos)].append(vcfLine)
                except: dico_line[int(pos)] = [vcfLine]
        sorted_dict = dict(sorted(dico_line.items(), key=operator.itemgetter(0)))
        for pos in sorted_dict:
            for vcfLine in sorted_dict[pos]: VCF.write(vcfLine)
        VCF.close()
        # Bgzip
        cmd_bgzip = "bgzip "+pathMitotip.replace(".gz","")
        process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitoTip bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Tabix
        cmd_tabix = "tabix -p vcf "+pathMitotip
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitoTip tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])
    #***** MitImpact *****#
    if dicoInit["update"]==True or not os.path.isfile(pathMitimpact) or os.path.getsize(pathMitimpact)==0:
        printcolor(" Download MitImpact\n","0","greenlighter",False,dicoInit["colorBool"])
        try : os.remove(pathMitimpact) ; os.remove(pathMitimpact+".tbi")
        except: pass
        # Check version
        try: mitimpactVersion = subprocess.check_output("curl -s -L https://mitimpact.css-mendel.it/ | grep -o -m 1 -P \"db_.+.txt.zip\"",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1]
        except: printError("KerMit download","Unable to check MitImpact version\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        mitimpactUrl = "https://mitimpact.css-mendel.it/cdn/MitImpact_"+mitimpactVersion
        # Download
        cmd_download = "wget -q "+mitimpactUrl+" -P "+dicoInit["pathTmp"]
        process = subprocess.Popen([cmd_download], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitImpact\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Unzip
        cmd_unzip = "unzip -p "+dicoInit["pathTmp"]+"/MitImpact_"+mitimpactVersion+" > "+dicoInit["pathTmp"]+"/MitImpact_"+mitimpactVersion.replace(".zip","")
        process = subprocess.Popen([cmd_unzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitImpact unzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])      
        # Convert to VCF
        filedate = datetime.datetime.strptime(time.ctime(os.path.getmtime(dicoInit["pathTmp"]+"/MitImpact_"+mitimpactVersion.replace(".zip",""))), "%a %b %d %H:%M:%S %Y").strftime("%Y%m%d")
        headerVCF = "##fileformat=VCFv4.2\n"\
                    "##fileDate="+filedate+"\n"\
                    "##source=https://mitimpact.css-mendel.it?format=vcf\n"\
                    "##reference=http://www.ncbi.nlm.nih.gov/nuccore/251831106\n"\
                    "##contig=<ID=MT,length=16569,assembly=B37>\n"\
                    "##INFO=<ID=P2,Number=.,Type=String,Description=\"Categorical prediction for PolyPhen ver. 2.2.2 (HumDiv dataset)\">\n"\
                    "##INFO=<ID=SIFT,Number=.,Type=String,Description=\"Categorical prediction for SIFT ver. 5.0.3.\">\n"\
                    "##INFO=<ID=FATW,Number=.,Type=String,Description=\"Categorical prediction for FatHmm ver. 2.3\">\n"\
                    "##INFO=<ID=FAT,Number=.,Type=String,Description=\"Categorical prediction for FatHmm ver. 2.2 (unweighted version).\">\n"\
                    "##INFO=<ID=PROV,Number=.,Type=String,Description=\"Categorical prediction for PROVEAN ver. 1.3.\">\n"\
                    "##INFO=<ID=MUTASS,Number=.,Type=String,Description=\"Categorical prediction for MutationAssessor ver. 2.0.\">\n"\
                    "##INFO=<ID=EFINSP,Number=.,Type=String,Description=\"Categorical prediction for EFIN SP (trained with SwissProt dataset).\">\n"\
                    "##INFO=<ID=EFINHD,Number=.,Type=String,Description=\"Categorical prediction for EFIN HD (trained with HumDiv dataset)\">\n"\
                    "##INFO=<ID=CADD,Number=.,Type=String,Description=\"Categorical prediction for CADD ver. 1.2\">\n"\
                    "##INFO=<ID=PANTHER,Number=.,Type=String,Description=\"Categorical prediction for PANTHER\">\n"\
                    "##INFO=<ID=PHDSNP,Number=.,Type=String,Description=\"Categorical prediction for PhD-SNP\">\n"\
                    "##INFO=<ID=SNAP,Number=.,Type=String,Description=\"Categorical prediction for SNAP\">\n"\
                    "##INFO=<ID=METASNP,Number=.,Type=String,Description=\"Categorical prediction for METASNP\">\n"\
                    "##INFO=<ID=CAROL,Number=.,Type=String,Description=\"Categorical prediction for CAROL consensus method.\">\n"\
                    "##INFO=<ID=CONDEL,Number=.,Type=String,Description=\"Categorical prediction for Condel consensus method.\">\n"\
                    "##INFO=<ID=COVEC,Number=.,Type=String,Description=\"Categorical prediction for COVEC ver. 0.4 Weighted Majority Rule consensus method\">\n"\
                    "##INFO=<ID=MTOOLBOX,Number=.,Type=String,Description=\"Categorical prediction for MtoolBox consensus method\">\n"\
                    "##INFO=<ID=APOGEE,Number=.,Type=String,Description=\"Consensus assessment of all categorical predictions obtained during the bootstrap loop (July 2015)\">\n"\
                    "##INFO=<ID=MTASTER,Number=.,Type=String,Description=\"Categorical prediction (disease_causing, disease_causing_automatic, polymorphism)\">\n"\
                    "##INFO=<ID=MTCLASS,Number=.,Type=String,Description=\"Categorical prediction (neutral, damaging) of the SVM classifier.\">\n"\
                    "#CHROM  POS ID  REF ALT QUAL    FILTER  INFO\n"
        TSV = open(dicoInit["pathTmp"]+"/MitImpact_"+mitimpactVersion.replace(".zip",""),'r')
        lst_lines = TSV.read().split("\n")
        TSV.close()
        VCF = open(pathMitimpact.replace(".gz",""),"w")
        VCF.write(headerVCF)
        dico_line = {}
        for line in lst_lines:
            if line!="" and not "Gene_symbol" in line:
                split_line = line.split("\t")
                # Format prediction string
                if split_line[22]=="benign": p2 = "N"
                elif split_line[22].replace(" ","_")=="possibly_damaging": p2 = "M"
                elif split_line[22]=="probably_damaging": p2 = "D"
                elif split_line[22]=="unknown" or split_line[22]==".": p2 = "."
                if split_line[24]=="neutral": sift = "N"
                elif split_line[24]=="deleterious": sift = "D"
                elif split_line[24]==".": sift = "."
                if split_line[26]=="neutral": fatw = "N"
                elif split_line[26]=="deleterious": fatw = "D"
                elif split_line[26]==".": fatw = "."
                if split_line[28]=="neutral": fat = "N"
                elif split_line[28]=="deleterious": fat = "D"
                elif split_line[28]==".": fat = "."
                if split_line[30]=="neutral": prov = "N"
                elif split_line[30]=="deleterious": prov = "D"
                elif split_line[30]==".": prov = "."
                if split_line[32]=="neutral_impact": mutass = "N"
                elif split_line[32]=="low_impact": mutass = "M"
                elif split_line[32]=="medium_impact": mutass = "M"
                elif split_line[32]=="high_impact": mutass = "D"
                elif split_line[32]==".": mutass = "."
                if split_line[35]=="neutral": efinsp = "N"
                elif split_line[35]=="damaging": efinsp = "D"
                elif split_line[35]==".": efinsp = "."
                if split_line[37]=="neutral": efinhd = "N"
                elif split_line[37]=="damaging": efinhd = "D"
                elif split_line[37]==".": efinhd = "."
                if split_line[40]=="neutral": cadd = "N"
                elif split_line[40]=="deleterious": cadd = "D"
                elif split_line[40]==".": cadd = "."
                if split_line[44]=="neutral": panther = "N"
                elif split_line[44]=="disease": panther = "D"
                elif split_line[44]==".": panther = "."
                if split_line[46]=="neutral": phdsnp = "N"
                elif split_line[46]=="disease": phdsnp = "D"
                elif split_line[46]==".": phdsnp = "."
                if split_line[48]=="neutral": snap = "N"
                elif split_line[48]=="disease": snap = "D"
                elif split_line[48]==".": snap = "."
                if split_line[50]=="neutral": metasnp = "N"
                elif split_line[50]=="disease": metasnp = "D"
                elif split_line[50]==".": metasnp = "."
                if split_line[53]=="neutral": carol = "N"
                elif split_line[53]=="deleterious": carol = "D"
                elif split_line[53]==".": carol = "."
                if split_line[55]=="neutral": condel = "N"
                elif split_line[55]=="deleterious": condel = "D"
                elif split_line[55]==".": condel = "."
                if split_line[57]=="neutral": covec = "N"
                elif split_line[57]=="deleterious": covec = "D"
                elif split_line[57]==".": covec = "."
                if split_line[59]=="neutral": mtoolbox = "N"
                elif split_line[59]=="deleterious": mtoolbox = "D"
                elif split_line[59]==".": mtoolbox = "."
                if split_line[69]=="N": apogee = "N"
                elif split_line[69]=="P": apogee = "D"
                elif split_line[69]==".": apogee = "."
                if split_line[71]=="polymorphism": mtaster = "N"
                elif split_line[71]=="disease_causing_automatic" or split_line[71]=="disease_causing": mtaster = "D"
                elif split_line[71]==".": mtaster = "."
                if split_line[91]=="neutral": mtclass = "N"
                elif split_line[91]=="damaging": mtclass = "D"
                elif split_line[91]==".": mtclass = "."
                vcfLine =  "MT\t"+split_line[1]+"\t.\t"+split_line[3]+"\t"+split_line[4]+"\t.\tPASS\t"+\
                           "P2="+p2+";SIFT="+sift+";FATW="+fatw+";FAT="+fat+";PROV="+prov+";MUTASS="+mutass+";EFINSP="+efinsp+\
                           ";EFINHD="+efinhd+";CADD="+cadd+";PANTHER="+panther+";PHDSNP="+phdsnp+";SNAP="+snap+";METASNP="+metasnp+";CAROL="+carol+\
                           ";CONDEL="+condel+";COVEC="+covec+";MTOOLBOX="+mtoolbox+";APOGEE="+apogee+";MTASTER="+mtaster+";MTCLASS="+mtclass+"\n"
                try: dico_line[int(split_line[1])].append(vcfLine)
                except: dico_line[int(split_line[1])] = [vcfLine]
        sorted_dict = dict(sorted(dico_line.items(), key=operator.itemgetter(0)))
        for pos in sorted_dict:
            for vcfLine in sorted_dict[pos]: VCF.write(vcfLine)
        VCF.close()
        # bgzip
        cmd_bgzip = "bgzip "+pathMitimpact.replace(".gz","")
        process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitImpact bgzip\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Touch to keep modif date
        cmd_touch = "touch -cr "+dicoInit["pathTmp"]+"/MitImpact_"+mitimpactVersion.replace(".zip","")+" "+pathMitimpact
        process = subprocess.Popen([cmd_touch], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitImpact touch\n    "+err.decode('utf-8'),dicoInit["colorBool"])
        # Tabix
        cmd_tabix = "tabix -p vcf "+pathMitimpact
        process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit download","MitImpact tabix\n    "+err.decode('utf-8'),dicoInit["colorBool"])



#***** RETRIEVE data files version *****#
def get_data_version(dicoInit):
    dicoDataVersion = {}
    # Required data
    pathClinvar = dicoInit["pathDirData"]+"/clinvar_MT.vcf.gz"
    pathHmtdbPatho = dicoInit["pathDirData"]+"/HMTDB_var_pa_tot.tsv.gz"
    pathHmtdbHealthy = dicoInit["pathDirData"]+"/HMTDB_var_tot.tsv.gz"
    pathMitomapPatho = dicoInit["pathDirData"]+"/MITOMAP_disease.vcf.gz"
    pathMitomapPoly = dicoInit["pathDirData"]+"/MITOMAP_polymorphisms.vcf.gz"
    pathMitotip = dicoInit["pathDirData"]+"/MITOTIP_scores.vcf.gz"
    pathMitimpact = dicoInit["pathDirData"]+"/MitImpact_db.vcf.gz"
    # For vcf get "fileDate="
    try: version = subprocess.check_output("zgrep -m 1 \"fileDate\" "+pathMitomapPatho+" | cut -d \"=\" -f 2",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1] ; dicoDataVersion["MITOMAP patho"] = version[:4]+"-"+version[4:6]+"-"+version[-2:]
    except: print("oioio") ; printError("KerMit download","Unable to check MITOMAP disease version",dicoInit["colorBool"])
    try: dicoDataVersion["Clinvar"] = subprocess.check_output("zgrep -m 1 \"fileDate\" "+pathClinvar+" | cut -d \"=\" -f 2",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1]
    except: printError("KerMit download","Unable to check Clinvar polymorphisms version",dicoInit["colorBool"])
    try: version = subprocess.check_output("zgrep -m 1 \"fileDate\" "+pathMitomapPoly+" | cut -d \"=\" -f 2",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1] ; dicoDataVersion["MITOMAP_poly"] = version[:4]+"-"+version[4:6]+"-"+version[-2:]
    except: printError("KerMit download","Unable to check MITOMAP polymorphisms version",dicoInit["colorBool"])
    try: version = subprocess.check_output("zgrep -m 1 \"fileDate\" "+pathMitimpact+" | cut -d \"=\" -f 2",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1] ; dicoDataVersion["MitImpact"] = version[:4]+"-"+version[4:6]+"-"+version[-2:]
    except: printError("KerMit download","Unable to check MitImpact version",dicoInit["colorBool"])
    try: version = subprocess.check_output("zgrep -m 1 \"fileDate\" "+pathMitotip+" | cut -d \"=\" -f 2",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1] ; dicoDataVersion["MitoTip"] = version[:4]+"-"+version[4:6]+"-"+version[-2:]
    except: printError("KerMit download","Unable to check MitoTip version",dicoInit["colorBool"])
    # for Hmtdb & get file creation date (keep with touch from zip)
    try: dicoDataVersion["HMTDB healthy"] = datetime.datetime.strptime(time.ctime(os.path.getmtime(pathHmtdbHealthy)), "%a %b %d %H:%M:%S %Y").strftime("%Y-%m-%d")
    except: printError("KerMit download","Unable to check Hmtdb Healthy date",dicoInit["colorBool"])
    try: dicoDataVersion["HMTDB patho"] = datetime.datetime.strptime(time.ctime(os.path.getmtime(pathHmtdbPatho)), "%a %b %d %H:%M:%S %Y").strftime("%Y-%m-%d")
    except: printError("KerMit download","Unable to check Hmtdb Patho date",dicoInit["colorBool"])   
    return dicoDataVersion



#***** FORMAT non_VCF vcfanno header *****#
def reheaderVCF(dicoInit):
    VCF = open(dicoInit["pathVCFout"].replace(".gz",""),'r')
    data = VCF.read()
    VCF.close()
    dicoReplace = {"calculated by self of overlapping values in column 4 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Eukaryota",\
                   "calculated by self of overlapping values in column 5 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Metazoa",\
                   "calculated by self of overlapping values in column 6 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Bilateria",\
                   "calculated by self of overlapping values in column 7 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Vertebrata",\
                   "calculated by self of overlapping values in column 8 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Tetrapoda",\
                   "calculated by self of overlapping values in column 9 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Amniota",\
                   "calculated by self of overlapping values in column 10 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Mammalia",\
                   "calculated by self of overlapping values in column 11 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Eutheria",\
                   "calculated by self of overlapping values in column 12 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Euarchontoglires",\
                   "calculated by self of overlapping values in column 13 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons conservation % for phylum Primates",\
                   "calculated by self of overlapping values in column 14 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Eukaryota",\
                   "calculated by self of overlapping values in column 15 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Metazoa",\
                   "calculated by self of overlapping values in column 16 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Bilateria",\
                   "calculated by self of overlapping values in column 17 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Vertebrata",\
                   "calculated by self of overlapping values in column 18 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Tetrapoda",\
                   "calculated by self of overlapping values in column 19 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Amniota",\
                   "calculated by self of overlapping values in column 20 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Mammalia",\
                   "calculated by self of overlapping values in column 21 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Eutheria",\
                   "calculated by self of overlapping values in column 22 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Euarchontoglires",\
                   "calculated by self of overlapping values in column 23 from /home/dooguy/Dev_prog/Dev_KerMit/data/MitoKons.bed.gz":"MitoKons jensen-shannon divergence score for Primates",\
                   "calculated by self of overlapping values in column 3 from /home/dooguy/Dev_prog/Dev_KerMit/data/HMTDB_var_tot.tsv.gz":"HMTDB site variability for Healthy",\
                   "calculated by self of overlapping values in column 3 from /home/dooguy/Dev_prog/Dev_KerMit/data/HMTDB_var_pa_tot.tsv.gz":"HMTDB site variability for Pathologic",\
                   "calculated by self of overlapping values in column 4 from /home/dooguy/Dev_prog/Dev_KerMit/data/chrM_HP.bed.gz":"Homopolymers count",\
                   "#CHROM":"##KerMit="+dicoInit["kermitVersion"]+" (https://github.com/dooguypapua/KerMit)\n#CHROM"}
    for oldStr in dicoReplace:
        data = data.replace(oldStr,dicoReplace[oldStr])
    VCF = open(dicoInit["pathVCFout"].replace(".gz",""),'w')
    VCF.write(data)
    VCF.close()
    # sed to remove "(from /home/dooguy/Dev_prog/Dev_KerMit/data/MitImpact_db.vcf.gz)"
    cmd_sed = "sed -i -E s/\"\(from.+\)\"/\"\"/g "+dicoInit["pathVCFout"].replace(".gz","")
    process = subprocess.Popen([cmd_sed], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: printError("KerMit postprocessing","sed format header\n    "+err.decode('utf-8'),dicoInit["colorBool"])



#***** CREATE EXCEL OUTPUT *****#
def makeXLSX(dicoInit,dicoHaplo):
    workbook = xlsxwriter.Workbook(dicoInit["pathXLSXout"])
    worksheet = workbook.add_worksheet()
    # Add a bold format to use to highlight cells.
    rowFormat = workbook.add_format({'align': 'center','valign': 'vcenter','border':1,'font_size':9,'text_wrap': True})
    mergeFormat1 = workbook.add_format({'bold': True,'align': 'center','valign': 'vcenter','border':1,'font_size':10,'bg_color':"#87deaa"})
    headerFormat1 = workbook.add_format({'bold': True,'align': 'center','valign': 'vcenter','border':1,'font_size':9,'bg_color':"#87deaa"})
    mergeFormat2 = workbook.add_format({'bold': True,'align': 'center','valign': 'vcenter','border':1,'font_size':10,'bg_color':"#beb7c8"})
    headerFormat2 = workbook.add_format({'bold': True,'align': 'center','valign': 'vcenter','border':1,'font_size':9,'bg_color':"#beb7c8"})
    # Write some data headers.
    worksheet.merge_range('A1:C1','VARIANT',mergeFormat1)
    worksheet.write('A2','Pos',headerFormat1) ; worksheet.write('A3','',headerFormat1)
    worksheet.write('B2','Ref',headerFormat1) ; worksheet.write('B3','',headerFormat1)
    worksheet.write('C2','Alt',headerFormat1) ; worksheet.write('C3','',headerFormat1)
    worksheet.merge_range('D1:H1','CALLING',mergeFormat2)
    worksheet.write('D2','Filter',headerFormat2) ; worksheet.write('D3','',headerFormat2)
    worksheet.write('E2','HP',headerFormat2) ; worksheet.write('E3','',headerFormat2)
    worksheet.write('F2','Qual',headerFormat2) ; worksheet.write('F3','',headerFormat2)
    worksheet.write('G2','AF',headerFormat2) ; worksheet.write('G3','',headerFormat2)
    worksheet.write('H2','DP',headerFormat2) ; worksheet.write('H3','',headerFormat2)
    worksheet.merge_range('I1:M1','ANNOTATION',mergeFormat1)
    worksheet.write('I2','Csq',headerFormat1) ; worksheet.write('I3','',headerFormat1)
    worksheet.write('J2','Impact',headerFormat1) ; worksheet.write('J3','',headerFormat1)
    worksheet.write('K2','Gene',headerFormat1) ; worksheet.write('K3','',headerFormat1)
    worksheet.write('L2','AA',headerFormat1) ; worksheet.write('L3','',headerFormat1)
    worksheet.write('M2',dicoHaplo['hg']+" ("+dicoHaplo['qual']+")",headerFormat1) ; worksheet.write('M3','',headerFormat1)
    worksheet.merge_range('N1:S1','MITOMAP',mergeFormat2)
    worksheet.write('N2','Count',headerFormat2) ; worksheet.write('N3','',headerFormat2)
    worksheet.write('O2','Freq',headerFormat2) ; worksheet.write('O3','',headerFormat2)
    worksheet.write('P2','Status',headerFormat2) ; worksheet.write('P3','',headerFormat2)
    worksheet.write('Q2','Disease',headerFormat2) ; worksheet.write('Q3','',headerFormat2)
    worksheet.write('R2','Hom',headerFormat2) ; worksheet.write('R3','',headerFormat2)
    worksheet.write('S2','Htz',headerFormat2) ; worksheet.write('S3','',headerFormat2)
    worksheet.merge_range('T1:W1','ClinVar',mergeFormat1)
    worksheet.write('T2','Id',headerFormat1) ; worksheet.write('T3','',headerFormat1)
    worksheet.write('U2','Disease',headerFormat1) ; worksheet.write('U3','',headerFormat1)
    worksheet.write('V2','Status',headerFormat1) ; worksheet.write('V3','',headerFormat1)
    worksheet.write('W2','Significance',headerFormat1) ; worksheet.write('W3','',headerFormat1)
    worksheet.merge_range('X1:Y1','HMTDB Site Var',mergeFormat2)
    worksheet.write('X2','Healthy',headerFormat2) ; worksheet.write('X3','',headerFormat2)
    worksheet.write('Y2','Disease',headerFormat2) ; worksheet.write('Y3','',headerFormat2)
    worksheet.merge_range('Z1:AA1','MitoTip',mergeFormat1)
    worksheet.write('Z2','Score',headerFormat1) ; worksheet.write('Z3','',headerFormat1)
    worksheet.write('AA2','Rank',headerFormat1) ; worksheet.write('AA3','',headerFormat1)
    worksheet.merge_range('AB1:AU1','MitImpact',mergeFormat2)
    worksheet.write('AB2','Poly2',headerFormat2) ; worksheet.write('AB3','',headerFormat2)
    worksheet.write('AC2','SIFT',headerFormat2) ; worksheet.write('AC3','',headerFormat2)
    worksheet.write('AD2','FatHmmW',headerFormat2) ; worksheet.write('AD3','',headerFormat2)
    worksheet.write('AE2','FatHmm',headerFormat2) ; worksheet.write('AE3','',headerFormat2)
    worksheet.write('AF2','PROVEAN',headerFormat2) ; worksheet.write('AF3','',headerFormat2)
    worksheet.write('AG2','MutAss',headerFormat2) ; worksheet.write('AG3','',headerFormat2)
    worksheet.write('AH2','EFIN-SP',headerFormat2) ; worksheet.write('AH3','',headerFormat2)
    worksheet.write('AI2','EFIN-HD',headerFormat2) ; worksheet.write('AI3','',headerFormat2)
    worksheet.write('AJ2','CADD',headerFormat2) ; worksheet.write('AJ3','',headerFormat2)
    worksheet.write('AK2','PANTHER',headerFormat2) ; worksheet.write('AK3','',headerFormat2)
    worksheet.write('AL2','PhD-SNP',headerFormat2) ; worksheet.write('AL3','',headerFormat2)
    worksheet.write('AM2','SNAP',headerFormat2) ; worksheet.write('AM3','',headerFormat2)
    worksheet.write('AN2','METASNP',headerFormat2) ; worksheet.write('AN3','',headerFormat2)
    worksheet.write('AO2','CAROL',headerFormat2) ; worksheet.write('AO3','',headerFormat2)
    worksheet.write('AP2','Condel',headerFormat2) ; worksheet.write('AP3','',headerFormat2)
    worksheet.write('AQ2','COVEC',headerFormat2) ; worksheet.write('AQ3','',headerFormat2)
    worksheet.write('AR2','MtoolBox',headerFormat2) ; worksheet.write('AR3','',headerFormat2)
    worksheet.write('AS2','APOGEE',headerFormat2) ; worksheet.write('AS3','',headerFormat2)
    worksheet.write('AT2','MutTaster',headerFormat2) ; worksheet.write('AT3','',headerFormat2)
    worksheet.write('AU2','Mitoclass',headerFormat2) ; worksheet.write('AU3','',headerFormat2)
    worksheet.merge_range('AV1:BE1','Conservation percent',mergeFormat1)
    worksheet.write('AV2','Euk',headerFormat1) ; worksheet.write('AV3','',headerFormat1)
    worksheet.write('AW2','Met',headerFormat1) ; worksheet.write('AW3','',headerFormat1)
    worksheet.write('AX2','Bil',headerFormat1) ; worksheet.write('AX3','',headerFormat1)
    worksheet.write('AY2','Ver',headerFormat1) ; worksheet.write('AY3','',headerFormat1)
    worksheet.write('AZ2','Tet',headerFormat1) ; worksheet.write('AZ3','',headerFormat1)
    worksheet.write('BA2','Amn',headerFormat1) ; worksheet.write('BA3','',headerFormat1)
    worksheet.write('BB2','Mam',headerFormat1) ; worksheet.write('BB3','',headerFormat1)
    worksheet.write('BC2','Eut',headerFormat1) ; worksheet.write('BC3','',headerFormat1)
    worksheet.write('BD2','Eua',headerFormat1) ; worksheet.write('BD3','',headerFormat1)
    worksheet.write('BE2','Pri',headerFormat1) ; worksheet.write('BE3','',headerFormat1)
    worksheet.merge_range('BF1:BO1','Jensen-shannon divergence score',mergeFormat2)
    worksheet.write('BF2','Euk',headerFormat2) ; worksheet.write('BF3','',headerFormat2)
    worksheet.write('BG2','Met',headerFormat2) ; worksheet.write('BG3','',headerFormat2)
    worksheet.write('BH2','Bil',headerFormat2) ; worksheet.write('BH3','',headerFormat2)
    worksheet.write('BI2','Ver',headerFormat2) ; worksheet.write('BI3','',headerFormat2)
    worksheet.write('BJ2','Tet',headerFormat2) ; worksheet.write('BJ3','',headerFormat2)
    worksheet.write('BK2','Amn',headerFormat2) ; worksheet.write('BK3','',headerFormat2)
    worksheet.write('BL2','Mam',headerFormat2) ; worksheet.write('BL3','',headerFormat2)
    worksheet.write('BM2','Eut',headerFormat2) ; worksheet.write('BM3','',headerFormat2)
    worksheet.write('BN2','Eua',headerFormat2) ; worksheet.write('BN3','',headerFormat2)
    worksheet.write('BO2','Pri',headerFormat2) ; worksheet.write('BO3','',headerFormat2)
    # Create data
    dicoData = {}
    lstHG = dicoHaplo['hg_found'].split(",")
    row = 3
    for variant in VCF(dicoInit["pathVCFout"].replace(".gz","")):
        lstVepAnno = variant.INFO.get('CSQ').split(",")
        for i in range(len(variant.ALT)):
            dicoData[row] = {}
            dicoData[row][0] = variant.start+1
            dicoData[row][1] = variant.REF
            dicoData[row][2] = variant.ALT[i]
            dicoData[row][3] = str(variant.FILTER).replace("None","PASS").capitalize()
            dicoData[row][4] = max(variant.INFO.get('HP').split(","))
            dicoData[row][5] = round(float(variant.QUAL))
            try: af = round(float(variant.format('AF')[0]),3)
            except: af = round(float(variant.format('AF')[0][i]),3)
            dicoData[row][6] = af
            dicoData[row][7] = variant.format('DP')[0][0]
            # VEP Annotation (Allele|Consequence|IMPACT|SYMBOL|Protein_position|Amino_acids")
            splitVep = lstVepAnno[i].split("|")
            dicoData[row][8] = splitVep[1]
            dicoData[row][9] = splitVep[2].capitalize()
            if splitVep[3]=="": dicoData[row][10] = "."
            else: dicoData[row][10] = splitVep[3]
            if splitVep[4]=="" or splitVep[1]=="FS": dicoData[row][11] = "."
            elif "/" in splitVep[5]: dicoData[row][11] = splitVep[5].split("/")[0]+splitVep[4]+splitVep[5].split("/")[1]
            else:dicoData[row][11] = splitVep[5]+splitVep[4]
            if str(variant.start+1)+variant.ALT[i] in lstHG: dicoData[row][12] = "Haplo"
            else: dicoData[row][12] = "."
            # MITOMAP Annotation
            if variant.INFO.get('MMDac')!=None:
                if type(variant.INFO.get('MMDac'))==tuple: count = variant.INFO.get('MMDac')[i]
                else: count = variant.INFO.get('MMDac')
            else:
                if type(variant.INFO.get('MMPac'))==tuple: count = variant.INFO.get('MMPac')[i]
                else: count = variant.INFO.get('MMPac')
            if count==None: count = 0
            try: dicoData[row][13] = int(count)
            except: dicoData[row][13] = 0
            if variant.INFO.get('MMDaf')!=None:
                if type(variant.INFO.get('MMDaf'))==tuple: freq = variant.INFO.get('MMDaf')[i]
                else: freq = variant.INFO.get('MMDaf')
            else:
                if type(variant.INFO.get('MMPaf'))==tuple: freq = variant.INFO.get('MMPaf')[i]
                else: variant.INFO.get('MMPaf')
            try: dicoData[row][14] = round(freq,1)
            except: dicoData[row][14] = 0.0
            if type(variant.INFO.get('MMDds'))==tuple: dicoData[row][15] = str(variant.INFO.get('MMDds')[i]).replace("None",".")
            else: dicoData[row][15] = str(variant.INFO.get('MMDds')).replace("None",".")
            if type(variant.INFO.get('MMDd'))==tuple: dicoData[row][16] = str(variant.INFO.get('MMDd')[i]).replace("None",".")          
            else: dicoData[row][16] = str(variant.INFO.get('MMDd')).replace("None",".")
            if type(variant.INFO.get('MMDhom'))==tuple: dicoData[row][17] = str(variant.INFO.get('MMDhom')[i]).replace("None",".")
            else: dicoData[row][17] = str(variant.INFO.get('MMDhom')).replace("None",".")
            if type(variant.INFO.get('MMDhtz'))==tuple: dicoData[row][18] = str(variant.INFO.get('MMDhtz')[i]).replace("None",".")
            else: dicoData[row][18] = str(variant.INFO.get('MMDhtz')).replace("None",".")
            # ClinVar Annotation
            dicoData[row][19] = str(variant.INFO.get('CVid')).replace("None",".")
            dicoData[row][20] = str(variant.INFO.get('CVdn')).replace("None",".").replace("|","\n").replace(",_",",").replace(",","\n")
            dicoData[row][21] = str(variant.INFO.get('CVrev')).replace("None",".").replace(",_",",").replace(",","\n")
            dicoData[row][22] = str(variant.INFO.get('CVsig')).replace("None",".")
            # HMTDB (for deletion make mean)
            hmtdbh = sum([float(i) for i in variant.INFO.get('HMTDBh').split(",")])/len(variant.INFO.get('HMTDBh').split(","))
            hmtdbp = sum([float(i) for i in variant.INFO.get('HMTDBp').split(",")])/len(variant.INFO.get('HMTDBp').split(","))
            dicoData[row][23] = round(hmtdbh,3)
            dicoData[row][24] = round(hmtdbp,3)
            # MitoTIP
            if variant.INFO.get('MTs')==None: dicoData[row][25] = "."
            else: dicoData[row][25] = round(float(variant.INFO.get('MTs')),2)
            if variant.INFO.get('MTq')=="Q1": rank = "D+"
            elif variant.INFO.get('MTq')=="Q2": rank = "D"
            elif variant.INFO.get('MTq')=="Q3": rank = "N"
            elif variant.INFO.get('MTq')=="Q4": rank = "N+"
            else: rank = "."
            dicoData[row][26] = rank
            # MitImpact
            col = 27
            for tool in ["P2","SIFT","FATW","FAT","PROV","MUTASS","EFINSP","EFINHD","CADD","PANTHER","PHDSNP","SNAP","METASNP","CAROL","CONDEL","COVEC","MTOOLBOX","APOGEE","MTASTER","MTCLASS"]:
                dicoData[row][col] = str(variant.INFO.get('P2')).replace("None",".")
                col+=1
            # Conservation & JS (for deletion make mean)
            for kermit in ["CONSEUK","CONSMET","CONSBIL","CONSVER","CONSTET","CONSAMN","CONSMAM","CONSEUT","CONSEUA","CONSPRI","JSEUK","JSMET","JSBIL","JSVER","JSTET","JSAMN","JSMAM","JSEUT","JSEUA","JSPRI"]:
                if variant.INFO.get(kermit)==None: value = 0.0
                else:
                    setValue = set(variant.INFO.get(kermit).split(","))
                    try: setValue.remove(None)
                    except: pass
                    value = max(setValue)
                dicoData[row][col] = value
                col+=1
            row+=1
    # Adjust column width
    dico_witdh={0:5,1:10,2:10,3:10,4:2,5:6,6:5,7:6,8:10,9:8,10:7,11:5,12:15,13:5,14:4,15:12,16:20,17:4,18:4,19:6,20:30,21:25,22:20,23:7,24:7,25:5,26:4,27:5,28:4,29:8,30:7,31:7,32:6,33:7,34:7,35:4,36:7,37:7,38:4,39:7,40:5,41:6,42:5,43:8,44:6,45:9,46:9,47:3,48:3,49:3,50:3,51:3,52:4,53:4,54:3,55:3,56:3,57:4,58:4,59:4,60:4,61:4,62:4,63:4,64:4,65:4,66:4}
    for col in dico_witdh: worksheet.set_column(col,col,dico_witdh[col])
    # Adjust row height
    for i in range(3,row,1):
        row_height = 0
        for j in dico_witdh:
            row_height = max(row_height,len(str(dicoData[i][j]))/dico_witdh[j]+(len(str(dicoData[i][j])) % dico_witdh[j] > 0))
        if row_height>1: worksheet.set_row(i,(row_height)*15)
    # Write data
    for i in range(3,row,1):
        for j in range(col+1): worksheet.write(i,j,dicoData[i][j],rowFormat)
    # Freeze pane on the top row.
    worksheet.freeze_panes(3,3)    
    # autofilter
    worksheet.autofilter(2, 0, row, col)
    # Save .xls
    workbook.close()    