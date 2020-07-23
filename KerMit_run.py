import os
import sys
import shutil
import subprocess
import datetime
import xlsxwriter
from math import ceil
from cyvcf2 import VCF
from KerMit_init import *
from KerMit_display import *
from KerMit_data import *




#***** Check command exists *****#
def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0



def auto_format_cell_width(ws):
    for letter in range(1,ws.max_column):
        maximum_value = 0
        for cell in ws[get_column_letter(letter)]:
            val_to_check = len(str(cell.value))
            if val_to_check > maximum_value:
               maximum_value = val_to_check
        ws.column_dimensions[get_column_letter(letter)].width = maximum_value + 1



#***** HAPLOGREP ANNOTATE *****#
def launchVep(dicoInit):
    # Disable downstream/upstream with "--distance 0"
    cmd_vep = dicoInit["pathVEP"]+" --offline --quiet --force_overwrite --use_transcript_ref --check_ref --no_stats "\
              "--cache --dir_cache "+dicoInit["pathDirData"]+" --species homo_sapiens --assembly GRCh38 --fasta "+dicoInit["pathFasta"]+" "\
              "--input_file "+dicoInit["pathVCFin"]+" --output_file "+dicoInit["pathVEPvcf"]+" --vcf --fork "+str(dicoInit["nbThread"])+" --distance 0"
    process = subprocess.Popen([cmd_vep], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: printError("KerMit annotate","VEP Annotation failed\n    "+err.decode('utf-8'),dicoInit["colorBool"])



#***** HAPLOGREP ANNOTATE *****#
def haplogrep_annotate(dicoInit):
    cmd_haplogrep = dicoInit["pathHaplogrep"]+" classify --in "+dicoInit["pathVCFin"]+" --out "+dicoInit["pathHaplogrepOut"]+" --format vcf --hetLevel 0.8 --phylotree 17 --extend-report"
    process = subprocess.Popen([cmd_haplogrep], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: printError("KerMit annotate","Haplogrep classify failed\n    "+err.decode('utf-8'),dicoInit["colorBool"])
    IN = open(dicoInit["pathHaplogrepOut"],'r')
    data = IN.read().replace("\"","")
    IN.close()
    if "ERROR" in data:
        printError("KerMit annotate","Haplogrep error\n"+data,dicoInit["colorBool"])
    else:
        split_line = data.split("\n")[1].split("\t")
        haplogroup = split_line[2]
        rank = split_line[3]
        quality = split_line[4]
        notfound_snps = split_line[5].replace(" ",",")
        found_snps = split_line[6].replace(" ",",")
        addLines = "##INFO=<ID=HAPLO,Number=0,Type=Flag,Description=\\\"Variant used to define haplogroup\\\">\\n"
        addLines+= "##haplogrep_version=\\\"v2.2.5\\\"\\n"
        addLines+= "##haplogrep_phylotree=\\\"17\\\"\\n"
        addLines+= "##haplogrep_haplogroup=\\\""+haplogroup+"\\\"\\n"
        addLines+= "##haplogrep_rank="+rank+"\\n"
        addLines+= "##haplogrep_quality="+quality+"\\n"
        addLines+= "##haplogrep_notfound=\\\""+notfound_snps+"\\\"\\n"
        addLines+= "##haplogrep_found=\\\""+found_snps+"\\\"\\n"
        cmd_sed = "sed -i s/\"#CHROM\"/\""+addLines+"#CHROM\"/g "+dicoInit["pathVEPformatedvcf"]
        process = subprocess.Popen([cmd_sed], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()
        if process.returncode!=0: printError("KerMit annotate","Haplogrep add to VCF failed\n    "+err.decode('utf-8'),dicoInit["colorBool"])



#***** CREATE toml & launch vcfanno *****#
def vcfanno(dicoInit):
    # Create conf toml
    pathToml = dicoInit["pathTmp"]+"/conf.toml"
    TOML = open(pathToml,'w')
    TOML.write("[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/clinvar_MT.vcf.gz\"\n"\
               "fields = [\"ALLELEID\",\"CLNDN\",\"CLNREVSTAT\",\"CLNSIG\"]\n"\
               "ops=[\"self\",\"self\",\"self\",\"self\"]\n"\
               "names=[\"CVid\",\"CVdn\",\"CVrev\",\"CVsig\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/MITOTIP_scores.vcf.gz\"\n"\
               "fields = [\"MTs\",\"MTq\"]\n"\
               "ops=[\"self\",\"self\"]\n"\
               "names=[\"MTs\",\"MTq\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/HMTDB_var_pa_tot.tsv.gz\"\n"\
               "columns = [3]\n"\
               "ops=[\"self\"]\n"\
               "names=[\"HMTDBp\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/HMTDB_var_tot.tsv.gz\"\n"\
               "columns = [3]\n"\
               "ops=[\"self\"]\n"\
               "names=[\"HMTDBh\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/MITOMAP_disease.vcf.gz\"\n"\
               "fields = [\"AC\",\"AF\",\"homoplasmy\",\"heteroplasmy\",\"Disease\",\"DiseaseStatus\",\"HGFL\"]\n"\
               "ops=[\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\"]\n"\
               "names=[\"MMDac\",\"MMDaf\",\"MMDhom\",\"MMDhtz\",\"MMDd\",\"MMDds\",\"MMDhg\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/MITOMAP_polymorphisms.vcf.gz\"\n"\
               "fields = [\"AC\",\"AF\"]\n"\
               "ops=[\"self\",\"self\"]\n"\
               "names=[\"MMPac\",\"MMPaf\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/MitImpact_db.vcf.gz\"\n"\
               "fields = [\"P2\",\"SIFT\",\"FATW\",\"FAT\",\"PROV\",\"MUTASS\",\"EFINSP\",\"EFINHD\",\"CADD\",\"PANTHER\",\"PHDSNP\",\"SNAP\",\"METASNP\",\"CAROL\",\"CONDEL\",\"COVEC\",\"MTOOLBOX\",\"APOGEE\",\"MTASTER\",\"MTCLASS\"]\n"\
               "ops=[\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\"]\n"\
               "names=[\"P2\",\"SIFT\",\"FATW\",\"FAT\",\"PROV\",\"MUTASS\",\"EFINSP\",\"EFINHD\",\"CADD\",\"PANTHER\",\"PHDSNP\",\"SNAP\",\"METASNP\",\"CAROL\",\"CONDEL\",\"COVEC\",\"MTOOLBOX\",\"APOGEE\",\"MTASTER\",\"MTCLASS\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/chrM_HP.bed.gz\"\n"\
               "columns = [4]\n"\
               "ops=[\"self\"]\n"\
               "names=[\"HP\"]\n"\
               "[[annotation]]\n"\
               "file=\""+dicoInit["pathDirData"]+"/MitoKons.bed.gz\"\n"\
               "columns = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]\n"\
               "ops=[\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\"]\n"\
               "names=[\"CONSEUK\",\"CONSMET\",\"CONSBIL\",\"CONSVER\",\"CONSTET\",\"CONSAMN\",\"CONSMAM\",\"CONSEUT\",\"CONSEUA\",\"CONSPRI\",\"JSEUK\",\"JSMET\",\"JSBIL\",\"JSVER\",\"JSTET\",\"JSAMN\",\"JSMAM\",\"JSEUT\",\"JSEUA\",\"JSPRI\"]"
               )
    TOML.close()
    # Launch vcfanno
    cmd_vcfanno = dicoInit["pathVcfanno"]+" "+pathToml+" "+dicoInit["pathVEPformatedvcf"]+" > "+dicoInit["pathVCFout"].replace(".gz","")
    process = subprocess.Popen([cmd_vcfanno], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: printError("KerMit annotate","vcfanno\n    "+err.decode('utf-8'),dicoInit["colorBool"])



#***** CREATE EXCEL OUTPUT *****#
def makeXLSX(dicoInit):
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
    worksheet.merge_range('I1:L1','ANNOTATION',mergeFormat1)
    worksheet.write('I2','Csq',headerFormat1) ; worksheet.write('I3','',headerFormat1)
    worksheet.write('J2','Impact',headerFormat1) ; worksheet.write('J3','',headerFormat1)
    worksheet.write('K2','Gene',headerFormat1) ; worksheet.write('K3','',headerFormat1)
    worksheet.write('L2','AA',headerFormat1) ; worksheet.write('L3','',headerFormat1)
    worksheet.merge_range('M1:R1','MITOMAP',mergeFormat2)
    worksheet.write('M2','Count',headerFormat2) ; worksheet.write('M3','',headerFormat2)
    worksheet.write('N2','Freq',headerFormat2) ; worksheet.write('N3','',headerFormat2)
    worksheet.write('O2','Status',headerFormat2) ; worksheet.write('O3','',headerFormat2)
    worksheet.write('P2','Disease',headerFormat2) ; worksheet.write('P3','',headerFormat2)
    worksheet.write('Q2','Hom',headerFormat2) ; worksheet.write('Q3','',headerFormat2)
    worksheet.write('R2','Htz',headerFormat2) ; worksheet.write('R3','',headerFormat2)
    worksheet.merge_range('S1:V1','ClinVar',mergeFormat1)
    worksheet.write('S2','Id',headerFormat1) ; worksheet.write('S3','',headerFormat1)
    worksheet.write('T2','Disease',headerFormat1) ; worksheet.write('T3','',headerFormat1)
    worksheet.write('U2','Status',headerFormat1) ; worksheet.write('U3','',headerFormat1)
    worksheet.write('V2','Significance',headerFormat1) ; worksheet.write('V3','',headerFormat1)
    worksheet.merge_range('W1:X1','HMTDB Site Var',mergeFormat2)
    worksheet.write('W2','Healthy',headerFormat2) ; worksheet.write('W3','',headerFormat2)
    worksheet.write('X2','Disease',headerFormat2) ; worksheet.write('X3','',headerFormat2)
    worksheet.merge_range('Y1:Z1','MitoTip',mergeFormat1)
    worksheet.write('Y2','Score',headerFormat1) ; worksheet.write('Y3','',headerFormat1)
    worksheet.write('Z2','Rank',headerFormat1) ; worksheet.write('Z3','',headerFormat1)
    worksheet.merge_range('AA1:AT1','MitImpact',mergeFormat2)
    worksheet.write('AA2','Poly2',headerFormat2) ; worksheet.write('AA3','',headerFormat2)
    worksheet.write('AB2','SIFT',headerFormat2) ; worksheet.write('AB3','',headerFormat2)
    worksheet.write('AC2','FatHmmW',headerFormat2) ; worksheet.write('AC3','',headerFormat2)
    worksheet.write('AD2','FatHmm',headerFormat2) ; worksheet.write('AD3','',headerFormat2)
    worksheet.write('AE2','PROVEAN',headerFormat2) ; worksheet.write('AE3','',headerFormat2)
    worksheet.write('AF2','MutAss',headerFormat2) ; worksheet.write('AF3','',headerFormat2)
    worksheet.write('AG2','EFIN-SP',headerFormat2) ; worksheet.write('AG3','',headerFormat2)
    worksheet.write('AH2','EFIN-HD',headerFormat2) ; worksheet.write('AH3','',headerFormat2)
    worksheet.write('AI2','CADD',headerFormat2) ; worksheet.write('AI3','',headerFormat2)
    worksheet.write('AJ2','PANTHER',headerFormat2) ; worksheet.write('AJ3','',headerFormat2)
    worksheet.write('AK2','PhD-SNP',headerFormat2) ; worksheet.write('AK3','',headerFormat2)
    worksheet.write('AL2','SNAP',headerFormat2) ; worksheet.write('AL3','',headerFormat2)
    worksheet.write('AM2','METASNP',headerFormat2) ; worksheet.write('AM3','',headerFormat2)
    worksheet.write('AN2','CAROL',headerFormat2) ; worksheet.write('AN3','',headerFormat2)
    worksheet.write('AO2','Condel',headerFormat2) ; worksheet.write('AO3','',headerFormat2)
    worksheet.write('AP2','COVEC',headerFormat2) ; worksheet.write('AP3','',headerFormat2)
    worksheet.write('AQ2','MtoolBox',headerFormat2) ; worksheet.write('AQ3','',headerFormat2)
    worksheet.write('AR2','APOGEE',headerFormat2) ; worksheet.write('AR3','',headerFormat2)
    worksheet.write('AS2','MutTaster',headerFormat2) ; worksheet.write('AS3','',headerFormat2)
    worksheet.write('AT2','Mitoclass',headerFormat2) ; worksheet.write('AT3','',headerFormat2)
    worksheet.merge_range('AU1:BD1','Conservation percent',mergeFormat1)
    worksheet.write('AU2','Euk',headerFormat1) ; worksheet.write('AU3','',headerFormat1)
    worksheet.write('AV2','Met',headerFormat1) ; worksheet.write('AV3','',headerFormat1)
    worksheet.write('AW2','Bil',headerFormat1) ; worksheet.write('AW3','',headerFormat1)
    worksheet.write('AX2','Ver',headerFormat1) ; worksheet.write('AX3','',headerFormat1)
    worksheet.write('AY2','Tet',headerFormat1) ; worksheet.write('AY3','',headerFormat1)
    worksheet.write('AZ2','Amn',headerFormat1) ; worksheet.write('AZ3','',headerFormat1)
    worksheet.write('BA2','Mam',headerFormat1) ; worksheet.write('BA3','',headerFormat1)
    worksheet.write('BB2','Eut',headerFormat1) ; worksheet.write('BB3','',headerFormat1)
    worksheet.write('BC2','Eua',headerFormat1) ; worksheet.write('BC3','',headerFormat1)
    worksheet.write('BD2','Pri',headerFormat1) ; worksheet.write('BD3','',headerFormat1)
    worksheet.merge_range('BE1:BN1','Jensen-shannon divergence score',mergeFormat2)
    worksheet.write('BE2','Euk',headerFormat2) ; worksheet.write('BE3','',headerFormat2)
    worksheet.write('BF2','Met',headerFormat2) ; worksheet.write('BF3','',headerFormat2)
    worksheet.write('BG2','Bil',headerFormat2) ; worksheet.write('BG3','',headerFormat2)
    worksheet.write('BH2','Ver',headerFormat2) ; worksheet.write('BH3','',headerFormat2)
    worksheet.write('BI2','Tet',headerFormat2) ; worksheet.write('BI3','',headerFormat2)
    worksheet.write('BJ2','Amn',headerFormat2) ; worksheet.write('BJ3','',headerFormat2)
    worksheet.write('BK2','Mam',headerFormat2) ; worksheet.write('BK3','',headerFormat2)
    worksheet.write('BL2','Eut',headerFormat2) ; worksheet.write('BL3','',headerFormat2)
    worksheet.write('BM2','Eua',headerFormat2) ; worksheet.write('BM3','',headerFormat2)
    worksheet.write('BN2','Pri',headerFormat2) ; worksheet.write('BN3','',headerFormat2)
    # Create data
    dicoData = {}
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
            # MITOMAP Annotation
            if variant.INFO.get('MMDac')!=None:
                if type(variant.INFO.get('MMDac'))==tuple: count = variant.INFO.get('MMDac')[i]
                else: count = variant.INFO.get('MMDac')
            else:
                if type(variant.INFO.get('MMPac'))==tuple: count = variant.INFO.get('MMPac')[i]
                else: count = variant.INFO.get('MMPac')
            if count==None: count = 0
            try: dicoData[row][12] = int(count)
            except: dicoData[row][12] = 0
            if variant.INFO.get('MMDaf')!=None:
                if type(variant.INFO.get('MMDaf'))==tuple: freq = variant.INFO.get('MMDaf')[i]
                else: freq = variant.INFO.get('MMDaf')
            else:
                if type(variant.INFO.get('MMPaf'))==tuple: freq = variant.INFO.get('MMPaf')[i]
                else: variant.INFO.get('MMPaf')
            try: dicoData[row][13] = round(freq,1)
            except: dicoData[row][13] = 0.0
            if type(variant.INFO.get('MMDds'))==tuple: dicoData[row][14] = str(variant.INFO.get('MMDds')[i]).replace("None",".")
            else: dicoData[row][14] = str(variant.INFO.get('MMDds')).replace("None",".")
            if type(variant.INFO.get('MMDd'))==tuple: dicoData[row][15] = str(variant.INFO.get('MMDd')[i]).replace("None",".")          
            else: dicoData[row][15] = str(variant.INFO.get('MMDd')).replace("None",".")
            if type(variant.INFO.get('MMDhom'))==tuple: dicoData[row][16] = str(variant.INFO.get('MMDhom')[i]).replace("None",".")
            else: dicoData[row][16] = str(variant.INFO.get('MMDhom')).replace("None",".")
            if type(variant.INFO.get('MMDhtz'))==tuple: dicoData[row][17] = str(variant.INFO.get('MMDhtz')[i]).replace("None",".")
            else: dicoData[row][17] = str(variant.INFO.get('MMDhtz')).replace("None",".")
            # ClinVar Annotation
            dicoData[row][18] = str(variant.INFO.get('CVid')).replace("None",".")
            dicoData[row][19] = str(variant.INFO.get('CVdn')).replace("None",".").replace("|","\n").replace(",_",",").replace(",","\n")
            dicoData[row][20] = str(variant.INFO.get('CVrev')).replace("None",".").replace(",_",",").replace(",","\n")
            dicoData[row][21] = str(variant.INFO.get('CVsig')).replace("None",".")
            # HMTDB (for deletion make mean)
            hmtdbh = sum([float(i) for i in variant.INFO.get('HMTDBh').split(",")])/len(variant.INFO.get('HMTDBh').split(","))
            hmtdbp = sum([float(i) for i in variant.INFO.get('HMTDBp').split(",")])/len(variant.INFO.get('HMTDBp').split(","))
            dicoData[row][22] = round(hmtdbh,3)
            dicoData[row][23] = round(hmtdbp,3)
            # MitoTIP
            if variant.INFO.get('MTs')==None: dicoData[row][24] = "."
            else: dicoData[row][24] = round(float(variant.INFO.get('MTs')),2)
            if variant.INFO.get('MTq')=="Q1": rank = "D+"
            elif variant.INFO.get('MTq')=="Q2": rank = "D"
            elif variant.INFO.get('MTq')=="Q3": rank = "N"
            elif variant.INFO.get('MTq')=="Q4": rank = "N+"
            else: rank = "."
            dicoData[row][25] = rank
            # MitImpact
            col = 26
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
    dico_witdh={0:5,1:10,2:10,3:10,4:2,5:6,6:5,7:6,8:10,9:8,10:7,11:5,12:5,13:4,14:12,15:20,16:4,17:4,18:6,19:30,20:25,21:20,22:7,23:7,24:5,25:4,26:5,27:4,28:8,29:7,30:7,31:6,32:7,33:7,34:4,35:7,36:7,37:4,38:7,39:5,40:6,41:5,42:8,43:6,44:9,45:9,46:3,47:3,48:3,49:3,50:3,51:4,52:4,53:3,54:3,55:3,56:4,57:4,58:4,59:4,60:4,61:4,62:4,63:4,64:4,65:4}
    for col in dico_witdh: worksheet.set_column(col,col,dico_witdh[col])
    # Adjust row height
    for i in range(3,row,1):
        row_height = 0
        for j in dico_witdh:
            row_height = max(row_height,ceil(len(str(dicoData[i][j]))/dico_witdh[j]))
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



#***** BGZIP & TABIX *****#
def bgzipTabix(dicoInit):
    cmd_bgzip = "bgzip -f "+dicoInit["pathVCFout"].replace(".gz","")
    process = subprocess.Popen([cmd_bgzip], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: printError("KerMit postprocessing","bgzip output\n    "+err.decode('utf-8'),dicoInit["colorBool"])
    # Tabix if required
    cmd_tabix = "tabix -f -p vcf "+dicoInit["pathVCFout"]
    process = subprocess.Popen([cmd_tabix], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if process.returncode!=0: printError("KerMit postprocessing","tabix output\n    "+err.decode('utf-8'),dicoInit["colorBool"])