import os
import sys
import shutil
import subprocess
import datetime
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
        # return results for xlsx output
        return {'hg':haplogroup, 'qual':quality,'hg_found':found_snps}


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