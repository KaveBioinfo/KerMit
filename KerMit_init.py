import os
import shutil
import subprocess
import datetime                                                                                                                                                                              
from KerMit_display import *
from KerMit_run import *



#***** DicoInit *****#
def init(lstArgv):
    dicoInit = { "startTime":datetime.datetime.now(), \
                 "pathDirData":os.path.dirname(os.path.realpath(__file__))+"/data", \
                 "pathDirSrc":os.path.dirname(os.path.realpath(__file__))+"/src", \
                 "pathFasta":os.path.dirname(os.path.realpath(__file__))+"/data/chrM.fasta", \
                 "pathVcfvalidator":os.path.dirname(os.path.realpath(__file__))+"/src/vcf_validator", \
                 "pathHaplogrep":os.path.dirname(os.path.realpath(__file__))+"/src/haplogrep", \
                 "pathVcfanno":os.path.dirname(os.path.realpath(__file__))+"/src/vcfanno", \
                 "pathTmp":"/tmp/kermit", \
                 "update":False, \
                 "bgzip":True, \
                 "nbThread":1 \
                }
    if "--color false" in str(lstArgv).lower(): DicoInit["colorBool"] = False
    else: dicoInit["colorBool"] = True
    return dicoInit



#***** Check files reuired by KerMit *****#
def checkKermitFiles(dicoInit):
    if not os.path.isfile(dicoInit["pathFasta"]): printError("KerMit files","`chrM.fasta` file not found", dicoInit["colorBool"])
    if not cmd_exists(dicoInit["pathVcfvalidator"]): printError("KerMit files","`vcf_validator` executable not found", dicoInit["colorBool"])
    if not cmd_exists(dicoInit["pathHaplogrep"]): printError("KerMit files","`haplogrep` executable not found", dicoInit["colorBool"])
    if not cmd_exists(dicoInit["pathVcfanno"]): printError("KerMit files","`vcfanno` executable not found", dicoInit["colorBool"])
    if dicoInit["bgzip"]==True:
        if not cmd_exists("bgzip"): printError("KerMit files","`bgzip` executable not found (please install)", dicoInit["colorBool"])
        if not cmd_exists("tabix"): printError("KerMit files","`tabix` executable not found (please install)", dicoInit["colorBool"])
    if not os.path.isdir(dicoInit["pathDirData"]+"/homo_sapiens/100_GRCh38/MT"): printError("KerMit files","`vep_cache` folder not found", dicoInit["colorBool"])


#***** Check files reuired by KerMit *****#
def getArgs(lstArgv,dicoInit):
    if "--help" in lstArgv or "-h" in lstArgv or len(lstArgv)==1: displayUsage(dicoInit["colorBool"])
    else:
        for i in range(1,len(lstArgv),2):
            # Check input VCF
            if lstArgv[i]=="--in":
                try: dicoInit["pathVCFin"] = lstArgv[i+1]
                except: printError("KerMit input","Missing input VCF value",dicoInit["colorBool"])
                if not os.path.isfile(dicoInit["pathVCFin"]): printError("KerMit input","Input VCF not found `"+dicoInit["pathVCFin"]+"`",dicoInit["colorBool"])
            # Check output VCF
            elif lstArgv[i]=="--out":
                try: dicoInit["pathVCFout"] = lstArgv[i+1]
                except: printError("KerMit input","Missing output VCF value",dicoInit["colorBool"])
                if not os.access(os.path.dirname(os.path.abspath(dicoInit["pathVCFout"])), os.W_OK): printError("KerMit input","Output VCF folder not writable `"+os.path.dirname(os.path.abspath(dicoInit["pathVCFout"]))+"`",dicoInit["colorBool"])
            elif lstArgv[i]=="--vep":
                try: dicoInit["pathVEP"] = lstArgv[i+1]
                except: printError("KerMit input","Missing VEP executable value",dicoInit["colorBool"])
                if not cmd_exists(dicoInit["pathVEP"]): printError("KerMit input","Invalid VEP executable `"+dicoInit["pathVEP"]+"`",dicoInit["colorBool"])
                try: dicoInit["VEPversion"] = subprocess.check_output(dicoInit["pathVEP"]+" | grep \"ensembl-vep\" | cut -d \":\" -f 2 | sed s/\" \"//g",stderr=subprocess.STDOUT, shell=True).decode().replace("\n",",")[:-1]
                except: printError("KerMit input","Unable to check VEP version `"+dicoInit["pathVEP"]+"`",dicoInit["colorBool"])
            # Check output XLS
            elif lstArgv[i]=="--xls":
                try: dicoInit["pathXLSXout"] = lstArgv[i+1]
                except: printError("KerMit input","Missing output XLSX value",dicoInit["colorBool"])
                if not os.access(os.path.dirname(os.path.abspath(dicoInit["pathXLSXout"])), os.W_OK): printError("KerMit input","Output XLSX folder not writable `"+os.path.dirname(os.path.abspath(dicoInit["pathXLSXout"]))+"`",dicoInit["colorBool"])
            elif lstArgv[i]=="--update":
                try: value = lstArgv[i+1]
                except: printError("KerMit input","Missing --update boolean value",dicoInit["colorBool"])
                if value.lower()=="true": dicoInit["update"] = True
                elif value.lower()=="false": dicoInit["update"] = False
                else: printError("KerMit input","Invalid `--update` argument (required true or false)",dicoInit["colorBool"])
            elif lstArgv[i]=="--bgzip":
                try: value = lstArgv[i+1]
                except: printError("KerMit input","Missing --bgzip boolean value",dicoInit["colorBool"])
                if value.lower()=="true": dicoInit["bgzip"] = True
                elif value.lower()=="false": dicoInit["bgzip"] = False
                else: printError("KerMit input","Invalid `--bgzip` argument (required true or false)",dicoInit["colorBool"])
            elif lstArgv[i]=="--tmp":
                try: dicoInit["pathTmp"] = lstArgv[i+1]+"/kermit"
                except: printError("KerMit input","Missing --tmp path value",dicoInit["colorBool"])
                if os.path.isdir(dicoInit["pathTmp"]): shutil.rmtree(dicoInit["pathTmp"])
                if not os.access(os.path.dirname(dicoInit["pathTmp"]), os.W_OK): printError("KerMit input","Temp folder not writable `"+os.path.dirname(os.path.abspath(dicoInit["pathTmp"]))+"`",dicoInit["colorBool"])
                os.makedirs(dicoInit["pathTmp"])
            elif lstArgv[i]=="--thread":
                try: value = lstArgv[i+1]
                except: printError("KerMit input","Missing thread value",dicoInit["colorBool"])
                try: dicoInit["nbThread"] = int(value)
                except: printError("KerMit input","Invalid `--thread` argument (required integer)",dicoInit["colorBool"])
    # Check input VCF
    if not "pathVCFin" in dicoInit: printError("KerMit input","Missing `--in` mandatory argument",dicoInit["colorBool"])
    if not "pathVCFout" in dicoInit: printError("KerMit input","Missing `--out` mandatory argument",dicoInit["colorBool"])
    if not "pathVEP" in dicoInit: printError("KerMit input","Missing `--vep` mandatory argument",dicoInit["colorBool"])
    # Check output vcf/xls extension  
    if dicoInit["bgzip"]==True and not dicoInit["pathVCFout"][-6:]=="vcf.gz": printError("KerMit input","Invalid output VCF extension with bgzip mode (required `.vcf.gz`)",dicoInit["colorBool"])
    if dicoInit["bgzip"]==False and not dicoInit["pathVCFout"][-4:]==".vcf": printError("KerMit input","Invalid output VCF extension (required `.vcf`)",dicoInit["colorBool"])
    if "pathXLSXout" in dicoInit and not dicoInit["pathXLSXout"][-5:]==".xlsx": printError("KerMit input","Invalid output XLSX extension (required `.xlsx`)",dicoInit["colorBool"])
    # Validate input VCF
    boolvalid,lst_errors = validateVCF(dicoInit,dicoInit["pathVCFin"])
    if boolvalid==False: printError("KerMit input","Invalid VCF input file `"+os.path.basename(dicoInit["pathVCFin"])+"`\n  "+"\n  ".join(lst_errors),dicoInit["colorBool"])
    # Temporary files
    dicoInit["pathVEPvcf"] = dicoInit["pathTmp"]+"/vep.vcf"
    dicoInit["pathVEPformatedvcf"] = dicoInit["pathTmp"]+"/vep_formated.vcf"
    dicoInit["pathHaplogrepOut"] = dicoInit["pathTmp"]+"/hablogrep.txt"
