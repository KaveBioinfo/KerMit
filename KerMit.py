#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : KerMit.py                        |    
# description     : Annotate mitochondrial VCF       |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200712                         |
# version         : 0.2                             |
# python_version  : 3.8.2                            |
#==========================================================================================================
# USAGE: KerMit.py pathVCF
# OUT  : 
#==========================================================================================================
import sys
import shutil
from KerMit_init import *
from KerMit_display import *
from KerMit_data import *
from KerMit_run import *



#***** Arguments *****#
# Init
dicoInit = init(sys.argv)
# HeaderAbout
displayAbout(dicoInit["colorBool"])
# Check required files
checkKermitFiles(dicoInit)
# Check input args
getArgs(sys.argv,dicoInit)
# Digest
displayDigest(dicoInit)



#***** DATA files *****#
printcolor("\n 🅳 🅰 🆃 🅰\n","0","greenlighter",False,dicoInit["colorBool"])
check_and_download(dicoInit)
dicoDataVersion = get_data_version(dicoInit)
# Data display
displayData(dicoDataVersion,dicoInit["colorBool"])


# Execution header
printcolor("\n 🅴 🆇 🅴 🅲 🆄 🆃 🅸 🅾 🅽\n","0","white",False,dicoInit["colorBool"])

#***** VEP Annotation *****#
printcolor(" VEP annot     :","0","white",False,dicoInit["colorBool"])
launchVep(dicoInit)
printcolor(" done\n","0","white",False,dicoInit["colorBool"])



#***** FORMAT VEP vcf *****#
printcolor(" VEP reformat  :","0","white",False,dicoInit["colorBool"])
vep_reformat(dicoInit)
printcolor(" done\n","0","white",False,dicoInit["colorBool"])



#***** ADD haplogroup *****#
printcolor(" Haplogrep     :","0","white",False,dicoInit["colorBool"])
dicoHaplo = haplogrep_annotate(dicoInit)
printcolor(" done\n","0","white",False,dicoInit["colorBool"])



#***** ADD external annotation (vcfanno) *****#
printcolor(" Vcfanno       :","0","white",False,dicoInit["colorBool"])
vcfanno(dicoInit)
printcolor(" done\n","0","white",False,dicoInit["colorBool"])



#***** CREATE excel  *****#
printcolor(" Output XLSX   :","0","white",False,dicoInit["colorBool"])
makeXLSX(dicoInit,dicoHaplo)
printcolor(" done\n","0","white",False,dicoInit["colorBool"])



#***** CREATE VCF *****# 
printcolor(" Output VCF    :","0","white",False,dicoInit["colorBool"])
reheaderVCF(dicoInit) # FORMAT VCF header
if dicoInit["bgzip"]==True: bgzipTabix(dicoInit) # Bgzip/tabix if required
printcolor(" done\n","0","white",False,dicoInit["colorBool"])



# Clean temp folder
printcolor(" Clean temp    :","0","white",False,dicoInit["colorBool"])
shutil.rmtree(dicoInit["pathTmp"])
printcolor(" done\n\n","0","white",False,dicoInit["colorBool"])
