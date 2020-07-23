import socket
import getpass
import datetime
import operator                                                                                                                                                                          


#***** RGB COLOR PRINT *****#
def printcolor(text,style,fgRGB,bgRGB,colorBool):
    # PREDEFINED color
    dicoColor = { "white":"211,215,207", "red":"255;85;85", "green":"113;180;120", "greenlight":"156;207;161", "greenlighter":"205;222;135" }
    if fgRGB in dicoColor: fgRGB = dicoColor[fgRGB]
    # STYLE : normal (0),bold (1),light (2),italic (3),underline (4),slow blink (5),rapid blink (6), crossed-out (9)
    if colorBool==True:
        if bgRGB: print("\x1b["+style+";38;2;"+fgRGB+";48;2;"+bgRGB+"m"+text+"\x1b[0m",end='')
        else: print("\x1b["+style+";38;2;"+fgRGB+"m"+text+"\x1b[0m",end='')
    else: print(text,end='')



#***** RGB COLOR PRINT *****#
def printError(errorTag,errorStr,colorBool):
    printcolor("\n 🅴 🆁 🆁 🅾 🆁\n ["+errorTag+"] "+errorStr+"\n","0","red",False,colorBool)
    exit("\n")


#***** USAGE DISPLAY *****#
def displayUsage(colorBool):
    printcolor("\n 🆄 🆂 🅰 🅶 🅴\n","0","white",False,colorBool)
    printcolor(" The typical command for running KerMit is as follows:\n\n","0","white",False,colorBool)
    printcolor(" python KerMit.py [options] --in path_input_vcf --out path_output_vcf --vep path_vep\n","0","white",False,colorBool)     
    printcolor("\n Mandatory arguments:\n","0","white",False,colorBool)
    printcolor("   --in                        Path to input VCF file\n","0","white",False,colorBool)
    printcolor("   --out                       Path to output VCF file\n","0","white",False,colorBool)
    printcolor("   --vep                       Path to VEP executable\n","0","white",False,colorBool)
    printcolor("\n Options arguments:\n","0","white",False,colorBool)
    printcolor("   --xlsx                      Path to output XLSX file\n","0","white",False,colorBool)    
    printcolor("   --update                    Update source files (MITOMAP, HmtDB, ClinVar, MitoTip) [false]\n","0","white",False,colorBool)
    printcolor("   --bgzip                     bgzip and tabix output vcf [true]\n","0","white",False,colorBool)
    printcolor("   --tmp                       Path to temporary folder [/tmp]\n","0","white",False,colorBool)
    printcolor("   --thread                    Number of threads to use [1]\n","0","white",False,colorBool)
    printcolor("   --help                      Print help\n","0","white",False,colorBool)
    printcolor("   --color                     Enable color output [true]\n","0","white",False,colorBool)
    exit("\n")


#***** Header *****#
def displayAbout(colorBool):
    # Header
    printcolor("\nＫｅｒＭｉｔ 🐸\n","1","green",False,colorBool)
    # About
    printcolor("\n 🅰 🅱 🅾 🆄 🆃\n","0","green",False,colorBool)
    printcolor(" Description   : Mitochondrial VCF annotation\n","0","green",False,colorBool)
    printcolor(" Version       : v0.1\n","0","green",False,colorBool)
    printcolor(" Copyright     : KAVE - CHU Angers\n","0","green",False,colorBool)
    printcolor(" Support       : https://github.com/dooguypapua/KerMit\n","0","green",False,colorBool)


#***** Digest *****#
def displayDigest(dicoInit):
    printcolor("\n 🅳 🅸 🅶 🅴 🆂 🆃\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" Date          : "+dicoInit["startTime"].ctime()+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" Hostname      : "+socket.gethostname()+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" User          : "+getpass.getuser()+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" VEPversion    : "+dicoInit["VEPversion"]+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" Threads       : "+str(dicoInit["nbThread"])+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" Input VCF     : "+dicoInit["pathVCFin"]+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" Output VCF    : "+dicoInit["pathVCFout"]+"\n","0","greenlight",False,dicoInit["colorBool"])
    if "pathXLSXout" in dicoInit: printcolor(" Output XLSX   : "+dicoInit["pathXLSXout"]+"\n","0","greenlight",False,dicoInit["colorBool"])
    printcolor(" Temp          : "+dicoInit["pathTmp"]+"\n","0","greenlight",False,dicoInit["colorBool"])



def displayData(dicoDataVersion,colorBool):
    sorted_dict = dict(sorted(dicoDataVersion.items(), key=operator.itemgetter(0)))
    for key in sorted_dict:
        printcolor(" "+key.ljust(14)+": "+dicoDataVersion[key]+"\n","0","greenlighter",False,colorBool)
    printcolor(" "+"MitoKons".ljust(14)+": 2019-07-16\n","0","greenlighter",False,colorBool)
