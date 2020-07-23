![kermit logo](data/Kermit_black.png)


## Description

<b>KerMit annotate your mitochondrial VCF with many dedicated annotation tools and sources:</b><br/>
- VEP (http://www.ensembl.org/info/docs/tools/vep/index.html)<br/>
- MITOMAP (https://mitomap.org/MITOMAP)<br/>
- MitoTIP (https://www.mitomap.org/foswiki/bin/view/MITOMAP/MitoTipInfo)<br/>
- ClinVar (https://www.ncbi.nlm.nih.gov/clinvar)<br/>
- HMTDB (https://www.hmtdb.uniba.it/)<br/>
- MitImpact (https://mitimpact.css-mendel.it/)<br/>
- Haplogrep (https://github.com/seppinho/haplogrep-cmd)<br/>
- MitoKons (tool to perform deep phylum mtDNA conservation analysis)<br/>
- Homopolymers regions
```diff
- Fast VCF annotation is perform using vcfanno (https://github.com/brentp/vcfanno/).
```

## Requirements
- python3<br/>
- vep (http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)<br/>
- cyvcf2 (https://pypi.org/project/cyvcf2)<br/>
- xlsxwriter (https://pypi.org/project/xlsxwriter)<br/><br/>


## Usage
```markdown
The typical command for running KerMit is as follows:

 python KerMit.py [options] --in path_input_vcf --out path_output_vcf --vep path_vep

 Mandatory arguments:
   --in                        Path to input VCF file
   --out                       Path to output VCF file
   --vep                       Path to VEP executable

 Options arguments:
   --xlsx                      Path to output XLSX file
   --update                    Update source files [false]
   --bgzip                     bgzip and tabix output vcf [true]
   --tmp                       Path to temporary folder [/tmp]
   --thread                    Number of threads to use [1]
   --help                      Print help
   --color                     Enable color output [true]
```


