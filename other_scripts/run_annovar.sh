#!/bin/bash
#

cd /pub6/Temp/Liaojl/Code/softwares/annovar
annotate_variation.pl was in annovar directory
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/
perl annotate_variation.pl -out /pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/output/SP116478_MGAT4C_anno -build hg19 -hgvs /pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/run_ANNOVAR/SP116478_MGAT4C_mut_dat.tsv humandb/

