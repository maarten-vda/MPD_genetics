#!/bin/bash

##Assumes bcftools is in $PATH and vcfs/diagnosed_cohort.vcf.gz and vcfs/undiagnosed_cohort.vcf.gz which both need to be tabix indexed

TSV=## Put filepath here
OUTPUT_TSV=## Put filepath here

skip_first_line=true
zcat "$TSV" | while IFS=$'\t' read -r line
do
    if $skip_first_line; then
        skip_first_line=false
        header=$(zgrep -m 1 '' $TSV)
        echo -e "$header\tundiagnosed_ac\tundiagnosed_an\tundiagnosed_homozygous_count\tdiagnosed_ac\tdiagnosed_an\tdiagnosed_homozygous_count" > "$OUTPUT_TSV"
        continue
    fi
    gnomad_id=$(echo $line | cut -d' ' -f1)
    chr="chr"$(echo $gnomad_id | cut -d'-' -f1)
    pos=$(echo $gnomad_id | cut -d'-' -f2)
    ref=$(echo $gnomad_id | cut -d'-' -f3)
    alt=$(echo $gnomad_id | cut -d'-' -f4)
    ### -r restricts the regions (requires .tbi file), -i chooses the right ID, -H suppresses header ###
    undiagnosed_string=$(bcftools query -r "${chr}:${pos}" -i 'REF=="'"${ref}"'" & ALT=="'"${alt}"'" & GT!="./." & GT!="0/0"' -f '%AC,%AN,[%GT,]' vcfs/undiagnosed_cohort.vcf.gz)
    diagnosed_string=$(bcftools query -r "${chr}:${pos}" -i 'REF=="'"${ref}"'" & ALT=="'"${alt}"'" & GT!="./." & GT!="0/0"' -f '%AC,%AN,[%GT,]' vcfs/diagnosed_cohort.vcf.gz)
    undiagnosed_ac=$(echo $undiagnosed_string | cut -d ',' -f1)
    undiagnosed_an=$(echo $undiagnosed_string |	cut -d ',' -f2)
    undiagnosed_samples=$(echo $undiagnosed_string | cut -d ',' -f3-)
    IFS=',' read -r -a undiagnosed_array <<< "$undiagnosed_samples"
    undiagnosed_homozygous_count=0
    for element in "${undiagnosed_array[@]}"; do
        allele1=$(echo $element | cut -d '/' -f1)
        allele2=$(echo $element | cut -d '/' -f2)
        if [[ $allele1 = $allele2 ]]; then
            ((undiagnosed_homozygous_count+=1))
        fi
    done
    diagnosed_ac=$(echo $diagnosed_string | cut -d ',' -f1)
    diagnosed_an=$(echo $diagnosed_string | cut -d ',' -f2)
    diagnosed_samples=$(echo $diagnosed_string | cut -d ',' -f3-)
    IFS=',' read -r -a diagnosed_array <<< "$diagnosed_samples"
    diagnosed_homozygous_count=0
    for element in "${diagnosed_array[@]}"; do
        allele1=$(echo $element | cut -d '/' -f1)
        allele2=$(echo $element | cut -d '/' -f2)
        if [[ $allele1 = $allele2 ]]; then
            ((diagnosed_homozygous_count+=1))
        fi
    done
## Set the variables to 0 if there is no variable, which happens when the record could not be found in a cohort ##
    if [[ -z "$undiagnosed_ac" ]]; then
        undiagnosed_ac=0
    fi
    if [[ -z "$undiagnosed_an" ]]; then
        undiagnosed_an=0
    fi
    if [[ -z "$diagnosed_ac" ]]; then
       	diagnosed_ac=0
    fi
    if [[ -z "$diagnosed_an" ]]; then
       	diagnosed_an=0
    fi
    echo -e "$line\t$undiagnosed_ac\t$undiagnosed_an\t$undiagnosed_homozygous_count\t$diagnosed_ac\t$diagnosed_an\t$diagnosed_homozygous_count" >> "$OUTPUT_TSV"
done
