#!/bin/bash

help() {
        echo ''
        echo "[USAGE] ./chip_ngs_compare.sh [chip_vcf] [chip_pos] [sample_ID] [out_dir] [gvcf]"
}

if [ $# -lt 5 ]; then
        help
else
	chip_vcf=$1
	chip_pos=$2
	sample_ID=$3
	out_dir=$4
	gvcf=$5

	python make_ref.py $chip_vcf $chip_pos $sample_ID $out_dir

	#python compare_array_chip_tabix.py result.list $gvcf $sample_ID

fi
