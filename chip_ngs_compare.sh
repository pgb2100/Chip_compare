#!/bin/bash

export PATH=/mnt/lustre2/BI_Tools/tools/anaconda2/bin/:/cm/shared/apps/sge/2011.11p1/bin/linux-x64/:/cm/shared/Tools/jre-1.8.0/bin:/cm/shared/Tools/ldc2-0.17.1-linux-x86_64/bin:/cm/shared/apps/gcc/4.8.2/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64:/cm/local/apps/ipmitool/1.8.13:/cm/shared/apps/slurm/14.03.0/sbin:/cm/shared/apps/slurm/14.03.0/bin:/cm/local/apps/cluster-tools/bin:/cm/local/apps/cmd/sbin:/cm/local/apps/cmd/bin:/cm/shared/apps/cmgui:/usr/lib64/qt-3.3/bin:/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/opt/ibutils/bin:/opt/dell/srvadmin/bin:/opt/dell/srvadmin/sbin:/root/bin:/usr/local/share/vcftools/vcftools_0.1.12b/bin:/cm/shared/Tools/bin:/cm/shared/CJP/bin:/cm/shared/apps/sge/2011.11p1/bin/linux-x64

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
