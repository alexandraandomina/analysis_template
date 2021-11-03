#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$(ls $lists_dir | sed "${job_num}q;d")

cd $output_dir
mkdir -p $job_num
cd $job_num

module load /cvmfs/vae.gsi.de/centos7/modules/linux-centos7-x86_64/gcc-8.1.0-gcc-4.8.5-oyp4lmr

echo "loading " $ownroot
source $ownroot

echo "executing $build_dir/analyse -i $filelist -t aTree -o PbPb13GeV.root -n -1"

$build_dir/analyse -i $filelist -t aTree -o PbPb13GeV.root -n -1

date $format
echo JOB FINISHED!