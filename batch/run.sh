#!/bin/bash

file_list=$1
output_dir=$2

ownroot=/lustre/nyx/hades/user/mmamaev/install/root-6.20.04-centos7-cxx17/bin/thisroot.sh

current_dir=$(pwd)
partition=main
time=8:00:00
build_dir=/lustre/cbm/users/aandomin/na61/analysis_template/build/src

lists_dir=${output_dir}/lists/
log_dir=${output_dir}/log/

mkdir -p $output_dir
mkdir -p $log_dir
mkdir -p $lists_dir

split -l 30 -d -a 3 --additional-suffix=.list "$file_list" $lists_dir

n_runs=$(ls $lists_dir/*.list | wc -l)

job_range=1-$n_runs

echo file list=$file_list
echo output_dir=$output_dir
echo log_dir=$log_dir
echo lists_dir=$lists_dir
echo n_runs=$n_runs
echo job_range=$job_range

sbatch -J analysis_template -p $partition -t $time -a $job_range -e ${log_dir}/%A_%a.e -o ${log_dir}/%A_%a.o --export=output_dir=$output_dir,file_list=$file_list,ownroot=$ownroot,lists_dir=$lists_dir,build_dir=$build_dir -- /cvmfs/vae.gsi.de/centos7/containers/centos7-user_container_20211027T1055.sif /lustre/cbm/users/aandomin/na61/analysis_template/batch/batch_run.sh