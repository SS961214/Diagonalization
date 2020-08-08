#!/bin/sh
#SBATCH --get-user-env
#SBATCH -p prawnew
#SBATCH --mem=2000M
#SBATCH --time=0-01:00:00

programDIR=$(cd $(dirname $0); pwd)
rootDIR=${programDIR%%/bin/*}

mkdir -p $rootDIR/data/Diagonalization
cd $rootDIR/data/Diagonalization

OUTPUT=zheev_cpu_$(hostname).txt
echo \#$(hostname)>$OUTPUT
set +x
for n in $(seq 1 100);do $rootDIR/bin/Diagonalization/CPU/testing_zheev_cpu.out $((500*$n)) 1>/dev/null; done 2>>$OUTPUT
set -x
