#!/bin/sh
#SBATCH --get-user-env
#SBATCH -p prawnew
#SBATCH --mem=10GB
#SBATCH --time=0-01:00:00

# ---------- Inputs ---------- #
DEVICE="GPU"
ROUTINE="zheevd"
# ---------- (END)Inputs ---------- #

programDIR=$(cd $(dirname $0); pwd)
rootDIR=${programDIR%%/bin/*}
program=$rootDIR/bin/Diagonalization/${DEVICE}/testing_${ROUTINE}_${DEVICE}.out
OUTPUT=${ROUTINE}_${DEVICE}_($(hostname)).txt

mkdir -p $rootDIR/data/Diagonalization
cd $rootDIR/data/Diagonalization

echo \#$(hostname)>$OUTPUT
echo \# 1.Dim of matrix 2.T_diagonalization 3.T_matrixProducts>>$OUTPUT
set +x
for n in $(seq 1 100);do $program $((500*$n)) 1>/dev/null; done 2>>$OUTPUT
set -x
echo Program: $program
echo Output file: $(pwd)/$OUTPUT
