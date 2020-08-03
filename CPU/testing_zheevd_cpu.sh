#!/bin/sh                                                                                                      
#SBATCH --get-user-env                                                                                         
#SBATCH -p prawnew                                                                                             
#SBATCH --mem=2000M                                                                                            
#SBATCH --time=0-01:00:00                                                                                      

mkdir -p $HOME/data/Diagonalization
cd $HOME/data/Diagonalization

OUTPUT=zheevd_cpu_$(hostname).txt
echo \#$(hostname) > $OUTPUT
echo 'for n in $(seq 1 100);do /home/sugimoto/bin/Diagonalization/CPU/testing_zheevd_cpu.out $((500*$n)) 1>/dev/null 2>$OUTPUT; done'
for n in $(seq 1 100);do /home/sugimoto/bin/Diagonalization/CPU/testing_zheevd_cpu.out $((500*$n)) 1>/dev/null;done 2>$OUTPUT
