echo "cd $PWD ; module add R ; Rscript precache.R --datadir fs " | qsub -l nodes=1:ppn=1,mem=8gb -N fs
echo "cd $PWD ; module add R ; Rscript precache.R --datadir rd " | qsub -l nodes=1:ppn=1,mem=8gb -N rd
echo "cd $PWD ; module add R ; Rscript precache.R --datadir rl " | qsub -l nodes=1:ppn=1,mem=8gb -N rl
echo "cd $PWD ; module add R ; Rscript precache.R --datadir na12878 " | qsub -l nodes=1:ppn=1,mem=24gb -N na12878
