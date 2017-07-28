echo "cd $PWD ; module add R ; Rscript precache.R --datadir fs " | qsub -l nodes=1:ppn=1,mem=8gb -N fs.precache
echo "cd $PWD ; module add R ; Rscript precache.R --datadir rd " | qsub -l nodes=1:ppn=1,mem=8gb -N rd.precache
echo "cd $PWD ; module add R ; Rscript precache.R --datadir rl " | qsub -l nodes=1:ppn=1,mem=8gb -N rl.precache
echo "cd $PWD ; module add R ; Rscript precache.R --datadir na12878 " | qsub -l nodes=1:ppn=1,mem=8gb -N na12878.precache
echo "cd $PWD ; module add R ; Rscript precache.R --datadir chm " | qsub -l nodes=1:ppn=1,mem=8gb -N chm.precache
