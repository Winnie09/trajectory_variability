
## tb xde
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep .sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tb_xde.txt
done

