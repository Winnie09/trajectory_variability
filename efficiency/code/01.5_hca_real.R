### hca real data ===============================================================
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/code/ | grep .EM.pm.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.pm_hca_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/code/ | grep run21.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tradeSeq_hca_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime/code/ | grep pm.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.pm_hca_tde.txt
done
