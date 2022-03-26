## tb xde
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep lamian.pm.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.pm_tb_xde.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep condiments.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> condiments_tb_xde.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep lamian.chisq.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.chisq_tb_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep monocle2_trajtest_corr.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> monocle2TrajtestCorr_tb_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep monocle2_trajtest.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> monocle2Trajtest_tb_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/tb/code/sex/run/ | grep tradeSeq.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tradeSeq_tb_xde.txt
done
