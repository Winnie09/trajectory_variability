for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/cellprop/code/ | grep .sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian_hca_simu_cellprop.txt
done

###
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/code/ | grep run01.EM_pm| grep .sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian_hca_simu_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/simu/testvar/addMultiSignalUsingExpr/code/ | grep run10.sh| grep .sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tradeSeq_hca_simu_xde.txt
done

###
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/code/ | grep .EM.pm.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian_hca_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testvar/code/ | grep run21.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tradeSeq_hca_xde.txt
done

###
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run06.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian_covid_Mod_Mi.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/compcode/tradeSeq/ | grep run.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tradeSeq_covid_Mod_Mi.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/real/testtime/code/ | grep pm.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian_hca_tde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run17.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian_covid_cellprop.txt
done



