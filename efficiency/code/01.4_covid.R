### covid
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run06.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.pm_covid_xde.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/compcode/tradeSeq/ | grep run.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> tradeSeq_covid_xde.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run17.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.pm_covid_cellprop.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run26.5.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> phenopath3_covid_xde.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run27.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> monocle2trajtest_covid_xde.txt
done


## ============== for some unknown reason, the following files do not exist anymore 
## new methods: check
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep "run35.2.sh | run34.2.sh | run36.2.sh" | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> phenopath100_covid_Mod_Mi.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep " run34.2.sh" | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> phenopath100_covid_Mod_Mi.txt
done

for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep "run36.2.sh" | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> phenopath100_covid_Mod_Mi.txt
done


### real data: covid
for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run28.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> condiments_covid_xde.txt ## rerunning this run28.sh 20220326
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep run29.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> monocle2TrajtestCorr_covid_xde.txt
done


for i in `ls /home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/code/ | grep Lamian.chisq.sh.e | sed 's/.*.sh.e//g'`
do
sacct -j $i --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss >> Lamian.chisq_covid_xde.txt
done



