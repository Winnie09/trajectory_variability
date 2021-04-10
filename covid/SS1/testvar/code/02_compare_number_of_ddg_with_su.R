Res.ss <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/SS1/testvar/result/EM_pm/Se_Mi/numeric_res.rds')
Res.su <- readRDS('covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Se_Mi/numeric_res.rds')
Res.su2 <- readRDS('covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/Mod_Mi/numeric_res.rds')
res.ss <- Res.ss$statistics
res.su <- Res.su$statistics
res.su2 <- Res.su2$statistics

dg.ss.se_mi <- rownames(res.ss)[res.ss[,1] < 0.05]
dg.su.se_mi <- rownames(res.su)[res.su[,1] < 0.05]
dg.su.mod_mi <- rownames(res.su2)[res.su1[,1] < 0.05]

venn.diagram(
  x = list(dg.su.se_mi, dg.ss.se_mi, dg.su.mod_mi),
  category.names = c("Su: Se vs. Mi" , "SS: Se vs. Mi " , "Su: Mod vs. Mi"),
  filename = '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/SS1/testvar/plot/venn_diagramm_su_ss2_number_of_ddg.png',
  output=TRUE
)
