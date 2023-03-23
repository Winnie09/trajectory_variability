d <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/tree_variability/res/infer_tree_structure_res.rds')
d <- d$order

ct <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hca/data/proc/ct/sc.rds')
table(ct[d[[1]]])
table(ct[d[[2]]])
table(ct[d[[3]]])

saveRDS(d[[1]],file='/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/traj/result/erythroid/pseudotime.rds')
saveRDS(d[[2]],file='/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/traj/result/monocyte/pseudotime.rds')
saveRDS(d[[3]],file='/scratch/users/whou10@jhu.edu/Wenpin/trajectory_variability/hcahar/real/traj/result/lymph/pseudotime.rds')
