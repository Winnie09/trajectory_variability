rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/covid/Su_2020_Cell/testvar/useSaverImp/result/EM_pm/'


Res = readRDS(paste0(rdir, 'Mi_MF/numeric_res.rds'))

s = Res$statistics
sum(s[,1] < 0.05)
# 0
head(s)


Res = readRDS(paste0(rdir, 'Mod_MF/numeric_res.rds'))
s = Res$statistics
sum(s[,1] < 0.05)
# 79
#  [1] "AC025164.1" "AC073111.5" "ACADVL"     "ALKBH7"     "APBB1"     
#  [6] "ARID1A"     "ATOX1"      "ATRX"       "CA5B"       "CAST"      
# [11] "CD68"       "CDC42"      "CHD4"       "CLIP1"      "CWC25"     
# [16] "DDX3X"      "DYNLRB1"    "EIF1AY"     "EIF1B"      "EIF2S3"    
# [21] "EIF3H"      "ENOSF1"     "FAM192A"    "FBXO6"      "FNBP4"     
# [26] "GLYR1"      "GRIPAP1"    "GTF2I"      "HADHB"      "HAUS1"     
# [31] "HP1BP3"     "KDM6A"      "KRTCAP2"    "LSM14A"     "MAP4K2"    
# [36] "MCFD2"      "MEF2C"      "MRPL47"     "MT1F"       "NDRG3"     
# [41] "NSA2"       "PAK1"       "PDLIM5"     "PLEKHB1"    "PRKX"      
# [46] "RALBP1"     "RARRES3"    "RB1"        "RNF181"     "RPL23"     
# [51] "RPL26L1"    "RPS4X"      "RRBP1"      "RSBN1L"     "SEC14L1"   
# [56] "SEC61B"     "SEC61G"     "SF3B6"      "SH3BGRL3"   "SH3GLB1"   
# [61] "SMIM15"     "SPAG9"      "STAG2"      "SUDS3"      "TBCD"      
# [66] "TIMM17A"    "TMED5"      "TMEM165"    "TXN"        "UBE2I"     
# [71] "UTP4"       "VILL"       "VIM-AS1"    "WBP4"       "ZBTB20"    
# [76] "ZDHHC20"    "ZFAND1"     "ZFX"        "ZRSR2"  
head(s)



Res = readRDS(paste0(rdir, 'Mod_Mi_F/numeric_res.rds'))
s = Res$statistics
sum(s[,1] < 0.05)
# 550
s <- s[order(s[,1], -s[,3]), ]
rownames(s)[1:100]
 #  [1] "AC245014.3" "FOSB"       "VIM"        "FAM207A"    "SFPQ"      
 #  [6] "RPL21"      "ZC3H15"     "TOP1"       "EEF1A1"     "DDX5"      
 # [11] "MATR3-1"    "MRPL57"     "AL357060.1" "ATRAID"     "TNFAIP3"   
 # [16] "ZNF593"     "Z93241.1"   "EEF2"       "RPL26"      "PSME1"     
 # [21] "NFKBIA"     "NDUFS8"     "RPL23A"     "NDUFB9"     "TRAPPC2L"  
 # [26] "PSME2"      "YME1L1"     "WSB1"       "AMD1"       "TSPYL2"    
 # [31] "RPS6"       "PABPC1"     "NDUFAF3"    "BTG1"       "HEXIM1"    
 # [36] "ABRACL"     "RPL9"       "RPS15A"     "IER2"       "NDUFA3"    
 # [41] "SSSCA1"     "COX5B"      "LGALS1"     "SIVA1"      "NKTR"      
 # [46] "PSMB1"      "NDUFB3"     "RPL4"       "SLC38A2"    "SNRNP35"   
 # [51] "NDUFA4"     "SARAF"      "ZFP36L2"    "SAMD9"      "FAU"       
 # [56] "MZT2B"      "PTMA"       "COX8A"      "EIF3D"      "MAGOH"     
 # [61] "JUNB"       "RAN"        "SNRPC"      "ILF3-DT"    "CMPK1"     
 # [66] "PIGBOS1"    "RPL35"      "WAPL"       "CCNI"       "MESD"      
 # [71] "DUT"        "UQCRFS1"    "RPL7A"      "EIF2S2"     "IFITM1"    
 # [76] "HELQ"       "EMP3"       "S100A6"     "RPA1"       "ZFP36L1"   
 # [81] "PSMB3"      "DNAJB12"    "COQ7"       "GADD45B"    "RHOG"      
 # [86] "NDUFB4"     "RPL14"      "MYL12A"     "CALM1"      "TMSB10"    
 # [91] "PRDX1"      "C12orf57"   "RPS27A"     "SNRPD3"     "H2AFV"     
 # [96] "TUBA1A"     "RPL10"      "ATP5F1D"    "SNRPD1"     "SF3A1"     


Res = readRDS(paste0(rdir, 'Mod_Mi_M/numeric_res.rds'))
s = Res$statistics
sum(s[,1] < 0.05)
# [1] 273
#   [1] "KLHDC7B"    "RAP1A"      "SMIM15"     "NKAP"       "CDK4"      
#   [6] "ECHDC1"     "STT3B"      "RBM17"      "AK2"        "RPL7L1"    
#  [11] "ICAM3"      "C11orf58"   "PPP2R2A"    "BAG4"       "TGOLN2"    
#  [16] "SLC25A46"   "DCAF5"      "SRP14"      "FAM192A"    "ANXA2R"    
#  [21] "SEC61B"     "GADD45B"    "RNF114"     "SSU72"      "SMAP2"     
#  [26] "FKBP11"     "AC245014.3" "EPRS"       "NOL7"       "ATP5MC1"   
#  [31] "CD69"       "S100A10"    "PWP1"       "PPIE"       "NUCKS1"    
#  [36] "S100A6"     "GPR174"     "RRBP1"      "EGR1"       "VIM"       
#  [41] "PTGDR"      "NKIRAS2"    "ACSS1"      "DYM"        "ATP1A1"    
#  [46] "ERICH1"     "GALT"       "ARRDC1"     "CCNK"       "PPIB"      
#  [51] "MAL"        "AC006369.1" "S1PR5"      "MRPL28"     "PLP2"      
#  [56] "YME1L1"     "S100A4"     "FCGRT"      "FAM207A"    "JUN"       
#  [61] "HDDC3"      "SNRPC"      "PA2G4"      "IL32"       "NSRP1"     
#  [66] "SERBP1"     "DLST"       "MPLKIP"     "ERCC5"      "CORO1B"    
#  [71] "FOXP1"      "MALSU1"     "DUSP1"      "RNPS1"      "TIMM17B"   
#  [76] "CTSZ"       "FOSB"       "NR4A1"      "INPP4B"     "CPNE1"     
#  [81] "PTPN6"      "S100A11"    "PDK3"       "RNF181"     "LAPTM4A"   
#  [86] "NUTM2B-AS1" "GADD45GIP1" "LAPTM5"     "C7orf26"    "SCCPDH"    
#  [91] "SNHG8"      "SPATA13"    "ISCA2"      "NFKBIA"     "SNRNP35"   
#  [96] "NOC3L"      "ZNF706"     "PPIF"       "ATOX1"      "ANXA5"    



