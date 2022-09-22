[Input ML Trees]
#treefile = ../01_iqtree-scaling/iqtree.contree
treefile = ../01_iqtree-scaling/iqtree.ufboot

[General Commands]
numsites = 1768348
nthreads = 16
thorough
log_pen

[Calibrations]
# 1 Haasiophis terrasanctus Tchernov et al., 2000
mrca = Alethinophidia_stem Epicrates_crassus Anilius_scytale
min = Alethinophidia_stem 93.9
max = Alethinophidia_stem 93.9
# 2 Titanoboa cerrejonensis Head et al., 2009
mrca = Boinae_stem Epicrates_crassus Eryx_conicus
min = Boinae_stem 58
# 3 Procerophis sahnii Rage et al., 2008
mrca = Colubroides_stem Acrochordus_javanicus Pareas_nuchalis
min = Colubroides_stem 54
# 4 Elapidae indet McCartney et al., 2014
mrca = Elapidae_stem Atractaspis_engaddensis Calliophis_maculiceps
min = Elapidae_stem 24.9
# 5 Natrix longivertebrata Rage and Szyndlar, 1986
mrca = Natricidae_crow Thamnophis_marcianus Xenochrophis_psicator
min = Natricidae_crow 13.8
# 6 Paleoheterodon tiheni Holman, 1964
mrca = NW_Dipsadidae_stem Thermophis_baileyi Elapomorphus_quinquelineatus
min = NW_Dipsadidae_stem 12.5
# 7 Vipera cf. V. antiqua Szyndlar and BÃ¶hme, 1993
mrca = Viperinae_stem Causus_maculatus Azemiops_feae
min = Viperinae_stem 22.1
# 8 Sistrurus sp. indet. Parmley and Holman, 2007
mrca = Crotalinae_stem Azemiops_feae Agkistrodon_contortrix
min = Crotalinae_stem 10.3

[Priming Command]
#Run only once on best ML tree to determine optimization parameters
#prime

[Best Optimization Parameters]
opt = 3
moredetail
optad = 3
moredetailad
optcvad = 1
moredetailcvad

[Cross Validation Analysis]
#Run only to figure out smoothing parameter, choose lowest chisq or when reaches stationarity
#randomcv
#cviter = 5
#cvsimaniter = 10000000000
#cvstart = 100000
#cvstop = 0.000000000001
#cvmultstep = 0.1
#cvoutfile = 01_treepl_cv.out

[Best Smoothing Parameter]
smooth = 1e-09

[Output File]
outfile = FINAL_treepl.tre
#mapspace
