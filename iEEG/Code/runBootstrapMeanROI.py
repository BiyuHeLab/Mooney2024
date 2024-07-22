from pathlib import Path
from pwmv_ecog.analysis.mean_roi import bootstrap_results, MeanROIAnalysis
from pwmv_ecog.manuscript.nature_figures import SUBJECTS

mra = MeanROIAnalysis(SUBJECTS).load()
saveDir = "../Analysis/MeanROI/data/bootstrap"
roiTests = {roi: ["Post>Pre"] for roi in ["EVC", "HLVC", "FPN"]}
bootstrap_results(mra, Path(saveDir), numBootstraps=2000, roiTests=roiTests)
