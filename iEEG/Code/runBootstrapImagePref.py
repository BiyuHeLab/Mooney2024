from pathlib import Path
from pwmv_ecog.analysis.image_preference import bootstrap_results, ImagePreferenceAnalysis
from pwmv_ecog.manuscript.nature_figures import SUBJECTS

ipa = ImagePreferenceAnalysis(SUBJECTS).load()
saveDir = "../Analysis/ImagePref/data/bootstrap"
roiTests = {"EVC": ["Post>0"], "HLVC": ["Post>0"]}
bootstrap_results(ipa, Path(saveDir), numBootstraps=2000, roiTests=roiTests)
