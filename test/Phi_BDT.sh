#!/bin/bash
python generate_histos.py phi 1 BDT unblind rootfiles/latest_productions/ZMesonGamma_Data.root histos/latest_productions/CR_Phi_BDT_Sidebands.root
python generate_histos.py phi 0 BDT blind rootfiles/latest_productions/ZMesonGamma_Data.root histos/latest_productions/SR_Phi_BDT_Data.root
python generate_histos.py phi 0 BDT unblind rootfiles/latest_productions/ZMesonGammaPhi_Signal.root histos/latest_productions/SR_Phi_BDT_Signal.root
python normalize_for_BkgEstimation.py BDT phi
python plot_histos_bkgEstimation.py 0.0001 1 1 histos/latest_productions/SR_Phi_BDT_Data.root histos/latest_productions/SR_Phi_BDT_Signal.root histos/latest_productions/CR_Phi_BDT_SidebandsNorm.root
