import ROOT
import os
import argparse

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile')
p.add_argument('meson_option', help='Type <<rho>> for rho, <<phi>> for phi') #flag for type of meson
args = p.parse_args()


#INPUT FILES
if args.meson_option == "rho" :
    fIn_bkg  = ROOT.TFile("../histos/latest_productions/CR_Rho_preselection_Sidebands.root")
    #tree_bkg = fIn_bkg.Get("tree_output_forMVA")
    tree_bkg = fIn_bkg.Get("tree_output")
    fIn_sig  = ROOT.TFile("../histos/latest_productions/SR_Rho_preselection_Signal.root")
    #tree_sig = fIn_sig.Get("tree_output_forMVA")
    tree_sig = fIn_sig.Get("tree_output")

elif args.meson_option == "phi" :
    fIn_bkg  = ROOT.TFile("../histos/latest_productions/CR_Phi_preselection_Sidebands.root")
    #tree_bkg = fIn_bkg.Get("tree_output_forMVA")
    tree_bkg = fIn_bkg.Get("tree_output")
    fIn_sig  = ROOT.TFile("../histos/latest_productions/SR_Phi_preselection_Signal.root")
    #tree_sig = fIn_sig.Get("tree_output_forMVA")
    tree_sig = fIn_sig.Get("tree_output")

#OUTPUT FILE
fOut = ROOT.TFile("outputs/Nominal_training.root","RECREATE")

#START MVA
ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","Transformations=I;D;P;G,D","AnalysisType=Classification"]))

dataloader = ROOT.TMVA.DataLoader()

#VARIABLES FROM THE TREE -------------------------------------------------------------------------------

#first track:  _firstCandPt    _firstCandIso    _firstCandEta
#second track: _secondCandPt   _secondCandIso   _secondCandEta
#best pair:   _bestCouplePt   _coupleIso   _bestCoupleEta   _bestCoupleDeltaR  mass_KK (using m_kk variable is not correct since you use it to define the bkg estimation, use it only to compute scatter plots to see correlation with other variables)
#photon:   _photonEt  _photonEta
#jet: _bestJetPt
#invariant mass: mass_KKg

#Pay attention to the order, it must be the same in the function_smuggler.py
dataloader.AddVariable("firstTrkIsoCh","F")
dataloader.AddVariable("pairIso0","F")
#dataloader.AddVariable("_bestCouplePt/_HpT","F")
#dataloader.AddVariable("_photonEt/_HpT","F")
#dataloader.AddVariable("bestPairPt/ZMass","F")
#dataloader.AddVariable("photonEt/ZMass","F")################remove
#dataloader.AddVariable("photonEta","F")
dataloader.AddVariable("bestPairEta","F")##################add
#dataloader.AddVariable("_JetNeutralEmEnergy","F")
#dataloader.AddVariable("_JetChargedHadEnergy","F")
#dataloader.AddVariable("_JetNeutralHadEnergy","F")
#dataloader.AddVariable("mesonGammaMass","F") #used just to check the correlation of Hmass and other input variables
#dataloader.AddVariable("_HpT","F") #used just to check the correlation of Hmass and other input variables
#dataloader.AddVariable("_bestCouplePt","F")
#dataloader.AddVariable("_photonEt","F")
#dataloader.AddVariable("_metPt","F")
#dataloader.AddVariable("_nJets","F")
#dataloader.AddVariable("_bestJetPt/mesonGammaMass","F")
#dataloader.AddVariable("_bestCoupleDeltaR","F")
dataloader.AddVariable("MesonGammaDeltaPhi","F")
#dataloader.AddVariable("_secondTrkIso","F")


##-------------------------------------------------------------------------------------------------------

sig_weight = 1.
bkg_weight = 1.

print "before AddSignalTree"
dataloader.AddSignalTree(tree_sig, sig_weight)
print "before AddBackgroundTree"
dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

###############dataloader.SetWeightExpression("_BDTweight") #_BDTweight is the weight variable of the tree

mycutSig = ROOT.TCut("")
mycutBkg = ROOT.TCut("")


dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"]))

#dataloader.PrepareTrainingAndTestTree(mycutBkg, ":".join(["!V","nTrain_Background=0:nTest_Background=0"]))

#method_btd1   = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT1", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))
########method_btd    = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining","SeparationType=GiniIndex","UseBaggedBoost:BaggedSampleFraction=0.5","H=True","VerbosityLevel=Verbose"]))
#method_bdt2    = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT2", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))
#method_Fisher  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kFisher, "BoostedFisher","H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.25:!Boost_DetailedMonitoring" )
#method_MLP     = factory.BookMethod(dataloader, ROOT.TMVA.Types.kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
#method_MLPBFGS = factory.BookMethod(dataloader, ROOT.TMVA.Types.kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
#method_btd3    = factory.BookMethod( dataloader, ROOT.TMVA.Types.kBDT, "BDT3","!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
method_btd    = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=600:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:UseRandomisedTrees:UseNvars=4:UseBaggedBoost:BaggedSampleFraction=0.5"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

#weightfile_dir = "default/weights/TMVAClassification_BDT.weights.xml"

#if evaluate_BDT_systematic:
    # weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_shifted.xml"
    # weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_shifted.xml"
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_secondPart.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_secondPart.xml"
#elif test_on_signal_shifted_up:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_up.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_up.xml"
#elif test_on_signal_shifted_down:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_down.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_down.xml"
#elif test_on_signal_sin2:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_sin2.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_sin2.xml"
#elif test_on_signal_cos:
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_cos_with_flat_theta.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_cos_with_flat_theta.xml"
#else:
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele.xml"
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Wmass.xml"
 #   weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_firstPart.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Wmass.xml"
  #  weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_firstPart.xml"

#if isMuon:
 #   rename_weightfile = "mv " + weightfile_dir + " " + weightfile_mu
  #  os.system(rename_weightfile)
#else:
 #   rename_weightfile = "mv " + weightfile_dir + " " + weightfile_ele
  #  os.system(rename_weightfile)
