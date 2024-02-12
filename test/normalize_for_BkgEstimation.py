import ROOT
import argparse
import numpy as np


p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Category', help='Type the category') #flag for bkg estimation
p.add_argument('Meson', help='Type phi or rho')
args = p.parse_args()
CAT = args.Category
Meson = args.Meson

#Take data rootfiles for control regions
if Meson == "phi" :
    DataSR = ROOT.TFile("histos/latest_productions/SR_Phi_"+CAT+"_Data.root")
    Sidebands = ROOT.TFile("histos/latest_productions/CR_Phi_"+CAT+"_Sidebands.root")
elif Meson == "rho" :
    DataSR = ROOT.TFile("histos/latest_productions/SR_Rho_"+CAT+"_Data.root")
    Sidebands = ROOT.TFile("histos/latest_productions/CR_Rho_"+CAT+"_Sidebands.root")
    #Sidebands = ROOT.TFile("histos/latest_production/histos_CR_"+CAT+"_RightSideband.root")

#Output file creation
if Meson == "phi" :
    fOut = ROOT.TFile("histos/latest_productions/CR_Phi_"+CAT+"_SidebandsNorm.root","RECREATE")
elif Meson == "rho" :
    fOut = ROOT.TFile("histos/latest_productions/CR_Rho_"+CAT+"_SidebandsNorm.root","RECREATE")

fOut.cd()

#Get the list of histograms
list_histos = []

keylist = DataSR.GetListOfKeys()
#print keylist
key = ROOT.TKey()
#print "prova"
for key in keylist :
    #print "prova"
    obj_class = ROOT.gROOT.GetClass(key.GetClassName())
    if not obj_class.InheritsFrom("TH1") :
        continue
    if not (key.ReadObj().GetName() == "h_efficiency" or key.ReadObj().GetName() == "h_cutOverflow"): #h_efficiency and h_cutOverflow is a plot plotted in other way   
        list_histos.append( key.ReadObj().GetName() )

    #for histo_name in list_histos:    
    #print "prova"

#list_histos = ["h_ZMass","h_mesonMass","h_firstTrk_pT","h_secondTrk_pT","h_firstTrk_Eta","h_secondTrk_Eta","h_firstTrk_Phi","h_secondTrk_Phi","h_bestCouplePt","h_bestCoupleEta","h_bestCoupleDeltaR","h_bestJetPt","h_bestJetEta","h_firstTrk_Iso","h_firstTrk_Iso_ch","h_secondTrk_Iso","h_secondTrk_Iso_ch","h_couple_Iso","h_couple_Iso_ch","h_photon_energy","h_photon_eta","h_nJets_25","h_nMuons","h_nElectrons","h_nPhotons38WP80","h_nPhotons20WP90","h_decayChannel","h_couple_Iso_neutral","h_met_pT","h_dPhiGammaTrk","h_pTOverHmass","h_eTOverHmass","h_JetChargedEmEnergy","h_JetNeutralEmEnergy","h_JetChargedHadEnergy","h_JetNeutralHadEnergy","h_massResolution","h_genPhotonEt","h_genMesonPt"]#,"h_BDT_out"]
for histo_name in list_histos:
    
    print "###############"
    print "histo_name = ",histo_name

    histoSR = DataSR.Get(histo_name)
    #print histoSR
    histoCR = Sidebands.Get(histo_name)
    #print histoCR

    if histo_name == "h_ZMass" or histo_name == "h_InvMass_TwoTrk_Photon_NoPhiMassCut":
        #print "CRIntegral =", histoCR.Integral()
        CRintegral = histoCR.Integral() - histoCR.Integral(histoCR.GetXaxis().FindBin(80.),histoCR.GetXaxis().FindBin(100.5)) #since in this plot there is the blind window for data in SR, this trick is to make the divide properly. Remember to bypass it for the unblinding
        #print "ZMass CR integral = ", CRintegral
    else:
        CRintegral = histoCR.Integral()

    SRintegral = histoSR.Integral()
    
    print "histo SR integral = ", SRintegral
    print "histo CR integral = ", CRintegral

    if not CRintegral == 0:
        histoCR.Scale(SRintegral/CRintegral)
    histoCR.Write()

    if histo_name == "h_ZMass":
        print SRintegral/CRintegral

    print "histo CR integral after the normalization = ", histoCR.Integral()
    print "###############"
'''
#Variables to go in the output tree ################################################################################
_ZMass = np.zeros(1, dtype=float)

tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('ZMass',_ZMass,'_ZMass/D') 

#EVENTS LOOP ########################################################################################################
for jentry in range(1,fOut.Get("h_ZMass").GetNbinsX()+1):
    #Retrieve variables from the histo
    ZMass  = fOut.Get("h_ZMass").GetBinCenter(jentry)
    print ZMass
    #FILL TREE ########################################################################################################
    _ZMass[0] = ZMass
    tree_output.Fill()

#Tree writing ##########################################################################################################
tree_output.Write()
'''
fOut.Close()