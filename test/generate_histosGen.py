import ROOT
import argparse
import math
import numpy as np
import sys
from array import array

#Following bools are given as input
#isDataBlind = False #Bool for blind analysis
isPhiAnalysis  = True # for Z -> Phi Gamma
isRhoAnalysis  = False # for Z -> Rho Gamma

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
#p.add_argument('CR_option', help='Type <<0>> for SR, <<1>> for CR') #flag for bkg estimation
#p.add_argument('isBlindAnalysis', help='Type <<blind>> or <<unblind>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("ZMesonGammaGen/mytree")


#HISTOS ###########################################################################################################
histo_map = dict()
list_histos = ["h_ZMass", "h_MesonMass", "h_firstTrkPt", "h_secondTrkPt", "h_firstTrkEta", "h_secondTrkEta", "h_firstTrkPhi", "h_secondTrkPhi", "h_firstTrkEnergy", "h_secondTrkEnergy","h_mesonPt", "h_mesonEta", "h_mesonPhi", "h_mesonEnergy", "h_photonPt", "h_photonEta", "h_photonPhi","h_photonEnergy", "h_ZpT", "h_ZEta", "h_ZPhi", "h_ZEnergy", "h_theta_pol", "h_bigTrkPt", "h_smallTrkPt"] 

histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{Z}", 300, 50., 200.) 
if   isPhiAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 1., 1.05) 
elif isRhoAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.5, 1.) 
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 100, 0.,70.)
if   isPhiAnalysis: histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0., 55.)
elif isRhoAnalysis: histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0., 50.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"#eta of the 1st track", 100, -2.5, 2.5)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 2nd track", 100, -2.5, 2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi of the 1st track", 100, -math.pi, math.pi)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 2nd track", 100, -math.pi, math.pi)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"energy of the 1st track", 100, 0., 140.)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"energy of the 2nd track", 100, 0., 140.)
histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"p_{T} of the meson", 100, 0., 140.)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"#eta_{meson}", 100, -2.5,2.5)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"#phi_{meson}", 100, -math.pi, math.pi)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"E_{T} of the meson", 100, 0., 250.)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"p_{T} of the #gamma", 100, 0., 140.)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"#phi_{gamma}", 100, -math.pi, math.pi)
histo_map[list_histos[17]] = ROOT.TH1F(list_histos[17],"E_{T} of the #gamma", 100, 0., 250.)
histo_map[list_histos[18]] = ROOT.TH1F(list_histos[18],"p_{T} of Z", 100, 0., 250.)
histo_map[list_histos[19]] = ROOT.TH1F(list_histos[19],"#eta_{Z}", 100, -2.5,2.5)
histo_map[list_histos[20]] = ROOT.TH1F(list_histos[20],"#phi_{Z}", 100, -math.pi, math.pi)
histo_map[list_histos[21]] = ROOT.TH1F(list_histos[21],"E_{T} of Z", 100, 0., 250.)
histo_map[list_histos[22]] = ROOT.TH1F(list_histos[22],"#theta_{pol}", 100, 0., 3.3)
histo_map[list_histos[23]] = ROOT.TH1F(list_histos[23],"p_{T} of the bigger track", 100, 0., 70.)
histo_map[list_histos[24]] = ROOT.TH1F(list_histos[24],"p_{T} of the smaller track", 100, 0., 70.)






#CREATE OUTPUT ROOTFILE ############################################################################################
fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()

#Variables to go in the output tree ################################################################################
_ZMass          = np.zeros(1, dtype=float)  
_mesonMass      = np.zeros(1, dtype=float)
_firstTrkPt     = np.zeros(1, dtype=float)
_secondTrkPt    = np.zeros(1, dtype=float)
_firstTrkEta    = np.zeros(1, dtype=float)  
_secondTrkEta   = np.zeros(1, dtype=float)
_bestPairPt     = np.zeros(1, dtype=float)
_bestPairEta    = np.zeros(1, dtype=float)  
_photonEt       = np.zeros(1, dtype=float)
_photonEta      = np.zeros(1, dtype=float) 
_theta          = np.zeros(1, dtype=float) 




tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('ZMass',_ZMass,'ZMass/D')
tree_output.Branch('mesonMass',_mesonMass,'MesonMass/D')
tree_output.Branch('firstTrkPt',_firstTrkPt,'firstTrkPt/D')
tree_output.Branch('secondTrkPt',_secondTrkPt,'secondTrkPt/D')
tree_output.Branch('firstTrkEta',_firstTrkEta,'firstTrkEta/D')
tree_output.Branch('secondTrkEta',_secondTrkEta,'secondTrkEta/D')
tree_output.Branch('mesonPt',_bestPairPt,'mesonPt/D')
tree_output.Branch('mesonEta',_bestPairEta,'mesonEta/D')
tree_output.Branch('photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('photonEta',_photonEta,'_photonEta/D')
tree_output.Branch('theta',_theta,'_theta/D')



#EVENTS LOOP ########################################################################################################
nentries = mytree.GetEntriesFast()
for jentry in xrange(nentries):
    ientry = mytree.LoadTree( jentry )
    if ientry < 0:
        break
    nb = mytree.GetEntry(jentry)
    if nb <= 0:
        print "nb < 0"
        continue

    #Retrieve variables from the tree 
    ZMass          = mytree.genZ_mass
    mesonMass      = mytree.genMeson_mass
    firstTrkPt     = mytree.genTrackminus_pT      
    secondTrkPt    = mytree.genTrackplus_pT
    firstTrkEta    = mytree.genTrackminus_eta
    secondTrkEta   = mytree.genTrackplus_eta
    firstTrkPhi    = mytree.genTrackminus_phi
    secondTrkPhi   = mytree.genTrackplus_phi
    firstTrkEnergy = mytree.genTrackminus_E
    secondTrkEnergy= mytree.genTrackplus_E
    mesonPt        = mytree.genMeson_pT      
    mesonEta       = mytree.genMeson_eta
    mesonPhi       = mytree.genMeson_phi
    mesonEnergy    = mytree.genMeson_E
    photonPt       = mytree.genGamma_pT
    photonEta      = mytree.genGamma_eta
    photonPhi      = mytree.genGamma_phi
    photonEt       = mytree.genGamma_E
    ZpT            = mytree.genZ_pT
    ZEta           = mytree.genZ_eta
    ZPhi           = mytree.genZ_phi
    ZEnergy        = mytree.genZ_E
    theta          = mytree.theta_pol
    bigTrkPt       = mytree.genTrackBig_pT
    smallTrkPt     = mytree.genTrackSmall_pT



    eventWeight = 1.
    #FILL HISTOS #####################################################################################################
    #if DATA -> Blind Analysis on Z inv mass plot
    
    histo_map["h_ZMass"].Fill(ZMass, eventWeight)          
    histo_map["h_MesonMass"].Fill(mesonMass, eventWeight)
    histo_map["h_firstTrkPt"].Fill(firstTrkPt, eventWeight)
    histo_map["h_secondTrkPt"].Fill(secondTrkPt, eventWeight)
    histo_map["h_firstTrkEta"].Fill(firstTrkEta, eventWeight)    
    histo_map["h_secondTrkEta"].Fill(secondTrkEta, eventWeight)   
    histo_map["h_firstTrkPhi"].Fill(firstTrkPhi, eventWeight)    
    histo_map["h_secondTrkPhi"].Fill(secondTrkPhi, eventWeight)
    histo_map["h_firstTrkEnergy"].Fill(firstTrkEnergy, eventWeight)    
    histo_map["h_secondTrkEnergy"].Fill(secondTrkEnergy, eventWeight)
    histo_map["h_mesonPt"].Fill(mesonPt, eventWeight)
    histo_map["h_mesonEta"].Fill(mesonEta, eventWeight)
    histo_map["h_mesonPhi"].Fill(mesonPhi, eventWeight)
    histo_map["h_mesonEnergy"].Fill(mesonEnergy, eventWeight)
    histo_map["h_photonPt"].Fill(photonPt, eventWeight)
    histo_map["h_photonEta"].Fill(photonEta, eventWeight)
    histo_map["h_photonPhi"].Fill(photonPhi, eventWeight)
    histo_map["h_photonEnergy"].Fill(photonEt, eventWeight)
    histo_map["h_ZpT"].Fill(ZpT, eventWeight)
    histo_map["h_ZEta"].Fill(ZEta, eventWeight)
    histo_map["h_ZPhi"].Fill(ZPhi, eventWeight)
    histo_map["h_ZEnergy"].Fill(ZEnergy, eventWeight)
    histo_map["h_theta_pol"].Fill(theta, eventWeight)
    histo_map["h_bigTrkPt"].Fill(bigTrkPt, eventWeight)
    histo_map["h_smallTrkPt"].Fill(smallTrkPt, eventWeight)




    

    #FILL TREE ########################################################################################################
    _ZMass[0]          = ZMass
    _mesonMass[0]      = mesonMass
    _firstTrkPt[0]     = firstTrkPt
    _secondTrkPt[0]    = secondTrkPt
    _firstTrkEta[0]    = firstTrkEta
    _secondTrkEta[0]   = secondTrkEta
    _bestPairPt[0]     = mesonPt
    _bestPairEta[0]    = mesonEta     
    _photonEt[0]       = photonEt
    _photonEta[0]      = photonEta 
    _theta[0]          = theta
    
    tree_output.Fill()

#HISTO LABELS #########################################################################################################
histo_map["h_ZMass"].GetXaxis().SetTitle("m_{Z} [GeV/c^2]")
histo_map["h_ZMass"].SetTitle("Z invariant mass")
histo_map["h_MesonMass"].GetXaxis().SetTitle("m_{meson } [GeV/c^2]")
histo_map["h_MesonMass"].SetTitle("Meson mass")
histo_map["h_firstTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_firstTrkPt"].SetTitle("Transverse momentum of the first charged particle ")
histo_map["h_secondTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_secondTrkPt"].SetTitle("Transverse momentum of the second charged particle")
histo_map["h_firstTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_firstTrkEta"].SetTitle("Pseudorapidity of the first charged particle")
histo_map["h_secondTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_secondTrkEta"].SetTitle("Pseudorapidity of the second charged particle")
histo_map["h_firstTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_firstTrkPhi"].SetTitle("Azimuthal angle of the first charged particle")
histo_map["h_secondTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_secondTrkPhi"].SetTitle("Azimuthal angle of the second charged particle")
histo_map["h_firstTrkEnergy"].GetXaxis().SetTitle("E_{trk} [GeV]")
histo_map["h_firstTrkEnergy"].SetTitle("Energy of the first track")
histo_map["h_secondTrkEnergy"].GetXaxis().SetTitle("E_{trk} [GeV]")
histo_map["h_secondTrkEnergy"].SetTitle("Energy of the second track")
histo_map["h_mesonPt"].GetXaxis().SetTitle("pT_{meson} [GeV/c]") 
histo_map["h_mesonPt"].SetTitle("Transverse momentum of the meson")
histo_map["h_mesonEta"].GetXaxis().SetTitle("#eta")
histo_map["h_mesonEta"].SetTitle("Pseudorapidity of the meson")
histo_map["h_mesonPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_mesonPhi"].SetTitle("Azimuthal angle of the meson")
histo_map["h_mesonEnergy"].GetXaxis().SetTitle("E_{meson} [GeV]")
histo_map["h_mesonEnergy"].SetTitle("Energy of the meson")
histo_map["h_photonPt"].GetXaxis().SetTitle("p_{T}^{#gamma} [GeV/c]")
histo_map["h_photonPt"].SetTitle("Transverse momentum of the photon")
histo_map["h_photonEta"].GetXaxis().SetTitle("#eta_{#gamma}")
histo_map["h_photonEta"].SetTitle("Eta of the photon")
histo_map["h_photonPhi"].GetXaxis().SetTitle("#phi_{#gamma} [rad]")
histo_map["h_photonPhi"].SetTitle("phi of the photon")
histo_map["h_photonEnergy"].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
histo_map["h_photonEnergy"].SetTitle("Energy of the photon")
histo_map["h_ZpT"].GetXaxis().SetTitle("pT_{Z} [GeV/c]")
histo_map["h_ZpT"].SetTitle("Transverse momentum of the Z")
histo_map["h_ZEta"].GetXaxis().SetTitle("#eta_{Z}")
histo_map["h_ZEta"].SetTitle("Eta of the Z")
histo_map["h_ZPhi"].GetXaxis().SetTitle("#phi_{Z} [rad]")
histo_map["h_ZPhi"].SetTitle("phi of the Z")
histo_map["h_ZEnergy"].GetXaxis().SetTitle("E_{Z} [GeV]")
histo_map["h_ZEnergy"].SetTitle("Energy of the Z")
histo_map["h_theta_pol"].SetTitle("angle of polarization")
histo_map["h_theta_pol"].GetXaxis().SetTitle("#theta [rad]")
histo_map["h_theta_pol"].SetTitle("angle of polarization")
histo_map["h_bigTrkPt"].GetXaxis().SetTitle("p_{T}^{1} [GeV]")
histo_map["h_bigTrkPt"].SetTitle("p_{T} of the leading track")
histo_map["h_smallTrkPt"].GetXaxis().SetTitle("p_{T}^{2} [GeV]")
histo_map["h_smallTrkPt"].SetTitle("p_{T} of the subleading track")





#Tree writing ##########################################################################################################
tree_output.Write()



#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()

fOut.Close()

    