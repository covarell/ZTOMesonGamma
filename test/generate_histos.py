import ROOT
import argparse
import math
import numpy as np
import sys
from array import array

#Following bools are given as input
isDataBlind   = False #Bool for blind analysis
isPhiAnalysis = False # for Z -> Phi Gamma
isRhoAnalysis = True # for Z -> Rho Gamma

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

# PARSER and INPUT #############################################################################################
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('CR_option', help='Type <<0>> for SR, <<1>> for CR') #flag for bkg estimation
p.add_argument('isBlindAnalysis', help='Type <<blind>> or <<unblind>>') #flag for loose selection or tight selection (from BDT output)
p.add_argument('rootfile_name', help='Type rootfile name')
p.add_argument('outputfile_option', help='Provide output file name')

args = p.parse_args()
fInput = ROOT.TFile(args.rootfile_name)
output_filename = args.outputfile_option
mytree = fInput.Get("ZMesonGamma/mytree")
h_Events = fInput.Get("ZMesonGamma/hEvents")

#SPLIT: I can split a string of chars before and after a split code (that could be a string or a symbol)
#then I can take the string standing before or after with [0] or [1], respectively. 
#[:-n] is not an emoji, it deletes last n symbols from the string
samplename =(args.rootfile_name.split("ZMesonGamma_")[1])[:-5] 
print "samplename = ", samplename

if args.isBlindAnalysis == "blind":
    isDataBlind = True
    print "isDataBlind = ",isDataBlind

CRflag = int(args.CR_option)
if CRflag > 0 :
    print "Processing the control region ", CRflag
else :
    print "Processing the signal region"


#Normalization for MC dataset ################################################################################

#normalization_weight = 

#Combine luminosity
#luminosity2018A = 14.00 #fb^-1
#luminosity2018B = 3.41 #fb^-1   
#luminosity2018C = 6.94 #fb^-1
#luminosity2018D = 31.93 #fb^-1
luminosity = 39.54 #total lumi delivered during the trigger activity: 39.54 #fb^-1

'''if not samplename == "Data":
    weightSum = 0.
else: weightSum = 1.'''

#HISTOS #########################################################################################################
histo_map = dict()
list_histos = ["h_ZMass", "h_MesonMass", "h_firstTrkPt", "h_secondTrkPt", "h_firstTrkEta", "h_secondTrkEta", "h_firstTrkPhi", "h_secondTrkPhi", "h_bestPairPt", "h_bestPairEta", "h_bestPairDeltaR", "h_bestJetPt", "h_bestJetEta", "h_firstTrkIso","h_firstTrkIsoCh","h_secondTrkIso","h_secondTrkIsoCh","h_pairIso","h_pairIsoCh", "h_photonEnergy", "h_photonEta", "h_nJets25","h_nMuons","h_nElectrons", "h_pairIsoNeutral", "h_efficiency"] 

if samplename == "Data":  histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{Z}", 300, 0., 700.)    
else: histo_map[list_histos[0]]  = ROOT.TH1F(list_histos[0],"M_{Z}", 300, 40., 140.) 
if isPhiAnalysis : histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 1., 1.05) 
elif isRhoAnalysis: histo_map[list_histos[1]]  = ROOT.TH1F(list_histos[1],"M_{meson}", 100, 0.5, 1.) 
histo_map[list_histos[2]]  = ROOT.TH1F(list_histos[2],"p_{T} of the 1st track", 100, 0.,70.)
if   isPhiAnalysis: histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0.,55.)
elif isRhoAnalysis: histo_map[list_histos[3]]  = ROOT.TH1F(list_histos[3],"p_{T} of the 2nd track", 100, 0.,50.)
histo_map[list_histos[4]]  = ROOT.TH1F(list_histos[4],"#eta of the 1st track", 100, -2.5,2.5)
histo_map[list_histos[5]]  = ROOT.TH1F(list_histos[5],"#eta of the 2nd track", 100, -2.5,2.5)
histo_map[list_histos[6]]  = ROOT.TH1F(list_histos[6],"#phi of the 1st track", 100, -math.pi, math.pi)
histo_map[list_histos[7]]  = ROOT.TH1F(list_histos[7],"#phi of the 2nd track", 100, -math.pi, math.pi)
histo_map[list_histos[8]]  = ROOT.TH1F(list_histos[8],"p_{T} of the meson", 100, 0.,140.)
histo_map[list_histos[9]]  = ROOT.TH1F(list_histos[9],"#eta_{meson}", 100, -2.5,2.5)
if   isPhiAnalysis: histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#Delta R_{meson}", 100, 0.,0.026)
elif isRhoAnalysis: histo_map[list_histos[10]] = ROOT.TH1F(list_histos[10],"#Delta R_{meson}", 100, 0.,0.07)
histo_map[list_histos[11]] = ROOT.TH1F(list_histos[11],"p_{T} of the jet", 100, 0.,190.)
histo_map[list_histos[12]] = ROOT.TH1F(list_histos[12],"#eta of the jet", 100, -2.5,2.5)
histo_map[list_histos[13]] = ROOT.TH1F(list_histos[13],"Iso of the 1st track", 100, 0.,1.)
histo_map[list_histos[14]] = ROOT.TH1F(list_histos[14],"Iso_ch of the 1st track", 100, 0.,1.)
histo_map[list_histos[15]] = ROOT.TH1F(list_histos[15],"Iso of the 2nd track", 100, 0.,1.)
histo_map[list_histos[16]] = ROOT.TH1F(list_histos[16],"Iso_ch of the 2nd track", 100, 0.,1.)
histo_map[list_histos[17]] = ROOT.TH1F(list_histos[17],"Iso of the meson", 100, 0.,1.)
histo_map[list_histos[18]] = ROOT.TH1F(list_histos[18],"Iso_ch of the meson", 100, 0., 1.)
histo_map[list_histos[19]] = ROOT.TH1F(list_histos[19],"E_{T} of the #gamma", 100, 0., 250.)
histo_map[list_histos[20]] = ROOT.TH1F(list_histos[20],"#eta_{#gamma}", 100, -2.5,2.5)
histo_map[list_histos[21]] = ROOT.TH1F(list_histos[21],"n. of jets over pre-filters",  7, -0.5,6.5)
histo_map[list_histos[22]] = ROOT.TH1F(list_histos[22],"n. of muons", 6, -0.5, 5.5)
histo_map[list_histos[23]] = ROOT.TH1F(list_histos[23],"n. of electrons", 5, -0.5, 5.5)
histo_map[list_histos[24]] = ROOT.TH1F(list_histos[24],"Iso_neutral of the meson", 100, 0.,1.)
histo_map[list_histos[25]] = ROOT.TH1F(list_histos[25],"Efficiency steps", 5, 0.,5.)


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
_pairIsoCh      = np.zeros(1, dtype=float)
_pairIso0       = np.zeros(1, dtype=float)
_pairIso        = np.zeros(1, dtype=float)
_bestJetPt      = np.zeros(1, dtype=float)
_photonEt       = np.zeros(1, dtype=float)
_photonEta      = np.zeros(1, dtype=float)  
_eventWeight    = np.zeros(1, dtype=float)
_nJets          = np.zeros(1, dtype=float)
_photonEt       = np.zeros(1, dtype=float)
_photonEta      = np.zeros(1, dtype=float) 
_firstTrkIso    = np.zeros(1, dtype=float)  
_secondTrkIso   = np.zeros(1, dtype=float)  
_firstTrkIsoCh  = np.zeros(1, dtype=float)  
_secondTrkIsoCh = np.zeros(1, dtype=float) 
_bestPairDeltaR = np.zeros(1, dtype=float) 


tree_output = ROOT.TTree('tree_output','tree_output')
tree_output.Branch('ZMass',_ZMass,'ZMass/D')
tree_output.Branch('mesonMass',_mesonMass,'MesonMass/D')
tree_output.Branch('firstTrkPt',_firstTrkPt,'firstTrkPt/D')
tree_output.Branch('secondTrkPt',_secondTrkPt,'secondTrkPt/D')
tree_output.Branch('firstTrkEta',_firstTrkEta,'firstTrkEta/D')
tree_output.Branch('secondTrkEta',_secondTrkEta,'secondTrkEta/D')
tree_output.Branch('bestPairPt',_bestPairPt,'bestPairPt/D')
tree_output.Branch('bestPairEta',_bestPairEta,'bestPairEta/D')
tree_output.Branch('pairIsoCh',_pairIsoCh,'_pairIsoCh/D')
tree_output.Branch('pairIso',_pairIso,'_pairIso/D')
tree_output.Branch('pairIso0',_pairIso0,'_pairso0/D')
tree_output.Branch('bestJetPt',_bestJetPt,'_bestJetPt/D')
tree_output.Branch('photonEt',_photonEt,'_photonEt/D')
tree_output.Branch('photonEta',_photonEta,'_photonEta/D')
tree_output.Branch('firstTrkIso',_firstTrkIso,'_firstTrkIso/D')
tree_output.Branch('firstTrkIsoCh',_firstTrkIsoCh,'_firstTrkIsoCh/D')
tree_output.Branch('secondTrkIso',_secondTrkIso,'_secondTrkIso/D')
tree_output.Branch('secondTrkIsoCh',_secondTrkIsoCh,'_secondTrkIsoCh/D')
tree_output.Branch('bestPairDeltaR',_bestPairDeltaR,'_bestPairDeltaR/D')
##tree_output.Branch('eventWeight',_eventWeight,'_eventWeight/D')
tree_output.Branch('nJets',_nJets,'_nJets/D')

#------------- counters -----------------
nEventsMesonAnalysis  = 0
nEventsOverLeptonVeto = 0


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
    ZMass          = mytree.ZMassFrom2KPhoton
    mesonMass      = mytree.MesonMass
    firstTrkPt     = mytree.firstTrkPt         
    secondTrkPt    = mytree.secondTrkPt 
    firstTrkEta    = mytree.firstTrkEta
    secondTrkEta   = mytree.secondTrkEta
    firstTrkPhi    = mytree.firstTrkPhi
    secondTrkPhi   = mytree.secondTrkPhi
    mesonPt        = mytree.bestCouplePt      
    mesonEta       = mytree.bestCoupleEta
    photonEt       = mytree.photon_eT
    photonEta      = mytree.photon_eta
    jetPt          = mytree.bestJet_pT   
    jetEta         = mytree.bestJet_eta
    MesonIso       = mytree.iso_couple
    MesonIsoCh     = mytree.iso_couple_ch #bestIso_couple_ch
    MesonIso0      = mesonPt/(mesonPt + (mytree.pairSumPt05 - mytree.pairSumPt05Ch))
    firstTrkiso    = mytree.iso_K1
    firstTrkisoCh  = mytree.iso_K1_ch
    secondTrkiso   = mytree.iso_K2
    secondTrkisoCh = mytree.iso_K2_ch
    nJets          = mytree.nJets25
    nElectrons     = mytree.nElectrons20
    nMuons         = mytree.nMuons10
    isPhiEvent     = mytree.isPhi
    isRhoEvent     = mytree.isRho
    #phi angle folding
    coupleDeltaPhi = math.fabs(mytree.firstTrkPhi - mytree.secondTrkPhi)
    if coupleDeltaPhi > 3.14:
        coupleDeltaPhi = 6.28 - coupleDeltaPhi
    deltaR = math.sqrt((mytree.firstTrkEta - mytree.secondTrkEta)**2 + (coupleDeltaPhi)**2)    


    #If I'm performing a PhiGamma analysis I don't want to choose those events tagged as a RhoGamma events, and viceversa
    if isPhiAnalysis and not isPhiEvent: 
        continue
    if isRhoAnalysis and not isRhoEvent: 
        continue
    nEventsMesonAnalysis+=1


    #Define Control and Signal regions: ------------------------------------------
    if isPhiAnalysis: #for Phi meson
        if CRflag == 0 and not (mesonMass > 1.008 and mesonMass < 1.032) :
            continue

        if CRflag == 1 and (mesonMass > 1.008 and mesonMass < 1.032) :
            continue

    if isRhoAnalysis: #for Rho meson
        if CRflag == 0 and not (mesonMass > 0.62 and mesonMass < 0.92) :
            continue

        if CRflag == 1 and (mesonMass > 0.62 and mesonMass < 0.92) :
            continue

    #NORMALIZATION -------------------------------------------------------------------
    #normalization for MC
    if not samplename == "Data" :
        PUWeight    = mytree.PU_Weight
        weight_sign = mytree.MC_Weight/abs(mytree.MC_Weight) #just take the sign of the MC gen weight
        #eventWeight =  luminosity * normalization_weight * weight_sign * PUWeight
    else:
        eventWeight = 1.

    #Lepton veto
    if nElectrons > 0: continue
    if nMuons     > 0: continue
    nEventsOverLeptonVeto += 1

    #FILL HISTOS #####################################################################################################
    #if DATA -> Blind Analysis on Z inv mass plot
    if samplename == "Data":
        if isDataBlind:
            if ZMass < 80. or ZMass > 110.:
                histo_map["h_ZMass"].Fill(ZMass)#, eventWeight)
        else:
            histo_map["h_ZMass"].Fill(ZMass)#, eventWeight)
    else:
        histo_map["h_ZMass"].Fill(ZMass)#, eventWeight)
            
    histo_map["h_MesonMass"].Fill(mesonMass)#	, eventWeight)
    histo_map["h_firstTrkPt"].Fill(firstTrkPt)#, eventWeight)
    histo_map["h_secondTrkPt"].Fill(secondTrkPt)#, eventWeight)
    histo_map["h_firstTrkEta"].Fill(firstTrkEta,)# eventWeight)    
    histo_map["h_secondTrkEta"].Fill(secondTrkEta)#, eventWeight)   
    histo_map["h_firstTrkPhi"].Fill(firstTrkPhi)#, eventWeight)    
    histo_map["h_secondTrkPhi"].Fill(secondTrkPhi)#, eventWeight)  
    histo_map["h_bestPairPt"].Fill(mesonPt)#, eventWeight)
    histo_map["h_bestPairEta"].Fill(mesonEta)#, eventWeight)
    histo_map["h_pairIsoCh"].Fill(MesonIsoCh)#, eventWeight)
    histo_map["h_bestJetPt"].Fill(jetPt)#, eventWeight)
    histo_map["h_bestJetEta"].Fill(jetEta)#, eventWeight)
    histo_map["h_photonEnergy"].Fill(photonEt)#, eventWeight)
    histo_map["h_photonEta"].Fill(photonEta)#, eventWeight)
    histo_map["h_bestPairDeltaR"].Fill(deltaR)#/mesonMass)#, eventWeight)
    histo_map["h_pairIso"].Fill(MesonIso)#, eventWeight)
    histo_map["h_firstTrkIso"].Fill(firstTrkiso)#, eventWeight)
    histo_map["h_firstTrkIsoCh"].Fill(firstTrkisoCh)#, eventWeight)
    histo_map["h_secondTrkIso"].Fill(secondTrkiso)#, eventWeight)
    histo_map["h_secondTrkIsoCh"].Fill(secondTrkisoCh)#, eventWeight)
    histo_map["h_pairIsoNeutral"].Fill(MesonIso0)#, eventWeight)
    histo_map["h_nJets25"].Fill(nJets)#, eventWeight)
    histo_map["h_nMuons"].Fill(nMuons)#, eventWeight)
    histo_map["h_nElectrons"].Fill(nElectrons)#, eventWeight)

    #FILL TREE ########################################################################################################
    _ZMass[0]          = ZMass
    _mesonMass[0]      = mesonMass
    _firstTrkPt[0]     = firstTrkPt
    _secondTrkPt[0]    = secondTrkPt
    _firstTrkEta[0]    = firstTrkEta
    _secondTrkEta[0]   = secondTrkEta
    _bestPairPt[0]     = mesonPt
    _bestPairEta[0]    = mesonEta
    _pairIsoCh[0]      = MesonIsoCh
    _pairIso[0]        = MesonIso
    _pairIso0[0]       = MesonIso0
    _bestJetPt[0]      = jetPt     
    _photonEt[0]       = photonEt
    _photonEta[0]      = photonEta 
    _firstTrkIso[0]    = firstTrkiso        
    _secondTrkIso[0]   = secondTrkiso        
    _firstTrkIsoCh[0]  = firstTrkisoCh       
    _secondTrkIsoCh[0] = secondTrkisoCh
    _bestPairDeltaR[0] = deltaR
    _nJets[0]          = nJets

#HISTO LABELS #########################################################################################################
histo_map["h_ZMass"].GetXaxis().SetTitle("m_{trk^{+}trk^{-}#gamma} [GeV/c^2]")
histo_map["h_ZMass"].SetTitle("Tracks+Photon invariant mass (Cut on phi inv. mass)")
histo_map["h_MesonMass"].GetXaxis().SetTitle("m_{trk^{+}trk^{-} } [GeV/c^2]")
histo_map["h_MesonMass"].SetTitle("Tracks invariant mass")
histo_map["h_firstTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_firstTrkPt"].SetTitle("Transverse momentum of the first charged particle (pT_{max} of the couple)")
histo_map["h_secondTrkPt"].GetXaxis().SetTitle("pT_{trk} [GeV/c]")
histo_map["h_secondTrkPt"].SetTitle("Transverse momentum of the second charged particle (pT_{min} of the couple)")
histo_map["h_firstTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_firstTrkEta"].SetTitle("Pseudorapidity of the first charged particle (pT_{max} of the couple)")
histo_map["h_secondTrkEta"].GetXaxis().SetTitle("#eta")
histo_map["h_secondTrkEta"].SetTitle("Pseudorapidity of the second charged particle (pT_{min} of the couple)")
histo_map["h_firstTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_firstTrkPhi"].SetTitle("Azimuthal angle of the first charged particle (pT_{max} of the couple)")
histo_map["h_secondTrkPhi"].GetXaxis().SetTitle("#phi [rad]")
histo_map["h_secondTrkPhi"].SetTitle("Azimuthal angle of the second charged particle (pT_{min} of the couple)")
histo_map["h_bestPairPt"].GetXaxis().SetTitle("pT_{trk^{+}trk^{-}} [GeV]") 
histo_map["h_bestPairPt"].SetTitle("Transverse momentum of the couple")
histo_map["h_bestPairEta"].GetXaxis().SetTitle("#eta")
histo_map["h_bestPairEta"].SetTitle("Pseudorapidity of the couple")
histo_map["h_bestJetPt"].GetXaxis().SetTitle("pT_{jet} [GeV]")
histo_map["h_bestJetPt"].SetTitle("Transverse momentum of the jet")
histo_map["h_bestJetEta"].GetXaxis().SetTitle("#eta_{jet}")
histo_map["h_bestJetEta"].SetTitle("Pseudorapidity of the jet")
histo_map["h_firstTrkIso"].GetXaxis().SetTitle("sum pT_{trk_{1}}/pT_{trk_{1}}")
histo_map["h_firstTrkIso"].SetTitle("Isolation of the trk_{1} candidate")
histo_map["h_firstTrkIsoCh"].GetXaxis().SetTitle("sum pT_{trk_{1}}/pT_{trk_{1}}")
histo_map["h_firstTrkIsoCh"].SetTitle("Isolation of the trk_{1} candidate")
histo_map["h_secondTrkIso"].GetXaxis().SetTitle("sum pT_{trk_{2}}/pT_{trk_{2}}")
histo_map["h_secondTrkIso"].SetTitle("Isolation of the trk_{2} candidate")
histo_map["h_secondTrkIsoCh"].GetXaxis().SetTitle("sum pT_{trk_{2}}/pT_{trk_{2}}")
histo_map["h_secondTrkIsoCh"].SetTitle("Isolation of the trk_{1} candidate")
histo_map["h_pairIso"].GetXaxis().SetTitle("sum pT_{2trk}/pT_{2trk}")
histo_map["h_pairIso"].SetTitle("Isolation of the couple candidate")
histo_map["h_pairIsoCh"].GetXaxis().SetTitle("sum pT_{2trk}/pT_{2trk}")
histo_map["h_pairIsoCh"].SetTitle("Isolation of the couple candidate")
histo_map["h_bestPairDeltaR"].GetXaxis().SetTitle("#DeltaR_{trk^{+}trk^{-} }")
histo_map["h_bestPairDeltaR"].SetTitle("Delta R of the couple")
histo_map["h_photonEnergy"].GetXaxis().SetTitle("E_{T}^{#gamma}[GeV]")
histo_map["h_photonEnergy"].SetTitle("Energy of the photon")
histo_map["h_photonEta"].GetXaxis().SetTitle("#eta_{#gamma}")
histo_map["h_photonEta"].SetTitle("Eta of the photon")
histo_map["h_nJets25"].GetXaxis().SetTitle("nJets over 25 GeV")
histo_map["h_nJets25"].SetTitle("# of jets with pT > 25 GeV")
histo_map["h_nMuons"].GetXaxis().SetTitle("nMuons over selection")
histo_map["h_nMuons"].SetTitle("# of muons")
histo_map["h_nElectrons"].GetXaxis().SetTitle("nElectrons over selection")
histo_map["h_nElectrons"].SetTitle("# of electrons")
histo_map["h_efficiency"].GetXaxis().SetTitle("")
histo_map["h_efficiency"].GetYaxis().SetTitle("#epsilon (%)")


#Tree writing ##########################################################################################################
tree_output.Write()

bin1content  = h_Events.GetBinContent(1)
bin2content  = h_Events.GetBinContent(2)
bin3content  = h_Events.GetBinContent(3)
bin4content  = h_Events.GetBinContent(4)
bin5content  = h_Events.GetBinContent(5)


nSignal      = bin1content
scale_factor = 100/nSignal

histo_map["h_efficiency"].Fill(0.5,bin1content*scale_factor)
histo_map["h_efficiency"].Fill(1.5,bin2content*scale_factor)
histo_map["h_efficiency"].Fill(2.5,bin3content*scale_factor)
histo_map["h_efficiency"].Fill(3.5,bin4content*scale_factor)
histo_map["h_efficiency"].Fill(4.5,bin5content*scale_factor)

histo_map["h_efficiency"].GetXaxis().SetBinLabel(1,"Events processed")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(2,"Events triggered")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(3,"Photon requested")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(4,"Best couple found")
histo_map["h_efficiency"].GetXaxis().SetBinLabel(5,"trk-cand pT selection")


c11 = ROOT.TCanvas()
c11.cd()
histo_map["h_efficiency"].SetFillColor(1) 
histo_map["h_efficiency"].SetFillStyle(3003)
ROOT.gStyle.SetPaintTextFormat("4.2f %")
ROOT.gStyle.SetOptStat(0)
histo_map["h_efficiency"].SetMarkerSize(1.4)
histo_map["h_efficiency"].GetXaxis().SetRangeUser(0.,7.1)
#histo_map["h_efficiency"].GetYaxis().SetRangeUser(0.,30.)
#histo_map["h_efficiency"].SetMaximum(max(histo_map["h_efficiency"].GetHistogram().GetMaximum(),30.))
histo_map["h_efficiency"].Draw("HIST TEXT0")
if not samplename == "Data" :
    c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Signal/h_efficiency.pdf")
    c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Signal/h_efficiency.png")
else:
    if isPhiAnalysis and CRflag == 0:
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Phi/SR/h_efficiency.pdf")
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Phi/SR/h_efficiency.png")
    elif isPhiAnalysis and CRflag == 1:
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Phi/CR/h_efficiency.pdf")
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Phi/CR/h_efficiency.png")
    elif isRhoAnalysis and CRflag == 0:
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Rho/SR/h_efficiency.pdf")
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Rho/SR/h_efficiency.png")
    elif isRhoAnalysis and CRflag == 1: 
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Rho/CR/h_efficiency.pdf")
        c11.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Data/Rho/CR/h_efficiency.png")


#HISTOS WRITING ########################################################################################################
fOut.cd()
for hist_name in list_histos:
    histo_map[hist_name].Write()
fOut.Close()

    