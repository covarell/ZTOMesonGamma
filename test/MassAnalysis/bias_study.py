import ROOT
import argparse
import tdrstyle, CMS_lumi

#BOOL
doBiasStudy = True
isPhiGammaAnalysis = False
isRhoGammaAnalysis = False


#This script does the fit of a combined PDF (signal pdf + bkg pdf) to a generated dataset

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi, <<Rho>> for Rho') #flag for bkg estimation
p.add_argument('CAT_option', help='Type <<BDT>>') #flag for bkg estimation

args = p.parse_args()

if args.Decay_channel_option == "Phi":
    isPhiGammaAnalysis = True
    CHANNEL = "Phi"
    print "Z -> PhiGamma analysis"
if args.Decay_channel_option == "Rho":
    isRhoGammaAnalysis = True
    CHANNEL = "Rho"
    print "Z -> RhoGamma analysis"

CATEGORY = args.CAT_option


#CMS-style plotting ---------------------------------------------------------------
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}"


#Some useful things before starting --------------------------------------------------------------------
ROOT.gROOT.SetBatch(True) #Supress the opening of many Canvas's


#Define the observable --------------------------
mass = ROOT.RooRealVar("ZMass","ZMass",75.,200.,"GeV")######same range as in the workspace (75 or 84)

#Data input -------------------------------------------------------------
if CATEGORY == "BDT" and isPhiGammaAnalysis: fileInputData = ROOT.TFile("histos/latest_productions/CR_Phi_BDT_SidebandsNorm.root")
if CATEGORY == "BDT" and isRhoGammaAnalysis: fileInputData = ROOT.TFile("histos/latest_productions/CR_Rho_BDT_SidebandsNorm.root")
fileInputData.cd()
#mytreeData = fileInputData.Get("tree_output")
h_ZMass = fileInputData.Get("h_ZMass")

#Signal input -------------------------------------------------------------
if CATEGORY == "BDT" and isPhiGammaAnalysis: fileInputSignal = ROOT.TFile("histos/latest_productions/SR_Phi_BDT_Signal.root")
if CATEGORY == "BDT" and isRhoGammaAnalysis: fileInputSignal = ROOT.TFile("histos/latest_productions/SR_Rho_BDT_Signal.root")
fileInputSignal.cd()
mytreeSignal = fileInputSignal.Get("tree_output")

'''
#Initialize a Voigtian pdf
if isPhiGammaAnalysis:
    mean  = ROOT.RooRealVar('mean', 'mean', 91.381, 80., 110.) 
    width = ROOT.RooRealVar('width', 'width', 0.0001, 0., 10.)
    sigma = ROOT.RooRealVar('sigma', 'sigma', 3.312, 0., 10.) 

else:
    mean  = ROOT.RooRealVar('mean', 'mean', 91.3480, 80., 110.) 
    width = ROOT.RooRealVar('width', 'width', 3.452, 0., 10.)
    sigma = ROOT.RooRealVar('sigma', 'sigma', 0.591, 0., 10.) 

mean.setConstant(1)
width.setConstant(1)
sigma.setConstant(1)

signalPDF = ROOT.RooVoigtian("voigtian_signal", "signalPDF", mass, mean, width, sigma)
'''

inputWS = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root")  
inputWS.cd()
workspace = inputWS.Get("workspace_"+CHANNEL+"")

signalPDF = workspace.pdf("voigtian_signal") #SIGNAL PDF from the workspace


'''
#Initialize a Landau pdf
if isPhiGammaAnalysis:
    land_mean = ROOT.RooRealVar('land_mean', 'land_mean', 102.6, 95.,110.) 
    land_sigma = ROOT.RooRealVar('land_sigma', 'land_sigma', 13.8, 10., 20.) 
else:
    land_mean = ROOT.RooRealVar('land_mean', 'land_mean', 103.8, 95.,110.) 
    land_sigma = ROOT.RooRealVar('land_sigma', 'land_sigma', 14.3, 10., 20.)

land_mean.setConstant(1)
land_sigma.setConstant(1)

bkgPDF_landau = ROOT.RooLandau("landau_bkg", "bkgPDF", mass, land_mean, land_sigma)
'''
bkgPDF_landau = workspace.pdf("landau_bkg") #BKG PDF from the workspace

'''
#Initialize a Chebichev
if isPhiGammaAnalysis:
    a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",-0.622,-5,5.)
    b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",-0.434,-5.,5.)
    c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.49,-5.,1)
    d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",-0.310,-0.5,0.)

else:
    a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",-0.5955,-5,5.)
    b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",-0.438,-5.,5.)
    c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.49,-5.,1)
    d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",-0.2598,-0.5,0.)

a_bkg.setConstant(1)
b_bkg.setConstant(1)
c_bkg.setConstant(1)
d_bkg.setConstant(1)

bkgPDF_chebichev = ROOT.RooChebychev("chebychev","bkgPDF",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg,d_bkg))
'''
bkgPDF_chebichev = workspace.pdf("chebychev") #BKG PDF from the workspace


#Define the POI -------------------------------------------------------------------------------------
B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",0.000001,-0.001,0.001) #parameter of interest

#Number of events in m_KKg
fInput_sig = ROOT.TFile("rootfiles/latest_productions/ZMesonGamma"+CHANNEL+"_Signal.root")
#mytree_sig     = fInput_sig.Get("ZMesonGamma/mytree")
h_Events       = fInput_sig.Get("ZMesonGamma/hEvents")
NsigPassed     = mytreeSignal.GetEntriesFast()
totalSigEvents = h_Events.GetBinContent(1)
#sigEfficiency  = float(NsigPassed/totalSigEvents)
if isPhiGammaAnalysis: sigEfficiency = 0.0272
else : sigEfficiency = 0.0197

#Ndata = mytreeData.GetEntriesFast()
Ndata = h_ZMass.Integral()
#Ndata = 3580.
#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 0.5*Ndata, 2*Ndata)

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig = ROOT.RooRealVar("cross_sig","The signal cross section",1928./0.036)#1928pb 
lumi = ROOT.RooRealVar("lumi","The luminosity",39540.)#pb^-1
eff  = ROOT.RooRealVar("eff","efficiency ",sigEfficiency) 
if CHANNEL == "Phi": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",0.489) # = BR(phi -> K+K-)
if CHANNEL == "Rho": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> PiPi)
Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig,eff,B_R_meson,B_R_))
#Nsig = ROOT.RooRealVar("Nsig", "Nsig", 0., 0., 50.)

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig,Nbkg)
print "arg list N size = ",argListN.getSize()


totPDF_landau = ROOT.RooAddPdf("totPDF_landau","Land",ROOT.RooArgList(signalPDF,bkgPDF_landau), argListN)
totPDF_chebi = ROOT.RooAddPdf("totPDF_chebi","Chebi",ROOT.RooArgList(signalPDF,bkgPDF_chebichev), argListN)

print "PDFs added!"


#Final useful information for debugging
print ""
print "######################################################"
print "Ndata          = ",Ndata
print "totalSigEvents = ",totalSigEvents
print "NsigPassed     = ",NsigPassed
print "######################################################"
print ""


#Construct the Toy-MC machinery
genPDF = totPDF_chebi
fitPDF = totPDF_landau

mcstudy = ROOT.RooMCStudy(genPDF, ROOT.RooArgSet(mass),ROOT.RooFit.FitModel(fitPDF), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(1), ROOT.RooFit.FitOptions(ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))
mcstudy.generateAndFit(500)

#Plot the distributions of the fitted parameter, the error and the pull
POIval_frame = mcstudy.plotParam(B_R_, ROOT.RooFit.Bins(40))
POIerr_frame = mcstudy.plotError(B_R_, ROOT.RooFit.Bins(40))
POIpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Bins(40), ROOT.RooFit.FitGauss(1))

#Plot distribution of minimized likelihood
NLLframe = mcstudy.plotNLL(ROOT.RooFit.Bins(40))

#Actually plot
canvas = ROOT.TCanvas()

canvas.Divide(2,2)
canvas.cd(1)
POIval_frame.Draw()
canvas.cd(2)
POIerr_frame.Draw()
canvas.cd(3)
POIpull_frame.Draw()
canvas.cd(4)
NLLframe.Draw()

'''
canvas.cd()
POIpull_frame.Draw()
'''
canvas.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/"+CHANNEL+"/Bias/BRpull.png")
canvas.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/"+CHANNEL+"/Bias/BRpull.pdf")

