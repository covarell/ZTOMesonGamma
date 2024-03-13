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
mass = ROOT.RooRealVar("ZMass","ZMass",75.,200.,"GeV")###########same range of the workspace (75 or 84)

#Signal input -------------------------------------------------------------
if CATEGORY == "BDT" and isPhiGammaAnalysis: fileInputSignal = ROOT.TFile("histos/latest_productions/SR_Phi_BDT_Signal.root")
if CATEGORY == "BDT" and isRhoGammaAnalysis: fileInputSignal = ROOT.TFile("histos/latest_productions/SR_Rho_BDT_Signal.root")
fileInputSignal.cd()
mytreeSignal = fileInputSignal.Get("tree_output")

#Input for unblinded fit
if CATEGORY == "BDT" and isPhiGammaAnalysis: fileInputFit = ROOT.TFile("histos/latest_productions/SR_Phi_BDT_Data.root")
if CATEGORY == "BDT" and isRhoGammaAnalysis: fileInputFit = ROOT.TFile("histos/latest_productions/SR_Rho_BDT_Data.root")
fileInputFit.cd()
tree_fit = fileInputFit.Get("tree_output")


inputWS = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root")  
inputWS.cd()
workspace = inputWS.Get("workspace_"+CHANNEL+"")

signalPDF = workspace.pdf("voigtian_signal") #SIGNAL PDF from the workspace
#workspace.var("mean").setConstant(1) 
#workspace.var("width").setConstant(1)
#workspace.var("sigma").setConstant(1)


bkgPDF = workspace.pdf("landau_bkg") #BKG PDF from the workspace
#workspace.var("land_mean").setConstant(1)
#workspace.var("land_sigma").setConstant(1)

'''
#Initialize a Landau pdf
land_mean = ROOT.RooRealVar('land_mean', 'land_mean', 103., 102.,103.) 
land_sigma = ROOT.RooRealVar('land_sigma', 'land_sigma', 14., 13., 15.) 
bkgPDF = ROOT.RooLandau("landau_bkg", "bkgPDF", mass, land_mean, land_sigma)
'''

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
#Ndata = h_ZMass.Integral()
if isPhiGammaAnalysis:
    Ndata = 3580.
else:
    Ndata = 52000.
#N bkg estimation
Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",Ndata, 0.5*Ndata, 2*Ndata)

#Define PDFs for the fit --------------------------------------------------------------------------------------------------
cross_sig = ROOT.RooRealVar("cross_sig","The signal cross section",1928./0.036)#1928pb 
lumi = ROOT.RooRealVar("lumi","The luminosity",39540.)#pb^-1
eff  = ROOT.RooRealVar("eff","efficiency ",sigEfficiency) 
if CHANNEL == "Phi": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",0.489) # = BR(phi -> K+K-)
if CHANNEL == "Rho": B_R_meson = ROOT.RooRealVar("B_R_meson","b_r_meson",1.) # = BR(rho -> PiPi)
#Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3*@4",ROOT.RooArgList(lumi,cross_sig,eff,B_R_meson,B_R_))
Nsig = ROOT.RooRealVar("Nsig", "Nsig", 50., 0., 1500.)

#ADD PDFs --------------------------------------------------------------------------------------------------
argListN = ROOT.RooArgList(Nsig,Nbkg)
print "arg list N size = ",argListN.getSize()


totPDF_landau = ROOT.RooAddPdf("totPDF_landau","Land",ROOT.RooArgList(signalPDF,bkgPDF), argListN)
#totPDF_landau=bkgPDF


print "PDFs added!"


#Final useful information for debugging
print ""
print "######################################################"
print "Ndata          = ",Ndata
print "totalSigEvents = ",totalSigEvents
print "NsigPassed     = ",NsigPassed
print "######################################################"
print ""

#prepare a rebinned TH1 for Z mass ---------------------------------------------------------------
xLowRange  = 75.#84
xHighRange = 200.
if isPhiGammaAnalysis:
    h_mZ = ROOT.TH1F("h_mZ","h_mZ", int((xHighRange - xLowRange)/2), xLowRange, xHighRange)
else:
    h_mZ = ROOT.TH1F("h_mZ","h_mZ", int((xHighRange - xLowRange)/2), xLowRange, xHighRange)

nentries = tree_fit.GetEntriesFast()
#print "nEntries_sig = ",nentries_sig


for jentry in xrange(nentries):
    ientry = tree_fit.LoadTree( jentry )
    if ientry < 0:
        print "break"
        break
    nb = tree_fit.GetEntry(jentry)
    if nb <= 0:
        print "nb < 0"
        continue
    
    h_mZ.Fill(tree_fit.ZMass)


observed_data = ROOT.RooDataHist("observed_data", "observed_data", ROOT.RooArgList(mass), h_mZ)
#nEntries = observed_data.numEntries() 
#print "nEntries = ",nEntries


#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_voig = totPDF_landau.fitTo(observed_data,ROOT.RooFit.Save())

#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_voig = ROOT.TCanvas()
canvas_voig.cd()

#Landau frame
if isPhiGammaAnalysis:
    xframe_voig = mass.frame(60)
else:
    xframe_voig = mass.frame(120)

observed_data.plotOn(xframe_voig)
totPDF_landau.plotOn(xframe_voig,ROOT.RooFit.NormRange("full"),ROOT.RooFit.Range("full"),ROOT.RooFit.Name("signal_pdf"), ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_voig.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_voig.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_voig.SetMaximum(1.3*xframe_voig.GetMaximum())
Set_par = ROOT.RooArgSet(Nsig,Nbkg)
totPDF_landau.paramOn(xframe_voig,ROOT.RooFit.Layout(0.65,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1)))
xframe_voig.getAttText().SetTextSize(0.02)
xframe_voig.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly


#Calculate Chi square and parameters 
nParam_voig = fitResult_voig.floatParsFinal().getSize()
chi2_voig = xframe_voig.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_voig = "{:.2f}".format(chi2_voig) #Crop the chi2 to 2 decimal digits
print "Chi square voightian = ",chi2_voig
print "n param voightian = ",nParam_voig
print ""
#print "prob = ",xframe_voig.Prob()


leg1 = ROOT.TLegend(0.6,0.42,0.82,0.80) #right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry(cut_chi2_voig,"#chi^{2}/ndof = " + cut_chi2_voig + " / " + str(nParam_voig),"brNDC")

leg1.Draw()

CMS_lumi.CMS_lumi(canvas_voig, iPeriod, iPos) #Print integrated lumi and energy information

if isPhiGammaAnalysis:
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_unblind.pdf")
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_unblind.png")
else:
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_unblind.pdf")
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_unblind.png")

'''
#Now save the data and the PDF into a Workspace, for later use for statistical analysis
#and save the workspace into a ROOT file
workspace = ROOT.RooWorkspace("w")
getattr(workspace,'import')(totPDF_landau)
getattr(workspace,'import')(observed_data)

fOut = ROOT.TFile("unblinded_fit","RECREATE")
fOut.cd()
workspace.Write()
fOut.Close()


#Open the rootfile and get the workspace from the exercise_6
fInput = ROOT.TFile("unblinded_fit")
fInput.cd()
w = fInput.Get("w")

w.var("land_mean").setConstant(1)
w.var("land_sigma").setConstant(1)
w.var("mean").setConstant(1)
w.var("width").setConstant(1)
w.var("sigma").setConstant(1)
w.var("Nbkg").setConstant(1)


#Set the RooModelConfig and let it know what the content of the workspace is about
model = ROOT.RooStats.ModelConfig()
model.SetWorkspace(w)
model.SetPdf("totPDF_landau")

#Here we explicitly set the value of the parameters for the null.  
#We want no signal contribution, so Nsig = 0
Nsig = w.var("Nsig")
poi = ROOT.RooArgSet(Nsig)
nullParams = poi.snapshot()
nullParams.setRealValue("Nsig",0.)

#Build the profile likelihood calculator
plc = ROOT.RooStats.ProfileLikelihoodCalculator()

plc.SetData(w.data("observed_data"))
plc.SetModel(model)
plc.SetParameters(poi)
plc.SetNullParameters(nullParams)

#We get a HypoTestResult out of the calculator, and we can query it.
hypo_test_result = plc.GetHypoTest()
print "-------------------------------------------------"
print "The p-value for the null is ", hypo_test_result.NullPValue()
print "Corresponding to a signifcance of ", hypo_test_result.Significance()
print "-------------------------------------------------"
'''