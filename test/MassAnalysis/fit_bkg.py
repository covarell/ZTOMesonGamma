import ROOT
import tdrstyle, CMS_lumi
import argparse

#####################################################
# This script takes the fit_Signal.py workspace and
# update it with a RooMultiPDF containing bkg models
# with their normalizations
#####################################################

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

#bool initialization
isRhoGammaAnalysis = False
isPhiGammaAnalysis = False

#INPUT and OUTPUT #############################################################################################
#Input
p = argparse.ArgumentParser(description='Select rootfile to plot')
p.add_argument('Decay_channel_option', help='Type <<Phi>> for Phi, <<Rho>> for Rho') #flag for bkg estimation
args = p.parse_args()

if args.Decay_channel_option == "Phi":
    isPhiGammaAnalysis = True
    CHANNEL = "Phi"
    print "Z -> PhiGamma analysis"
if args.Decay_channel_option == "Rho":
    isRhoGammaAnalysis = True
    CHANNEL = "Rho"
    print "Z -> RhoGamma analysis"


#CMS-style plotting ---------------------------------------------------------------
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.lumiTextOffset = 0.25
CMS_lumi.cmsTextSize = 0.8
#CMS_lumi.cmsTextOffset = 0.4
CMS_lumi.lumi_13TeV = "39.54 fb^{-1}" 

#Parameters of the PDF ---------------------------------------------------------------
mass = ROOT.RooRealVar("ZMass","ZMass",70.,200.,"GeV")
mass.setRange("full",70.,200.)

#Initialize a Landau pdf
land_mean = ROOT.RooRealVar('land_mean', 'land_mean', 103., 95.,110.) 
land_sigma = ROOT.RooRealVar('land_sigma', 'land_sigma', 13., 10., 20.) 
bkgPDF_landau = ROOT.RooLandau("landau_bkg", "bkgPDF", mass, land_mean, land_sigma)

#Input file and tree ---------------------------------------------------------------
if isPhiGammaAnalysis:
    fileInput = ROOT.TFile("histos/latest_productions/CR_Phi_BDT_SidebandsNorm.root")###background normalizzato alla signal region di m_KK blinded (vedi normalize_for_Bkg_Estimation.py)
else :
    fileInput = ROOT.TFile("histos/latest_productions/CR_Rho_BDT_SidebandsNorm.root")
fileInput.cd()
#tree = fileInput.Get("tree_output")

h_mZ = fileInput.Get("h_ZMass")
if isPhiGammaAnalysis:
    h_mZ.Rebin(2)
#else : h_mZ.Rebin(3)

#Retrieve observed_data from the tree, insert the variable also ---------------------------------------------------------------
#observed_data = ROOT.RooDataSet("observed_data","observed_data",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
observed_data = ROOT.RooDataHist("observed_data", "observed_data", ROOT.RooArgList(mass), h_mZ)############
#nEntries = observed_data.numEntries()
#nEntries =h_mZ.Integral(h_mZ.FindBin(70.), h_mZ.FindBin(200.))
nEntries=observed_data.sumEntries() 
print "nEntries in bkg distribution = ",nEntries

#recupero l'istogramma della SR e della CR e il tree della SR
if isPhiGammaAnalysis:
    SR_input = ROOT.TFile("histos/latest_productions/SR_Phi_BDT_Data.root")
else :
    SR_input = ROOT.TFile("histos/latest_productions/SR_Rho_BDT_Data.root")

#h_mZ_SR = SR_input.Get("h_ZMass")# mi serve per calcolare l'integrale sulle sidebands di mkkg
#h_mZ_CR = fileInput.Get("h_ZMass")# mi serve per calcolare l'integrale sulle sidebands di mkkg
mZ_SR_tree = SR_input.Get("tree_output")
True_data = ROOT.RooDataSet("True_data","True_data",ROOT.RooArgSet(mass),ROOT.RooFit.Import(mZ_SR_tree))#Dati su cui fare il fit con Combine
#True_data = ROOT.RooDataHist("True_data", "True_data", ROOT.RooArgList(mass), h_mZ_SR)


#Do the fit ------------------------------------------------------------------------------------------------------------------
fitResult_landau = bkgPDF_landau.fitTo(observed_data,ROOT.RooFit.Save())

#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_landau = ROOT.TCanvas()
canvas_landau.cd()

#Landau frame
if isPhiGammaAnalysis:
    xframe_landau = mass.frame(60)
else:
    xframe_landau = mass.frame(120)

observed_data.plotOn(xframe_landau)
bkgPDF_landau.plotOn(xframe_landau,ROOT.RooFit.NormRange("full"),ROOT.RooFit.Range("full"),ROOT.RooFit.Name("bkgPDF_landau"), ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_landau.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_landau.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_landau.SetMaximum(1.3*xframe_landau.GetMaximum())
bkgPDF_landau.paramOn(xframe_landau, ROOT.RooFit.Layout(0.65,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_landau.getAttText().SetTextSize(0.02)
xframe_landau.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam_landau = fitResult_landau.floatParsFinal().getSize()
chi2_landau = xframe_landau.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_landau = "{:.2f}".format(chi2_landau) #Crop the chi2 to 2 decimal digits
print "Chi square convolution = ",chi2_landau
print "n param convolution = ",nParam_landau
print ""

leg1 = ROOT.TLegend(0.5,0.52,0.72,0.90) #right positioning
leg1.SetHeader(" ")
leg1.SetNColumns(1)
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.AddEntry(cut_chi2_landau,"#chi^{2}/ndof = " + cut_chi2_landau + " / " + str(nParam_landau),"brNDC")

leg1.Draw()

CMS_lumi.CMS_lumi(canvas_landau, iPeriod, iPos) #Print integrated lumi and energy information

if isPhiGammaAnalysis:
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_bkg_BDT.pdf")
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_bkg_BDT.png")
else:
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_BDT.pdf")
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_BDT.png")


#create Workspace ------------------------------------------------------------------------------------------------------------------------------
norm  = nEntries#*(h_mZ_SR.Integral(h_mZ_SR.FindBin(70.),h_mZ_SR.FindBin(200.)))/(h_mZ_CR.Integral(h_mZ_CR.FindBin(70.),h_mZ_CR.FindBin(200.))-h_mZ_CR.Integral(h_mZ_CR.FindBin(80.),h_mZ_CR.FindBin(100.)))
print "bkg normalization:", norm
bkg_norm = ROOT.RooRealVar(bkgPDF_landau.GetName()+ "_norm", bkgPDF_landau.GetName()+ "_norm", norm, 0.5*norm, 2*norm)

inputWS = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root")  
inputWS.cd()
workspace = inputWS.Get("workspace_"+CHANNEL+"")
getattr(workspace,'import')(bkgPDF_landau)
getattr(workspace,'import')(True_data)
getattr(workspace,'import')(bkg_norm)
print "integral BKG :",bkg_norm.Print()

fOut = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root","UPDATE")
fOut.cd()
workspace.Write()
print "-------------------------------------------"
print "Final print to check the workspace update:"
workspace.Print()

fOut.Close()
