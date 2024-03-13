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
mass = ROOT.RooRealVar("ZMass","ZMass",70.,180.,"GeV")
mass.setRange("full",70.,180.)

#Initialize a Landau pdf
land_mean = ROOT.RooRealVar('land_mean', 'land_mean', 103., 95.,110.) 
land_sigma = ROOT.RooRealVar('land_sigma', 'land_sigma', 13., 10., 20.) 
bkgPDF_landau = ROOT.RooLandau("landau_bkg", "bkgPDF", mass, land_mean, land_sigma)
##################################################################################################
#Initialize a Gaussian
gaus_mean = ROOT.RooRealVar('gaus_mean', 'gaus_mean', 103., 95.,105.)
gaus_sigma = ROOT.RooRealVar('gaus_sigma', 'gaus_sigma', 10., 0., 50.)
#gaus_sigma_r = ROOT.RooRealVar('gaus_sigma', 'gaus_sigma', 10., 0., 50.)
#gaus = ROOT.RooGaussian("gaus", "gaus", mass, gaus_mean, gaus_sigma)

#Initialize an exponential
#e_par = ROOT.RooRealVar("exp_par", "exp_par", 10, -100, 100.)
#exp = ROOT.RooExponential("exp", "exp", mass, e_par)


#Initialize convolution
#if not isPhiGammaAnalysis:
#    n_1 = ROOT.RooRealVar("n1","number1 of events",20000,0,320000) 
#    n_2 = ROOT.RooRealVar("n2","number2 of events",15000,0,320000)
#else: 
#    n_1 = ROOT.RooRealVar("n1","number1 of events",25000,0,370000) 
#    n_2 = ROOT.RooRealVar("n2","number2 of events",50000,0,370000)

#n=ROOT.RooArgList()
#n.add(n_1)
#n.add(n_2) 
#pdfs=ROOT.RooArgList()
#pdfs.add(gaus)
#pdfs.add(exp)
#convolution_pdf = ROOT.RooAddPdf("convolution", "convolution", pdfs, n)
#convolution_pdf = ROOT.RooFFTConvPdf("convolution", "convolution", mass, gaus, exp)
#convolution_pdf = gaus
#convolution_pdf = ROOT.RooGExpModel("conv", "conv", mass, gaus_mean, gaus_sigma, e_par, gaus_mean_SF, gaus_sigma_SF, e_par_SF)
#convolution_pdf = ROOT.RooNovosibirsk("convolution", "convolution", mass, gaus_mean, gaus_sigma, e_par)
###############
#Initialize CB
#alpha_l = ROOT.RooRealVar("alphal", "alphal", -1, -10, 10.)
#n_l = ROOT.RooRealVar("nl", "nl", 1., 0., 10.)
alpha_r = ROOT.RooRealVar("alphar", "alphar", 1., -10., 0.)#1,-10,0
n_r = ROOT.RooRealVar("nr", "nr", 1., 0., 15.)#1,0,15
#convolution_pdf = ROOT.RooCrystalBall("cb", "cb", mass, gaus_mean, gaus_sigma, gaus_sigma_r, alpha_l, n_l, alpha_r, n_r)
convolution_pdf = ROOT.RooCBShape("cb", "cb", mass, gaus_mean, gaus_sigma, alpha_r, n_r)

#Input file and tree ---------------------------------------------------------------
if isPhiGammaAnalysis:
    fileInput = ROOT.TFile("histos/latest_productions/CR_Phi_preselection_SidebandsNorm.root")###background normalizzato alla signal region di mkk blindata
else :
    fileInput = ROOT.TFile("histos/latest_productions/CR_Rho_preselection_SidebandsNorm.root")
fileInput.cd()
#tree = fileInput.Get("tree_output")

h_mZ = fileInput.Get("h_ZMass")
if isPhiGammaAnalysis:
    h_mZ.Rebin(2)

#Retrieve observed_data from the tree, insert the variable also ---------------------------------------------------------------
#observed_data = ROOT.RooDataSet("observed_data","observed_data",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
observed_data = ROOT.RooDataHist("observed_data", "observed_data", ROOT.RooArgList(mass), h_mZ)############
#nEntries = observed_data.numEntries()
#nEntries =h_mZ.Integral(h_mZ.FindBin(70.), h_mZ.FindBin(200.))
nEntries=observed_data.sumEntries() 
print "nEntries in bkg distribution = ",nEntries

#recupero l'istogramma della SR e della CR e il tree della SR
if isPhiGammaAnalysis:
    SR_input = ROOT.TFile("histos/latest_productions/SR_Phi_preselection_Data.root")
else :
    SR_input = ROOT.TFile("histos/latest_productions/SR_Rho_preselection_Data.root")

h_mZ_SR = SR_input.Get("h_ZMass")# mi serve per calcolare l'integrale sulle sidebands di mkkg
h_mZ_CR = fileInput.Get("h_ZMass")# mi serve per calcolare l'integrale sulle sidebands di mkkg
mZ_SR_tree = SR_input.Get("tree_output")
True_data = ROOT.RooDataSet("True_data","True_data",ROOT.RooArgSet(mass),ROOT.RooFit.Import(mZ_SR_tree))#Dati su cui fare il fit con Combine
#True_data = ROOT.RooDataHist("True_data", "True_data", ROOT.RooArgList(mass), h_mZ_SR)


#Do the fit ------------------------------------------------------------------------------------------------------------------
fitResult_landau = bkgPDF_landau.fitTo(observed_data,ROOT.RooFit.Save())
fitResult_conv = convolution_pdf.fitTo(observed_data,ROOT.RooFit.Save())


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
print "Chi square landau = ",chi2_landau
print "n param landa = ",nParam_landau
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
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_bkg_preselection.pdf")
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_bkg_preselection.png")
else:
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_preselection.pdf")
    canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_preselection.png")


canvas_conv = ROOT.TCanvas()
canvas_conv.cd()

#Landau frame
if isPhiGammaAnalysis:
    xframe_conv = mass.frame(60)
else:
    xframe_conv = mass.frame(120)

observed_data.plotOn(xframe_conv)
convolution_pdf.plotOn(xframe_conv,ROOT.RooFit.NormRange("full"),ROOT.RooFit.Range("full"),ROOT.RooFit.Name("conv"), ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_conv.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_conv.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_conv.SetMaximum(1.3*xframe_conv.GetMaximum())
convolution_pdf.paramOn(xframe_conv, ROOT.RooFit.Layout(0.65,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_conv.getAttText().SetTextSize(0.02)
xframe_conv.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam_conv = fitResult_conv.floatParsFinal().getSize()
chi2_conv = xframe_conv.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_conv = "{:.2f}".format(chi2_conv) #Crop the chi2 to 2 decimal digits
print "Chi square convolution = ",chi2_conv
print "n param convolution = ",nParam_conv
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
leg1.AddEntry(cut_chi2_conv,"#chi^{2}/ndof = " + cut_chi2_conv + " / " + str(nParam_conv),"brNDC")

leg1.Draw()

CMS_lumi.CMS_lumi(canvas_conv, iPeriod, iPos) #Print integrated lumi and energy information

if isPhiGammaAnalysis:
    canvas_conv.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_bkg_conv_preselection.pdf")
    canvas_conv.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_bkg_conv_preselection.png")
else:
    canvas_conv.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_conv_preselection.pdf")
    canvas_conv.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_conv_preselection.png")

'''
#create Workspace ------------------------------------------------------------------------------------------------------------------------------
norm  = nEntries#*(h_mZ_SR.Integral(h_mZ_SR.FindBin(70.),h_mZ_SR.FindBin(200.)))/(h_mZ_CR.Integral(h_mZ_CR.FindBin(70.),h_mZ_CR.FindBin(200.))-h_mZ_CR.Integral(h_mZ_CR.FindBin(80.),h_mZ_CR.FindBin(100.)))
print "bkg normalization:", norm
bkg_norm = ROOT.RooRealVar(bkgPDF_landau.GetName()+ "_norm", bkgPDF_landau.GetName()+ "_norm", norm, 0.5*norm, 2*norm)


land_mean.setConstant(1)
land_sigma.setConstant(1)
bkg_norm.setConstant(1)

inputWS = ROOT.TFile("workspaces/workspace_"+CHANNEL+"_preselection.root")  
inputWS.cd()
workspace = inputWS.Get("workspace_"+CHANNEL+"_preselection")
getattr(workspace,'import')(bkgPDF_landau)
getattr(workspace,'import')(True_data)
getattr(workspace,'import')(bkg_norm)
print "integral BKG :",bkg_norm.Print()

fOut = ROOT.TFile("workspaces/workspace_"+CHANNEL+"_preselection.root","UPDATE")
fOut.cd()
workspace.Write()
print "-------------------------------------------"
print "Final print to check the workspace update:"
workspace.Print()

fOut.Close()
'''