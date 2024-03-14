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
mass = ROOT.RooRealVar("ZMass","ZMass",75.,130.,"GeV")#######84-200
mass.setRange("full",75.,130.)#####84-200

#Initialize a Landau pdf
land_mean = ROOT.RooRealVar('land_mean', 'land_mean', 103., 95.,110.) 
land_sigma = ROOT.RooRealVar('land_sigma', 'land_sigma', 13., 10., 20.) 
landau_pdf = ROOT.RooLandau("landau_pdf", "landau_pdf", mass, land_mean, land_sigma)
'''
#Initialize a CrystallBall
gaus_mean = ROOT.RooRealVar('gaus_mean', 'gaus_mean', 103., 95.,105.)
gaus_sigma = ROOT.RooRealVar('gaus_sigma', 'gaus_sigma', 10., 0., 50.)
alpha_r = ROOT.RooRealVar("alphar", "alphar", 1., -10., 10.)#1,-10,0
n_r = ROOT.RooRealVar("nr", "nr", 1., 0., 20.)#1,0,15
convolution_pdf = ROOT.RooCBShape("cb", "cb", mass, gaus_mean, gaus_sigma, alpha_r, n_r)
'''
#Initialize a Chebychev pdf
a_bkg = ROOT.RooRealVar("a_bkg","a_bkg",0.01,-5,5.)
b_bkg = ROOT.RooRealVar("b_bkg","b_bkg",0.3,-5.,5.)
c_bkg = ROOT.RooRealVar("c_bkg","c_bkg",0.01,-5.,1)
d_bkg = ROOT.RooRealVar("d_bkg","d_bkg",-0.05,-0.5,0.)
#e_bkg = ROOT.RooRealVar("e_bkg_","e_bkg",-0.05,-0.1,0.)

chebychev_pdf = ROOT.RooChebychev("chebychev_pdf","chebychev_pdf",mass,ROOT.RooArgList(a_bkg,b_bkg,c_bkg,d_bkg))

#landau_pdf=chebychev_pdf################

#Input file and tree ---------------------------------------------------------------
if isPhiGammaAnalysis:
    fileInput = ROOT.TFile("histos/latest_productions/CR_Phi_BDT_Sidebands.root")
else :
    fileInput = ROOT.TFile("histos/latest_productions/CR_Rho_BDT_Sidebands.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

h_mZ = fileInput.Get("h_ZMass")
#if isPhiGammaAnalysis:
#   h_mZ.Rebin(2)
#else : h_mZ.Rebin(2)

#Retrieve observed_data from the tree, insert the variable also ---------------------------------------------------------------
observed_data = ROOT.RooDataSet("observed_data","observed_data",ROOT.RooArgSet(mass),ROOT.RooFit.Import(tree))
#observed_data = ROOT.RooDataHist("observed_data", "observed_data", ROOT.RooArgList(mass), h_mZ)############
nEntries = observed_data.numEntries()
#nEntries =h_mZ.Integral(h_mZ.FindBin(70.), h_mZ.FindBin(130.))
#nEntries=observed_data.sumEntries() 
print "nEntries in bkg distribution = ",nEntries

#recupero l'istogramma della SR e della CR e il tree della SR
if isPhiGammaAnalysis:
    SR_input = ROOT.TFile("histos/latest_productions/SR_Phi_BDT_Data.root")
else :
    SR_input = ROOT.TFile("histos/latest_productions/SR_Rho_BDT_Data.root")

h_mZ_SR = SR_input.Get("h_ZMass")# mi serve per calcolare l'integrale sulle sidebands di mkkg
h_mZ_CR = fileInput.Get("h_ZMass")# mi serve per calcolare l'integrale sulle sidebands di mkkg
mZ_SR_tree = SR_input.Get("tree_output")
True_data = ROOT.RooDataSet("True_data","True_data",ROOT.RooArgSet(mass),ROOT.RooFit.Import(mZ_SR_tree))#Dati su cui fare il fit con Combine
#True_data = ROOT.RooDataHist("True_data", "True_data", ROOT.RooArgList(mass), h_mZ_SR)


#Do the fit ------------------------------------------------------------------------------------------------------------------
fitResult_landau = landau_pdf.fitTo(observed_data,ROOT.RooFit.Save())
fitResult_cheby = chebychev_pdf.fitTo(observed_data,ROOT.RooFit.Save())

#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_landau = ROOT.TCanvas()
canvas_landau.cd()

#Landau frame
if isPhiGammaAnalysis:
    xframe_landau = mass.frame(60)
else:
    xframe_landau = mass.frame(120)

observed_data.plotOn(xframe_landau)
landau_pdf.plotOn(xframe_landau,ROOT.RooFit.NormRange("full"),ROOT.RooFit.Range("full"),ROOT.RooFit.Name("landau_pdf"), ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_landau.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_landau.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_landau.SetMaximum(1.3*xframe_landau.GetMaximum())
landau_pdf.paramOn(xframe_landau, ROOT.RooFit.Layout(0.65,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_landau.getAttText().SetTextSize(0.02)
xframe_landau.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam_landau = fitResult_landau.floatParsFinal().getSize()
chi2_landau = xframe_landau.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_landau = "{:.2f}".format(chi2_landau) #Crop the chi2 to 2 decimal digits
print "Chi square chebychev = ",chi2_landau
print "n param chebychev = ",nParam_landau
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
#    canvas_landau.SaveAs("/afs/cern.ch/user/c/covarell/www/zphigamma/fit_bkg_BDT.pdf")
   canvas_landau.SaveAs("/afs/cern.ch/user/c/covarell/www/zphigamma/fit_bkg_BDT_phi.png")
else:
 #   canvas_landau.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_BDT.pdf")
    canvas_landau.SaveAs("/afs/cern.ch/user/c/covarell/www/zphigamma/fit_bkg_BDT_rho.png")


canvas_cheby = ROOT.TCanvas()
canvas_cheby.cd()

#Landau frame
if isPhiGammaAnalysis:
    xframe_cheby = mass.frame(60)
else:
    xframe_cheby = mass.frame(120)

observed_data.plotOn(xframe_cheby)
chebychev_pdf.plotOn(xframe_cheby,ROOT.RooFit.NormRange("full"),ROOT.RooFit.Range("full"),ROOT.RooFit.Name("chebychev_pdf"), ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_cheby.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_cheby.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_cheby.SetMaximum(1.3*xframe_cheby.GetMaximum())
chebychev_pdf.paramOn(xframe_cheby, ROOT.RooFit.Layout(0.65,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1))) #,ROOT.RooFit.Layout(0.65,0.90,0.90)
xframe_cheby.getAttText().SetTextSize(0.02)
xframe_cheby.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly

#Calculate Chi square and parameters 
nParam_cheby = fitResult_cheby.floatParsFinal().getSize()
chi2_cheby = xframe_cheby.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_cheby = "{:.2f}".format(chi2_cheby) #Crop the chi2 to 2 decimal digits
print "Chi square chebychev = ",chi2_cheby
print "n param chebychev = ",nParam_cheby
print ""

leg2 = ROOT.TLegend(0.5,0.52,0.72,0.90) #right positioning
leg2.SetHeader(" ")
leg2.SetNColumns(1)
leg2.SetFillColorAlpha(0,0.)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)
leg2.AddEntry(cut_chi2_cheby,"#chi^{2}/ndof = " + cut_chi2_cheby + " / " + str(nParam_cheby),"brNDC")

leg2.Draw()

CMS_lumi.CMS_lumi(canvas_landau, iPeriod, iPos) #Print integrated lumi and energy information

if isPhiGammaAnalysis:
    #canvas_cheby.SaveAs("/fit_bkg_BDT_2.pdf")
    canvas_cheby.SaveAs("/afs/cern.ch/user/c/covarell/www/zphigamma/fit_bkg_BDT_2_phi.png")
else:
    #canvas_cheby.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_bkg_BDT_2.pdf")
    canvas_cheby.SaveAs("/afs/cern.ch/user/c/covarell/www/zphigamma/fit_bkg_BDT_2_rho.png")


# Multipdf ------------------------------------------------------------------------------------------------------------------------------
cat = ROOT.RooCategory("pdf_index","Index of Pdf which is active")
mypdfs = ROOT.RooArgList()
mypdfs.add(landau_pdf)
mypdfs.add(chebychev_pdf)


multipdf = ROOT.RooMultiPdf("multipdf","All Pdfs",cat,mypdfs)




#create Workspace ------------------------------------------------------------------------------------------------------------------------------
norm  = nEntries*(h_mZ_SR.Integral(h_mZ_SR.FindBin(75.),h_mZ_SR.FindBin(130.))-h_mZ_CR.Integral(h_mZ_SR.FindBin(80.),h_mZ_SR.FindBin(100.)))/(h_mZ_CR.Integral(h_mZ_CR.FindBin(75.),h_mZ_CR.FindBin(130.))-h_mZ_CR.Integral(h_mZ_CR.FindBin(80.),h_mZ_CR.FindBin(100.)))
print "bkg normalization:", norm
#bkg_norm = ROOT.RooRealVar(landau_pdf.GetName()+ "_norm", landau_pdf.GetName()+ "_norm", norm, 0.5*norm, 2*norm)
bkg_norm = ROOT.RooRealVar(multipdf.GetName()+ "_norm", multipdf.GetName()+ "_norm", norm,0.5*norm, 2*norm)

#land_mean.setConstant(1)
#land_sigma.setConstant(1)
#a_bkg.setConstant(1)
#b_bkg.setConstant(1)
#c_bkg.setConstant(1)
#d_bkg.setConstant(1)
#e_bkg.setConstant(1)

inputWS = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root")  
inputWS.cd()
workspace = inputWS.Get("workspace_"+CHANNEL+"")
getattr(workspace,'import')(landau_pdf)
getattr(workspace,'import')(chebychev_pdf)
getattr(workspace,'import')(True_data)
getattr(workspace,'import')(cat)
#getattr(workspace,'import')(multipdf)
getattr(workspace,'import')(bkg_norm)
print "integral BKG :",bkg_norm.Print()

fOut = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root","UPDATE")
fOut.cd()
workspace.Write()
print "-------------------------------------------"
print "Final print to check the workspace update:"
workspace.Print()

fOut.Close()


