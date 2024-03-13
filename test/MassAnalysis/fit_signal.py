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
mass = ROOT.RooRealVar("ZMass","ZMass",91.,84.,100.,"GeV")
mass.setRange("full",84.,100.)


#Initialize a Voigtian pdf
mean  = ROOT.RooRealVar('mean', 'mean', 90., 80., 110.) 
width = ROOT.RooRealVar('width', 'width', 5., 0., 10.)
sigma = ROOT.RooRealVar('sigma', 'sigma', 5., 0., 10.) 

sigPDF_voig = ROOT.RooVoigtian("voigtian_signal", "signalPDF", mass, mean, width, sigma)


#Input file and tree ---------------------------------------------------------------
if isPhiGammaAnalysis:
    fileInput = ROOT.TFile("histos/latest_productions/SR_Phi_BDT_Signal.root")
else :
    fileInput = ROOT.TFile("histos/latest_productions/SR_Rho_BDT_Signal.root")
fileInput.cd()
tree = fileInput.Get("tree_output")

#prepare a rebinned TH1 for Z mass ---------------------------------------------------------------
xLowRange  = 84.
xHighRange = 100.

h_mZ = ROOT.TH1F("h_mZ","h_mZ", int(xHighRange - xLowRange), xLowRange, xHighRange)

nentries_sig = tree.GetEntriesFast()
#print "nEntries_sig = ",nentries_sig

tot=0.
for jentry in xrange(nentries_sig):
    ientry = tree.LoadTree( jentry )
    if ientry < 0:
        print "break"
        break
    nb = tree.GetEntry(jentry)
    if nb <= 0:
        print "nb < 0"
        continue
    #tot+=tree.eventWeight
    #print "eventWeight =", tree.eventWeight
    #print "somma =", tot

    h_mZ.Fill(tree.ZMass, tree.eventWeight)

#Retrieve observed_data from the tree, insert the variable also ---------------------------------------------------------------
observed_data = ROOT.RooDataHist("observed_data", "observed_data", ROOT.RooArgList(mass), h_mZ)
#nEntries = observed_data.numEntries() 
#print "nEntries = ",nEntries


#Do the fit ------------------------------------------------------------------------------------------------------------------------------
fitResult_voig = sigPDF_voig.fitTo(observed_data,ROOT.RooFit.Save())

#Plot ------------------------------------------------------------------------------------------------------------------------
canvas_voig = ROOT.TCanvas()
canvas_voig.cd()

#Landau frame
if isPhiGammaAnalysis:
    xframe_voig = mass.frame(60)
else:
    xframe_voig = mass.frame(120)

observed_data.plotOn(xframe_voig)
sigPDF_voig.plotOn(xframe_voig,ROOT.RooFit.NormRange("full"),ROOT.RooFit.Range("full"),ROOT.RooFit.Name("signal_pdf"), ROOT.RooFit.LineColor(ROOT.kBlue))
xframe_voig.SetTitle("#sqrt{s} = 13 TeV       lumi = 39.54/fb")
xframe_voig.GetXaxis().SetTitle("m_{ditrk,#gamma} [GeV]")
xframe_voig.SetMaximum(1.3*xframe_voig.GetMaximum())
sigPDF_voig.paramOn(xframe_voig,ROOT.RooFit.Layout(0.65,0.94,0.91),ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(1)))
xframe_voig.getAttText().SetTextSize(0.02)
xframe_voig.Draw() #remember to draw the frame before the legend initialization to fill the latter correctly


#Calculate Chi square and parameters 
nParam_voig = fitResult_voig.floatParsFinal().getSize()
chi2_voig = xframe_voig.chiSquare()#Returns chi2. Remember to remove the option XErrorSize(0) from data.PlotOn
cut_chi2_voig = "{:.2f}".format(chi2_voig) #Crop the chi2 to 2 decimal digits
print "Chi square voightian = ",chi2_voig
print "n param voightian = ",nParam_voig
print ""


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
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_signal.pdf")
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Phi/Fit/fit_signal.png")
else:
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_signal.pdf")
    canvas_voig.SaveAs("/eos/user/e/eferrand/ZMesonGamma/CMSSW_10_6_27/src/ZMesonGammaAnalysis/ZTOMesonGamma/plots/Rho/Fit/fit_signal.png")


#create Workspace ------------------------------------------------------------------------------------------------------------------------------
#norm = fileInput.Get("h_ZMass").Integral()  # get the normalization of MC signal (area under MC signal)
norm = h_mZ.Integral(h_mZ.FindBin(84.), h_mZ.FindBin(100.))
print "norm = ", norm
sig_norm = ROOT.RooRealVar(sigPDF_voig.GetName() + "_norm", sigPDF_voig.GetName() + "_norm", norm)

mean.setConstant(1)
width.setConstant(1)
sigma.setConstant(1)
sig_norm.setConstant(1)


workspace = ROOT.RooWorkspace("workspace_"+CHANNEL+"")
getattr(workspace,'import')(sigPDF_voig)
getattr(workspace,'import')(sig_norm)

print "-------------------------------------------"
print "Final print to check the workspace update:"
workspace.Print()


fOutput = ROOT.TFile("workspaces/workspace_"+CHANNEL+".root","RECREATE")
fOutput.cd()
workspace.Write()
#fOutput.Write()
fOutput.Close()


