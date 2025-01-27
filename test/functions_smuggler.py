#######################################################################
#                                                                     #
# <<  Hi! I smuggle bad coded functions. But, ehy, they're free! >>   #
#                                                                     #
#######################################################################                                                                                                                                                       

import ROOT
import math
import os
from array import array
import time


#------- Arrays and reader for the BDT -------#

trk1isoCh_array            = array('f', [0.])
pairIso_array              = array('f', [0.])
#mesonPt_array              = array('f', [0.])
#photonEt_array             = array('f', [0.])
bestPairEta_array          = array('f', [0.])
#photonEta_array            = array('f', [0.])
#nJet_array                 = array('f', [0.])
#JetNeutralEmEnergy_array   = array('f', [0.])
#JetChargedHadEnergy_array  = array('f', [0.])
#JetNeutralHadEnergy_array  = array('f', [0.])
#metPt_array      = array('f', [0.])
#bestJetPt_array  = array('f', [0.])
mesonGammaDeltaPhi_array  = array('f', [0.])

reader = ROOT.TMVA.Reader("!Color")


class Simplified_Workflow_Handler:

    def __init__(self,signalname,dataname,isBDT):

                
        ###################################################################################
        #                                                                                 #
        #------------------------ Add BDT variables to the reader ------------------------#
        #                                                                                 #
        ###################################################################################

        reader.AddVariable("firstTrkIsoCh",trk1isoCh_array)
        reader.AddVariable("pairIso0",pairIso_array)
        #reader.AddVariable("bestPairPt/ZMass",mesonPt_array)
        #reader.AddVariable("photonEt/ZMass",photonEt_array)
        #reader.AddVariable("photonEta",photonEta_array)
        reader.AddVariable("bestPairEta",bestPairEta_array)
        #reader.AddVariable("_photonEt/_HpT",photonEt_array)
        #reader.AddVariable("_bestCouplePt/_HpT",mesonPt_array)
        #reader.AddVariable("_photonEt/_HpT",photonEt_array)
        #reader.AddVariable("_JetNeutralEmEnergy",JetNeutralEmEnergy_array)
        #reader.AddVariable("_JetChargedHadEnergy",JetChargedHadEnergy_array)   
        #reader.AddVariable("_JetNeutralHadEnergy",JetNeutralHadEnergy_array)
        #reader.AddVariable("_metPt",metPt_array)
        #reader.AddVariable("_bestJetPt/mesonGammaMass",bestJetPt_array)
        reader.AddVariable("MesonGammaDeltaPhi",mesonGammaDeltaPhi_array)

        if isBDT:
            reader.BookMVA("BDT","MVA/default/weights/TMVAClassification_BDT.weights.xml")

    #Get BDT output function ###########################################################################################################

    def get_BDT_output(self,trk1isoCh,pairIso,mesonGammaMass,bestPairEta,mesonGammaDeltaPhi):#photonEt,photonEta
        trk1isoCh_array[0] = trk1isoCh
        pairIso_array[0]   = pairIso
        #mesonPt_array[0]   = mesonPt/mesonGammaMass
        #photonEt_array[0]  = photonEt/mesonGammaMass
        #photonEta_array[0]  = photonEta
        bestPairEta_array[0] = bestPairEta
        #JetNeutralEmEnergy_array[0]  = JetNeutralEmEn
        #JetChargedHadEnergy_array[0] = JetChargedHadEn
        #JetNeutralHadEnergy_array[0] = JetNeutralHadEn
        #metPt_array[0]     = metPt
        #bestJetPt_array[0] = bestJetPt/mesonGammaMass
        mesonGammaDeltaPhi_array[0] = mesonGammaDeltaPhi

        return reader.EvaluateMVA("BDT")

    '''#Photon scaling function ###########################################################################################################
    def get_photon_scale(self,ph_pt, ph_eta):

        local_ph_pt = ph_pt
        if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
            local_ph_pt = 499.
        
        local_ph_eta = ph_eta
        if local_ph_eta >= 2.5:
            local_ph_eta = 2.49
        if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
            local_ph_eta = -2.5

        scale_factor_ID = ph_ID_scale_histo_2018.GetBinContent( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )
        ph_ID_err       = ph_ID_scale_histo_2018.GetBinError( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )

        etaBin = 1 if abs(local_ph_eta) < 1.48 else 4
        scale_factor_pixVeto = ph_pixVeto_scale_histo_2018.GetBinContent(etaBin)
        ph_pixVeto_err       = ph_pixVeto_scale_histo_2018.GetBinError(etaBin)
        
        scale_factor = scale_factor_ID * scale_factor_pixVeto
        tot_err      = math.sqrt( scale_factor_pixVeto * scale_factor_pixVeto * ph_ID_err * ph_ID_err + scale_factor_ID * scale_factor_ID * ph_pixVeto_err * ph_pixVeto_err )

        return scale_factor, tot_err
        '''


