import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("USER")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #Could be not the coprrect one, but should contain the one without "condDBv2"
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.parseArguments()

################################################################################################################
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,
#                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
#                       era='2018-Prompt')  
################################################################################################################



#Input source
if options.runningOnData:
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v33') # OLD ONE : 102X_dataRun2_Sep2018ABC_v2
    inputFiles = { '/store/data/Run2018B/Tau/MINIAOD/UL2018_MiniAODv2-v2/70000/09FD3540-D36B-6249-9817-D3BAC2F02E74.root'}
    
else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1')  # OLD ONE : 102X_upgrade2018_realistic_v18
   inputFiles={'file:/eos/user/p/pellicci/MesonGamma_root/2023/Zrhogamma_miniAOD/Zrhogamma_2018UL_12.root','file:/eos/user/p/pellicci/MesonGamma_root/2023/Zrhogamma_miniAOD/Zrhogamma_2018UL_11.root'} 
   #inputFiles={'file:/eos/user/p/pellicci/MesonGamma_root/2023/Zphigamma_miniAOD/Zphigamma_2018UL_12.root','file:/eos/user/p/pellicci/MesonGamma_root/2023/Zphigamma_miniAOD/Zphigamma_2018UL_11.root'}
   input_path = '/eos/user/p/pellicci/MesonGamma_root/2023/Zrhogamma_miniAOD/'
   #input_path = '/eos/user/p/pellicci/MesonGamma_root/2023/Zphigamma_miniAOD/'



#INPUT FILE LIST
'''                                                                                                                                                                                                   
    For the given path, get the List of all files in the directory tree                                                                                                                               
'''                                                                                                                                                                                                   
def getListOfFiles(dirName):                                                                                                                                                                          
    # create a list of file and sub directories                                                                                                                                                       
    # names in the given directory                                                                                                                                                                    
    listOfFile = os.listdir(dirName)                                                                                                                                                                  
    allFiles = list()                                                                                                                                                                                 
    # Iterate over all the entries                                                                                                                                                                    
    for entry in listOfFile:                                                                                                                                                                          
        # Create full path                                                                                                                                                                            
        fullPath = os.path.join(dirName, entry)                                                                                                                                                       
        # If entry is a directory then get the list of files in this directory                                                                                                                        
        if os.path.isdir(fullPath):                                                                                                                                                                   
            allFiles = allFiles + getListOfFiles(fullPath)                                                                                                                                            
        else:                                                                                                                                                                                         
            allFiles.append("file:" + fullPath)                                                                                                                                                                 
                                                                                                                                                                                                      
    return allFiles     

# Get the list of all files in directory tree at given path                                                                                                                                           
listOfFiles = getListOfFiles(input_path) ######                                                                                                                                                          
print(listOfFiles)  

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles), #inputFiles or listOfFiles
                             duplicateCheckMode = cms.untracked.string ('noDuplicateCheck')
                             )


# Output file
if options.runningOnData:
    process.TFileService = cms.Service("TFileService", fileName = cms.string("ZMesonGamma_Data.root"))
else:
    process.TFileService = cms.Service("TFileService", fileName = cms.string("ZMesonGamma_Signal.root"))



###############################################################################################################################
#                                                                                                                               #
#-------------------------------------------------------- Jet corrections ----------------------------------------------------#
#                                                                                                                             #
#                                      https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC                                     #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?redirectedfrom=CMS.WorkBookJetEnergyCorrections #
#                                                                                                                             #
###############################################################################################################################

#Use the latest JECs
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if not options.runningOnData:      #This loop is for MC
   print "Jet Energy Corrections on Monte Carlo will be applied"
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')# The right corrections are taken through the chosen global tag
else:
   print "Jet Energy Corrections on Data will be applied "
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')# The right corrections are taken through the chosen global tag

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = jetCorrectionsList
)




process.load("ZMesonGammaAnalysis.ZTOMesonGamma.ZMesonGammaGen_cfi")
process.ZMesonGamma.runningOnData = options.runningOnData

# Apply JEC to both MC and data
process.ZMesonGamma.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

###############################################
#                                             #
#------------------ Trigger ------------------#
#                                             #
###############################################


#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for muons
#process.trigger_filter = hlt.triggerResultsFilter.clone()
#process.trigger_filter.triggerConditions = cms.vstring('HLT_Photon35_TwoProngs35_v*')
#process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
#process.trigger_filter.l1tResults = cms.InputTag("")
#process.trigger_filter.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

process.seq = cms.Path(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.ZMesonGamma)

process.schedule = cms.Schedule(process.seq)