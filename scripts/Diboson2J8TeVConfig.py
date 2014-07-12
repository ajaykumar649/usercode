from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray, kYellow

def theConfig(Nj, mH, isElectron = False, initFile = [], includeSignal = True,
        btagged = False, includeExtras = False):
    pars = Wjj2DFitterPars()
    #pars.MCDirectory = "root://cmseos:1094//eos/uscms/store/user/lnujj/RDtrees_co_11Feb14/"
    pars.MCDirectory ="/uscms_data/d3/ajay/test/forPhil/CMSSW_5_3_2_patch4/src/ElectroWeakAnalysis/VPlusJets/test/OutData/"
    pars.SecondaryDirectory = pars.MCDirectory
    pars.isElectron = isElectron
    pars.btagSelection = btagged
    pars.boostedSelection = False
    pars.useTopSideband = True
    pars.initialParametersFile = initFile
    pars.extras = []
    pars.Njets = Nj
    pars.mHiggs = 126.
    if isElectron:
        flavorString = 'el'
    else:
        flavorString = 'mu'
    pars.cuts ='(hvbf_event && hvbf_wjj_m > 43.0 && hvbf_wjj_m < 127.0 &&  hvbf_waj_pt>30. && hvbf_wbj_pt >30. && hvbf_aj_pt >30. && hvbf_bj_pt >30.)'+\
            '&& (hvbf_event_met_pfmet >25. && hvbf_wjj_deta < 1.5 )'#+\
            #'&& (mvavbfWjetsmu >0.3 && mvavbfTopmu >0.3)'

    if isElectron:
	pars.cuts ='(hvbf_event && hvbf_wjj_m > 43.0 && hvbf_wjj_m < 127.0 && hvbf_waj_pt>30. && hvbf_wbj_pt >30. && hvbf_aj_pt >30. &&  hvbf_bj_pt>30.)'+\
                '&& (hvbf_event_met_pfmet >30. && hvbf_wjj_deta < 1.5  )'#+\
                #'&& (mvavbfWjetsel >0.3 && mvavbfTopel >0.3)'

        #implement topselection, btagged or anti-btagged cuts
    pars.btagVeto = False
    if pars.useTopSideband:
        pars.cuts += '&&(hvbf_topWm>0.0  )'
    else:
        #pars.cuts +='&&(hvbf_wjet1_btagCSV  <0.679 && hvbf_wjet2_btagCSV < 0.679 )'
        #pars.cuts +='&&(hvbf_wjet1_btagCSV  <0.244 && hvbf_wjet2_btagCSV < 0.244 )'
        pars.cuts +='&&(hvbf_event==1 && bjet_veto==0 )' 
	#if isElectron:
	#	pars.cuts += '&&( hvbf_topWm > 43.&&  hvbf_topWm < 127.)'
    if pars.useTopSideband:
        pars.backgrounds = ['top', 'WpJ']
        pars.yieldConstraints = {}
        #pars.yieldConstraints = {'top' : 0.0001}
    else:
        pars.backgrounds = ['VBF_WW', 'top', 'WpJ']
        pars.yieldConstraints = {'top' : 0.000007, 'WpJ' : 0.2}
	#pars.yieldConstraints = {'top' : 0.20, 'WpJ' : 0.20} ##Relaxed Constraint Cross-Check
        pars.constrainShapes = ['WpJ']
    if pars.btagSelection:
        if includeExtras:
            pars.extras = ['WZ']
        pars.backgrounds = ['VBF_WW', 'WHbb', 'top', 'WpJ']
        pars.yieldConstraints = {'WHbb' : 0.000001, 'top' : 0.07, 'WpJ' : 0.05}
	##        pars.yieldConstraints = {'WHbb' : 0.000001, 'top' : 0.20, 'WpJ' : 0.20} ##Relaxed Constraint Cross-Check
        pars.constrainShapes = ['WpJ']
    # you need a files entry and a models entry for each of the fit
    # compoents in backgrounds and signals
    # the files should a list with entries like (filename, Ngen, xsec)
	#################### Global Convolution Models ####################
    #pars.GlobalConvModels=[27]
    pars.GlobalConvModels=[-1]
    pars.GlobalConvModelsAlt=pars.GlobalConvModels
    #####################  VBF_WW: #######################################
    pars.VBF_WWFiles = [
        (pars.MCDirectory + 'RD_%s_phantom_CMSSW532.root' % (flavorString),
         333061, 0.0776*(3.0/2.0)),
        ]
    pars.VBF_WWFracOfData = -1
    pars.VBF_WWModels = [5]
    #pars.VBF_WWModels = [-1]
    pars.VBF_WWModelsAlt = pars.VBF_WWModels
    pars.VBF_WWConvModels = pars.GlobalConvModels
    pars.VBF_WWConvModelsAlt = pars.VBF_WWConvModels
    ### WZ separately: ###
    pars.WZFiles = [(pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
                     10000267, 22.88),
                    ]
    pars.WZFracOfData = -1
    pars.WZModels = [13]
    pars.WZModelsAlt = pars.WZModels
    pars.WZConvModels = pars.GlobalConvModels
    pars.WZConvModelsAlt = pars.WZConvModels
    #####################  WpJ: #######################################
    wpj_kfactor = 1.16
    if pars.btagSelection:
        wpj_kfactor = 2.05 # from arXiv:1011.6647 [hep-ph]
    pars.WpJFiles = [
        # (pars.MCDirectory + 'RD_%s_WpJ_CMSSW532.root' % (flavorString),
        #  18353019+50768992, 36257.2),
         (pars.MCDirectory + 'RD_%s_W1Jets_CMSSW532.root' % (flavorString),
          19871598, 5400.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_W2Jets_CMSSW532.root' % (flavorString),
         33004921, 1750.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_W3Jets_CMSSW532.root' % (flavorString),
         15059503, 519.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_W4Jets_CMSSW532_old.root' % (flavorString),
         4369420, 214.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_ZpJ_CMSSW532.root' % (flavorString),
         30209426, 3503.71),
        ]
    #### Cross check using ScaleUp and ScaleDown samples as default ####
    ##set the scale factors such that the expected yield matches the deafult
    ##     crosscheck_wpjSU_kfactor=117078.9/93123.7;
    ##     pars.WpJFiles = [(pars.MCDirectory + 'RD_%s_WpJscaleup_CMSSW532.root' % (flavorString), 20784694, crosscheck_wpjSU_kfactor*36257.2), ]
    ##     crosscheck_wpjSD_kfactor=117078.9/199303.4;
    ##     pars.WpJFiles = [(pars.MCDirectory + 'RD_%s_WpJscaledown_CMSSW532.root' % (flavorString), 20760830, crosscheck_wpjSD_kfactor*36257.2), ]

    # To implement Template Morphing set pars.WpJModels=[-2]
    # and be sure to edit the WpJ*InputParameters.txt file so that the
    # naming scheme corresponds to the correct components/subcomponents.
    # E.g. the parameters from the shape fit to WpJ default MC should
    # contain the suffix Nom, while the overall yield shouldn't, and the
    # lines
    # fMU_WpJ = 0.0 +/- 100.0 L(-1 - 1)
    # fSU_WpJ = 0.0 +/- 100.0 L(-1 - 1)
    # should be added to the .txt file

    # pars.WpJModels = [-2]
    if pars.btagSelection:
        pars.WpJModels = [23]
        #pars.WpJModelsAlt = [10]
        #pars.WpJModels = [14]
        #pars.WpJModelsAlt = [10]
        pars.WpJAuxModels = [3]
        if isElectron:
            #pars.WpJFracOfData = 0.2098 #medium
            pars.WpJFracOfData = 0.51145 #loose
            #pars.WpJFracOfData = 0.4831 #loose with QCD
        else:
            #pars.WpJFracOfData = 0.1935 #medium
            pars.WpJFracOfData = 0.48379 #loose
        #pars.WpJFracOfData = -1
    else:
        #pars.WpJModels = [10]#with Wjet mva
        pars.WpJModels = [8] #without mva and with top mva
        #pars.WpJModelsAlt = [24]
        #pars.WpJAuxModelsAlt = [3]
        #pars.WpJModels = [-2]
        #pars.WpJModels = [23]
        pars.WpJModelsAlt = [10]
        #pars.WpJModels = [14]
        #pars.WpJModelsAlt = [10]
        pars.WpJAuxModels = [3]
        pars.WpJFracOfData = 0.5719#0.7268#0.55275#0.68986#Wjetsmva #0.7471#topmva    0.6974#nomva 0.564#0.58
        if isElectron:
            pars.WpJFracOfData =0.5877#0.71748#Wjetsmva  #0.7852#topmva  0.7478#nomva 0.5714#0.569
        #pars.WpJFracOfData = -1
    pars.WpJ_NomFiles = pars.WpJFiles
    #pars.WpJ_NomModels = [23]
    pars.WpJ_NomModels = [-1]
    pars.WpJ_NomAuxModels = [5]
    pars.WpJ_MUFiles = [ (pars.MCDirectory + 'RD_%s_WpJmatchingup_CMSSW532.root' % (flavorString), 20976007, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJmatchingup_CMSSW532_UserGen.root' % (flavorString), 20976007, 36257.2), ]##the centrally produced and user generated samples are added with equal weights
    pars.WpJ_MUModels = [-1]
    pars.WpJ_MDFiles = [ (pars.MCDirectory + 'RD_%s_WpJmatchingdown_CMSSW532.root' % (flavorString), 21364575, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJmatchingdown_CMSSW532_UserGen.root' % (flavorString), 21364575, 36257.2), ]
    pars.WpJ_MDModels = [-1]
    pars.WpJ_SUFiles = [ (pars.MCDirectory + 'RD_%s_WpJscaleup_CMSSW532.root' % (flavorString), 20784694, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJscaleup_CMSSW532_UserGen.root' % (flavorString), 20784694, 36257.2), ]
    pars.WpJ_SUModels = [-1]
    pars.WpJ_SDFiles = [ (pars.MCDirectory + 'RD_%s_WpJscaledown_CMSSW532.root' % (flavorString), 20760830, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJscaledown_CMSSW532_UserGen.root' % (flavorString), 20760830, 36257.2), ]
    pars.WpJ_SDModels = [-1]
    ##     #Fitting the alternate samples as a cross check (& to potentially inflate the fit errors)
    ##     pars.WpJ_MUModels = [10]
    ##     pars.WpJ_MDModels = [10]
    ##     pars.WpJ_SUModels = [10]
    ##     pars.WpJ_SDModels = [10]
    pars.WpJ_NomModelsAlt = pars.WpJ_NomModels
    pars.WpJ_NomAuxModelsAlt = pars.WpJ_NomAuxModels
    pars.WpJ_MUModelsAlt = pars.WpJ_MUModels
    pars.WpJ_MDModelsAlt = pars.WpJ_MDModels
    pars.WpJ_SUModelsAlt = pars.WpJ_SUModels
    pars.WpJ_SDModelsAlt = pars.WpJ_SDModels
    pars.WpJConvModels = pars.GlobalConvModels
    pars.WpJConvModelsAlt = pars.WpJConvModels
    if pars.useTopSideband:
        pars.WpJModels = [-1]
        pars.WpJFracOfData = -1
        #pars.WpJConvModels = [-1]
        #####################  ZpJ: #######################################
    pars.ZpJFiles = [
        (pars.MCDirectory + 'RD_%s_ZpJ_CMSSW532.root' % (flavorString),
         30209426, 3503.71),
        ]
    pars.ZpJFracOfData = -1
    #pars.ZpJFracOfData = 0.045
    pars.ZpJModels = [14]
    #pars.ZpJModels = [-1]
    #pars.ZpJAuxModels = [4]
    if pars.btagSelection:
        pars.ZpJModels = [0]
        pars.ZpJAuxModels = [3]
        pars.ZpJFracOfData = -1

    pars.ZpJModelsAlt = pars.ZpJModels
    pars.ZpJConvModels = pars.GlobalConvModels
    pars.ZpJConvModelsAlt = pars.ZpJConvModels
    #####################  top: #######################################
    if pars.useTopSideband:
	ttbar_kfactor=1.0
    else:
        if isElectron:
		ttbar_kfactor=1.1226#1.0503#Wjetsmva #1.0945 topmva # 1.1226 nomva # 0.9045;  #1.3841#1.1105 # update it after calculating from top region
        else:
		ttbar_kfactor=0.96521#0.93686#Wjetsmva    #0.9494 topmva # 0.96521 nomva  1.0686;  #1.4862#1.1180 #update it after calculating from top region
    pars.topFiles = [
        (pars.MCDirectory + 'RD_%s_TTbar_CMSSW532.root' % (flavorString),
         6893735, 225.197*ttbar_kfactor),
        (pars.MCDirectory + 'RD_%s_STopS_Tbar_CMSSW532.root' % (flavorString),
         139974, 1.75776*ttbar_kfactor),
        (pars.MCDirectory + 'RD_%s_STopS_T_CMSSW532.root' % (flavorString),
         259960, 3.89394*ttbar_kfactor),
        (pars.MCDirectory + 'RD_%s_STopT_Tbar_CMSSW532.root' % (flavorString),
         1935066, 30.0042*ttbar_kfactor),
        (pars.MCDirectory + 'RD_%s_STopT_T_CMSSW532.root' % (flavorString),
         3758221, 55.531*ttbar_kfactor),
        (pars.MCDirectory + 'RD_%s_STopTW_Tbar_CMSSW532.root' % (flavorString),
         493458, 11.1773*ttbar_kfactor),
        (pars.MCDirectory + 'RD_%s_STopTW_T_CMSSW532.root' % (flavorString),
         497657, 11.1773*ttbar_kfactor),
        ]

    if pars.btagSelection:
        if isElectron:
            #pars.topFracOfData = 0.7722 #medium
            pars.topFracOfData = 0.47349 #loose
            #pars.topFracOfData = 0.4472 #loose with QCD
        else:
            #pars.topFracOfData = 0.7914 #medium
            pars.topFracOfData = 0.50073 #loose
        #pars.topFracOfData = -1
        pars.topModels = [23]
        pars.topAuxModels = [3]
    else:
        pars.topFracOfData = -1
        #pars.topFracOfData = 0.41
        pars.topModels = [2724] #anti-btag selection
        #pars.topModels = [23] #anti-btag selection
        pars.topAuxModels = [3]
    #pars.topModels = [-1]
    pars.topModelsAlt = pars.topModels
    pars.topAuxModelsAlt = pars.topAuxModels
    pars.topConvModels = pars.GlobalConvModels
    pars.topConvModelsAlt = pars.topConvModels

    if pars.useTopSideband:
        pars.topModels = [-1]
        pars.topFracOfData = -1
        #####################  WHbb: #######################################
    pars.WHbbFiles = [
        (pars.MCDirectory + 'RD_%s_WH_WToLNu_HToBB_M-125_CMSSW532.root' % (flavorString), 999998, 0.6966*0.577*(0.1075+0.1057+0.1125)),
        ]

    pars.WHbbFracOfData = -1
    pars.WHbbModels = [5]

    pars.WHbbModelsAlt = pars.WHbbModels
    pars.WHbbConvModels = pars.GlobalConvModels
    pars.WHbbConvModelsAlt = pars.WHbbConvModels
    ###################################################################
    pars.WZPlotting = {'color' : kGreen+3, 'title' : 'WZ'}
    pars.VBF_WWPlotting = {'color' : kAzure+8, 'title' : 'VBF WW'}
    pars.WpJPlotting = { 'color' : kRed, 'title' : 'V+jets'}
    pars.ZpJPlotting = { 'color' : kYellow, 'title' : 'Z+jets'}
    pars.topPlotting = {'color' : kGreen+2, 'title' : 'top'}
    pars.QCDPlotting = {'color' : kGray, 'title' : 'multijet'}
    pars.WHbbPlotting = {'color' : kBlue, 'title' : 'WHbb'}

    pars.var = ['hvbf_wjj_m']
    pars.varRanges = {'hvbf_wjj_m': (14, 43., 127., []),}
    pars.sigRegionMin = 65.0
    pars.sigRegionMax = 95.0

    if pars.btagSelection:
        pars.varRanges = {'hvbf_wjj_m': (15, 40., 160., []),
                          }
    pars.varTitles = {'hvbf_wjj_m': 'm_{jj}',
                      }
    pars.varNames = {'hvbf_wjj_m': 'hvbf_wjj_m' }


    if pars.useTopSideband:
        pars.var = ['hvbf_topWm']
        pars.varRanges = {'hvbf_topWm': (14, 43., 127., []),}
        pars.varTitles = {'hvbf_topWm': 'Top m_{jj}',}
        pars.varNames = {'hvbf_topWm': 'hvbf_topWm' }
    pars.exclude = {}
    pars.blind = False
    # pars.v1binEdges = [50, 55.,60.,65.,70.,75.,80.,85.,95.,
    #                    105.,115.,125.,135.,150.,165.,180.,200.]
    # pars.v1nbins = len(pars.v1binEdges)-1

    pars.integratedLumi = 19300.

    # pars.binData = False
    pars.binData = False

    #Standard vs QCD cuts:
    pars.QCDcuts = pars.cuts
    #pars.cuts += '&&(event_met_pfmet>25)'
    pars.QCDcuts += '&&(event_met_pfmet>20)'


    pars.includeSignal = includeSignal
    pars.signals = []

    return customizeElectrons(pars) if isElectron else \
        customizeMuons(pars)

def customizeElectrons(pars):
    pars.DataFile = pars.MCDirectory + 'RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root'
    #if ((not pars.useTopSideband) and (not pars.btagSelection)):
    #    pars.QCDFiles = [(pars.SecondaryDirectory + 'RDQCD_WenuJets_Isog0p3NoElMVA_19p2invfb.root',
    #                      1,1), #The events come from the data sideband
    #                     ]
    #pars.backgrounds.append('QCD')
    ##         if pars.btagSelection:
    ##             pars.QCDFracOfData = 0.054
    ##             pars.QCDModels = [0]
    ##             pars.yieldConstraints['QCD'] = 0.5
    #     pars.QCDFracOfData = 0.07
    #     pars.QCDModels = [17]
    #pars.QCDModels = [-1]
    #     pars.yieldConstraints['QCD'] = 0.5
    #     pars.QCDModelsAlt = pars.QCDModels
    #     pars.QCDConvModels = pars.GlobalConvModels
    #     pars.QCDConvModelsAlt = pars.QCDConvModels


    pars.doEffCorrections = True
    pars.effToDo = ['lepton']
    pars.leptonEffFiles = {
        'id': ["EffTable2012/scaleFactor-Run2012ABC-GsfElectronToId.txt"],
        'reco': ["EffTable2012/scaleFactor-Run2012ABC-SCToElectron.txt"],
        'HLT': ["EffTable2012/efficiency-Run2012ABC-WP80ToHLTEle.txt"]
        }
    pars.lumiPerEpoch = [pars.integratedLumi]

    pars.integratedLumi = 19200.
    #pars.QCDcuts += '&&(W_electron_pt>30)'
    pars.cuts += '&&(W_electron_pt>30)'
    return pars

def customizeMuons(pars):

    pars.DataFile = pars.MCDirectory + 'RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root'
    pars.doEffCorrections = True
    pars.effToDo = ['lepton']
    pars.leptonEffFiles = {
        'id': ["EffTable2012/scaleFactor-Run2012ABC-RecoToIso.txt"],
        'HLT': ["EffTable2012/efficiency-Run2012ABC-IsoToIsoMuHLT.txt"]
        }
    pars.lumiPerEpoch = [pars.integratedLumi]

    #pars.QCDcuts += '&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)'
    pars.cuts += '&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)'
    return pars
