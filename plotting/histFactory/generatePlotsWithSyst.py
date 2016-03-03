from cp3_llbb.ZATools.ZACnC import *

def printInPyWithSyst(f, g, name = '', variable = '', cut = '', weight = '', binning = '', norm = "pdf_nominal", writeInPlotIt = 0):
    f.write( "        {\n")
    f.write( "        'name': '"+name+"',\n")
    f.write( "        'variable': '"+variable+"',\n")
    f.write( "        'plot_cut': '"+cut+"',\n")
    f.write( "        'weight': '"+weight+"',\n")
    if (norm != "pdf_nominal" and norm != ""):
        f.write( "        'normalize-to': '"+norm+"',\n")
    f.write( "        'binning': '"+binning+"'\n")
    f.write( "        },\n")
    if writeInPlotIt == 1 :
      g.write("'"+name+"':\n")
      g.write("  x-axis: '"+name+"'\n")
      g.write("  y-axis: 'Evt'\n")
      g.write("  y-axis-format: '%1% / %2$.0f GeV'\n")
      g.write("  normalized: false\n")
      g.write("  log-y: both\n")
      g.write("  save-extensions: ['png','pdf']\n")
      g.write("  show-ratio: true\n")
      #g.write("  no-data: true\n")

# binnings

nPV_binning = "(35,0,35)"
mll_binning = "(10,60,120)"
mbb_binning = "(20,0,600)"
mllbb_binning = "(20,0,1000)"

# mll

mll = "za_diLeptons[0].p4.M()"
mll_SYST = "za_SYST_diLeptons[0].p4.M()"
mllName = "mll"

# mbb

mbb = "za_diJets[0].p4.M()"
mbb_SYST = "za_SYST_diJets[0].p4.M()"
mbbName = "mbb"

# mllbb

mllbb = "za_diLepDiJets[0].p4.M()"
mllbb_SYST = "za_SYST_diLepDiJets[0].p4.M()"
mllbbName = "mllbb"

# PV N

nPV = "za_vertex_ndof.size()"
nPV_SYST = "za_SYST_vertex_ndof.size()"
nPVName = "nVX"

# Cuts

ll_weights = "(event_is_data !=1 ?( za_diLeptons[0].triggerSF * event_pu_weight * event_weight) : 1.0)"


twoLCond = []
twoLCondName = []
twoLCond.append("(za_mumu_Mll_cut && (za_mumu_fire_trigger_Mu17_Mu8_cut || za_mumu_fire_trigger_Mu17_TkMu8_cut) && za_diLeptons[0].isMM && za_diLeptons[0].triggerMatched)")
twoLCond.append("(za_elel_Mll_cut && za_elel_fire_trigger_Ele17_Ele12_cut && za_diLeptons[0].isTT && za_diLeptons[0].triggerMatched)")
twoLCond.append("(za_mumu_Mll_cut ||  za_elel_Mll_cut)")
twoLCond.append("(za_muel_Mll_cut ||  za_elmu_Mll_cut)")
twoLCondName.append("mm")
twoLCondName.append("ee")
twoLCondName.append("ll")
twoLCondName.append("me")

twoLtwoJCond = []
twoLtwoJCondName = []
twoLtwoBCond = []
twoLtwoBCondName = []
twoLthreeBCond = []
twoLthreeBCondName = []
twoLtwoSubJCond = []
twoLtwoSubJCondName = []
twoLOneBFatJetTCond = []
twoLOneBFatJetTCondName = []
twoLTwoBSubJetsLLCond = []
twoLTwoBSubJetsLLCondName = []
twoLTwoBSubJetsMMCond = []
twoLTwoBSubJetsMMCondName = []
twoLTwoBHighMassCond = []
twoLTwoBHighMassCondName = []

basicJcond="(elel_TwoJets_cut || mumu_TwoJets_cut || muel_TwoJets_cut || elmu_TwoJets_cut)"
basicTwoBcond="(Length$(za_diJets) > 0) && ((za_elel_TwoBjets_cut || za_mumu_TwoBjets_cut || za_muel_TwoBjets_cut || za_elmu_TwoBjets_cut) "
basicThreeBcond="(elel_ThreeBjets_cut || mumu_ThreeBjets_cut || muel_ThreeBjets_cut || elmu_ThreeBjets_cut)"
basicSubJcond="(elel_TwoSubJets_cut || mumu_TwoSubJets_cut || muel_TwoSubJets_cut || elmu_TwoSubJets_cut)"
basicOneBFatJetTcond="((elel_OneBFatJetT_cut || mumu_OneBFatJetT_cut || muel_OneBFatJetT_cut || elmu_OneBFatJetT_cut) && (Length$(za_diLepFatJets) > 0))"
#basicOneBFatJetTcond="(Length$(za_diLepFatJets) > 0)"
basicTwoBSubJetsLLcond="(elel_TwoBSubJetsLL_cut || mumu_TwoBSubJetsLL_cut )"
basicTwoBSubJetsMMcond="(elel_TwoSubJetsMM_cut || mumu_TwoSubJetsMM_cut || muel_TwoSubJetsMM_cut || elmu_TwoSubJetsMM_cut)"
basictwoLTwoBHighMasscond="(Length$(za_diJets) > 0 && Length$(za_diLeptons) > 0 && za_diJets[0].p4.Pt() > 200 && za_diLeptons[0].p4.Pt() > 200)"
basicJcondName="jj"
basicTwoBcondName="bb"
basicThreeBcondName="bbb"
basicSubJcondName="fj"
basicOneBFatJetTcondName="bfj"
basicTwoBSubJetsLLcondName="sbjsbj"
basictwoLTwoBHighMasscondName="highmass"

for x in range(0,4):
	twoLtwoJCond.append("("+twoLCond[x]+" && "+basicJcond+")")
	twoLtwoBCond.append("("+twoLCond[x]+" && "+basicTwoBcond+")")
	twoLthreeBCond.append("("+twoLCond[x]+" && "+basicThreeBcond+")")
	twoLtwoSubJCond.append("("+twoLCond[x]+" && "+basicSubJcond+")")
	twoLOneBFatJetTCond.append("("+twoLCond[x]+" && "+basicOneBFatJetTcond+")")
	twoLTwoBSubJetsLLCond.append("("+twoLCond[x]+" && "+basicTwoBSubJetsLLcond+")")
	twoLTwoBHighMassCond.append("("+twoLCond[x]+" && "+basictwoLTwoBHighMasscond+")")
	twoLtwoJCondName.append(twoLCondName[x]+basicJcondName)
	twoLtwoBCondName.append(twoLCondName[x]+basicTwoBcondName)
        twoLthreeBCondName.append(twoLCondName[x]+basicThreeBcondName)
	twoLtwoSubJCondName.append(twoLCondName[x]+basicSubJcondName)
	twoLOneBFatJetTCondName.append(twoLCondName[x]+basicOneBFatJetTcondName)
        twoLTwoBSubJetsLLCondName.append(twoLCondName[x]+basicTwoBcondName+basicTwoBSubJetsLLcondName)
	twoLTwoBHighMassCondName.append(twoLCondName[x]+basicTwoBcondName+basictwoLTwoBHighMasscondName)
	print twoLtwoJCond[x]


#test_highMass_cond = "(Length$(za_diJets) > 0 && Length$(za_diLeptons) > 0 && za_diJets[0].p4.Pt() > 200 && za_diLeptons[0].p4.Pt() > 200)"
#test_highMass =  twoL_cond + " * " + twoB_cond + " * " + test_highMass_cond + weights
#test_highMassName = "highPt"


# Writing the JSON

## 2 Muons 2 Jets :

systematics = {'__jecup':'jecup',
               '__jecdown':'jecdown',
               '__jerup':'jerup',
               '__jerdown':'jerdown'}

cutBtagsMM = ["(za_mumu_DiJetBWP_MM_cut && za_mumu_LooseZCandidate_cut )","(za_elel_DiJetBWP_MM_cut && za_elel_LooseZCandidate_cut)", "((za_mumu_DiJetBWP_MM_cut && za_mumu_LooseZCandidate_cut ) || (za_elel_DiJetBWP_MM_cut && za_elel_LooseZCandidate_cut))"]
cutBtagsMM_SYST = ["za_SYST_mumu_DiJetBWP_MM_cut && za_SYST_mumu_LooseZCandidate_cut","za_SYST_elel_DiJetBWP_MM_cut && za_SYST_elel_LooseZCandidate_cut", "((za_SYST_mumu_DiJetBWP_MM_cut && za_SYST_mumu_LooseZCandidate_cut) || (za_SYST_elel_DiJetBWP_MM_cut && za_SYST_elel_LooseZCandidate_cut))"]

fjson = open('plots_syst.py', 'w')
fjson.write( "plots = [\n")
fyml = open('plots_syst.yml', 'w')

weights = "event_pu_weight * event_weight"
weights_puup = "event_pu_weight_up * event_weight"
weights_pudown = "event_pu_weight_down * event_weight"

llTrigSF = "(event_is_data !=1 ?( za_diLeptons[0].triggerSF * za_diLeptons[0].triggerMatched) : 1.0)"
llTrigSF_SYST = "(event_is_data !=1 ?( za_SYST_diLeptons[0].triggerSF * za_SYST_diLeptons[0].triggerMatched) : 1.0)"

#btagSF = "(event_is_data !=1 ? (jet_sf_csvv2_medium[za_diJets[0].idxJet1][0] * jet_sf_csvv2_medium[za_diJets[0].idxJet2][0] ) : 1.0)"
#btagSFup = "(event_is_data !=1 ? (jet_sf_csvv2_medium[za_diJets[0].idxJet1][1] * jet_sf_csvv2_medium[za_diJets[0].idxJet2][1] ) : 1.0)"
#btagSFdown = "(event_is_data !=1 ? (jet_sf_csvv2_medium[za_diJets[0].idxJet1][2] * jet_sf_csvv2_medium[za_diJets[0].idxJet2][2] ) : 1.0)"

btagSF = "(event_is_data !=1 ?  ( common::combineScaleFactors<2>({{{ jet_sf_csvv2_medium[za_diJets[0].idxJet1][0] , 1 }, { jet_sf_csvv2_medium[za_diJets[0].idxJet2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::NOMINAL) ) : 1.0) ";

btagSF_SYST = "(event_is_data !=1 ?  ( common::combineScaleFactors<2>({{{ jet_SYST_sf_csvv2_medium[za_SYST_diJets[0].idxJet1][0] , 1 }, { jet_SYST_sf_csvv2_medium[za_SYST_diJets[0].idxJet2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::NOMINAL) ) : 1.0) ";


btagSFup = "(event_is_data !=1 ? ( common::combineScaleFactors<2>({{{ jet_sf_csvv2_medium[za_diJets[0].idxJet1][0] , jet_sf_csvv2_medium[za_diJets[0].idxJet1][2] }, { jet_sf_csvv2_medium[za_diJets[0].idxJet2][0] , jet_sf_csvv2_medium[za_diJets[0].idxJet2][2] }}}, {{1, 1}, {1, 1}}, common::Variation::UP) ) : 1.0) ";

btagSFdown = "(event_is_data !=1 ? ( common::combineScaleFactors<2>({{{ jet_sf_csvv2_medium[za_diJets[0].idxJet1][0] , jet_sf_csvv2_medium[za_diJets[0].idxJet1][1] }, { jet_sf_csvv2_medium[za_diJets[0].idxJet2][0] , jet_sf_csvv2_medium[za_diJets[0].idxJet2][1] }}}, {{1, 1}, {1, 1}}, common::Variation::DOWN) ) : 1.0) ";

llIdIsoSF = "(za_diLeptons[0].isElEl ? common::combineScaleFactors<2>({{{ electron_sf_hww_wp[za_diLeptons[0].idxLep1][0] , 1 }, { electron_sf_hww_wp[za_diLeptons[0].idxLep2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::NOMINAL) : common::combineScaleFactors<2>({{{ muon_sf_hww_wp[za_diLeptons[0].idxLep1][0] , 1 }, { muon_sf_hww_wp[za_diLeptons[0].idxLep2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::NOMINAL))"

llIdIsoSF_elup = "(za_diLeptons[0].isElEl ? common::combineScaleFactors<2>({{{ electron_sf_hww_wp[za_diLeptons[0].idxLep1][0] , electron_sf_hww_wp[za_diLeptons[0].idxLep1][2] }, { electron_sf_hww_wp[za_diLeptons[0].idxLep2][0] , electron_sf_hww_wp[za_diLeptons[0].idxLep2][2] }}}, {{1, 1}, {1, 1}}, common::Variation::UP) : common::combineScaleFactors<2>({{{ muon_sf_hww_wp[za_diLeptons[0].idxLep1][0] , 1 }, { muon_sf_hww_wp[za_diLeptons[0].idxLep2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::UP))"

llIdIsoSF_eldown = "(za_diLeptons[0].isElEl ? common::combineScaleFactors<2>({{{ electron_sf_hww_wp[za_diLeptons[0].idxLep1][0] , electron_sf_hww_wp[za_diLeptons[0].idxLep1][1] }, { electron_sf_hww_wp[za_diLeptons[0].idxLep2][0] , electron_sf_hww_wp[za_diLeptons[0].idxLep2][1] }}}, {{1, 1}, {1, 1}}, common::Variation::DOWN) : common::combineScaleFactors<2>({{{ muon_sf_hww_wp[za_diLeptons[0].idxLep1][0] , 1 }, { muon_sf_hww_wp[za_diLeptons[0].idxLep2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::DOWN))"

llIdIsoSF_muup = "(za_diLeptons[0].isElEl ? common::combineScaleFactors<2>({{{ electron_sf_hww_wp[za_diLeptons[0].idxLep1][0] , 1 }, { electron_sf_hww_wp[za_diLeptons[0].idxLep2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::NOMINAL) : common::combineScaleFactors<2>({{{ muon_sf_hww_wp[za_diLeptons[0].idxLep1][0] ,  muon_sf_hww_wp[za_diLeptons[0].idxLep1][2] }, { muon_sf_hww_wp[za_diLeptons[0].idxLep2][0] , muon_sf_hww_wp[za_diLeptons[0].idxLep2][2] }}}, {{1, 1}, {1, 1}}, common::Variation::UP))"

llIdIsoSF_mudown = "(za_diLeptons[0].isElEl ? common::combineScaleFactors<2>({{{ electron_sf_hww_wp[za_diLeptons[0].idxLep1][0] , 1 }, { electron_sf_hww_wp[za_diLeptons[0].idxLep2][0] , 1 }}}, {{1, 1}, {1, 1}}, common::Variation::NOMINAL) : common::combineScaleFactors<2>({{{ muon_sf_hww_wp[za_diLeptons[0].idxLep1][0] ,  muon_sf_hww_wp[za_diLeptons[0].idxLep1][1] }, { muon_sf_hww_wp[za_diLeptons[0].idxLep2][0] , muon_sf_hww_wp[za_diLeptons[0].idxLep2][1] }}}, {{1, 1}, {1, 1}}, common::Variation::DOWN))"



pdf_up = "event_pdf_weight_up"
pdf_down = "event_pdf_weight_down"
scale_0 = "event_scale_weights[0]"
scale_1 = "event_scale_weights[1]"
scale_2 = "event_scale_weights[2]"
scale_3 = "event_scale_weights[3]"
scale_4 = "event_scale_weights[4]"
scale_5 = "event_scale_weights[5]"

llweights = {'': llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__puup': llIdIsoSF+'*'+llTrigSF+'*'+weights_puup,
             '__pudown': llIdIsoSF+'*'+llTrigSF+'*'+weights_pudown,
             '__pdfup': pdf_up+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__pdfdown': pdf_down+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__scale0': scale_0+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__scale1': scale_1+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__scale2': scale_2+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__scale3': scale_3+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__scale4': scale_4+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights,
             '__scale5': scale_5+'*'+llIdIsoSF+'*'+llTrigSF+'*'+weights             
    }

llbbweights = {'': [llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_nominal"],
             '__puup': [llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights_puup, "pdf_nominal"],
             '__pudown': [llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights_pudown, "pdf_nominal"],
             '__btagup': [llIdIsoSF+'*'+llTrigSF+'*'+btagSFup+'*'+weights, "pdf_nominal"],
	     '__btagdown': [llIdIsoSF+'*'+llTrigSF+'*'+btagSFdown+'*'+weights, "pdf_nominal"],
             '__elIDup': [llIdIsoSF_elup+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_nominal"],
             '__elIDdown': [llIdIsoSF_eldown+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_nominal"],
             '__muIDup': [llIdIsoSF_muup+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_nominal"],
             '__muIDdown': [llIdIsoSF_mudown+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_nominal"],
             '__pdfup': [pdf_up+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_up"],
             '__pdfdown': [pdf_down+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "pdf_down"],
             '__scale0': [scale_0+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "scale_0"],
             '__scale1': [scale_1+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "scale_1"],
             '__scale2': [scale_2+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "scale_2"],
             '__scale3': [scale_3+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "scale_3"],
             '__scale4': [scale_4+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "scale_4"],
             '__scale5': [scale_5+'*'+llIdIsoSF+'*'+llTrigSF+'*'+btagSF+'*'+weights, "scale_5"]
    }


## Plots for combine


options = options_()

for x in range(0,3):
    for cutkey in options.cut :
        for s,wpair in llbbweights.iteritems() :
            w = wpair[0]
            n = wpair[1]
            print 'cutkey : ', cutkey
            ### get M_A and M_H ###
            #mH[0] = float(options.mH_list[cutkey])
            #mA[0] = float(options.mA_list[cutkey])

            ### SIGNAL Region ###
            printInPyWithSyst(fjson, fyml, 
                    name = twoLtwoBCondName[x]+"SR"+cutkey+s, 
                    variable = '0.5', 
                    cut = options.cut[cutkey]+" && "+cutBtagsMM[x], 
                    weight = w, 
                    binning = "(1,0,1)",
                    norm=n, 
                    writeInPlotIt = (1 if s==''  else 0)
                    )
            ### BACKGROUND Region ###
            printInPyWithSyst(fjson, fyml, 
                    name = mllName+'_'+twoLtwoBCondName[x]+"BR"+cutkey+s, 
                    variable = mll, 
                    cut = "!"+options.cut[cutkey]+" && "+cutBtagsMM[x], 
                    weight = w, 
                    binning = mll_binning, 
                    norm=n,
                    writeInPlotIt = (1 if s==''  else 0)
                    )
        for s1,s2 in systematics.iteritems() :
            w = llTrigSF_SYST.replace('SYST',s2)+'*'+btagSF_SYST.replace('SYST',s2)+'*'+weights
            ### SIGNAL Region ###
            printInPyWithSyst(fjson, fyml, 
                    name=twoLtwoBCondName[x]+"SR"+cutkey+s1, 
                    variable="0.5", 
                    cut= " !event_is_data && "+options.cut[cutkey]+" && "+cutBtagsMM_SYST[x].replace('SYST',s2), 
                    weight=w, 
                    binning="(1,0,1)", 
                    writeInPlotIt= 0)
            ### BACKGROUND Region ###
            printInPyWithSyst(fjson, fyml,
                    name=mllName+'_'+twoLtwoBCondName[x]+"BR"+cutkey+s1,
                    variable=mll_SYST.replace('SYST',s2),
                    cut= " !event_is_data && !"+options.cut[cutkey]+" && "+cutBtagsMM_SYST[x].replace('SYST',s2),
                    weight=w, 
                    binning=mll_binning, 
                    writeInPlotIt= 0)

        



printInPyWithSyst(fjson, fyml,
            name = 'jet_sf_csvv2_medium',
            variable = 'jet_sf_csvv2_medium.size()',
            cut = cutBtagsMM[x],
            weight = w,
            binning = '(10,0,2)',
            writeInPlotIt = 0
            )




## Control Plots in ll category:
### TO be modified to include the normalization

'''
for x in range(0,2):
    for s,w in llweights.iteritems() : 
	print x, s, w 
        # N of vertices
        #printInPyWithSyst(fjson, fyml, name=nPVName+'_'+twoLCondName[x]+s, variable=nPV, cut=twoLCond[x], weight=w, binning=nPV_binning, writeInPlotIt= (1 if s==''  else 0))
        # M_ll
        printInPyWithSyst(fjson, fyml, name=mllName+'_'+twoLCondName[x]+s, variable=mll, cut=twoLCond[x], weight=w, binning=mll_binning, writeInPlotIt= (1 if s==''  else 0))
'''


printInPyWithSyst(fjson, fyml, name="pdf_weight_down", variable="event_pdf_weight_down", cut=cutBtagsMM[0], weight=llIdIsoSF+'*'+llTrigSF+'*'+weights, binning="(40,0.5,1.5)", norm="", writeInPlotIt= 1)
printInPyWithSyst(fjson, fyml, name="pdf_weight_up", variable="event_pdf_weight_up", cut=cutBtagsMM[0], weight=llIdIsoSF+'*'+llTrigSF+'*'+weights, binning="(40,0.5,1.5)", norm="", writeInPlotIt= 1)

## Control Plots in llbb category:

for x in range(0,2):
    for s,wpair in llbbweights.iteritems() :
        w = wpair[0]
        n = wpair[1]
        print x, s, w, cutBtagsMM[x]
        # M_ll
        printInPyWithSyst(fjson, fyml, name=mllName+'_'+twoLtwoBCondName[x]+s, variable=mll, cut=cutBtagsMM[x], weight=w, binning=mll_binning, norm=n, writeInPlotIt= (1 if s==''  else 0))
        # M_bb
        printInPyWithSyst(fjson, fyml, name=mbbName+'_'+twoLtwoBCondName[x]+s, variable=mbb, cut=cutBtagsMM[x], weight=w, binning=mbb_binning, norm=n, writeInPlotIt= (1 if s==''  else 0))
        # M_llbb
        printInPyWithSyst(fjson, fyml, name=mllbbName+'_'+twoLtwoBCondName[x]+s, variable=mllbb, cut=cutBtagsMM[x], weight=w, binning=mllbb_binning, norm=n, writeInPlotIt= (1 if s==''  else 0))

    for s1,s2 in systematics.iteritems() :
        w = llTrigSF_SYST.replace('SYST',s2)+'*'+btagSF_SYST.replace('SYST',s2)+'*'+weights
        # M_ll
        printInPyWithSyst(fjson, fyml, name=mllName+'_'+twoLtwoBCondName[x]+s1, variable=mll_SYST.replace('SYST',s2), cut= " !event_is_data && "+ cutBtagsMM_SYST[x].replace('SYST',s2), weight=w, binning=mll_binning, writeInPlotIt= (1 if s==''  else 0))
        # M_bb
        printInPyWithSyst(fjson, fyml, name=mbbName+'_'+twoLtwoBCondName[x]+s1, variable=mbb_SYST.replace('SYST',s2), cut= " !event_is_data && "+ cutBtagsMM_SYST[x].replace('SYST',s2), weight=w, binning=mbb_binning, writeInPlotIt= (1 if s==''  else 0))
        # M_llbb
        printInPyWithSyst(fjson, fyml, name=mllbbName+'_'+twoLtwoBCondName[x]+s1, variable=mllbb_SYST.replace('SYST',s2), cut= " !event_is_data && "+ cutBtagsMM_SYST[x].replace('SYST',s2), weight=w, binning=mllbb_binning, writeInPlotIt= (1 if s==''  else 0))


fjson.write( "        ]\n")

