#!/usr/bin/env python3


import glob
import ROOT
import numpy as np

basePath = "/storage/af/group/btl-production/upload/QAQC_SM/qaqc-gui_output/SM_QAQC_Production/"

def get_TFiles(minRunNum = 0, maxRunNum=10000, moduleNum = 32110020008497):
    '''
    return list of TFiles from desired modules, with minimum run number and maximum run number specified
    '''
    files_list = glob.glob(basePath+"*/module_"+str(moduleNum)+"_analysis_both_calibs.root")
    files_list_filtered = []
    for file in files_list:
        file_str = str(file)
        runNum = int(file[file.find("run")+3:file.rfind('/')])
        if runNum>=minRunNum and runNum<=maxRunNum:
            files_list_filtered.append(ROOT.TFile(file_str))
    return files_list_filtered

def get_graph_list(TFileList = [], histtype = 'g_L_light_yield_vs_bar'):
    '''
    return list where each entry is the designated histogram corresponding to the input TFile
    '''
    return [TFile.Get(histtype) for TFile in TFileList]


def average_tests_inverseVariance(TFileList, histtype = 'g_L_light_yield_vs_bar', invertCh = False):
    '''
    return TGraph where the ith point is the weighted average of the ith point of each graph in the list.
    Weighting is done via inverse variance
    All TGraphs from TFiles in list must have the same number of points
    '''
    TGraphList = get_graph_list(TFileList, histtype)
    NPoints = TGraphList[0].GetN()
    if not invertCh:
        iter = range(NPoints)
    else:
        iter = range(NPoints-1, -1, -1)
    inverseVariance_TGraph  = ROOT.TGraphErrors()
    for i, ch in enumerate(iter):
        weighted_vals_num = []
        errors = []
        for TGraph in TGraphList:
            weighted_vals_num.append(TGraph.GetPointY(ch)/TGraph.GetErrorY(ch)**2)
            errors.append(1/TGraph.GetErrorY(ch)**2)
            #print(TGraph.GetPointY(ch))
            #weighted_vals_num.append(TGraph.GetPointY(ch))
        inverseVariance_TGraph.SetPointX(i,i)
        #inverseVariance_TGraph.SetPointY(i, np.average(np.array(weighted_vals_num)))
        inverseVariance_TGraph.SetPointY(i, np.sum(np.array(weighted_vals_num))/np.sum(np.array(errors)))
        inverseVariance_TGraph.SetPointError(i, 0, (1/np.sum(np.array(errors)))**0.5)
    return inverseVariance_TGraph

def graphByBar(TFileList, leftSideGraphName, rightSideGraphName, averagedGraphName):
    
    graph = average_tests_inverseVariance(TFileList, leftSideGraphName)
    graph.SetLineColor(ROOT.kRed)
    graph_right = average_tests_inverseVariance(TFileList, rightSideGraphName, True)
    graph_right.SetLineColor(ROOT.kBlue)
    graph_both = average_tests_inverseVariance(TFileList, averagedGraphName)
    graph_both.SetLineColor(ROOT.kBlack)

    graph.SetMarkerStyle(20)
    graph_right.SetMarkerStyle(20)
    graph_both.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kRed)
    graph_right.SetMarkerColor(ROOT.kBlue)
    graph_both.SetMarkerColor(ROOT.kBlack)
    graph.SetMarkerSize(0.5)
    graph_right.SetMarkerSize(0.5)
    graph_both.SetMarkerSize(0.5)

    leg = ROOT.TLegend(0.15, 0.15, 0.3, 0.3)
    leg.AddEntry(graph, "left side", "lep")
    leg.AddEntry(graph_right, "right side", "lep")
    leg.AddEntry(graph_both, "combined", "lep")

    
    return graph,graph_right, graph_both, leg

   
def uncertainty_by_channel(TFileList, histtype='g_lyso_pc_per_kev_vs_ch'):
    '''
    return TGraphErrors where the data point is the standard deviation in the measurements from the TFileList
    All TGraphs from TFiles in list must have the same number of points 
    '''
    TGraphList = get_graph_list(TFileList, histtype)
    std_TGraph  = ROOT.TGraphErrors()
    NPoints = TGraphList[0].GetN()
    for i in range(NPoints):
        vals_byChannel = []
        for TGraph in TGraphList:
            vals_byChannel.append(TGraph.GetPointY(i))
        std_TGraph.SetPointX(i,i)
        std_val = np.std(np.array(vals_byChannel))/np.average(np.array(vals_byChannel))*100
        #print(std_val)
        std_TGraph.SetPointY(i, std_val)
        std_TGraph.SetPointError(i, 0, (std_val/np.sqrt(len(TGraphList))))
    return std_TGraph

def get_LYSO_pedestal_width(TFileList):
    '''
    return trial-averaged pedestal values for each channel from LYSO events
    '''
    #TGraphList = get_graph_list(TFileList, histtype)
    pedestalWidth_TGraph  = ROOT.TGraphErrors()
    #NPoints = TGraphList[0].GetN()
    for ch in range(32):
        vals_byChannel = []
        errors = []
        for TFile in TFileList:
            if ch<10:
                ch_str = '0'+str(ch)
            else:
                ch_str = ch
            fit = TFile.Get(f"lyso_ch{ch_str}_pedestal_offset_fit")
            vals_byChannel.append(fit.GetParameter(2)/fit.GetParError(2)**2)
            errors.append(1/fit.GetParError(2)**2)
        pedestalWidth_TGraph.SetPointX(ch,ch)
        std_val = np.average(np.array(vals_byChannel))/np.average(np.array(vals_byChannel))*100
        #print(std_val)
        pedestalWidth_TGraph.SetPointY(ch, np.sum(np.array(vals_byChannel))/np.sum(np.array(errors)))
        pedestalWidth_TGraph.SetPointError(ch, 0, (1/np.sum(np.array(errors)))**0.5)
    return pedestalWidth_TGraph