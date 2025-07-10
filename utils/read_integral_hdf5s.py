#!/usr/bin/env python3

import sys
sys.path.append(".")


import glob
import ROOT
import numpy as np
import h5py as h5


basePath = "/storage/af/group/btl-production/upload/QAQC_SM/qaqc-gui_output/SM_QAQC_Production/"
num_channels = 32
triggerGroup1_channelRange = (0,8)
triggerGroup2_channelRange = (8,16)
triggerGroup3_channelRange = (16,24)
triggerGroup4_channelRange = (24,32)
triggerGroup_Ranges = [triggerGroup1_channelRange, triggerGroup2_channelRange, triggerGroup3_channelRange, triggerGroup4_channelRange]
testingSlot = 0
calib_file_path = '/storage/af/group/btl-production/upload/QAQC_SM/qaqc-gui_output/calibs_sodium/master_calib.root'

def get_hdf5_list(minRunNum = 0, maxRunNum=10000, moduleNum = 32110020008497):
    '''
    return list of hdf5 file from desired modules, with minimum run number and maximum run number specified
    '''
    files_list = glob.glob(basePath+"*/module_"+str(moduleNum)+"_integrals.hdf5")
    files_list_filtered = []
    for file in files_list:
        file_str = str(file)
        runNum = int(file[file.find("run")+3:file.rfind('/')])
        if runNum>=minRunNum and runNum<=maxRunNum:
            files_list_filtered.append(file_str)
    return files_list_filtered


def get_total_calibrations(type="lyso"):
    '''
    return dictionary with calibration value (based on slot and channel) for a given testing channel
    '''
    calibs = {}
    calib_graph = ROOT.TFile(calib_file_path).Get(f"g_{type}_slot{testingSlot}")
    for ch in range(num_channels):
        calibs[ch] = calib_graph.GetPointY(ch)
    return calibs

def get_LYSO_pedestals(hdf5_file_name):
    '''
    return pedestal values for each channel from the run associated with the provided HDF5 file
    '''
    ROOT_file = ROOT.TFile(hdf5_file_name.replace("_integrals.hdf5", "_analysis_both_calibs.root"))
    pedestal_vals = {}
    for ch in range(num_channels):
        if ch<10:
                ch_str = '0'+str(ch)
        else:
            ch_str = str(ch)
        fit = ROOT_file.Get(f"lyso_ch{ch_str}_pedestal_offset_fit")
        pedestal_vals[ch] = fit.GetParameter(1)
    return pedestal_vals

def getChargeArrays_PedestalAdded(hdf5_file_name, pedestal_dict, calib_dict, type='lyso'):
    '''
    return array with integrals with pedestals added (subtacted, but pedestals are typically negative) and calibrated based on testing slot and channel
    there are four subarray (one for each trigger group), and within each subarray there is the integrated charge by channel per event
    '''
    charge_files = f = h5.File(hdf5_file_name, 'r')
    triggerGroupArrs = []
    
    for triggerGroupNum, triggerGroupRange in enumerate(triggerGroup_Ranges):
        triggerGroupCharges = []
        for ch in range(triggerGroupRange[0], triggerGroupRange[1]):
            triggerGroupCharges.append((np.array(charge_files[f'ch{ch}'][f'{type}_charge'])-pedestal_dict[ch])*calib_dict[ch])
        triggerGroupArrs.append(np.stack(triggerGroupCharges, axis=1))
    return triggerGroupArrs

def getCrosstalkMatrix(hdf5_file_name, type='lyso'):
    '''
    return crosstalk matrix associated with the data from a given hdf5 file - array of 4 will be returned, one for each trigger group
    '''
    if type=='lyso':
        pedestals = get_LYSO_pedestals(hdf5_file_name)
    calib_dict = get_total_calibrations()
    data_arrays = getChargeArrays_PedestalAdded(hdf5_file_name, pedestals, calib_dict, type)
    crosstalk_arrays = []
    for triggerGroup, charges in enumerate(data_arrays):
        triggerChannel = np.argmax(charges, axis=1)
        arr_for_stack = []
        for i in range(8):
            triggeredEventCharges = charges[triggerChannel==i]
            triggeredCharge = triggeredEventCharges[:,i]
            triggeredCharge = np.expand_dims(triggeredCharge, axis=1)
            crosstalk_denom = np.broadcast_to(triggeredCharge, (len(triggeredCharge),8))
            crosstalk_vals = triggeredEventCharges/crosstalk_denom
            crosstalk_array = np.average(crosstalk_vals, axis=0)
            arr_for_stack.append(crosstalk_array)
        stacked = np.stack(arr_for_stack, axis=1)
        crosstalk_arrays.append(stacked)

    return crosstalk_arrays


def getCrosstalkMatrcies_allFiles(hdf5_list, type='lyso'):
    '''
    return crosstalk matrices with matrix elements averaged over the trials, four arrays will be returned, one for each trigger group
    '''
    full_matrices = []
    for hdf5file in hdf5_list:
        full_matrices.append(np.array(getCrosstalkMatrix(hdf5file, type)))
    full_matrices = np.array(full_matrices)
    return np.average(full_matrices, axis=0)