#!/usr/bin/env python3

import sys
sys.path.append(".")

import h5py as h5
import glob
import numpy as np

basePath = "/storage/af/group/btl-production/upload/QAQC_SM/qaqc-gui_output/SM_QAQC_Production/"

triggerGroup1_channelRange = (0,8)
triggerGroup2_channelRange = (8,16)
triggerGroup3_channelRange = (16,24)
triggerGroup4_channelRange = (24,32)

atten = 5.85 #waveform attentuation from PiPad

def get_hdf5_file(runNum, moduleNum = 32110020008497):
    '''
    return list of TFiles from desired modules, with minimum run number and maximum run number specified
    '''
    file = glob.glob(basePath+f"run0{runNum}/module_"+str(moduleNum)+".hdf5")
    return str(file[0])

def get_trigger_mask(hdf5File, ch, type='lyso'):
    '''
    return mask with events for which set channel was the trigger channel and max entry index considered
    '''
    integrals_file = hdf5File.replace(".hdf5", "_integrals.hdf5")
    file = h5.File(integrals_file, 'r')

    if ch<8: ch_range = triggerGroup1_channelRange
    elif ch<16: ch_range = triggerGroup2_channelRange
    elif ch<24: ch_range = triggerGroup3_channelRange
    else: ch_range = triggerGroup4_channelRange

    triggerGroupWaveforms = []
    numEvents = []
    for groupCh in range(ch_range[0], ch_range[1]):
        numEvents.append(np.array(file[f'ch{groupCh}'][f"{type}_charge"]).shape[0])
    maxEntry = min(numEvents)
    for groupCh in range(ch_range[0], ch_range[1]):
        triggerGroupWaveforms.append(file[f'ch{groupCh}'][f"{type}_charge"][0:maxEntry])
    triggerGroupWaveforms = np.stack(triggerGroupWaveforms, axis=1)
    triggerChannel = np.argmax(triggerGroupWaveforms, axis=1)
    # if type=='spe': #rough approx to make sure at least one channel in event has an spe
    #     print(np.max(triggerGroupWaveforms, axis=1))
    #     spe_event = np.max(triggerGroupWaveforms, axis=1)>2.5 #rough lower bound to be confident that we have an spe event
    #     #return np.logical_and(triggerChannel==ch, spe_event), maxEntry
    #     return triggerChannel==ch, maxEntry
    return triggerChannel==ch, maxEntry

def convert_values(hdf5File, averaged_waveform, type):
    '''
    return the averaged waveform with time on the x axis and voltage on the y axis
    input is the averaged waveform in digital units. The number of points in the array is the number of times the digitizer sampled
    this code is borrowed from Tony's code for integrating the waveforms
    '''
    xinc = 1/(hdf5File[type].attrs['drs4_frequency'] * 10**6)
    points = hdf5File[type].attrs['record_length']
    x = np.linspace(0, xinc * points, int(points))
    y = averaged_waveform/2**12
    return x*1e9, y


def average_waveforms(hdf5File, ch=0, type='lyso'):
    '''
    return sampled time stamps and event-averaged voltage values. The waveforms that were averaged must have the largest integrated charge for that event within their trigger group
    the digitizer sampled values are converted to voltages and time
    '''
    file = h5.File(hdf5File, 'r')
    waveform_data = np.array(file[type][f'ch{ch}'])
    triggerMask, maxEntry = get_trigger_mask(hdf5File, ch, type)
    waveforms_masked = waveform_data[0:maxEntry][triggerMask]
    averaged_waveform = np.average(waveforms_masked, axis=0)
    sampled_points, voltage_vals = convert_values(file, averaged_waveform, type)
    #if type=='lyso':
    #    voltage_vals=voltage_vals*atten
    return sampled_points, voltage_vals
