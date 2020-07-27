import argparse
import sys
import keras
from keras.utils.np_utils import to_categorical
from keras.models import Sequential, load_model
from keras.layers import Dense, Activation, Dropout
from keras import backend as K
from keras.optimizers import SGD
from keras import optimizers
from keras.callbacks import EarlyStopping, ModelCheckpoint
#import sklearn
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
from root_numpy import root2array, tree2array, array2root
from common_function import dataset, AMS, read_data, prepare_data, drawfigure, calc_sig, calc_sig_new, f1, f1_loss
import config_OPT_NN as conf
import ROOT 
from pathlib import Path
import os

#===============================================================================
#-------------------------------------------------------------------------------
#
# Shows the number of events of EW and QCD background as well as the individual
# mass' number of events, within mass window and cut-based. 
# Calculates the significance before the training.
#
#-------------------------------------------------------------------------------
#===============================================================================

parser = argparse.ArgumentParser()
parser.add_argument("--m",  "--model",     help="Model to use",                     default="GM", type=str)
parser.add_argument("--sd", "--subdir",    help="Subdirectory in OutputRoot",       default="",   type=str)
parser.add_argument("--rID", "--resultID", help="Identifier for the results' file", default="0",  type=str)
parser.add_argument("--w",  "--weight",    help="WeightNormalized (0), abs(WN) (1), or else (2)", default=0, type=int)
args = parser.parse_args()

#-------------------------------------------------------------------------------

def getSamples(mass):
    """
    Similar to config_OPT_NN.input_samples, except this class uses the output
    of a training for an individual mass
    """
    class input_samples_NN:
        # Prefixe added to the ntuples of the training output
        prestr = f"new_{args.m}_main"

        # Directory where ntuples are located
        filedir = f"OutputRoot/{args.sd}/"+prestr

        # Background samples
        bckgr = {
        'name'    : ['MVA.361292_MGaMcAtNloPy8EG_NNPDF30LO_A14NNPDF23LO_WZ_lvll_FxFx_ntuples.root',
                     'MVA.364284_Sherpa_222_NNPDF30NNLO_lllvjj_EW6_ntuples.root'],
        'xs'      : [1704., 47.],
        'nevents' : [3890000, 7325000]}

        # Signal samples
        sigGM = {
        'name'   : ['MVA.450765_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m200_ntuples.root',
                    'MVA.450766_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m250_ntuples.root',
                    'MVA.450767_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m300_ntuples.root', 
                    'MVA.450768_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m350_ntuples.root', 
                    'MVA.450769_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m400_ntuples.root',
                    'MVA.450770_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m450_ntuples.root',
                    'MVA.450771_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m500_ntuples.root',
                    'MVA.450772_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m600_ntuples.root',
                    'MVA.450773_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m700_ntuples.root',
                    'MVA.450774_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m800_ntuples.root',
                    'MVA.305035_MGPy8_A14NNPDF30NLO_VBS_H5p_lvll_900_qcd0_ntuples.root'],
        'xs'      : [7.0596,  7.710,  3.9238,   4.582, 2.4428,    2.95, 1.6113, 1.1005, 0.77398, 0.55433, 0.40394 ],
        'nevents' : [70000,   70000,   70000,  190000,  70000,   70000,  70000,  70000,   70000,   70000,   50000 ],
        'filtEff' : [1,     0.77156,       1, 0.77507,      1, 0.77891,      1,      1,       1,       1,       1 ]}
        sigHVT = {
        'name'   : ['MVA.307730_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0250_ntuples.root',
                    'MVA.307731_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0300_ntuples.root',
                    'MVA.309528_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0350_ntuples.root',
                    'MVA.307732_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0400_ntuples.root',
                    'MVA.309529_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0450_ntuples.root',
                    'MVA.307733_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0500_ntuples.root',
                    'MVA.307734_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0600_ntuples.root',
                    'MVA.307735_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0700_ntuples.root',
                    'MVA.307736_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0800_ntuples.root',
                    'MVA.307737_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0900_ntuples.root',
                    'MVA.307738_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1000_ntuples.root'],
        'xs'      : [24.42 , 10.54 , 2.299 , 0.7975, 0.3408, 0.1663, 0.087984, 0.049882, 0.02961, 1.3050, 4.4843],
        'nevents' : [190000, 190000, 190000, 190000, 190000, 190000,   190000,   175000,  190000, 160000, 160000],
        'filtEff' : [     1,      1,      1,      1,      1,      1,        1,        1,       1,      1,      1]}

    return input_samples_NN

def getCV():
    cv = 0
    for model in os.listdir(f"OutputModel/{args.sd}/"):
        if model.endswith("4_NN.h5"):
            print(model)            
            temp = float(model[int(model.find("CV")+3):int(model.find("CV")+8)])
            cv += temp
    return cv/4


def AMS(s, b):
    """ Approximate Median Significance """

    br = 0.00001
    sigma = math.sqrt(b+br)
    n = s+b+br
    radicand = 2 *( n * math.log (n*(b+br+sigma)/(b**2+n*sigma+br))-b**2/sigma*math.log(1+sigma*(n-b)/(b*(b+br+sigma))))
    significance = math.sqrt(radicand)
    
    return significance

def bkgEvents(cins, mass, ucut=0, uwei=0):
    """
    Returns the cut-based, normalize-weighted number of background events within
    the mass window - total, QCD and EW

    cins - config.input_samples or similar object
    mass - 
    ucut - user selection to use instead of cut-based selection
    uwei - method of normalizing weight
    """

    bkg = np.array([0.,0.])    # Number of events for QCD and EW

    for i in range(len(cins.bckgr['name'])):
        # Accessing the data
        bkg_file = ROOT.TFile(cins.filedir+cins.bckgr['name'][i])
        tree = bkg_file.Get('nominal')

        # Choosing the cut
        if ucut: cut = ucut
        else: cut = "M_jj>500 && Deta_jj>3.5"

        # Selecting the data
        cuts = 'Jet1Pt>0 && Jet2Pt>0 && M_WZ>({0}*0.88) && M_WZ<({0}*1.12) && {1}'.format(mass, cut)
        bkg_DF = pd.DataFrame(tree2array(tree, selection=cuts))

        # Calculating the number of events
        if   uwei==0: bkg[i] = sum(bkg_DF['WeightNormalized'])
        elif uwei==1: bkg[i] = sum(abs(bkg_DF['WeightNormalized']))
        elif uwei==2: bkg[i] = sum(bkg_DF['WeightNormalized'])

    return sum(bkg), bkg[0], bkg[1]

def sigEvents(cins, mass, model, ucut=0, uwei=0):
    """
    Returns the cut-based, normalize-weighted number of events of the signal
    for an individual mass within the mass window

    cins - config.input_samples or similar object
    mass - mass at which the numbe rof events is evaluated
    ucut - user selection to use instead of cut-based selection
    uwei - user weight to use instead of WeightNormalized
    """
    # Accessing the signal for the appropriate model
    if model == "GM":
        sigMO = cins.sigGM
    elif model == "HVT":
        sigMO = cins.sigHVT

    isf = int(np.where([f"{mass}" in sigMO['name'][i] for i in range(len(sigMO['name']))])[0])
    sig_file = ROOT.TFile(cins.filedir+sigMO['name'][isf])
    tree = sig_file.Get('nominal')
    
    # Choosing the cut
    if ucut: cut = ucut
    else: cut = "M_jj>500 && Deta_jj>3.5"

    # Selecting the data
    cuts = 'M_jj>100 && M_WZ>({0}*0.88) && M_WZ<({0}*1.12) && {1}'.format(mass, cut)
    sig_DF = pd.DataFrame(tree2array(tree, selection=cuts))

    # Calculating the number of events
    sig_DF = pd.DataFrame(tree2array(tree, selection=cuts))
    if   uwei==0: sig_nEv = sum(sig_DF['WeightNormalized'])
    elif uwei==1: sig_nEv = sum(abs(sig_DF['WeightNormalized']))
    elif uwei==2: sig_nEv = sum((abs(sig_DF['Weight'])<10)*sig_DF['WeightNormalized'])
    return sig_nEv

#-------------------------------------------------------------------------------

# Initializing arrays
res_arr  = np.zeros((14,5))
mass_arr = np.array([250,300,400,500,600,700,800])

# Input samples for cut-based analysis
cins_CB = conf.input_samples

cv = getCV()
print(f"Mean cut value : {cv}")

# Optimal cut values for each mass at tmass 250
ocv_250 = [0.47, 0.77, 0.97, 0.97, 0.97, 0.97, 0.97]
ocv_600 = [0.91, 0.91, 0.91, 0.97, 0.93, 0.97, 0.97]

for c in range(len(mass_arr)):
    mass = mass_arr[c]
    cv = ocv_600[c]

    # Number of events for cut-based analysis
    n_bkg, QCD, EW = bkgEvents(cins_CB, mass, uwei=args.w)
    n_sig = sigEvents(cins_CB, mass, args.m, uwei=args.w)

    # Number of events for NN analysis
    NN_bkg, NN_QCD, NN_EW = bkgEvents(getSamples(mass), mass, f"pSignal>{cv}", args.w)
    NN_sig = sigEvents(getSamples(mass), mass, args.m, f"pSignal>{cv}", args.w)

    # Printing number of events
    print(f"\nMass : {mass} \t Cut-based \t NN")
    print(   "Signal        \t {:.2f}    \t {:.2f}\
            \nQCD           \t {:.2f}    \t {:.2f}\
            \nEW            \t {:.2f}    \t {:.2f}\
            \nCV            \t           \t {:.2f}".format(n_sig, NN_sig, QCD, NN_QCD, EW, NN_EW,cv))
    
    # Calculating and printing significance
    ams = 0.
    NN_ams = 0.
    if np.isfinite(n_bkg):
        if n_bkg>0 and n_sig>0:
            ams = "{:.2f}".format(AMS(n_sig, n_bkg))
    if np.isfinite(NN_bkg):
        if NN_bkg>0 and NN_sig>0:
            NN_ams = "{:.2f}".format(AMS(NN_sig, NN_bkg))
    print(f"Significance \t {ams}     \t {NN_ams}")    

    # Updating results' array
    res_arr[c*2] = None, n_sig, QCD, EW, ams
    res_arr[c*2+1] = cv, NN_sig, NN_QCD, NN_EW, NN_ams

if args.rID != "0":
    # Creating a DataFrame of the results
    mass_rl = [m for m in mass_arr for _ in (0,1)]
    cols_rl = ["Cut-based","NN"]*7
    indx_rl = ["Cut value","Signal","QCD","EW","Significance"]
    pd.options.display.float_format = '{:.2f}'.format
    res_DF  = pd.DataFrame(res_arr.T, columns=pd.MultiIndex.from_tuples(list(zip(mass_rl,cols_rl))),index=indx_rl)

    # Saving the results
    Path(f"Results/{args.sd}").mkdir(parents=True, exist_ok=True)
    res_name = f"Results/{args.sd}/analysis_CBvsNN_{args.m}_cv{cv:.2f}"
    if args.rID != "0": res_name += f"_{args.rID}"
    res_name += ".csv"
    res_DF.to_csv(res_name)
    print(f"Saving results as {res_name}")

