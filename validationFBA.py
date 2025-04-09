import subprocess
import os
import shutil
import configparser
import datetime
import json
from pathlib import Path
import cobra
from collections import defaultdict
from six import iteritems, itervalues 
import pandas as pd
import sklearn.metrics  as sk
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from cobra import Model, Reaction, Metabolite
from itertools import compress

## Options

analyzeResults = True
plotHeatmaps = True #True when you have a selected model

plotAccTables = True  #True when you have a selected model
plotH37Ra = False
model_type = 'min_gf'
uptake = 0.1
cutoff = 0.001
## Specify working directory

dir = os.getcwd()

## Specify date
date = '280125'

## import data
#if uptake == 0.1:
#    media = pd.read_csv('/home/icancino/Documents/PhD/Reconstruction MTB analysis/validation_media_compounds_DR_IC_0.1.csv')

#else:   
#    media = pd.read_csv('/home/icancino/Documents/PhD/Reconstruction MTB analysis/validation_media_compounds_DR_IC.csv')

exp_results = pd.read_csv('/home/icancino/Documents/PhD/Reconstruction MTB analysis/phenotypes_binary_luiz.csv')
#exp_results.to_latex('phenotypes_binary_luiz.txt')
#exp_results = pd.read_csv('/home/icancino/Documents/PhD/Reconstruction MTB analysis/phenotypes_binary_acely.csv')
reference_model = 'iEK1011_2.0.json'
iEK= cobra.io.load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/'+ reference_model)
#exp_data_H37Rv = pd.read_csv('/home/icancino/Documents/PhD/Reconstruction MTB analysis/carbon_nitro_H37rv_metrics.csv')

## Create dict to names
species_dir = ['abscessus','marinum','smegmatis','aromaticivorans','H37Rv_new']


if model_type == 'min_gf':
    selected_models={'abscessus':'abscessusmore-sensitive651e-50',
                    'marinum':'marinummore-sensitive551e-50',
                    'smegmatis':'gf_smegmatissensitive651',
                    'aromaticivorans':'aromaticivoransfast651e-10',#aromaticivoransfast651e-10

                    'H37Rv':'new_H37Rvfast551e-10'} #minimal models

elif model_type == 'max_gf':
    selected_models={'abscessus':'abscessusfast301',
                'marinum':'marinummore-sensitive301e-50',
                'smegmatis':'gf_smegmatissensitive401e-10',
                'aromaticivorans': 'aromaticivoranssensitive301',#old 'aromaticivorans':'aromaticivoransmore-sensitive551e-10',
                'H37Rv':'new_H37Rvmore-sensitive551e-10'} #maximal models
     
elif model_type == 'max_genes_perc_gf':
    selected_models = {'abscessus':'abscessusfast301e-10',
            'marinum':'marinumfast301',
            'smegmatis':'gf_smegmatissensitive401e-10',
            'aromaticivorans':'aromaticivoranssensitive301',#'aromaticivoranssensitive301',#'aromaticivoransmore-sensitive301', # old 'aromaticivorans':'aromaticivoranssensitive551e-10', aromaticivoransmore-sensitive301
            'H37Rv':'new_H37Rvfast501'}
    

species_exp =['abscessus','marinum','smegmatis','aromaticivorans','H37Rv']

## Mass balance function
def is_mass_balanced(reaction):
    balance = defaultdict(int)
    for metabolite, coefficient in iteritems(reaction.metabolites):
        if metabolite.elements is None or len(metabolite.elements) == 0:
            return False
        for element, amount in iteritems(metabolite.elements):
            balance[element] += coefficient * amount
    return all(amount == 0 for amount in itervalues(balance))

media_bounds_dict= {}

#for media_no in media.keys(): 
#    if media_no.endswith('_bounds') == False:
#        media_bounds_dict[media_no] ={}
#        media_bounds_dict[media_no] = dict(zip(media[media_no],media[str(media_no)+'_bounds']))
#print(media_bounds_dict)


media_bounds_dict = {'7H10': {#'EX_co2_e':-1000, 
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glu__L_e': -uptake, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_cit_e': -uptake, 'EX_btn_e': -1000.0, 'EX_ pydxn_e': -1000.0},
            '3G': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glu__L_e': -uptake, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake}, 
            'R-C': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_arg__L_e': -uptake}, 
            'K-C': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_lys__L_e': -uptake}, 
            'Py-C': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_pyr_e': -uptake},
            'Ma-C': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_mal__L_e': -uptake},
            'V-C': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_ac_e': -uptake}, 
            'Rh-C': {#'EX_co2_e':-1000,
            'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_rmn_e': -uptake},
                '0-C': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_co2_e': -uptake}, 
                'HCO3-C': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_hco3_e': -uptake}, 
                'L-C': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_nh4_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_leu__L_e': -uptake}, 
                'K-N': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_lys__L_e': -uptake},
                'R-N': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_arg__L_e': -uptake},
                'W-R': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_trp__L_e': -uptake}, 
                'L-N': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_leu__L_e': -uptake}, 
                'A-N': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_ala__L_e': -uptake}, 
                'G-N': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_gly_e': -uptake},
                '0-N': {#'EX_co2_e':-1000,
                'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_n2_e': -uptake},
                    'NH4-N': {#'EX_co2_e':-1000,
                    'EX_k_e': -1000.0, 'EX_h_e': -1000.0, 'EX_pi_e': -1000.0, 'EX_o2_e': -1000.0, 'EX_na1_e': -1000.0, 'EX_so4_e': -1000.0, 'EX_cl_e': -1000.0, 'EX_zn2_e': -1000.0, 'EX_fe3_e': -1000.0, 'EX_fe2_e': -1000.0, 'EX_cu2_e': -1000.0, 'EX_mn2_e': -1000.0, 'EX_mobd_e': -1000.0, 'EX_ca2_e': -1000.0, 'EX_mg2_e': -1000.0, 'EX_h2o_e': -1000.0, 'EX_h2_e': -1000.0, 'EX_cobalt2_e': -1000.0, 'EX_glc__D_e': -uptake, 'EX_glyc_e': -uptake, 'EX_nh4_e': -uptake}}


## Perform FBA on models inside the different species_dir folders 
if  analyzeResults: 
    for species in species_dir:
        
        globals()[f'{species}_results'] = {()}

        globals()[f'{species}_FBA'] =pd.DataFrame({})
        globals()[f'{species}_results'] = pd.DataFrame({}) 
        globals()[f'{species}_flux_reactions'] = pd.DataFrame({})
        globals()[f'{species}_missing_rxn'] = pd.DataFrame({})

        #for media_no in media_bounds_dict.keys():
        #    globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes'] = pd.DataFrame({'Index': range(2000)})

        if species != 'H37Rv_new':
            for f in os.listdir(dir+'/'+species):
                filename = os.fsdecode(dir+'/'+species+'/'+f)
                print(f[:-4])
                if filename.endswith('.xml'):# and f[:-4] == selected_models[species]:
                    
                    model = cobra.io.read_sbml_model(filename)
                    model.solver = 'gurobi'
                    model.objective = 'Growth'
                    reaction_keys = model.reactions.__dict__['_dict'].keys()        
                    globals()[f"{filename[len(dir+'/'+species+'/'):-4]}_missing_rxn"] = set()

                    for media_no in media_bounds_dict.keys():
                        
                        model.reactions.BIOMASS_2.bounds = (0,0)    
                        model.reactions.ATPM.bounds = (3.15,1000)        
    
                        for exchange in model.exchanges:
                            model.reactions.get_by_id(exchange.id).lower_bound = 0

                        for rxn in model.reactions:
                            if rxn.id.__contains__('DM_'):
                                model.reactions.get_by_id(rxn.id).bounds = (0,0)
                        for sink in model.sinks:
                            model.reactions.get_by_id(sink.id).bounds = (0,1000)
            
                        for reaction in media_bounds_dict[media_no].keys():
                            if reaction in reaction_keys:
                                        #print('bounds',reaction,media_bounds[reaction])
                                model.reactions.get_by_id(reaction).lower_bound = media_bounds_dict[media_no][str(reaction)]
                                print((media_bounds_dict[media_no][str(reaction)]))
                            elif not pd.isna(reaction): ## add missing reaction in another file
                                print('Reaction ',reaction,' in media but not in the model')
                                globals()[f"{filename[len(dir+'/'+species+'/'):-4]}_missing_rxn"].add(reaction)
                                
                                continue 
                        globals()[f'{species}_missing_rxn'].loc[str(filename[len(dir+'/'+species+'/'):-4]),'Missing rxns'] =  str(globals()[f"{filename[len(dir+'/'+species+'/'):-4]}_missing_rxn"]) 
                        fluxes = model.optimize().fluxes
                    
                        #active_fluxes = fluxes[fluxes > 0][fluxes < 0]
                        #reactions_w_genes =  0
                        #for reaction in active_fluxes.keys():
                        #    if model.reactions.get_by_id(reaction).genes:
                        #        reactions_w_genes +=1
                        #for reaction in model.exchanges:
                        #    if reaction.id in active_fluxes.keys():
                        #        flux_df = pd.DataFrame({str(reaction.id): [active_fluxes[reaction.id]]})
                        #    else:
                        #        continue
                        #globals()[f'{species}_flux_reactions'].loc['Reactions with genes and active flux',str(filename[len(dir+'/'+species+'/'):-4])]  = reactions_w_genes/len(model.reactions)
                        # if not active_neg_fluxes.empty:
                        #     globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes']['ReactionsID_'+str(filename[len(dir+'/'+species+'/'):-4])]  = active_neg_fluxes.keys()
                        #     globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes'][str(filename[len(dir+'/'+species+'/'):-4])]  = np.array(active_neg_fluxes)
                        #     #globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes']= pd.concat([globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes'],pd.DataFrame({'ReactionsID_'+str(filename[len(dir+'/'+species+'/'):-4]): active_fluxes.keys(),str(filename[len(dir+'/'+species+'/'):-4]):np.array(active_fluxes)})])
                        # else:
                        #     globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes']['ReactionsID_'+str(filename[len(dir+'/'+species+'/'):-4])]  = 'All reactions are inactive'
                        #     globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes'][str(filename[len(dir+'/'+species+'/'):-4])]  = 0
                        #     #globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes']= pd.concat([globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes'],pd.DataFrame({'ReactionsID_'+str(filename[len(dir+'/'+species+'/'):-4]): "All inactive",str(filename[len(dir+'/'+species+'/'):-4]):0})])

                        if model.optimize().objective_value <= cutoff:
                            print(media_no,model.optimize().objective_value)
                            globals()[f'{species}_FBA'].loc[str(filename[len(dir+'/'+species+'/'):-4]),str(media_no)] =  0 #model.optimize().objective_value
    
                        elif model.optimize().objective_value > cutoff:
                            print(media_no,model.optimize().objective_value)
                            globals()[f'{species}_FBA'].loc[str(filename[len(dir+'/'+species+'/'):-4]),str(media_no)] = 1
                        else:
                            print('Error',media_no)
                    
                    genes_to_rxn = 0
                    iek_reactions = 0 
                    unbal_rxns = 0
                    for rxn in model.reactions:
                        if rxn.genes:
                            genes_to_rxn +=1
                        if rxn in iEK.reactions:
                            iek_reactions +=1
                        if is_mass_balanced(rxn) == False and str(rxn.id)[:2] != 'EX' and str(rxn.id)[:4] != 'sink' and str(rxn.id)[:2] != 'DM':
                            unbal_rxns +=1

                    globals()[f'{species}_results'].loc[filename[len(dir+'/'+species+'/'):-4],'Total Reactions'] = len(model.reactions) 
                    globals()[f'{species}_results'].loc[filename[len(dir+'/'+species+'/'):-4],'Total Metabolites'] = len(model.metabolites) 
                    globals()[f'{species}_results'].loc[filename[len(dir+'/'+species+'/'):-4],'Total Genes'] = len(model.genes)
                    globals()[f'{species}_results'].loc[filename[len(dir+'/'+species+'/'):-4],'Rxn in iEK1011_2 (%)'] = (iek_reactions*100/len(model.reactions))
                    globals()[f'{species}_results'].loc[filename[len(dir+'/'+species+'/'):-4],'Rxn associated to genes (%)'] = (genes_to_rxn*100/len(model.reactions))
                    globals()[f'{species}_results'].loc[filename[len(dir+'/'+species+'/'):-4],'Unbalanced reactions'] = unbal_rxns
            

      
        ## Compare results to experimental values in exp_results
            ACC = pd.DataFrame({})
            MCC = pd.DataFrame({})
            print(globals()[f'{species}_FBA'])

            for idx, res_key in enumerate(globals()[f'{species}_FBA'].index):
                print(idx, res_key)

                out = exp_results[species][[not np.isnan(x).any() for x in exp_results[species]]]
            
                cropped_FBA = globals()[f'{species}_FBA'].loc[res_key][[not np.isnan(x).any() for x in exp_results[species]]]
                ACC.loc[res_key,'ACC'] = sk.accuracy_score(out,cropped_FBA.dropna())
                print(ACC)
                try:
                    MCC.loc[res_key,'MCC'] = sk.matthews_corrcoef(out,cropped_FBA.dropna())

                except ZeroDivisionError:
                    MCC[res_key] = -1
                globals()[f'{species}_results'].loc[res_key,'Total Accuracy'] = ACC.loc[res_key,'ACC']
                globals()[f'{species}_results'].loc[res_key,'Total MCC'] =  MCC.loc[res_key,'MCC']
            

                globals()[f'{species}_FBA'].to_csv(species+f'_FBA.csv')
                globals()[f'{species}_results'].to_csv(species+f'_results.csv')
                globals()[f'{species}_missing_rxn'].to_csv(species+f'_missing_rxn.csv')

                print(globals()[f'{species}_results'])
        
if plotHeatmaps:
    #Plotting Heatmap
    name_list = ['7H10', '3G', 'L-Arginine (C)','L-Lysine (C)','Pyruvate (C)', 'Malate (C)', 'Acetate (C)', 'L-Rhamnose (C)','No carbon', 'Bicarbonate (C)', 'L-Leucine (C)', 'L-Lysine (N)','L-Arginine (N)','L-Tryptophane (N)', 'L-Leucine (N)','L-Alanine (N)', 'L-Glycine (N)', 'No nitrogen', 'Ammonium (N)']

    name_source = {id_source: name for id_source, name in zip(exp_results['Unnamed: 0'], name_list)}

    species_name  = {'abscessus':'M. abscessus', 'smegmatis':'M. smegmatis' , 'marinum':'M. marinum', 'aromaticivorans':'M. aromaticivorans'}

    cmap = {
        1: 'black',
        0: 'grey',
        np.nan: 'yellow'
    }
    
    
    def color_map(val):
        if pd.isna(val):
            return cmap[np.nan]
        return cmap[val]
    
    exp_results_bin = pd.concat([exp_results['abscessus'], exp_results['smegmatis'], exp_results['marinum'], exp_results['aromaticivorans']], axis =1)
    color_matrix = exp_results_bin.applymap(color_map)
    
    model_results_bin = pd.concat([abscessus_FBA.loc[selected_models['abscessus']], smegmatis_FBA.loc[selected_models['smegmatis']], marinum_FBA.loc[selected_models['marinum']], aromaticivorans_FBA.loc[selected_models['aromaticivorans']]], axis =1) #add gf
    color_matrix_models = model_results_bin.applymap(color_map)

    model_results_bin.reset_index(drop=True, inplace=True)
    model_results_new_name = model_results_bin.rename(columns={selected_models['abscessus']: "abscessus",selected_models['smegmatis']:"smegmatis", selected_models['marinum']: "marinum", selected_models['aromaticivorans'] #add gf
    :"aromaticivorans"})  
    comparison = model_results_new_name == exp_results_bin
    ## Change colors
    custom_colors = ['red', 'green']  # Fuchsia for mismatch, LimeGreen for match
    custom_cmap = ListedColormap(custom_colors)
    sec_colors = ['grey','black']
    custom_sec_col= ListedColormap(sec_colors)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(10,5))
    g1 =sns.heatmap(exp_results_bin,annot=False, fmt="", cmap=[cmap[0], cmap[np.nan], cmap[1]], ax= ax1, cbar=False,linewidths=1, linecolor='white')
    g2 =sns.heatmap(model_results_bin,annot=False, fmt="", cmap=[cmap[0], cmap[np.nan], cmap[1]],ax= ax2, cbar=False,linewidths=1, linecolor='white',yticklabels=False)
    g3 = sns.heatmap(comparison, cmap=custom_colors, ax= ax3,cbar=False,linewidths=1, linecolor='white',yticklabels=False)

    #ax3.title('Experiments vs Models', fontsize=16)
    #plt.xlabel('Species', fontsize=14)
    ax1.set_xticks(np.arange(len(list(species_name.keys()))) + 0.5, list(species_name.values()) ,style = 'italic', rotation = 45)
    ax2.set_xticks(np.arange(len(list(species_name.keys()))) + 0.5, list(species_name.values()) ,style = 'italic', rotation = 45)
    ax3.set_xticks(np.arange(len(list(species_name.keys()))) + 0.5, list(species_name.values()) ,style = 'italic', rotation = 45)

    #ax3.set_yticks(np.arange(len(list(model_results_bin.iloc[:,0]))) + 0.5, list(name_source.values()), rotation=0)

    ax1.set_ylabel('Media', fontsize=14)
    #ax1.set_xlabel('Experiments', rotation = 45)
    #ax2.set_xlabel('Model prediction', rotation = 45)
    #ax3.set_xlabel('Experiments vs Model', rotation = 45)
    ax2.tick_params(left=False, bottom=True) ## other options are right and top
    ax3.tick_params(left=False, bottom=True) ## other options are right and top

    #ax1.set_xticks(5.5, 'Experiments' ,style = 'italic', rotation =45)
    #ax1.set_yticks(np.arange(len(media_name))+0.5, media_name, rotation=0)
    ax1.set_yticks(np.arange(len(exp_results['Unnamed: 0'])) + 0.5, list(name_source.values()), rotation=0)

    

    plt.savefig(f'all_four_models_exp_heatmap_upt_{uptake}_co_{cutoff}_{date}_{model_type}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    #g1.title('Experimental Results', fontsize=16)
    #plt.xlabel('Species', fontsize=14)
    #g1.ylabel('Media', fontsize=14)
    #g1.xticks(np.arange(len(list(species_name.keys()))) + 0.5, list(species_name.values()) ,style = 'italic',rotation = 45)
    #g1.yticks(np.arange(len(exp_results['Unnamed: 0'])) + 0.5, list(name_source.values()), rotation=0)
    
    #fig, ax = plt.subplots(figsize=(5, 7))

    #for text, color in zip(ax2.texts, color_matrix_models.values.flatten()):
    #    text.set_color(color)
    #    text.set_fontsize(14)

    #g2.title('Model predictions', fontsize=16)
    #plt.xlabel('Species', fontsize=14)
    #g2.ylabel('Media', fontsize=14)
    #g2.xticks(np.arange(len(list(species_name.keys()))) + 0.5, list(species_name.values()) ,style = 'italic', rotation = 45)
    #g2.yticks(np.arange(len(list(model_results_bin.iloc[:,0]))) + 0.5, list(name_source.values()), rotation=0)
    #plt.savefig(f'gf_model_results_luiz_heatmap_{date}.png', dpi=300, bbox_inches='tight')
    #plt.show()
    #plt.close()


    #plt.figure(figsize=(5, 7))

    # Display the plot
    #plt.close()




if plotAccTables:
    species_name  = {'abscessus':'M. abscessus', 'smegmatis':'M. smegmatis' , 'marinum':'M. marinum', 'aromaticivorans':'M. aromaticivorans'}#, 'H37Ra':'M. tuberculosis H37Ra'}

    #Ra_cm = sk.confusion_matrix(exp_results['H37Ra'], H37Ra_FBA.loc[selected_models['H37Ra']])
    abs_cm = sk.confusion_matrix(exp_results['abscessus'], abscessus_FBA.loc[selected_models['abscessus']]) #add gf_
    smeg_cm = sk.confusion_matrix(exp_results['smegmatis'], smegmatis_FBA.loc[selected_models['smegmatis']]) #add gf_
    mari_cm = sk.confusion_matrix(exp_results['marinum'], marinum_FBA.loc[selected_models['marinum']])
    aroma_cm = sk.confusion_matrix(exp_results['aromaticivorans'].dropna(), list(compress(aromaticivorans_FBA.loc[selected_models['aromaticivorans']],exp_results['aromaticivorans'].notna())))
    cm_list = [abs_cm, smeg_cm, mari_cm, aroma_cm]#, Ra_cm]
    acc_list = ['%.2f'%(abscessus_results.loc[selected_models['abscessus'],'Total Accuracy']*100),'%.2f'%(smegmatis_results.loc[selected_models['smegmatis'],'Total Accuracy']*100),'%.2f'%(marinum_results.loc[selected_models['marinum'],'Total Accuracy']*100),'%.2f'%(aromaticivorans_results.loc[selected_models['aromaticivorans'],'Total Accuracy']*100)]#,'%.2f'%(H37Ra_results.loc[selected_models['H37Ra'],'Total Accuracy']*100)] #add gf_
    for idx, cm in enumerate(cm_list):
        plt.figure(figsize = (6,6))
        sns.heatmap(cm, 
                    annot=True,
                    fmt='g', 
                    xticklabels=['No Growth','Growth'],
                    yticklabels=['No Growth','Growth'],
                    cmap = sns.diverging_palette(20, 145, s=100, as_cmap=True))
        plt.ylabel('Experiments', fontsize=13)
        plt.title('Confusion Matrix for $\mathit{'+f'{list(species_name.values())[idx]}'+'}$', fontsize=16, pad=20)
        plt.gca().xaxis.set_label_position('top') 
        plt.xlabel('Model predictions', fontsize=13)
        plt.gca().xaxis.tick_top()

        plt.gca().figure.subplots_adjust(bottom=0.2)
        plt.gca().figure.text(0.43, 0.13, 'Accuracy of '+ acc_list[idx]+'%', ha='center', fontsize=13)
        plt.savefig(f'{list(species_name.keys())[idx]}_{date}_{model_type}.png', dpi=300, bbox_inches='tight')

        plt.show()
        #plt.close()
         
if plotH37Ra:
    name_list = ['7H10', 'm7H9', 'Arginine (Carbon)','Lysine (Carbon)','Pyruvate (Carbon)', 'Malate (Carbon)', 'Acetate (Carbon)', 'Rhamnose (Carbon)','No carbon', 'Bicarbonate (Carbon)', 'Leucine (Carbon)', 'Lysine (Nitrogen)','Arginine (Nitrogen)','Tryptophane (Nitrogen)', 'Leucine (Nitrogen)','Alanine (Nitrogen)', 'Alanine (Nitrogen)', 'Glycine (Nitrogen)', 'No nitrogen', 'Ammonium (Nitrogen)']

    name_source = {id_source: name for id_source, name in zip(exp_results['Unnamed: 0'], name_list)}

    cmap = {
        1: 'green',
        0: 'red',
        np.nan: 'yellow'
    }
    label =   ['$M. tuberculosis H37Ra$']
    def color_map(val):
        if pd.isna(val):
            return cmap[np.nan]
        return cmap[val]

    exp_results_bin = pd.DataFrame(exp_results['H37Ra'])
    color_matrix = exp_results_bin.applymap(color_map)

    fig, ax = plt.subplots(figsize=(2, 7))
    sns.heatmap(exp_results_bin,annot=False, fmt="", cmap=[cmap[0], cmap[np.nan], cmap[1]], ax=ax, cbar=False,linewidths=1, linecolor='white')

    for text, color in zip(ax.texts, color_matrix.values.flatten()):
        text.set_color(color)
        text.set_fontsize(14)

    #plt.title('Experimental Results', fontsize=16)
    #plt.xlabel('Species', fontsize=14)
    plt.ylabel('Media', fontsize=14)
    plt.xticks(np.arange(len(label)) + 0.5, label ,style = 'italic',rotation = 0)
    plt.yticks(np.arange(len(exp_results['Unnamed: 0'])) + 0.5, list(name_source.values()), rotation=0)
    plt.savefig(f'exp_results_H37Ra_luiz_heatmap_{date}_{model_type}.png', dpi=300, bbox_inches='tight')
    plt.show()
    #plt.close()
    print('map')

#for species in species_dir:
#    globals()[f'{species}_flux_reactions'].to_csv(species+'_flux_reactions.csv')
    #for media_no in media_bounds_dict.keys():
    #    globals()[f'{species}_{str(media_no).replace('-','_')}_fluxes'].to_csv(species+media_no+'_fluxes.csv')                  
        
#abscessus_FBA.to_csv('abscessus_FBA.csv') 
#aromaticivorans_FBA.to_csv('aromaticivorans_FBA.csv') 

#marinum_FBA.to_csv('marinum_FBA.csv') 
#smegmatis_FBA.to_csv('smegmatis_FBA.csv') 
#H37Ra_FBA.to_csv('H37Ra.csv') 
#H37Rv_FBA.to_csv('H37Rv.csv') 

#abscessus_results.to_csv('results_abscessus.csv')
#aromaticivorans_results.to_csv('results_aromaticivorans.csv')
#marinum_results.to_csv('results_marinum.csv')
#smegmatis_results.to_csv('results_smegmatis.csv')
#H37Rv_results.to_csv('results_H37Rv.csv')
#H37Ra_results.to_csv('results_H37Ra.csv')

