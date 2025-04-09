import subprocess
import os
import shutil
import configparser
import datetime
import json
from pathlib import Path
import cobra
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model, validate_sbml_model
from cobra import Model, Reaction, Metabolite, flux_analysis
from collections import defaultdict
from six import iteritems, itervalues 
import pandas as pd
import sklearn.metrics  as sk
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from cobra import Model, Reaction, Metabolite
from optlang import Model, Variable, Constraint, Objective
from cobra.flux_analysis import flux_variability_analysis


CompareUptakeRates = False # compares different uptake rates
fluxCholSecretion = True
fluxGlycUptake = True
fluxTrypUptake= True
fluxGalactose = True
fluxCO2Secretion = True
CO2seed = True
model_type= 'max_genes_perc_gf' #max_genes_perc_gf, max_gf or min_gf
date = '280225'
source = 'cbn' #'cbn' or 'ntr'
arom = 'arom'
uptake = 0.2 # 5, 1, 0.2 0.1
cutoff= 0.001 

if model_type == 'min_gf':
    selected_models={'abscessus':'abscessusmore-sensitive651e-50',
                    'marinum':'marinummore-sensitive551e-50',
                    'smegmatis':'gf_smegmatissensitive651',
                    'aromaticivorans':'aromaticivoransfast651e-10',  ##aromaticivoranssensitive651e-50
                    'H37Rv':'new_H37Rvfast551e-10'} #minimal models
    old_models = {'smegmatis':'smegmatissensitive651'} #min
    figsize = (10,15)
    figsize_2 = (5,6)


elif model_type == 'max_gf':
    selected_models={'abscessus':'abscessusfast301',
                'marinum':'marinummore-sensitive301e-50',
                'smegmatis':'gf_smegmatissensitive401e-10',
                'aromaticivorans':'aromaticivoransmore-sensitive301', #    old max model 'aromaticivorans':'aromaticivoransmore-sensitive551e-10',
                'H37Rv':'new_H37Rvmore-sensitive551e-10'} #maximal models
    old_models = {'smegmatis':'smegmatissensitive401e-10'} #max
    figsize = (10, 35)
elif model_type == 'max_genes_perc_gf':
    selected_models = {'abscessus':'abscessusfast301e-10',
            'marinum':'marinumfast301',
            'smegmatis':'gf_smegmatissensitive401e-10',
            'aromaticivorans':'aromaticivoranssensitive301',#'aromaticivoransmore-sensitive301', #  # old best model aromaticivoranssensitive551e-10
            'H37Rv':'new_H37Rvfast501'}
    old_models = {'smegmatis':'smegmatissensitive401e-10'}

if source == 'cbn':
    figsize = (10, 40)
    figsize_2 = (10,6)
    figsize_3 = (16,18)
else:
    figsize = (10, 40)
    figsize_2 = (10,5)
    figsize_3 = (16,14)


smegmatis = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/smegmatis/{selected_models["smegmatis"]}.xml')
marinum = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/marinum/{selected_models["marinum"]}.xml')
aromaticivorans = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/aromaticivorans/{selected_models["aromaticivorans"]}.xml')
abscessus = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/abscessus/{selected_models["abscessus"]}.xml')
Rv = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/H37Rv_new/{selected_models["H37Rv"]}.xml')
universe = read_sbml_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/universe/universe_bact_arch_3.xml')
if arom == 'arom':
    model_list  = [abscessus,marinum,smegmatis,aromaticivorans, Rv]
else:
    model_list = [abscessus,marinum,smegmatis, Rv]

if arom == 'arom':
    species_name  = ['abscessus','marinum','smegmatis','aromaticivorans', 'Rv']
else:
    species_name  = ['abscessus','marinum','smegmatis', 'Rv']

def checkCbnNtr(model):
    cbn_source = []
    ntr_source = []
    for x in model.exchanges:
        for i in x.metabolites:
            print(i.formula)
            if str(i.formula).__contains__('C') and str(i.formula).__contains__('H') and str(i.formula).__contains__('Ca') == False and str(i.formula).__contains__('Co') == False and str(i.formula).__contains__('Cl') == False and str(i.formula).__contains__('Cu') == False:
                cbn_source.append(x.id)
            if str(i.formula).__contains__('N') and str(i.formula).__contains__('Na') == False:
                ntr_source.append(x.id)
    return cbn_source, ntr_source
mari_cbn, mari_ntr = checkCbnNtr(marinum)
absc_cbn, absc_ntr = checkCbnNtr(abscessus)
if arom == 'arom':
    arom_cbn, arom_ntr = checkCbnNtr(aromaticivorans)
smeg_cbn, smeg_ntr = checkCbnNtr(smegmatis)
tb_cbn, tb_ntr = checkCbnNtr(Rv)
if arom == 'arom':
    cbn_source = set(mari_cbn + absc_cbn 
                 +arom_cbn 
                 +smeg_cbn+ tb_cbn)
    ntr_source = set(mari_ntr + absc_ntr 
                 +arom_ntr 
                 +smeg_ntr + tb_ntr)
else:
    cbn_source = set(mari_cbn + absc_cbn 
                 +smeg_cbn+ tb_cbn)
    ntr_source = set(mari_ntr + absc_ntr 
                 +smeg_ntr + tb_ntr)

categories = {
    "Amino acids": ["alanine", "glutamine", "glycine", "valine"],
    "Sugars": ["glucose", "fructose", "mannose", "galactose"],
    "Lipids": ["phosphatidyl", "glycerophospho", "ceramide"],
    "Short Chain Fatty Acids": ["acetate", "butyrate", "propionate"],
    "Siderophores": ["siderophore", "enterobactin"],
    "Cofactors": ["ATP", "NADH", "FAD", "CoA"],
}

# classifying metabolites
def categorize_metabolite(metabolite):
    for category, keywords in categories.items():
        if any(keyword in metabolite.lower() for keyword in keywords):
            return category
    return "Other"

# old carbon source and nitrogen, not used
carbon = ["EX_ac_e",
          "EX_ocdcea_e",
          "EX_hdca_e",
          "EX_val__L_e",
          "EX_thr__L_e",
          "EX_ser__L_e",
          "EX_pro__L_e",
          "EX_phe__L_e",
          "EX_orn_e",
          "EX_met__L_e",
          "EX_lys__L_e",
          "EX_leu__L_e",
          "EX_ile__L_e",
          "EX_his__L_e",
          "EX_glu__L_e",
          "EX_gly_e",
          "EX_asp__L_e",
          "EX_arg__L_e",
          "EX_asn__L_e",
          "EX_ala__L_e",
          "EX_chsterol_e",
          "EX_succ_e",
          "EX_pyr_e",
          "EX_ppa_e",
          "EX_mal__L_e",
          "EX_lac__L_e",
          "EX_glyc_e",
          "EX_glyc3p_e",
          "EX_glc__D_e",
          "EX_cit_e",
          "EX_trp__L_e",
          "EX_rmn_e"
          ]

nitro = ["EX_val__L_e",
        "EX_thr__L_e",
        "EX_ser__L_e",
        "EX_pro__L_e", 
        "EX_phe__L_e",
        "EX_orn_e", 
        "EX_met__L_e", 
        "EX_lys__L_e",
        "EX_leu__L_e",
        "EX_ile__L_e",
        "EX_his__L_e",
        "EX_glu__L_e",
        "EX_gly_e",
        "EX_asp__L_e",
        "EX_arg__L_e",
        "EX_asn__L_e",
        "EX_ala__L_e",
        "EX_trp__L_e"]

FBA = pd.DataFrame({})
FBA_nit = pd.DataFrame({})
#FBA['Carbon source'] = np.zeros(len(carbon))
#FBA_nit['Nitro source'] = np.zeros(len(nitro))

## Roisin minimal mediay
    
reactions_lower_bnd = {"EX_k_e":   -1000,
                       "EX_cl_e":  -1000,
                       "EX_fe3_e": -1000,
                       "EX_fe2_e": -1000,
                       "EX_na1_e":  -1000,
                       "EX_nh4_e": -1000,
                       #"EX_co2_e":   -1000,
                       "EX_so4_e": -1000,
                       "EX_cu2_e": -1000,
                       "EX_mn2_e": -1000,
                       "EX_mobd_e":-1000,
                       "EX_pi_e":  -1000,
                       "EX_o2_e":  -1000,
                       "EX_h_e":   -1000,
                       "EX_h2_e" : -1000,
                       "EX_h2o_e": -1000,
                       "EX_zn2_e": -1000,
                       "EX_cobalt2_e": -1000,
                       "EX_mg2_e" :    -1000,
                       "EX_ca2_e":     -1000
                        }



for idx1, model in enumerate(model_list):
    print(idx1,model)
    model.solver = "gurobi"
    model.objective = 'Growth' #Growth or BIOMASS__2
    model.reactions.BIOMASS_2.bounds= (0,0)
    model.reactions.ATPM.bounds = (3.15,1000)        

    reaction_keys = model.reactions.__dict__['_dict'].keys()
    for exchange in model.exchanges:
        model.reactions.get_by_id(exchange.id).lower_bound= 0
    for rxn in model.reactions:
        if rxn.id.__contains__('DM_'):
            model.reactions.get_by_id(rxn.id).bounds = (0,0)
    for sink in model.sinks:
        model.reactions.get_by_id(sink.id).bounds = (0,1000)
            
    for reaction_key in reactions_lower_bnd.keys():
        if reaction_key in reaction_keys:
            model.reactions.get_by_id(reaction_key).lower_bound = reactions_lower_bnd[reaction_key]
        else:
            print(reaction_key,' not in model ', model)
            continue
        ## Growth rate when growing on single carbon source
    if source == 'cbn':
        for idx, cbn in enumerate(cbn_source):
            try:  
                model.reactions.get_by_id(cbn).lower_bound = - uptake
                
                if model.optimize().objective_value <= cutoff: 
                    FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name), species_name[idx1]] = 0
                elif model.optimize().objective_value > cutoff: 
                    FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name), species_name[idx1]] = 1

                #FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name),species_name[idx1]] = model.optimize().objective_value
                #FBA_results[cbn] = model.optimize().fluxes
                model.reactions.get_by_id(cbn).lower_bound = 0
            except KeyError:
                FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name),species_name[idx1]] = np.nan

            ## Growth rate when growing in single nitrogen source
    elif source == 'ntr':
        for idx, ntr in enumerate(ntr_source):
            model.reactions.get_by_id("EX_nh4_e").lower_bound = 0
            model.reactions.get_by_id("EX_glyc_e").lower_bound = - uptake
            model.reactions.get_by_id("EX_glc__D_e").lower_bound = - uptake

            try:
                model.reactions.get_by_id(ntr).lower_bound = - uptake
                
                if model.optimize().objective_value <= cutoff: 
                    FBA.loc[str(universe.metabolites.get_by_id(ntr[3:]).name), species_name[idx1]] = 0
                elif model.optimize().objective_value > cutoff: 
                    FBA.loc[str(universe.metabolites.get_by_id(ntr[3:]).name), species_name[idx1]] = 1

                model.reactions.get_by_id(ntr).lower_bound = 0
            except KeyError:
                print(KeyError)
                FBA.loc[str(universe.metabolites.get_by_id(ntr[3:]).name),species_name[idx1]] = np.nan
                
total_nan_count_FBA = FBA.isna().sum().sum()
#total_nan_count_FBA_nit = FBA_nit.isna().sum().sum()


FBA.iloc[:, 0:] = (FBA.iloc[:, 0:] > cutoff).astype(int) # add 0 to np.nan values
#FBA_nit.iloc[:, 0:] = (FBA_nit.iloc[:, 0:] > cutoff).astype(int)

#FBA_results = pd.concat([FBA.iloc[:,0:],FBA_nit.iloc[:,0:]])
FBA_results = FBA
FBA_results.loc["Count"] = FBA_results.sum()
FBA_results.to_csv(f'Proposed_media_{model_type}_{source}_FBA_{date}_upt_{uptake}_{arom}.csv')
FBA_results.index = FBA_results.index.str.replace(' R R  2 3 Butanediol C4H10O2', 'Butanediol')
FBA_results.index = FBA_results.index.str.replace('Butyrate (n-C4:0)', 'Butyrate')
FBA_results.index = FBA_results.index.str.replace('N-Acetyl-D-glucosamine(anhydrous)N-Acetylmuramic acid', 'N-Acetyl-D-glucosamine N-Acetylmuramic acid')
FBA_results.index = FBA_results.index.str.replace('Beta alanylL leucine', '\u03B2-alanyl L-leucine')
FBA_results.index = FBA_results.index.str.replace('Beta alanylL glycine', '\u03B2-alanyl L-glycine')
FBA_results.index = FBA_results.index.str.replace('Beta alanylL alanine', '\u03B2-alanyl L-alanine')
FBA_results.index = FBA_results.index.str.replace('L-Carnosine', '\u03B2-alanyl L-histidine')
FBA_results.index = FBA_results.index.str.replace('Beta alanyl beta alanine	', '\u03B2-alanyl \u03B2 alanine')
FBA_results.index = FBA_results.index.str.replace('Beta Alaninamide', '\u03B2-alaninamide')
FBA_results.index = FBA_results.index.str.replace('Orotate C5H3N2O4', 'Orotate')
FBA_results.index = FBA_results.index.str.replace('Guanosine 3 phosphate C10H12N5O8P', 'Guanosine-3-phosphate')
FBA_results.index = FBA_results.index.str.replace('Lactose C12H22O11', 'Lactose')
FBA_results.index = FBA_results.index.str.replace('Thymine C5H6N2O2', 'Thymine')
FBA_results.index = FBA_results.index.str.replace('Linoleic acid (all cis C18:2) n-6', 'Linoleic acid')
FBA_results.index = FBA_results.index.str.replace('Propionate (n-C3:0)', 'Propionate')
FBA_results.index = FBA_results.index.str.replace('Lactose C12H22O11', 'Lactose')
FBA_results.index = FBA_results.index.str.replace('3  AMP C10H12N5O7P', '3-AMP')
FBA_results.index = FBA_results.index.str.replace('N Ribosylnicotinamide C11H15N2O5', 'N Ribosylnicotinamide')
FBA_results.index = FBA_results.index.str.replace('Phthiocerol dimycocerosate A (Mtb) [extracellular]', 'Phthiocerol dimycocerosate A')
FBA_results.index = FBA_results.index.str.replace('UMP C9H11N2O9P', 'UMP')
FBA_results.index = FBA_results.index.str.replace('Maltose C12H22O11', 'Maltose')
FBA_results.index = FBA_results.index.str.replace('Ferrypyoverdine  P putida KT2440 specific', 'Ferrypyoverdine')
FBA_results.index = FBA_results.index.str.replace('Nonanoate C9H17O2', 'Nonanoate')
FBA_results.index = FBA_results.index.str.replace('Cys Gly C5H10N2O3S', 'Cysteinylglycine')
FBA_results.index = FBA_results.index.str.replace('Coprogen unloaded (no Fe(III))', 'Coprogen')
FBA_results.index = FBA_results.index.str.replace('Pyoverdine  P putida specific', 'Pyoverdine')
FBA_results.index = FBA_results.index.str.replace('Glycylglycine C4H8N2O3', 'Glycylglycine')
FBA_results.index = FBA_results.index.str.replace('Beta-Alanine', '\u03B2-alanine')
FBA_results.index = FBA_results.index.str.replace('CO2 CO2', 'CO2')
FBA_results.index = FBA_results.index.str.replace('Alpha-L-Arabinan (3 subunits)', '\u03B1-L-Arabinan')
FBA_results.index = FBA_results.index.str.replace('Iron(III) chelated carboxymycobactin T (R=8 carbon, final carbon is carboxyl group)', 'Iron(III) chelated carboxymycobactin T')
FBA_results.index = FBA_results.index.str.replace('Tetradecanoate (n-C14:0)', 'Tetradecanoate')
FBA_results.index = FBA_results.index.str.replace('Ferroxamine minus Fe(3)', 'Ferroxamine')
FBA_results.index = FBA_results.index.str.replace('Urea CH4N2O', 'Urea')
FBA_results.index = FBA_results.index.str.replace('Octanoate (n-C8:0)', 'Octanoate')
FBA_results.index = FBA_results.index.str.replace('Glycolate C2H3O3', 'Glycolate')
FBA_results.index = FBA_results.index.str.replace('AMP C10H12N5O7P', 'AMP')
FBA_results.index = FBA_results.index.str.replace('Choline C5H14NO', 'Choline')
FBA_results.index = FBA_results.index.str.replace('CMP', 'CMP')
FBA_results.index = FBA_results.index.str.replace('Carboxymycobactin T (R=8 carbon, final carbon is carboxyl group)', 'Carboxymycobactin T')
FBA_results.index = FBA_results.index.str.replace('Beta 1 4 glycosidic  4  beta D mannuronic  1L guluronic acid', '\u03B2-1,4-glycosidic-4-\u03B2-D-mannuronic-1-L-guluronic acid')
FBA_results.index = FBA_results.index.str.replace('L histidinylhistidine', 'L-histidinylhistidine')
FBA_results.index = FBA_results.index.str.replace('D Arabinose C5H10O5', 'D-Arabinose')
FBA_results.index = FBA_results.index.str.replace('Hexanoate (n-C6:0)', 'Hexanoate')
FBA_results.index = FBA_results.index.str.replace('Beta alanyl beta alanine', '\u03B2-alanyl-\u03B2-alanine')
FBA_results.index = FBA_results.index.str.replace('NMN C11H14N2O8P', 'NMN')
FBA_results.index = FBA_results.index.str.replace('Iron bound extracellular staphyloferrin A', 'Iron bound staphyloferrin A')
FBA_results.index = FBA_results.index.str.replace('Taurine C2H7NO3S', 'Taurine')
FBA_results.index = FBA_results.index.str.replace('L glycinylglutamate', 'L-glycinylglutamate')
FBA_results.index = FBA_results.index.str.replace('L glycinylglutamate', 'L-glycinylglutamate')
FBA_results.index = FBA_results.index.str.replace('L glycinylglutamate', 'L-glycinylglutamate')

FBA_results.sort_index(inplace =True)

#FBA_results = FBA_results.rename_axis("Metabolite", axis="index")
#FBA_results["Category"] = FBA_results.index.apply(categorize_metabolite)


# Cluster results
#plt.fivgure(figsize=(10,30))
custom_xticks = ['$M. abscessus$','$M. marinum$','$M. smegmatis$','$M. aromaticivorans$' ,'$M. tuberculosis$']  

clustermap = sns.clustermap(FBA_results, cmap='Blues', annot=False, cbar=None, linewidths=0.5, linecolor='lightblue',col_cluster=False, row_cluster=True, xticklabels = True, yticklabels= True, figsize=figsize) #10,30 max_gene
clustermap.ax_row_dendrogram.set_visible(False)
clustermap.ax_col_dendrogram.set_visible(False)
clustermap.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
clustermap.cax.set_visible(False)
plt.title(f'Uptake rate of {uptake}')
plt.savefig(f'./heatmap_proposed_exp_{source}_{date}_{model_type}_upt_{uptake}_{arom}.png', bbox_inches = 'tight')
    
# Cluster by different media
sns.set(font_scale=1.5)

filtered_FBA_results  = FBA_results[FBA_results.iloc[:, 1:].nunique(axis=1) > 1]
clustermap = sns.clustermap(filtered_FBA_results, cmap='Blues', annot=False, cbar=None, linewidths=0.5, linecolor='lightblue',col_cluster=False, row_cluster=True, xticklabels = True,yticklabels= True, figsize=figsize) #10,30 max_gene
clustermap.ax_row_dendrogram.set_visible(False)
#clustermap.ax_col_dendrogram.set_visible(False)
clustermap.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
clustermap.cax.set_visible(False)
#plt.title(f'Uptake rate of {uptake}')
plt.savefig(f'./heatmap_proposed_dif_exp_{source}_{date}_{model_type}_upt_{uptake}_{arom}.png', bbox_inches = 'tight')
    

#
if arom == 'arom':
    test_media_fast = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 0) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1)& (FBA_results['aromaticivorans'] == 1)]
    if source == 'ntr':
        test_media_fast.index = test_media_fast.index + ' (N)'   
         #test_media_fast.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_fast.index.values.tolist()]
    else:
        test_media_fast.index = test_media_fast.index + ' (C)'   
    
    test_media_non_path = FBA_results[(FBA_results['abscessus'] == 0) & (FBA_results['marinum'] == 0) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1) & (FBA_results['aromaticivorans'] == 1)]
    #test_media_non_path.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_non_path.index.values.tolist()]
    if source == 'ntr':
        test_media_non_path.index = test_media_non_path.index + ' (N)'   
    else:
        test_media_non_path.index = test_media_non_path.index + ' (C)'   
    
    test_media_slow = FBA_results[(FBA_results['abscessus'] == 0) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 1) & (FBA_results['smegmatis'] == 0) & (FBA_results['aromaticivorans'] == 0)]
    #test_media_slow.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_slow.index.values.tolist()]

    test_media_path = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 1) & (FBA_results['smegmatis'] == 0) & (FBA_results['aromaticivorans'] == 0)]
    #test_media_path.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_path.index.values.tolist()]

    test_media_ntm = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1) & (FBA_results['aromaticivorans'] == 1)]
    #test_media_ntm_no_arom = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1)]

else:    
    test_media_fast = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 0) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1)] #& (FBA_results['aromaticivorans'] == 1)]
    #test_media_fast.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_fast.index.values.tolist()]

    test_media_non_path = FBA_results[(FBA_results['abscessus'] == 0) & (FBA_results['marinum'] == 0) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1)]# & (FBA_results['aromaticivorans'] == 1)]
    #test_media_non_path.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_non_path.index.values.tolist()]

    test_media_slow = FBA_results[(FBA_results['abscessus'] == 0) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 1) & (FBA_results['smegmatis'] == 0)] #& (FBA_results['aromaticivorans'] == 0)]
    #test_media_slow.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_slow.index.values.tolist()]

    test_media_path = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 1) & (FBA_results['smegmatis'] == 0)] #& (FBA_results['aromaticivorans'] == 0)]
    #test_media_path.loc[:,'ReactionName'] = [universe.metabolites.get_by_id(x[3:]).name for x in test_media_path.index.values.tolist()]

    test_media_ntm = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1)] #& (FBA_results['aromaticivorans'] == 1)]
    #test_media_ntm_no_arom = FBA_results[(FBA_results['abscessus'] == 1) & (FBA_results['marinum'] == 1) & (FBA_results['Rv'] == 0) & (FBA_results['smegmatis'] == 1)]

test_media_fast.to_csv(f'Test_media_fast_growers_{model_type}_{source}_{date}_upt_{uptake}_{arom}.csv')
test_media_non_path.to_csv(f'Test_media_non_patho_{model_type}_{source}_{date}_upt_{uptake}_{arom}.csv')
test_media_slow.to_csv(f'Test_media_slow_{model_type}_{source}_{date}_{arom}.csv')
test_media_path.to_csv(f'Test_media_path_{model_type}_{source}_{date}_upt_{uptake}_{arom}.csv')

test_media_fast_ntr = pd.read_csv('~/Documents/Carveme/Genomes/ModelsRef1241024Egg/Test_media_fast_growers_max_genes_perc_gf_ntr_280225_upt_0.2_arom.csv')
test_media_fast_cbn = pd.read_csv('~/Documents/Carveme/Genomes/ModelsRef1241024Egg/Test_media_fast_growers_max_genes_perc_gf_cbn_280225_upt_0.2_arom.csv')
test_media_non_path_ntr = pd.read_csv('~/Documents/Carveme/Genomes/ModelsRef1241024Egg/Test_media_non_patho_max_genes_perc_gf_ntr_280225_upt_0.2_arom.csv')
test_media_non_path_cbn = pd.read_csv('~/Documents/Carveme/Genomes/ModelsRef1241024Egg/Test_media_non_patho_max_genes_perc_gf_cbn_280225_upt_0.2_arom.csv')
all_results_fast = pd.concat([test_media_fast_cbn, test_media_fast_ntr], axis =0)
all_results_non_path = pd.concat([test_media_non_path_cbn, test_media_non_path_ntr], axis =0)
all_results_fast.set_index('Unnamed: 0', inplace =True)
all_results_non_path.set_index('Unnamed: 0', inplace =True)
all_results_non_path.index.name = None
all_results_fast.index.name = None

if len(test_media_fast) > 0:

    sns.set(font_scale=1.5)

    clustermap_fast = sns.clustermap(all_results_fast, cmap='Greens', annot=False, cbar=False, linewidths=0.5, linecolor='lightgreen',col_cluster=False, row_cluster=True, xticklabels = True, yticklabels= True, figsize=figsize_2) #10,30 max_gene
    clustermap_fast.ax_row_dendrogram.set_visible(False)
    clustermap_fast.ax_col_dendrogram.set_visible(False)
    clustermap_fast.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
    clustermap_fast.ax_heatmap.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, colors='black', length=  5, width=1)
    clustermap_fast.cax.set_visible(False)
    plt.setp(clustermap_fast.ax_heatmap.get_yticklabels(), rotation=0)
    plt.savefig(f'./heatmap_proposed_exp_fast_{source}_{date}_{model_type}_upt_{uptake}_{arom}.png', bbox_inches = 'tight')

if len(test_media_non_path) > 0:
    sns.set(font_scale=2)
    clustermap_non_path = sns.clustermap(all_results_non_path, cmap='Reds', annot=False, cbar=False, linewidths=0.5, linecolor='pink',col_cluster=False, row_cluster=True,xticklabels = True,  yticklabels= True, figsize=figsize_3) #10,30 max_gene
    clustermap_non_path.ax_row_dendrogram.set_visible(False)
    clustermap_non_path.ax_col_dendrogram.set_visible(False)
    clustermap_non_path.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
    clustermap_non_path.ax_heatmap.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, colors='black', length=7, width=2)
    clustermap_non_path.cax.set_visible(False)
    plt.setp(clustermap_non_path.ax_heatmap.get_yticklabels(), rotation=0)
    plt.savefig(f'./heatmap_proposed_exp_non_patho_{source}_{date}_upt_{uptake}_{model_type}_{arom}.png', bbox_inches = 'tight')


if len(test_media_slow) > 0:
    clustermap_slow = sns.clustermap(test_media_slow, cmap='Blues', annot=False, cbar=True, linewidths=0.5, linecolor='lightblue', col_cluster=False, row_cluster=True, yticklabels= True, figsize=figsize_2) #10,30 max_gene
    clustermap_slow.ax_row_dendrogram.set_visible(False)
    clustermap_slow.ax_col_dendrogram.set_visible(False)
    clustermap_slow.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
    clustermap_slow.ax_heatmap.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, colors='black', length=5, width=1)

    clustermap_slow.cax.set_visible(False)

    plt.savefig(f'./heatmap_proposed_exp_fast_{source}_{date}_upt_{uptake}_{model_type}_{arom}.png', bbox_inches = 'tight')

if len(test_media_path) > 0:
    clustermap_path = sns.clustermap(test_media_path, cmap='Blues', annot=False, cbar=True, linewidths=0.5, linecolor='lightblue', col_cluster=False, row_cluster=True, yticklabels= True, figsize=figsize_2) #10,30 max_gene
    clustermap_path.ax_row_dendrogram.set_visible(False)
    clustermap_path.ax_col_dendrogram.set_visible(False)
    clustermap_path.cax.set_visible(False)
    clustermap_path.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
    clustermap_path.ax_heatmap.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True, colors='black', length=5, width=1)

    plt.savefig(f'./heatmap_proposed_exp_path_{source}_{date}_{model_type}_{arom}.png', bbox_inches = 'tight')

#if len(test_media_ntm) > 0:
#    test_media_ntm.to_csv(f'Test_media_ntm_{model_type}_{source}_{date}_{arom}.csv')
if fluxGalactose:
    fva = pd.DataFrame({})
    fva_results = pd.DataFrame({})

    for idx1, model in enumerate(model_list):
        print(idx1,model)
        model.solver = "gurobi"
        model.objective = 'Growth' #Growth or BIOMASS__2
        model.reactions.BIOMASS_2.bounds= (0,0)
        model.reactions.ATPM.bounds = (3.15,1000)        
        #model.reactions.PGMT.bounds= (0,0.1)

        reaction_keys = model.reactions.__dict__['_dict'].keys()
        for exchange in model.exchanges:
            model.reactions.get_by_id(exchange.id).lower_bound= 0
        for rxn in model.reactions:
            if rxn.id.__contains__('DM_'):
                model.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in model.sinks:
            model.reactions.get_by_id(sink.id).bounds = (0,1000)
                
        for reaction_key in reactions_lower_bnd.keys():
            if reaction_key in reaction_keys:
                model.reactions.get_by_id(reaction_key).lower_bound = reactions_lower_bnd[reaction_key]
            else:
                print(reaction_key,' not in model ', model)
                continue
            ## Growth rate when growing on single carbon source
        #if model ==selected_models['smegmatis'] or model == selected_models['aromaticivorans']:  
            #print('galactitol')
        model.reactions.get_by_id('EX_gal_e').lower_bound = - uptake
        #FBA.loc[str(universe.metabolites.get_by_id('gal_c').name), species_name[idx1]] = model.optimize().objective_value

        #FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name),species_name[idx1]] = model.optimize().objective_value
        fva_results = flux_variability_analysis(model, model.reactions)
        model.reactions.get_by_id('EX_gal_e').lower_bound = 0
    #else:
        #    continue
        fva_results.loc[:,'maximum'].to_csv(f'FVA_galactose_{model}_{date}_{model_type}.csv')

        #except KeyError:
        #    FBA.loc[str(universe.metabolites.get_by_id('galt_c').name),species_name[idx1]] = np.nan
        

if fluxCholSecretion:
    fva = pd.DataFrame({})
    fva_results = pd.DataFrame({})

    for idx1, model in enumerate(model_list):
        print(idx1,model)
        model.solver = "gurobi"
        model.objective = 'Growth' #Growth or BIOMASS__2
        model.reactions.BIOMASS_2.bounds= (0,0)
        model.reactions.ATPM.bounds = (3.15,1000)        
        #model.reactions.PGMT.bounds= (0,0.1)

        reaction_keys = model.reactions.__dict__['_dict'].keys()
        for exchange in model.exchanges:
            model.reactions.get_by_id(exchange.id).lower_bound= 0
        for rxn in model.reactions:
            if rxn.id.__contains__('DM_'):
                model.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in model.sinks:
            model.reactions.get_by_id(sink.id).bounds = (0,1000)
                
        for reaction_key in reactions_lower_bnd.keys():
            if reaction_key in reaction_keys:
                model.reactions.get_by_id(reaction_key).lower_bound = reactions_lower_bnd[reaction_key]
            else:
                print(reaction_key,' not in model ', model)
                continue
            ## Growth rate when growing on single carbon source
        #if model ==selected_models['smegmatis'] or model == selected_models['aromaticivorans']:  
            #print('galactitol')
        model.reactions.get_by_id('EX_glyc_e').lower_bound =  - uptake
        model.reactions.get_by_id('EX_glc__D_e').lower_bound =  - uptake
        model.reactions.get_by_id('EX_pyr_e').lower_bound =  - uptake


        model.reactions.get_by_id('EX_chsterol_e').lower_bound =  0
        model.reactions.get_by_id('EX_chsterol_e').upper_bound =  1000
        #FBA.loc[str(universe.metabolites.get_by_id('gal_c').name), species_name[idx1]] = model.optimize().objective_value

        #FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name),species_name[idx1]] = model.optimize().objective_value
        fva_results = flux_variability_analysis(model, model.reactions)
        model.reactions.get_by_id('EX_chsterol_e').lower_bound = 0
    #else:
        #    continue
        fva_results.loc[:,'maximum'].to_csv(f'FVA_cholesterol_{model}_{date}_{model_type}.csv')

        #except KeyError:
        #    FBA.loc[str(universe.metabolites.get_by_id('galt_c').name),species_name[idx1]] = np.nan
if fluxCO2Secretion:
    fva = pd.DataFrame({})
    fva_results = pd.DataFrame({})
    fba_results = pd.DataFrame({})
    pfba_results = pd.DataFrame({})

    for idx1, model in enumerate(model_list):
        print(idx1,model)
        model.solver = "gurobi"
        model.objective = 'Growth' #Growth or BIOMASS__2
        model.reactions.BIOMASS_2.bounds= (0,0)
        model.reactions.ATPM.bounds = (3.15,1000)        
        #model.reactions.PGMT.bounds= (0,0.1)

        reaction_keys = model.reactions.__dict__['_dict'].keys()
        for exchange in model.exchanges:
            model.reactions.get_by_id(exchange.id).lower_bound= 0
        for rxn in model.reactions:
            if rxn.id.__contains__('DM_'):
                model.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in model.sinks:
            model.reactions.get_by_id(sink.id).bounds = (0,1000)
                
        for reaction_key in reactions_lower_bnd.keys():
            if reaction_key in reaction_keys:
                model.reactions.get_by_id(reaction_key).lower_bound = reactions_lower_bnd[reaction_key]
            else:
                print(reaction_key,' not in model ', model)
                continue
            ## Growth rate when growing on single carbon source
        #if model ==selected_models['smegmatis'] or model == selected_models['aromaticivorans']:  
            #print('galactitol')
        model.reactions.get_by_id('FDH').bounds = (-1000, 1000)
        model.reactions.get_by_id('EX_co2_e').lower_bound =  - uptake


        #FBA.loc[str(universe.metabolites.get_by_id('gal_c').name), species_name[idx1]] = model.optimize().objective_value

        #FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name),species_name[idx1]] = model.optimize().objective_value
        ## FBA normal
        #fba_results = model.optimize().fluxes

        ## pFBA
        pfba_results = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model).fluxes

        #fva_results = flux_variability_analysis(model, model.reactions)
        model.reactions.get_by_id('EX_co2_e').lower_bound = 0
    #else:
        #    continue
        pfba_results.to_csv(f'pFBA_FDH_co2_{model}_{date}_{model_type}.csv')

        #except KeyError:
        #    FBA.loc[str(universe.metabolites.get_by_id('galt_c').name),species_name[idx1]] = np.nan
#if CO2seed:

    ## implement seed        
if fluxTrypUptake: 
    pfba_results = pd.DataFrame({})

    for idx1, model in enumerate(model_list):
        print(idx1,model)
        model.solver = "gurobi"
        model.objective = 'Growth' #Growth or BIOMASS__2
        model.reactions.BIOMASS_2.bounds= (0,0)
        model.reactions.ATPM.bounds = (3.15,1000)        
        #model.reactions.PGMT.bounds= (0,0.1)

        reaction_keys = model.reactions.__dict__['_dict'].keys()
        for exchange in model.exchanges:
            model.reactions.get_by_id(exchange.id).lower_bound= 0
        for rxn in model.reactions:
            if rxn.id.__contains__('DM_'):
                model.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in model.sinks:
            model.reactions.get_by_id(sink.id).bounds = (0,1000)
                
        for reaction_key in reactions_lower_bnd.keys():
            if reaction_key in reaction_keys:
                model.reactions.get_by_id(reaction_key).lower_bound = reactions_lower_bnd[reaction_key]
            else:
                print(reaction_key,' not in model ', model)
                continue
            ## Growth rate when growing on single carbon source
        #if model ==selected_models['smegmatis'] or model == selected_models['aromaticivorans']:  
            #print('galactitol')
        model.reactions.get_by_id('EX_trp__L_e').lower_bound =  - uptake


        
        ## pFBA
        pfba_results = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model).fluxes

        #fva_results = flux_variability_analysis(model, model.reactions)
        model.reactions.get_by_id('EX_trp__L_e').lower_bound = 0
    #else:
        #    continue
        pfba_results.to_csv(f'pFBA_tryp_{model}_{date}_{model_type}.csv')



if fluxGlycUptake:
    fva = pd.DataFrame({})
    fva_results = pd.DataFrame({})

    for idx1, model in enumerate(model_list):
        print(idx1,model)
        model.solver = "gurobi"
        model.objective = 'Growth' #Growth or BIOMASS__2
        model.reactions.BIOMASS_2.bounds= (0,0)
        model.reactions.ATPM.bounds = (3.15,1000)        
        #model.reactions.PGMT.bounds= (0,0.1)

        reaction_keys = model.reactions.__dict__['_dict'].keys()
        for exchange in model.exchanges:
            model.reactions.get_by_id(exchange.id).lower_bound= 0
        for rxn in model.reactions:
            if rxn.id.__contains__('DM_'):
                model.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in model.sinks:
            model.reactions.get_by_id(sink.id).bounds = (0,1000)
                
        for reaction_key in reactions_lower_bnd.keys():
            if reaction_key in reaction_keys:
                model.reactions.get_by_id(reaction_key).lower_bound = reactions_lower_bnd[reaction_key]
            else:
                print(reaction_key,' not in model ', model)
                continue
            ## Growth rate when growing on single carbon source
        #if model ==selected_models['smegmatis'] or model == selected_models['aromaticivorans']:  
            #print('galactitol')
        model.reactions.get_by_id('EX_glyc_e').lower_bound =  - uptake

        #FBA.loc[str(universe.metabolites.get_by_id('gal_c').name), species_name[idx1]] = model.optimize().objective_value

        #FBA.loc[str(universe.metabolites.get_by_id(cbn[3:]).name),species_name[idx1]] = model.optimize().objective_value
        fva_results = flux_variability_analysis(model, model.reactions)
        model.reactions.get_by_id('EX_chsterol_e').lower_bound = 0
    #else:
        #    continue
        fva_results.loc[:,'maximum'].to_csv(f'FVA_glycerol_{model}_{date}_{model_type}.csv')

        #except KeyError:
        #    FBA.loc[str(universe.metabolites.get_by_id('galt_c').name),species_name[idx1]] = np.nan
        


if CompareUptakeRates: 
    # Compare uptake rates

    cbn_max_5  = pd.read_csv(f'Proposed_media_{model_type}_cbn_FBA_031224_upt_5_arom.csv')
    cbn_max_1= pd.read_csv(f'Proposed_media_{model_type}_cbn_FBA_031224_upt_1_arom.csv')
    cbn_max_02= pd.read_csv(f'Proposed_media_{model_type}_cbn_FBA_031224_upt_0.2_arom.csv')
    cbn_max_01= pd.read_csv(f'Proposed_media_{model_type}_cbn_FBA_031224_upt_0.1_arom.csv')
    cbn_max_5.sort_values(by= 'Unnamed: 0', inplace = True)
    cbn_max_5.set_index('Unnamed: 0', inplace=True)
    cbn_max_1.sort_values(by= 'Unnamed: 0', inplace = True)
    cbn_max_1.set_index('Unnamed: 0', inplace=True)
    cbn_max_02.sort_values(by= 'Unnamed: 0', inplace = True)
    cbn_max_02.set_index('Unnamed: 0', inplace=True)
    cbn_max_01.sort_values(by= 'Unnamed: 0', inplace = True)
    cbn_max_01.set_index('Unnamed: 0', inplace=True)


    ntr_max_5= pd.read_csv(f'Proposed_media_{model_type}_ntr_FBA_031224_upt_5_arom.csv')
    ntr_max_1= pd.read_csv(f'Proposed_media_{model_type}_ntr_FBA_031224_upt_1_arom.csv')
    ntr_max_02= pd.read_csv(f'Proposed_media_{model_type}_ntr_FBA_031224_upt_0.2_arom.csv')
    ntr_max_01= pd.read_csv(f'Proposed_media_{model_type}_ntr_FBA_031224_upt_0.1_arom.csv')
    ntr_max_5.sort_values(by= 'Unnamed: 0', inplace = True)
    ntr_max_5.set_index('Unnamed: 0', inplace=True)
    ntr_max_1.sort_values(by= 'Unnamed: 0', inplace = True)
    ntr_max_1.set_index('Unnamed: 0', inplace=True)
    ntr_max_02.sort_values(by= 'Unnamed: 0', inplace = True)
    ntr_max_02.set_index('Unnamed: 0', inplace=True)
    ntr_max_01.sort_values(by= 'Unnamed: 0', inplace = True)
    ntr_max_01.set_index('Unnamed: 0', inplace=True)


    results_cbn = [cbn_max_5,cbn_max_1,cbn_max_02,cbn_max_01] 
    results_ntr = [ntr_max_5,ntr_max_1,ntr_max_02,ntr_max_01] 

    def compare_dataframes(df_list):
        n = len(df_list)
        comparisons = {}
        for i in range(n):
            for j in range(i + 1, n):
                are_equal = df_list[i].equals(df_list[j])
                comparisons[f'DF{i+1} vs DF{j+1}'] = are_equal
        return comparisons
    # Example usage
    # Assuming df1 and df2 are your DataFrames
    result_cbn = compare_dataframes(results_cbn)
    result_ntr = compare_dataframes(results_ntr)
    #result_1_5 = compare_dataframes(cbn_max_5, cbn_max_1)

    for pair, is_equal in result_cbn.items():
        print(f"{pair}: {'CBN Equal' if is_equal else 'Not Equal'}")
    for pair, is_equal in result_ntr.items():
        print(f"{pair}: {'NTR Equal' if is_equal else 'Not Equal'}")

    comparison = cbn_max_01 == cbn_max_02
    rows_with_differences = cbn_max_01[~comparison.all(axis=1)]
    print('cbn_01 vs cbn_02 ', len(rows_with_differences))
    comparison = cbn_max_01 == cbn_max_1
    rows_with_differences = cbn_max_01[~comparison.all(axis=1)]

    print('cbn_01 vs cbn_1 ', len(rows_with_differences))
    comparison = cbn_max_01 == cbn_max_5
    rows_with_differences = cbn_max_01[~comparison.all(axis=1)]
    print('cbn_01 vs cbn_5 ', len(rows_with_differences))
    comparison = cbn_max_02 == cbn_max_1
    rows_with_differences = cbn_max_02[~comparison.all(axis=1)]

    print('cbn_02 vs cbn_1 ', len(rows_with_differences))
    comparison = cbn_max_02 == cbn_max_5
    rows_with_differences = cbn_max_02[~comparison.all(axis=1)]

    print('cbn_02 vs cbn_5 ', len(rows_with_differences))
    comparison = cbn_max_5 == cbn_max_1
    rows_with_differences = cbn_max_5[~comparison.all(axis=1)]

    print('cbn_5 vs cbn_1 ', len(rows_with_differences))

                    
