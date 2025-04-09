import subprocess
import os
import shutil
import datetime
from pathlib import Path
import cobra
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model, validate_sbml_model
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import moma
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from cobra import Model, Reaction, Metabolite
from optlang import Model, Variable, Constraint, Objective
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


plotFit =  False
escher = False
SmegGlyc = True
MariOcd = True
mtb = True
model_type= 'max_genes_perc_gf' #max_genes_perc_gf, max_gf or min_gf
date = '011224'

if model_type == 'min_gf':
    selected_models={'abscessus':'abscessusmore-sensitive651e-50',
                    'marinum':'marinummore-sensitive551e-50',
                    'smegmatis':'gf_smegmatissensitive651',
                    'aromaticivorans':'aromaticivoranssensitive651e-50',
                    'H37Rv':'new_H37Rvfast551e-10'} #minimal models
    old_models = {'smegmatis':'smegmatissensitive651'} #min

elif model_type == 'max_gf':
    selected_models={'abscessus':'abscessusfast301',
                'marinum':'marinummore-sensitive301e-50',
                'smegmatis':'gf_smegmatissensitive401e-10',
                'aromaticivorans':'aromaticivoransmore-sensitive301', # old best model 'aromaticivoransmore-sensitive551e-10'
                'H37Rv':'new_H37Rvmore-sensitive551e-10'} #maximal models
    old_models = {'smegmatis':'smegmatissensitive401e-10'} #max
elif model_type == 'max_genes_perc_gf':
    selected_models = {'abscessus':'abscessusfast301e-10',
            'marinum':'marinumfast301',
            'smegmatis':'gf_smegmatissensitive401e-10',
            'aromaticivorans':'aromaticivoransmore-sensitive301', #  # old best model aromaticivoranssensitive551e-10
            'H37Rv':'new_H37Rvfast501'}
    old_models = {'smegmatis':'smegmatissensitive401e-10'}

smegmatis = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/smegmatis/{selected_models["smegmatis"]}.xml')
marinum = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/marinum/{selected_models["marinum"]}.xml')
aromaticivorans = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/aromaticivorans/{selected_models["aromaticivorans"]}.xml')
abscessus = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/abscessus/{selected_models["abscessus"]}.xml')
Rv = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/H37Rv_new/{selected_models["H37Rv"]}.xml')

model_list  = [abscessus,marinum,smegmatis,aromaticivorans, Rv]
species_name  = ['abscessus','marinum','smegmatis','aromaticivorans', 'Rv']

expDataSmeg = pd.read_excel('./computed_flux_data_from_figures.xlsx', sheet_name = 'smegmatis', header = 1)
#expDataMari = pd.read_excel('./computed_flux_data_from_figures.xlsx', sheet_name = 'marinum', header = 1) #Excel values 


expDataMari = pd.read_excel('./computed_flux_data_from_figures.xlsx', sheet_name = 'marinum', header = 1)

time = expDataMari['time(h)'].values[0:4]
biomass = expDataMari['mg/ml'].values[0:4]
oleate = expDataMari['mmol/L'].values[0:4]
print(f"Fitting {len(time)} points")

def biomass_model(t, b_0, mu):
    return b_0 * np.exp(mu * t)

popt_biomass, _ = curve_fit(biomass_model, time, biomass, p0 = [0.0038,0.087])
b_0, mu = popt_biomass

print(f"Initial Biomass: {b_0:.3f}, Growth Rate: {mu:.3f}")

def oleate_model(t, ol_0, v_ol):
    b_t = biomass_model(t, b_0, mu)
    return ol_0 - v_ol * np.cumsum(b_t * np.gradient(t))

popt_ol, pcov_ol = curve_fit(oleate_model, time, oleate, p0 = [0.6,0.34])
ol_0, v_ol = popt_ol

print(f"Initial Oleate Concentration: {ol_0:.3f}, Uptake Rate: {v_ol:.3f}" )

# calculate r-squared oleate fit
ol_pred = oleate_model(time, *popt_ol)
ss_res_oleate = np.sum((oleate - ol_pred) ** 2)
ss_tot_oleate = np.sum((oleate - np.mean(oleate)) ** 2)
r_squared_oleate = 1 - (ss_res_oleate / ss_tot_oleate)

bm_pred = biomass_model(time, *popt_biomass)
ss_res_biomass = np.sum((biomass - bm_pred) ** 2)
ss_tot_biomass = np.sum((biomass - np.mean(biomass)) ** 2)
r_squared_biomass = 1 - (ss_res_biomass / ss_tot_biomass)

print(f"R-squared for oleate model: {r_squared_oleate:.3f}, R-squared for biomass: {r_squared_biomass:.3f}")


if plotFit: 
    # Plotting the data and fitted curves
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.scatter(time, biomass, label="Observed biomass")
    plt.plot(time, biomass_model(time, *popt_biomass), color="r", label="Fitted Biomass Model")
    plt.xlabel("Time")
    plt.ylabel("Biomass Concentration")
    plt.legend()
    plt.title("Biomass Growth Fit")

    plt.subplot(1, 2, 2)
    plt.scatter(time, oleate, label="Observed oleate")
    plt.plot(time, oleate_model(time, *popt_ol), color="r", label="Fitted Oleate Model")
    plt.xlabel("Time")
    plt.ylabel("Oleate Concentration")
    plt.legend()
    plt.title("Oleate Uptake Fit")

    plt.tight_layout()
    plt.show()



Sol = pd.DataFrame({})
flux_sol = pd.DataFrame({})
## Roisin minimal media
    
RMM_lower_bnd = {"EX_k_e":   -1000,
                       "EX_cl_e":  -1000,
                       "EX_fe3_e": -1000,
                       "EX_fe2_e": -1000,
                       "EX_na1_e":  -1000,
                       "EX_nh4_e": -1000,
                       "EX_co2_e":   -1000,
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

mSMM_lower_bnd = {"EX_k_e":   -1000,
                       "EX_cl_e":  -1000,
                       "EX_fe3_e": -1000,
                       "EX_fe2_e": -1000,
                       "EX_na1_e":  -1000,
                       "EX_nh4_e": -1000,
                       "EX_co2_e":   -1000,
                       "EX_so4_e": -1000,
                       "EX_cu2_e": -1000,
                       "EX_mn2_e": -1000,
                       #"EX_mobd_e":-1000,
                       "EX_pi_e":  -1000,
                       "EX_o2_e":  -1000,
                       "EX_h_e":   -1000,
                       "EX_h2_e" : -1000,
                       "EX_h2o_e": -1000,
                       "EX_zn2_e": -1000,
                       "EX_cobalt2_e": -1000,
                       "EX_mg2_e" :    -1000,
                       "EX_ca2_e":     -1000,
                       "EX_cit_e": -1000
                        }

# ---- Smeg -----

if SmegGlyc:
#    modelNorm= smegmatis
    for expIdx, GrowthRate in enumerate(expDataSmeg['GrowthRateGlyc']):
        smegmatis = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/smegmatis/{selected_models["smegmatis"]}.xml')

        modelExp = smegmatis
        modelExp.objective = 'Growth'
        modelExp.solver = "gurobi"
        modelExp.reactions.BIOMASS_2.bounds= (0,0)
        modelExp.reactions.ATPM.bounds = (3.15,1000)        

        reaction_keys = modelExp.reactions.__dict__['_dict'].keys()
        for exchange in modelExp.exchanges:
            modelExp.reactions.get_by_id(exchange.id).lower_bound = 0
            modelExp.reactions.get_by_id(exchange.id).upper_bound = 1000

        for rxn in modelExp.reactions:
            if rxn.id.__contains__('DM_'):
                modelExp.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in modelExp.sinks:
            modelExp.reactions.get_by_id(sink.id).bounds = (0,1000)
        for reaction_key in RMM_lower_bnd.keys():
            if reaction_key in reaction_keys:
                modelExp.reactions.get_by_id(reaction_key).lower_bound = RMM_lower_bnd[reaction_key]
        # -- add experimental values --
        error = [[0.2,0.2,0.2,0.2], [0.2,0.2,0.89,0.2] ,[0.2,0.2,0.25,0.2],[0.2,0.2,0.21,0.2]] 
        modelExp.reactions.get_by_id('EX_o2_e').bounds = (-1 * float(expDataSmeg['qO2Glyc'][expIdx]) - error[expIdx][0] * float(expDataSmeg['qO2Glyc'][expIdx]), -1 * float(expDataSmeg['qO2Glyc'][expIdx]) + error[expIdx][0] * float(expDataSmeg['qO2Glyc'][expIdx])) 
        modelExp.reactions.get_by_id('EX_co2_e').bounds = (1 * float(expDataSmeg['qCO2Glyc'][expIdx]) - error[expIdx][1] * float(expDataSmeg['qCO2Glyc'][expIdx]), 1 * float(expDataSmeg['qCO2Glyc'][expIdx]) + error[expIdx][1] * float(expDataSmeg['qCO2Glyc'][expIdx]))
        modelExp.reactions.get_by_id('EX_glyc_e').bounds = (-1 * float(expDataSmeg['UptakeRateGlyc'][expIdx]) - error[expIdx][2] * float(expDataSmeg['UptakeRateGlyc'][expIdx]), -1 * float(expDataSmeg['UptakeRateGlyc'][expIdx]) + error[expIdx][2] * float(expDataSmeg['UptakeRateGlyc'][expIdx]))
        modelExp.reactions.get_by_id('Growth').bounds = (1 * float(expDataSmeg['GrowthRateGlyc'][expIdx]) - error[expIdx][3] * float(expDataSmeg['GrowthRateGlyc'][expIdx]), 1 * float(expDataSmeg['GrowthRateGlyc'][expIdx]) + error[expIdx][3] * float(expDataSmeg['GrowthRateGlyc'][expIdx]))
        growthrate =  str(float(expDataSmeg['GrowthRateGlyc'][expIdx]))[0:5]
        Sol.loc[expIdx,'Experiment1'] = 'T. Hooper'
        Sol.loc[expIdx,'CarbonExp1'] = 'EX_glyc_e'
        Sol.loc[expIdx,'DilutionRateExp1'] = growthrate 
        Sol.loc[expIdx,'ErrorExp1'] = 'Error: ' + str(error)
        try:
            moma_solution = moma(modelExp)
            for idx, ex in enumerate(smegmatis.exchanges):
                flux_sol.loc[idx,f'Experiment1_smegGlyc_{growthrate}_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == ex.id][0]
                flux_sol.loc[idx,f'Experiment1_{growthrate}_flux'] = moma_solution.fluxes[ex.id]
            Sol.loc[expIdx,'SolutionStatusExp1'] = moma_solution
        except:
            Sol.loc[expIdx,'SolutionStatusExp1'] = 'Infeasible'
            
    for expIdx, GrowthRate in enumerate(expDataSmeg['GrowthRateSucc']):
        smegmatis = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/smegmatis/{selected_models["smegmatis"]}.xml')

        modelExp = smegmatis
        modelExp.objective = 'Growth'
        modelExp.solver = "gurobi"
        modelExp.reactions.BIOMASS_2.bounds= (0,0)
        modelExp.reactions.ATPM.bounds = (3.15,1000)        

        reaction_keys = modelExp.reactions.__dict__['_dict'].keys()
        for exchange in modelExp.exchanges:
            modelExp.reactions.get_by_id(exchange.id).lower_bound = 0
            modelExp.reactions.get_by_id(exchange.id).upper_bound = 1000

        for rxn in modelExp.reactions:
            if rxn.id.__contains__('DM_'):
                modelExp.reactions.get_by_id(rxn.id).bounds = (0,0)
        for sink in modelExp.sinks:
            modelExp.reactions.get_by_id(sink.id).bounds = (0,1000)
        for reaction_key in RMM_lower_bnd.keys():
            if reaction_key in reaction_keys:
                modelExp.reactions.get_by_id(reaction_key).lower_bound = RMM_lower_bnd[reaction_key]
        # -- add experimental values --
        error = [[1.5,1,0.4,0.4], [1.4,1,0.3,0.3] ,[1.4,1,0.3,0.3],[1.4,1,0.3,0.3]] #o2 errors when fitted to line (0.21,0.24,0.37 to 1.40); #co2 errors according to dil. rate, (0.01,0.73, 0.95, 0.15%) succ, growth
        #modelExp.reactions.get_by_id('EX_o2_e').bounds = (-1000, 1000) 
        #modelExp.reactions.get_by_id('EX_co2_e').bounds = (-1000, 1000) 

        modelExp.reactions.get_by_id('EX_o2_e').bounds = (-1 * float(expDataSmeg['qO2Succ'][expIdx]) - error[expIdx][0] * float(expDataSmeg['qO2Succ'][expIdx]), -1 * float(expDataSmeg['qO2Succ'][expIdx]) + error[expIdx][0] * float(expDataSmeg['qO2Succ'][expIdx])) 
        modelExp.reactions.get_by_id('EX_co2_e').bounds = (1 * float(expDataSmeg['qCO2Succ'][expIdx]) - error[expIdx][1] * float(expDataSmeg['qCO2Succ'][expIdx]), 1 * float(expDataSmeg['qCO2Succ'][expIdx]) + error[expIdx][1] * float(expDataSmeg['qCO2Succ'][expIdx]))
        modelExp.reactions.get_by_id('EX_succ_e').bounds = (-1 * float(expDataSmeg['UptakeRateSucc'][expIdx]) - error[expIdx][2] * float(expDataSmeg['UptakeRateSucc'][expIdx]), -1 * float(expDataSmeg['UptakeRateSucc'][expIdx]) + error[expIdx][2] * float(expDataSmeg['UptakeRateSucc'][expIdx]))
        modelExp.reactions.get_by_id('Growth').bounds = (1 * float(expDataSmeg['GrowthRateSucc'][expIdx]) - error[expIdx][3] * float(expDataSmeg['GrowthRateSucc'][expIdx]), 1 * float(expDataSmeg['GrowthRateSucc'][expIdx]) + error[expIdx][3] * float(expDataSmeg['GrowthRateSucc'][expIdx]))
        #print('dilution smeg succ: ',modelExp.reactions.get_by_id('Growth').bounds, 'uptake: ',modelExp.reactions.get_by_id('EX_succ_e').bounds)
        
        growthrate =  str(float(expDataSmeg['GrowthRateSucc'][expIdx]))[0:5]
        Sol.loc[expIdx,'ExperimentExp2'] = 'T. Hooper'
        Sol.loc[expIdx,'CarbonExp2'] = 'EX_succ_e'
        Sol.loc[expIdx,'DilutionRateExp2'] = growthrate
        Sol.loc[expIdx,'ErrorExp2'] = 'Error: ' + str(error[expIdx])
        try:
            moma_solution = moma(modelExp)
            for idx, ex in enumerate(smegmatis.exchanges):
                flux_sol.loc[idx,f'Experiment2_smegSucc_{growthrate}_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == ex.id][0]
                flux_sol.loc[idx,f'Experiment2_{growthrate}_flux'] = moma_solution.fluxes[ex.id]
        
            Sol.loc[expIdx,'SolutionStatusExp2'] = moma_solution
            print(moma_solution)
        except:
            Sol.loc[expIdx,'SolutionStatusExp2'] = 'Infeasible'
            

if MariOcd:
    marinum = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/marinum/{selected_models["marinum"]}.xml')
    modelExp = marinum
    modelExp.objective = 'Growth'
    modelExp.solver = "gurobi"
    modelExp.reactions.BIOMASS_2.bounds= (0,0)
    modelExp.reactions.ATPM.bounds = (3.15,1000)        

    reaction_keys = modelExp.reactions.__dict__['_dict'].keys()
    for exchange in modelExp.exchanges:
        modelExp.reactions.get_by_id(exchange.id).lower_bound = 0
        modelExp.reactions.get_by_id(exchange.id).upper_bound = 1000

    for rxn in modelExp.reactions:
        if rxn.id.__contains__('DM_'):
            modelExp.reactions.get_by_id(rxn.id).bounds = (0,0)
    for sink in modelExp.sinks:
        modelExp.reactions.get_by_id(sink.id).bounds = (0,1000)
    for reaction_key in mSMM_lower_bnd.keys():
        if reaction_key in reaction_keys:
            modelExp.reactions.get_by_id(reaction_key).lower_bound = mSMM_lower_bnd[reaction_key]
    # -- add experimental values --
    error = 0.1
    
    modelExp.reactions.get_by_id('EX_ocdcea_e').bounds = (-v_ol - error * v_ol , -v_ol + error * v_ol )
    modelExp.reactions.get_by_id('Growth').bounds = (mu - error * mu, 1 * mu + error * mu)
    Sol.loc[0,'Experiment3'] = 'Dong'
    Sol.loc[0,'CarbonExp3'] = 'EX_ocdcea_e'
    #Sol.loc[0,'DilutionRateExp3'] = str(float(expDataMari.loc[0,'GrowthRateOle']))[0:5] 
    Sol.loc[0,'DilutionRateExp3'] = str(mu)[0:5] 
    Sol.loc[0,'ErrorExp3'] = 'Error: ' + str(error)
    try:
        moma_solution = moma(modelExp)
        for idx, ex in enumerate(marinum.exchanges):
            flux_sol.loc[idx,'Experiment3_mariOleate_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == ex.id][0]
            flux_sol.loc[idx,'Experiment3_flux'] = moma_solution.fluxes[ex.id]
        Sol.loc[0,'SolutionStatusExp3'] = moma_solution

        print(moma_solution)
    except:
        Sol.loc[0,'SolutionStatusExp3'] = 'Infeasible'

if mtb:

    modelExp = Rv
    modelExp.objective = 'Growth'
    modelExp.solver = "gurobi"
    modelExp.reactions.BIOMASS_2.bounds= (0,0)
    modelExp.reactions.ATPM.bounds = (3.15,1000)        

    reaction_keys = modelExp.reactions.__dict__['_dict'].keys()
    for exchange in modelExp.exchanges:
        modelExp.reactions.get_by_id(exchange.id).lower_bound = 0
        modelExp.reactions.get_by_id(exchange.id).upper_bound = 1000

    for rxn in modelExp.reactions:
        if rxn.id.__contains__('DM_'):
            modelExp.reactions.get_by_id(rxn.id).bounds = (0,0)
    for sink in modelExp.sinks:
        modelExp.reactions.get_by_id(sink.id).bounds = (0,1000)
    for reaction_key in RMM_lower_bnd.keys():
        if reaction_key in reaction_keys:
            modelExp.reactions.get_by_id(reaction_key).lower_bound = RMM_lower_bnd[reaction_key]
    # -- add experimental values --
    error = 0.1
    
    modelExp.reactions.get_by_id('EX_glyc_e').bounds = (-0.39 - error * (0.39), -0.39 + error * (0.39))
    modelExp.reactions.get_by_id('EX_ocdcea_e').bounds = (-0.003 - error * (0.003), -0.003 + error * (0.003))
    modelExp.reactions.get_by_id('EX_co2_e').bounds = (0.23 - error * (0.23), 0.23 + error * (0.23))
    modelExp.reactions.get_by_id('Growth').bounds =  (0.01 - error * (0.01), 0.01 + error * (0.01))
    try:
        moma_solution = moma(modelExp)
        for idx, ex in enumerate(modelExp.exchanges):
            flux_sol.loc[idx,'Experiment4_MTBGlycBeste_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == ex.id][0]
            flux_sol.loc[idx,'Experiment4_flux'] = moma_solution.fluxes[ex.id]
        for idx, rxn in enumerate(modelExp.reactions):
            flux_sol.loc[idx,f'Experiment4_MTBGlycBeste_all_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == rxn.id][0]
            flux_sol.loc[idx,f'Experiment4_MTBGlycBeste_all_flux'] = moma_solution.fluxes[rxn.id]
            
        Sol.loc[0,'Experiment4'] = 'Beste'
        Sol.loc[0,'CarbonExp4'] = 'EX_glyc_e' + 'EX_ocdcea_e'
        Sol.loc[0,'DilutionRateExp4'] = str(0.01) 
        Sol.loc[0,'ErrorExp4'] = 'Error: ' + str(error)
        Sol.loc[0,'SolutionStatusExp4'] = moma_solution
        print(moma_solution)
    except:
        Sol.loc[0,'SolutionStatusExp4'] = 'Infeasible'
    # --- second set of data (from Borah)
    modelExp = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/H37Rv_new/{selected_models["H37Rv"]}.xml')

    modelExp.objective = 'Growth'
    modelExp.solver = "gurobi"
    modelExp.reactions.BIOMASS_2.bounds= (0,0)
    modelExp.reactions.ATPM.bounds = (3.15,1000)        

    reaction_keys = modelExp.reactions.__dict__['_dict'].keys()
    for exchange in modelExp.exchanges:
        modelExp.reactions.get_by_id(exchange.id).lower_bound = 0
        modelExp.reactions.get_by_id(exchange.id).upper_bound = 1000

    for rxn in modelExp.reactions:
        if rxn.id.__contains__('DM_'):
            modelExp.reactions.get_by_id(rxn.id).bounds = (0,0)
    for sink in modelExp.sinks:
        modelExp.reactions.get_by_id(sink.id).bounds = (0,1000)
    for reaction_key in RMM_lower_bnd.keys():
        if reaction_key in reaction_keys:
            modelExp.reactions.get_by_id(reaction_key).lower_bound = RMM_lower_bnd[reaction_key]
    # -- add experimental values --
    error = 0.1
    
    modelExp.reactions.get_by_id('EX_ac_e').bounds = (-0.26 - error * (0.26), -0.26 + error * (0.26))
    modelExp.reactions.get_by_id('EX_chsterol_e').bounds = (-0.0085 - error * (0.0085), -0.0085 + error * (0.0085))
    modelExp.reactions.get_by_id('EX_co2_e').bounds = (0.245 - error * (0.245), 0.245 + error * (0.245))
    modelExp.reactions.get_by_id('Growth').bounds =  (0.01 - error * (0.01), 0.01 + error * (0.01))
    try:
        moma_solution = moma(modelExp)
        for idx, ex in enumerate(modelExp.exchanges):
            flux_sol.loc[idx,'Experiment5_MTBace_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == ex.id][0]
            flux_sol.loc[idx,'Experiment5_flux'] = moma_solution.fluxes[ex.id]
        Sol.loc[0,'Experiment5'] = 'Borah'
        Sol.loc[0,'CarbonExp5'] = 'EX_ac_e' + 'EX_chsterol_e'
        Sol.loc[0,'DilutionRateExp5'] = str(0.01) 
        Sol.loc[0,'ErrorExp5'] = 'Error: ' + str(error)
        Sol.loc[0,'SolutionStatusExp5'] = moma_solution
        print(moma_solution)
    except:
        Sol.loc[0,'SolutionStatusExp5'] = 'Infeasible'
    # --- third set of data (from Borah)
    modelExp =  read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/H37Rv_new/{selected_models["H37Rv"]}.xml')

    modelExp.objective = 'Growth'
    modelExp.solver = "gurobi"
    modelExp.reactions.BIOMASS_2.bounds= (0,0)
    modelExp.reactions.ATPM.bounds = (3.15,1000)        

    reaction_keys = modelExp.reactions.__dict__['_dict'].keys()
    for exchange in modelExp.exchanges:
        modelExp.reactions.get_by_id(exchange.id).lower_bound = 0
        modelExp.reactions.get_by_id(exchange.id).upper_bound = 1000

    for rxn in modelExp.reactions:
        if rxn.id.__contains__('DM_'):
            modelExp.reactions.get_by_id(rxn.id).bounds = (0,0)
    for sink in modelExp.sinks:
        modelExp.reactions.get_by_id(sink.id).bounds = (0,1000)
    for reaction_key in RMM_lower_bnd.keys():
        if reaction_key in reaction_keys:
            modelExp.reactions.get_by_id(reaction_key).lower_bound = RMM_lower_bnd[reaction_key]
    # -- add experimental values --
    error = 0.1
    
    modelExp.reactions.get_by_id('EX_glyc_e').bounds = (-0.23 - error * (0.23), -0.23 + error * (0.23))
    modelExp.reactions.get_by_id('EX_ocdcea_e').bounds = (-0.002 - error * (0.002), -0.002 + error * (0.002))
    modelExp.reactions.get_by_id('EX_co2_e').bounds = (0.178 - error * (0.178), 0.178 + error * (0.178))
    modelExp.reactions.get_by_id('Growth').bounds =  (0.01 - error * (0.01), 0.01 + error * (0.01))
    try:
        moma_solution = moma(modelExp)
        for idx, ex in enumerate(modelExp.exchanges):
            flux_sol.loc[idx,'Experiment6_MTBGlyc_Bor_rxn'] = moma_solution.fluxes.index[moma_solution.fluxes.index == ex.id][0]
            flux_sol.loc[idx,'Experiment6_flux'] = moma_solution.fluxes[ex.id]
        Sol.loc[0,'Experiment6'] = 'Borah'
        Sol.loc[0,'CarbonExp6'] = 'EX_glyc_e' + 'EX_ocdcea_e'
        Sol.loc[0,'DilutionRateExp6'] = str(0.01) 
        Sol.loc[0,'ErrorExp6'] = 'Error: ' + str(error)
        Sol.loc[0,'SolutionStatusExp6'] = moma_solution
        print(moma_solution)
    except:
        Sol.loc[0,'SolutionStatusExp6'] = 'Infeasible'
     
    
      
flux_sol.to_csv(f'validation_flux_data_fluxes_{model_type}_{date}.csv')
Sol.to_csv(f'validation_flux_data_{model_type}_{date}.csv')
if escher: 
    Rv = read_sbml_model('./H37Rv_new/new_H37Rvfast501.xml')
    save_json_model(filename= 'new_H37Rvfast501.json', model =Rv)
