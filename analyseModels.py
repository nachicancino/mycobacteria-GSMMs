import numpy as np
import cobra
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model, validate_sbml_model
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
from cobra.util.solver import linear_reaction_coefficients
from supervenn import supervenn 
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from venn import venn
from collections import namedtuple
from scipy.stats import fisher_exact, false_discovery_control
from sklearn.linear_model import LinearRegression
from striprtf.striprtf import rtf_to_text
from graphviz import Digraph
import os
from Bio import SeqIO
import seaborn as sns


## OPTIONS
subsystem_plot_ratio  = True
addSubsystem =True
plotSuperVenn = True
plotGenesRxns =True
groupSubsystem = True
fisherTest = True
analyseGapfill = True
gapFillManual = False
modifyModels =False
checkMissingMet = False
ressemblanceMTB = True
heatmapDifMin = True
latexHierarchySubsystem =True
##
date = '280125'
level_subsystem = 'detailed' # or "detailed" 'high' or 'raw'
arom = 'arom' #or arom
model_type= 'max_genes_perc_gf' #max_genes_perc_gf, max_gf or min_gf
if model_type == 'min_gf':
    selected_models={'abscessus':'abscessusmore-sensitive651e-50',
                    'marinum':'marinummore-sensitive551e-50',
                    'smegmatis':'gf_smegmatissensitive651',
                    'aromaticivorans':'aromaticivoransfast651e-10', # aromaticivoranssensitive651e-10 aromaticivoranssensitive651e-50
                    'H37Rv':'new_H37Rvfast551e-10'} #minimal models
    old_models = {'smegmatis':'smegmatissensitive651'} #min
    plotSubsystemMaxgen  =False 
    plotSubsystemMin = True


elif model_type == 'max_gf':
    selected_models={'abscessus':'abscessusfast301',
                'marinum':'marinummore-sensitive301e-50',
                'smegmatis':'gf_smegmatissensitive401e-10',
                'aromaticivorans':'aromaticivoranssensitive301',#'aromaticivoransmore-sensitive301', # old best model 'aromaticivoransmore-sensitive551e-10'
                'H37Rv':'new_H37Rvmore-sensitive551e-10'} #maximal models
    old_models = {'smegmatis':'smegmatissensitive401e-10'} #max
elif model_type == 'max_genes_perc_gf':
    selected_models = {'abscessus':'abscessusfast301e-10',
            'marinum':'marinumfast301',
            'smegmatis':'gf_smegmatissensitive401e-10',
            'aromaticivorans':'aromaticivoranssensitive301', #'aromaticivoransmore-sensitive301', #  # old best model aromaticivoranssensitive551e-10
            'H37Rv':'new_H37Rvfast501'}
    old_models = {'smegmatis':'smegmatissensitive401e-10'}
    plotSubsystemMaxgen =True
    plotSubsystemMin = False

smegmatis = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/smegmatis/{selected_models["smegmatis"]}.xml')
marinum = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/marinum/{selected_models["marinum"]}.xml')
aromaticivorans = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/aromaticivorans/{selected_models["aromaticivorans"]}.xml')
abscessus = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/abscessus/{selected_models["abscessus"]}.xml')
Rv = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/H37Rv_new/{selected_models["H37Rv"]}.xml')
Rv_reference = read_sbml_model('/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/modified_H37Rvfast501e_10/modified_H37Rvfast501e_10.xml')
print(len(Rv_reference.reactions),len(Rv_reference.metabolites),len(Rv_reference.genes))
universe = read_sbml_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/universe/universe_bact_arch_3.xml')
## -- Smegmatis remove one of GLCRD_1  and GLCRAL_1 (already exist their equivalent GLCRD and GLCRAL )
## Replacing MCTC_1 by MCCC (same formula)
## for aromaticivorans remove 4CMLCL_kt (repeated: CMLDC)

#if model_type =='min_gf':
#    smegmatis.remove_reactions(['GLCRD_1','GLCRAL_1','MCTC_1'])
#    smegmatis.add_reactions([universe.reactions.get_by_id('MCCC')])
#    marinum.remove_reactions(['GUAtex','EX_gua_e'])
#    aromaticivorans.remove_reactions(['4CMLCL_kt','GUAtex','EX_gua_e']) 
#    Rv.remove_reactions(['GUAtex','EX_gua_e'])


genome_folder ='/home/icancino/Documents/Carveme/Genomes/eggnog_fasta_5'
genome_list  = ['/abscessus.fasta','/marinum.fasta','/smegmatis.fasta','/aromaticivorans.fasta', '/H37Rv.fasta']


abscessus_annotations = pd.read_excel(genome_folder+'/abscessus_annotations.xlsx', header= 2)
marinum_annotations  = pd.read_excel(genome_folder+'/marinum_annotations.xlsx', header= 2 )
smegmatis_annotations  = pd.read_csv(genome_folder+'/smegmatis_annotations.csv' , header= 4)
aromaticivorans_annotations  = pd.read_excel(genome_folder+'/aromaticivorans_annotations.xlsx' , header= 2)
Rv_annotations  = pd.read_excel(genome_folder+'/H37Rv_annotations.xlsx', header= 2 )


old_smegmatis = read_sbml_model(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/smegmatis/{old_models["smegmatis"]}.xml')
#old_abscessus = read_sbml_model(f'/home/icancino/Documents/Models/ModelsRef1201024/abscessus/{old_models["abscessus"]}.xml')
#old_marinum = read_sbml_model(f'/home/icancino/Documents/Models/ModelsRef1201024/marinum/{old_models["marinum"]}.xml')


if arom == 'arom':
    model_list  = [abscessus,marinum,smegmatis,aromaticivorans, Rv]
    species_name  = ['abscessus','marinum','smegmatis','aromaticivorans', 'Rv']
    species_folder = ['abscessus','marinum','smegmatis', 'aromaticivorans','H37rv_new']


elif arom == 'noarom':
    model_list = [abscessus,marinum,smegmatis,Rv]
    species_name = ['abscessus','marinum','smegmatis', 'Rv']
    species_folder = ['abscessus','marinum','smegmatis','H37rv_new']


## Removing RXN
#marinum.remove_reactions([marinum.reactions.MDH3])
## Adding SUBSYSTEMS
## REFERENCE MODEL

#Rv = load_json_model('/home/icancino/Documents/Carveme/Genomes/MTB H37Rv/iEK1011_2.0.json')
directory_path  = os.getcwd()


old_data_myco= {'M. tuberculosis H37Rv': {'Reactions': len(Rv.reactions),'Genes': 3906}, 
            'M. marinum':{'Reactions': len(marinum.reactions),'Genes': 5141},
            'M. smegmatis':{'Reactions': len(old_smegmatis.reactions),'Genes': 6416},
            'M. abscessus':{'Reactions': len(abscessus.reactions),'Genes': 4939},
            'M. aromaticivorans':{'Reactions': len(aromaticivorans.reactions),'Genes': 5866}}
gf_carve_data_myco= {'M. tuberculosis H37Rv': {'Reactions': len(Rv.reactions),'Genes': 3906}, 
            'M. marinum':{'Reactions': len(marinum.reactions),'Genes': 5141},
            'M. smegmatis':{'Reactions': len(smegmatis.reactions),'Genes': 6416},
            'M. abscessus':{'Reactions': len(abscessus.reactions),'Genes': 4939},
            'M. aromaticivorans':{'Reactions': len(aromaticivorans.reactions),'Genes': 5866}}
            
if analyseGapfill:
    smeg_id = []
    abs_id = []
    mari_id = []

    old_smeg_id = []
    old_abs_id = []
    old_mari_id = []

    for rxn in smegmatis.reactions:
        smeg_id.append(rxn.id)
    for rxn in abscessus.reactions:
        abs_id.append(rxn.id)
    for rxn in marinum.reactions:
        mari_id.append(rxn.id)
    
    for rxn in old_smegmatis.reactions:
        old_smeg_id.append(rxn.id)
    #for rxn in old_abscessus.reactions:
    #    old_abs_id.append(rxn.id)
    #for rxn in old_marinum.reactions:
    #    old_mari_id.append(rxn.id)
    
    differences = {'dif_smeg': set(smeg_id) - set(old_smeg_id)}#,
                #'dif_abs' : set(abs_id) - set(old_abs_id), 'dif_mari':set(mari_id) - set(old_mari_id)}


    print(differences['dif_smeg'])#, differences['dif_abs'], differences['dif_mari'])
if gapFillManual:
    media = pd.read_csv('/home/icancino/Documents/PhD/Reconstruction MTB analysis/validation_media_compounds_DR_IC.csv')
    media_bounds_dict= {}
    results_gf = pd.DataFrame({})
    for media_no in media.keys(): 
        if media_no.endswith('_bounds') == False:
            media_bounds_dict[media_no] ={}
            media_bounds_dict[media_no] = dict(zip(media[media_no],media[str(media_no)+'_bounds']))

    for species in old_models.keys():
        #os.mkdir(directory_path+'/man_gf_'+str(species))
        #os.chdir(directory_path+'/man_gf_'+str(species))
        filename = os.fsdecode('./'+species+'/'+old_models[species]+'.xml')        
        out = 'man_gf_'+ str(old_models[species])+'.xml'
        model = cobra.io.read_sbml_model(filename) 

        if species == 'smegmatis':
            reactions_to_add = list(differences['dif_smeg']) #['RMPA', 'RMI', 'RMK','RMNtex','RMNtpp']

        elif species == 'abscessus':
            reactions_to_add = list(differences['dif_abs']) #['DLYSOXGAT','MS_1','ALDD31',]
        
        #model.add_reactions(reactions_to_add)
        #cobra.io.write_sbml_model(model,str(out))

        model.objective = 'Growth'
        reaction_keys = model.reactions.__dict__['_dict'].keys()        
        for rxn_add in reactions_to_add:
            model.add_reactions([universe.reactions.get_by_id(rxn_add)])

            for media_no in media_bounds_dict.keys():

                model.reactions.BIOMASS_2.bounds = (0,0)        
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
                        model.reactions.get_by_id(reaction).lower_bound = int(media_bounds_dict[media_no][str(reaction)])
                    elif not pd.isna(reaction): ## add missing reaction in another file
                        print('Reaction ',reaction,' in media but not in the model')
                        #globals()[f"{filename[len(dir+'/'+species+'/'):-4]}_missing_rxn"].add(reaction)
                        
                        continue 
            
                if model.optimize().objective_value <= 0.0000001:
                    print(rxn_add,media_no,model.optimize().objective_value)
                    results_gf.loc[species+'_'+rxn_add,str(media_no)] = 0 
                elif model.optimize().objective_value > 0.0000001:
                    print( rxn_add,media_no,model.optimize().objective_value)
                    results_gf.loc[species+'_'+rxn_add,str(media_no)] = 1
                else:
                    print('Error',media_no)
            #if model.reactions.get_by_id(rxn) :
            model.remove_reactions([model.reactions.get_by_id(rxn_add)])
            #else:
        #    continue
    results_gf.to_csv('./manual_gapfilling_results.csv')

        #os.chdir('../')  
if ressemblanceMTB:
    tb_id = []

    smeg_id = []
    abs_id = []
    aro_id = []
    mari_id = []
    for rxn in smegmatis.reactions:
        smeg_id.append(rxn.id)
    for rxn in abscessus.reactions:
        abs_id.append(rxn.id)
    for rxn in aromaticivorans.reactions:
        aro_id.append(rxn.id)
    for rxn in marinum.reactions:
        mari_id.append(rxn.id)
    for rxn in Rv.reactions:
        tb_id.append(rxn.id)
    species = [smeg_id,abs_id,aro_id,mari_id,tb_id]

    tbper_ressemblance = pd.DataFrame({'smegmatis':(len(set.intersection(set(smeg_id),set(tb_id)))*100)/len(tb_id),
                                     'abscessus':(len(set.intersection(set(abs_id),set(tb_id)))*100)/len(tb_id),   
                                     'aromaticivorans':(len(set.intersection(set(aro_id),set(tb_id)))*100)/len(tb_id),   
                                     'marinum':(len(set.intersection(set(mari_id),set(tb_id)))*100)/len(tb_id)}, index =[0])
    
    per_ressemblance = pd.DataFrame({'smegmatis':(len(set.intersection(set(smeg_id),set(tb_id)))*100)/len(smeg_id),
                                     'abscessus':(len(set.intersection(set(abs_id),set(tb_id)))*100)/len(abs_id),   
                                     'aromaticivorans':(len(set.intersection(set(aro_id),set(tb_id)))*100)/len(aro_id),   
                                     'marinum':(len(set.intersection(set(mari_id),set(tb_id)))*100)/len(mari_id)}, index =[0])
    
if checkMissingMet:
    print('check Exchange Metabolites from data')
    def checkMet(met_missing_list,rxn_missing_list, rxn,model,universe):
        for idx_met, met in enumerate(universe.reactions.get_by_id(rxn).metabolites.keys()):

            try:
                if model.metabolites.get_by_id(str(met)):
                    #print(met, f'metabolite in {model}')
                    continue
            except KeyError:
                print(met, f'metabolite not in {model}')
                met_missing_list.append(met.id)
                m = 0
                for idx, x in enumerate(universe.metabolites.get_by_id(str(met)).reactions):
                    if x.id != rxn and idx == m:
                        m +=1
                        checkReaction(met_missing_list,rxn_missing_list, x.id, model,universe)
                    else:
                        continue
        return 
        
    def checkReaction(met_missing_list, rxn_missing_list,rxn, model, universe):
        try:
            if model.reactions.get_by_id(rxn):
                #print(rxn, f'Reaction in {model}')
                return 
        except KeyError:
            print(rxn, f'Reaction not in {model}')
            rxn_missing_list.append(rxn)
            checkMet(met_missing_list, rxn_missing_list, rxn, model, universe)
        return  

    rxn_search = ['EX_cit_e','CITtex','EX_ac_e','ACtex','EX_glu__L_e','EX_glyc_e','EX_glc__D_e','EX_glu__L_e','EX_acon_C_e','ACONCtex','EX_lac__L_e','L_LACD','L_LACD2','L_LACD3','L_LACtex','ACONTa','ACONTb','TAG','LIPY','TREH','GLCP', 'GLCS2']
    model = [smegmatis, abscessus,marinum,aromaticivorans,Rv]

    met_miss_abs = []
    met_miss_mari = []
    met_miss_smeg = []
    met_miss_Rv = []
    met_miss_aro = []
    rxn_miss_abs = []
    rxn_miss_mari = []
    rxn_miss_smeg = []
    rxn_miss_Rv = []
    rxn_miss_aro = []

    for i in range(len(rxn_search)):
        checkReaction(met_miss_abs,rxn_miss_abs, rxn_search[i],abscessus,universe)
    for i in range(len(rxn_search)):
        checkReaction(met_miss_mari,rxn_miss_mari, rxn_search[i],marinum,universe)
    for i in range(len(rxn_search)):
        checkReaction(met_miss_smeg,rxn_miss_smeg, rxn_search[i],smegmatis,universe)
    for i in range(len(rxn_search)):
        checkReaction(met_miss_Rv,rxn_miss_Rv, rxn_search[i],Rv,universe)
    for i in range(len(rxn_search)):
        checkReaction(met_miss_aro,rxn_miss_aro, rxn_search[i],aromaticivorans,universe)
#9 rxn
    print(set(met_miss_abs), '\n',
        set(met_miss_mari), '\n',
        set(met_miss_smeg),'\n',
        set(met_miss_Rv), '\n',
        set(met_miss_aro), '\n',
        (rxn_miss_abs), '\n',
        (rxn_miss_mari),'\n',
        (rxn_miss_smeg),'\n',
        (rxn_miss_Rv),'\n', 
        (rxn_miss_aro))


if modifyModels:

    for species in model_list:
        os.mkdir(directory_path+'/modified_'+str(species))
        os.chdir(directory_path+'/modified_'+str(species))
        model = species 
        out = 'modified_'+ str(species)+'.xml'

        if str(species) == selected_models['aromaticivorans']:
            reactions_to_add = [universe.reactions.get_by_id('ACONTa'),universe.reactions.get_by_id('ACONTb'), universe.reactions.get_by_id('EX_acon_C_e'), universe.reactions.get_by_id('ACONCtex'), universe.reactions.get_by_id('ACONCtupp')]
            reactions_to_rem = [model.reactions.get_by_id('ACONT')]
            met_to_add = [universe.metabolites.get_by_id('acon_C_c'),universe.metabolites.get_by_id('acon_C_p'),universe.metabolites.get_by_id('acon_C_e')]
        
        elif str(species) == selected_models['smegmatis']:
            reactions_to_add = [universe.reactions.get_by_id('EX_acon_C_e'), universe.reactions.get_by_id('ACONCtex'),universe.reactions.get_by_id('ACONCtex'), universe.reactions.get_by_id('ACONCtupp')]
            reactions_to_rem = [model.reactions.get_by_id('ACONT')]
            met_to_add = [universe.metabolites.get_by_id('acon_C_p'),universe.metabolites.get_by_id('acon_C_e')]

        elif str(species) == selected_models['abscessus']:
            #reactions_to_add = [universe.reactions.get_by_id('EX_acon_C_e'), universe.reactions.get_by_id('ACONCtex'), universe.reactions.get_by_id('ACONCtupp')]
            reactions_to_rem = [model.reactions.get_by_id('ACONT')]
            #met_to_add = [universe.metabolites.get_by_id('acon_C_p'),universe.metabolites.get_by_id('acon_C_e')]

        elif str(species) == selected_models['marinum']:
            reactions_to_add = [universe.reactions.get_by_id('EX_acon_C_e'), universe.reactions.get_by_id('ACONCtex'),universe.reactions.get_by_id('ACONCtex'), universe.reactions.get_by_id('ACONCtupp')]
            reactions_to_rem = [model.reactions.get_by_id('ACONT')]
            met_to_add = [universe.metabolites.get_by_id('acon_C_p'),universe.metabolites.get_by_id('acon_C_e')]

        elif str(species) == selected_models['H37Rv']:
            reactions_to_add = [universe.reactions.get_by_id('EX_acon_C_e'), universe.reactions.get_by_id('ACONCtex'), universe.reactions.get_by_id('ACONCtex'), universe.reactions.get_by_id('ACONCtupp')]
            reactions_to_rem = [model.reactions.get_by_id('ACONT')]
            met_to_add = [universe.metabolites.get_by_id('acon_C_p'),universe.metabolites.get_by_id('acon_C_e')]

        model.add_reactions(reactions_to_add)
        model.remove_reactions(reactions_to_rem)
        model.add_metabolites(met_to_add)
        cobra.io.write_sbml_model(model,out)
        os.chdir(directory_path)  
    
if addSubsystem:

    # Read the .rtf file
    with open("/home/icancino/Downloads/Subsystem-curation-aromaticivorans.rtf", "r", encoding="utf-8") as file:
    #with open("/home/icancino/Downloads/Subsystem-curation(2).rtf", "r", encoding="utf-8") as file:
        rtf_content = file.read()

    plain_text = rtf_to_text(rtf_content)

    print(plain_text)
    content = plain_text

    
    result_dict = {}
    for line in content.splitlines():
        if line.startswith('	*'):
        
            key = line[1:].split(':')[0][2:]
            if len(line[1:].split(':')) > 1:  
                values = list(line[1:].split(':')[1].split(','))
                values_clean = [str(v).strip() for v in values] 
                result_dict[key] = values_clean

    print(result_dict)
    bigg_ecoli_2 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iML1515.json')
    bigg_ecoli_4 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iAPECO1_1312.json')
    bigg_meth = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iAF692.json')
    bigg_pseud = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iJN1463.json')
    bigg_pseud_2 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iJN746.json')
    bigg_pneu = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iYL1228.json')
    bigg_pest = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iPC815.json')
    bigg_lacto = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iNF517.json')
    bigg_salmo =load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iYS1720.json')
    bigg_syn = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iSynCJ816.json')
    bigg_aur = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iYS854.json')
    bigg_aur2 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iSB619.json')
    bigg_bsub = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iYO844.json')
    bigg_scere = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iMM904.json')
    bigg_syn_2 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iJN678.json')
    bigg_abaum = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iCN718.json')
    bigg_iEK1008 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iEK1008.json')
    bigg_sson = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iSSON_1240.json')
    bigg_cdif = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iCN900.json')
    bigg_sson = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iSSON_1240.json')
    bigg_hpil = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iIT341.json')
    bigg_gmet = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iAF987.json')
    bigg_clju = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iHN637.json')
    bigg_selo = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iJB785.json')
    bigg_ecoli_5 = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/iAF1260.json')
    bigg_ecoli_core = load_json_model('/home/icancino/Documents/PhD/Reconstruction MTB analysis/e_coli_core.json')



    taylored_subsystem = {
                        'ABC transporters': ['MCBTFabcpp',  'FE3DCITabcpp',  'CU2abcpp',  'CU2t4pp',  'FE3DCIT2abcpp',  'FE3abcpp',  'FEENTERabcpp',  'CD2abcpp',  'ZN2abcpp',  'FEOXAMabcpp'],
                        'Alanine, aspartate and glutamate metabolism': ['GGDAPAH', 'GGSPMDS','GTHS','GG15DAPAH','GLNS_1', '3A2OA', 'GLUSx', 'ASPO5_2', 'ALAD_L2', 'ASPTA4', 'LASP2OA'],
                        'Amino sugar and nucleotide sugar metabolism': ['GALT','G1PACT',  'UAMRH',  'GDMANE',  'UAGCVT',  'UAGDP',  'GOFUCR',  'UDPGALM',  'UAPGR',  'UDPG4E',  'GF6PTA',  'PGAMT',  'UAG4Ei',  'UAGCVT_1',  'UAGDP_1',  'UDPACGLP'],
                        'Aminoacyl-tRNA biosynthesis':['GLUTRS','GLUTRS_3','GLUTRS_2'],
                        'Aminobenzoate degradation': ['ACOAT4','VNTDM', 'VNDH'],

                        'Arginine and proline metabolism': ['AGMPTRCtpp','ARGDI_1','PUTA3','HPROx','GGPTRCS','GGGABAH','GGPTRCO','GGSPMDO', 'GGDAPO', 'GGDAPS', 'GGDAPDxr', 'CBPAH', 'PRO1y','AGMDA','OCBT_1','AGMT',  'G5SD', 'G5SD2', 'P5CRx', '4ABUTD',  'NABTNO', 'ORNTA_1', 'APRTO2', 'AABTN', 'PRO1x', 'GG15DAPDxr', 'GGGABADxr', 'PTRCA2','DKMPPD','PHCHGS', 'PHCD', 'PY5CCR2',  'PTRCAT1',  'PY5CCR','GGGABADr'],
                        'Arginine biosynthesis': ['ACODA_1','OCBT2i','UREASE','AGPR',  'ACGS',  'ORNTAC',  'UREA',  'UREA_1',  'ACGK',  'ACOTA',  'ALPHNH','ARGSS2','OCBT_2', 'CBPS_1'],
                        'Arsenic resistance':['ASR2'],
                        'Ascorbate and aldarate metabolism' :['GLCRAL','GALCTD','D4DGCD', 'GLCRD_1','GLCRD','GLCRAL_1'], #'',"GLCRAL"
                        'Benzoate degradation': ['CMLDC','GLYBACT1r','SUCBZT1', 'SUCBZT2', 'BCOALIG', 'BCOALIG2','CO23OC', '3OXCOAT', 'OXOAEL', 'OP4ENH', 'PCADYOX', 'MUCCY_kt', 'HMSH', 'HOPNTAL',  '3OADPCOAT', 'ACALD',  '4CMLCL_kt',  'BZ12DOX',  'BZDIOLDH',  'DHACOAH' ],
                        'beta-Alanine metabolism': ['MCD','ABUTD','NBAHH_ir', 'MMSDHir'],
                        
                        'Biomass and maintenance functions':['Growth','BIOMASS_2'],
                        'Butanoate metabolism': ['MALDDH','PBUTT', 'BUTKr','HACD1','VACOAI','ACOAD1fr','SSCOARy', 'ACLS_a', 'MALEI',  'AACOAR_syn', 'BTDD_RR','ALCD2x','HXCT', 'ACACCT', 'ALCD4','ALCD4y'],
                        'Bacterial secretion system': ['CMCBTFtonex',  'CMCBTFtpp',  'CMCBTFtex',  'FE3PYOVD2tonex',  'FEOXAMtonex',  'FEENTERtonex',  'FEENTERtex'],
                        'Biosynthesis of amino sugars and nucleotide sugars': ['UAG2EMA','GALKr', 'UGLT', 'AGDC', 'ACM6PH', 'G6PDA', 'INOSTO'],
                        'Biosynthesis of siderophore group nonribosomal peptides': ['ASP3H','DHSKDH','MCBTS3','MCBYS2_1','MCBTS_1','PYOVDM1',  'PYOVDM2',  'FBACS',  'ORN5O'],
                        'Biosynthesis of steroids':['GGTT','HEXTT',''],
                        'Biosynthesis of unsaturated fatty acids': ['HATBH','DESAT18_1','DESAT16', 'HMR_0260',  'HMR_0188'],
                        'Biosynthesis of various alkaloids': ['FCOAHA'],
                        'Biotin metabolism' :['BTS3r', 'BTS_1'],
                        'Cellular community - prokaryotes': ['AI2tpp'],
                        'Citrate cycle':['OOR3r', 'SDH1'],
                        'Cysteine and methionine metabolism': ['METSR_S2','METSR_R2','METSR_R1', 'METS_R2', 'MHPGLUT2','MTAN','CYS', 'SAMU3','DM_scys__L_c', 'METSOXR1', 'METSOXR2', 'CYTOM', 'RHCYS', 'SHSL3', 'MDRPD', 'SLDx', 'UNK5', 'TRPAS1', 'TRPS3_1', 'ACSERHS', 'MTRK', 'CYSS_2', '5DOAN_1', 'SLDy'],
                        'Degradation of polycyclic aromatic hydrocarbon': ['SALCOD', 'BZSS'],
                        'D-Amino acids metabolism':['DAAD','ARGR'],
                        'Exchange': ['EX_araban__L_e', 'EX_madg_e', 'EX_2ameph_e, EX_fe3dhbzs3_e, EX_salchs4_e', 'EX_salchs4fe_e', 'EX_salchs2fe_e', 'EX_salchs2_e', 'EX_23dhbzs3_e'],
                        'Electrochemical potential-driven transport': ['MGt2pp', '4HOXPACt2pp', '4HBALDt2pp', 'Kt2pp', 'ASO3t2pp', 'CRO4t3pp', 'NH4tpp_1', 'THYMt4pp', 'CLt3_2pp', '3OXOADPt_pp','IDONt2rpp','GUAtpp', 'HYXNtpp', 'DXYLUDtpp', '5DGLCNt2rpp','NACt1pp',  '2DHGLCNkt_tpp',  '4HPTNtpp'],
                        'Fatty acid biosynthesis': ['ACCOAC',  'TDMS3',  'TMHAS3',  'TMHAS1',  'HDECH',  'FACOAL181_2',  'FACOAL80',  'FASm2402',  'FASm2002',  'FAS161',  'FAS160',  'FAS200',  'FAS180',  'FACOAL160',  'FAS260',  'FACOAL180',  'MCOATA',  'TDMS1',  'FASm180',  'FAS140',  'TDMS4',  'FACOAL161',  'FAS181',  'FACOAL180t2pp',  'TMHAS4',  'FAS240_L',  'FAS120',  'TDMS2',  'FAS80_L',  'FASm2202',  'PHDCATA',  'FACOAL140_1',  'FACOAL200',  'FAS100',  'TMHAS2',  'FA120ACPHi',  'FACOAL50i',  'FACOAL80t2pp',  'FACOAL140t2pp',  'FACOALP90t2pp',  'FACOAL90i',  'FACOAL70t2pp',  'FACOALP60t2pp',  'FACOAL120t2pp',  'FACOALP100t2pp',  'FACOAL181d6t2pp',  'FACOAL50t2pp',  'FACOAL160t2pp',  'FACOAL141t2pp',  'FA141ACPHi',  'FACOAL1821',  'FACOAL90t2pp',  'FACOAL141',  'FACOAL60i',  'FACOAL182t2pp',  'FACOAL100t2pp',  'FACOAL60t2pp',  'FACOAL_70_c',  'FACOAL181t2pp',  'FACOAL1812',  'FA80ACPHi',  'FACOAL161t2pp',  'FASC141ACP',  'FACOAL181d11tpp','KAS17', 'KAS7', 'KAS8', 'KAS13', 'KAS2','FASm240',  'FASm2802','FACOALPHDCA','FASm220','FAMPL4','FAMPL2','FASm1601','FASm300','FASm280','FASm2401','FASPHDCA','FASm1801','FASm2801','FASm320','FASm260','FASm2602','BIRA','FASm340','ACCD3','FASm2601','FASm2201','FASm2001', 'AACPS4'],
                        'Fatty acid elongation':['ECOAH5_1','ACChex_1', 'ACOAD7','HACD7i',  'ACACT8r',  'FACOAE140',  'AACPS10',  'HACD10i',  'FACOAE160',  'HACD23i',  'HACD6i',  'HACD3i',  'HACD8i',  'HACD11i', 'HACD19i',  'FACOAE180',  'HACD5i',  'HACD1i',  'HACD4i',  'FACOAE60',  'HACD26i', 'CTECOAI7', 'ACACT12',  'FACOAE120',  'HACD29i', 'HACD30i',  'ACACT10',  'FACOAE70',  'FACOAE90',  'FACOAE50', 'FACOAE1829Z12Z', 'FACOAE100',  'ACACT13',  'FACOAE161', 'FACOAE141',  'ACACT9',  'ACACT11',  'ACACT8',  'FACOAE181',  'FACOAE80', 'ACACT2r', 'ACACT4r', 'ACACT3r', 'ACACT5r', 'ACACT7r', 'arachACP', 'ACACT6r', 'HACD5', 'HACD8', 'HACD7', 'HACD3', 'HACD2', 'HACD4', 'HACD6'],
                        'Fatty acid degradation': ['FACOALP70t2pp','FACOAL50It2pp', 'FACOALP50t2pp','FACOAL40t2pp','FACOALP80t2pp','FACOALT60t2pp','KAT19', 'ECOAH24', 'ACOAD21f', 'ACOAD20f','AACPS11', 'AACPS3', 'AACPS2', 'AACPS9', 'AACPS7','ACOAD3', 'ACOAD6', 'ACOAD2','ECOAH27',  'HACD32i',  '24DECOAR',  'ECOAH15',  'FAD_19',  'FAD_11',  'FAD_9',  'ECOAH31',  'FAO11',  'ECOAH36',  'ECOAH2',  'FAD_18',  'VCOAD',  'HACD31i',  'HACD22i',  'ECOAH5',  'ECOAH4',  'ECOAH8',  'ECOAH7',  'FAD_10',  'HVCD',  'ACOAD6f',  'ACOAD4f',  'FAO10',  'ECOAH1',  'ECOAH37',  'FAD_12',  'FAD_14',  'FAD_8',  'ECOAH23',  'FAD_17',  'FAD_13',  'FAD_27',  'FAD_15',  'FAO3',  'FAO1',  'VECOAH',  'FAO2',  'ACOAD8f',  'FAD_16',  'ECOAH28',  'ECOAH14', 'ECOAH3',  'ACOAD5f',  'KAT31',  'ECOAH33',  'ACOAD32f',  'ECOAH21',  'ACOAD16f',  'ACOAD27f',  'KAT26', 'ACOAD30f',  'KAT1',  'CTI1',  '23CTI1',  'ACOAD19f',  'KAT27',  'HBCH',  'HACD28i',  'KAT14',  'KAT20',  'HACD25i',  'ACOAD10f' , 'ACOAD15f',  'KAT17',  'ECOAH25',  'ECOAH19',  'HACD20i',  'ACOAD11f',  'HACD15i',  'ACOAD22f',  'HACD14i',  'KAT18',  'KAT10',  'ACOAD3f',  'HACD12i',  'KAT29',  'ACOAD31f',  'KAT28',  'KAT33',  'ACOAD28f',  'KAT23',  'ACOAD29f',  'ACOAD18f','ACOAD26f',  'ECOAH16',  'ACOAD23f',  '23CTI3',  'ACOAD33f',  'HACD21i',  'ECOAH34',  'ACOAD24f',  'ECOAH22',  'KAT11',  'KAT25',  'KAT12',  'ECOAH18',  'KAT16',  'KAT30',  'HACD18i',  'KAT15',  'ECOAH17',  'ACOAD7f',  'ECOAH26',  'ACOAD2f',  'ACOAD25f', 'ACOAD34f',  'ECOAH32',  'ECOAH20',  'ACOAD1f',  'ACOAD12f',  'HACD33i',  'CTECOAI6',  'KAT32',  'KAT13',  'HACD27i',  'ECOAH35',  'HACD17i',  'ECOAH30',  '3HPTNCOAR',  'KAT21',  'HACD24i',  'ACOAD13f',  'ACOAD17f',  'ECOAH29',  'CTECOAI8',  'KAT24',  'ACOAD14f',  'KAT22',  'ECOAH38','HACD16i',  '23CTI2',  'PLCD', 'ACOAD4_1', 'ACOAD5_1', 'HACD1_1', 'ACACT5r_1', 'ACACT6r_1', 'ECOAH1_1'],
                        'Fatty Acid Metabolism':['HACD2i','ECOAH6','FA161ACPHi','ARACHTA'],
                        'Methane Metabolism' :['ALDD1'],
                        'Fructose and mannose metabolism':['GAPP', 'ALKP', 'HEX4', 'FRUK', 'HEX7','FCLK', 'MN6PP', 'FCI', 'RMI', 'RMPA', 'F1PP', 'F6PP', 'RMK','MANPGH'],
                        'Glycerophospholipid Metabolism': ['DASYN_HP','CDAPPA_HP','G3PD7','G3PD2','PGPP160','ETHAAL','CHOLK', 'PLIPA1E180', 'PLIPA1E160','GPDDA3_1', 'GPDDA1_1', 'GPDDA2_1', 'GLYCK_1', 'MI1PP', 'PGPP160190', 'PGPP190'],
                        'Glycolipid biosynthesis': ['MAS3',  'PKS2','TRESULT', 'TATS','PATS'],
                        'Glycerolipid metabolism': ['13PPDH','G3PT', 'DHAPT','GLYCDH','GLYK1','DHAK','GLYK','ALCD19'],
                        'Glutathione metabolism': ['AMPTASEGGLN', 'AMPTASEGS',  'AMPTASEAT',  'AMPTASEGGLU',  'AMPTASEATRP',  'AMPTASEAH', 'AMPTASEAL', 'AMPTASEGM', 'AMPTASEHH', 'AMPTASEHG',  'GLYPHEHYc', 'LEULEULAPc', 'GGTAe2','OPAH', 'GG15DAPO', 'GG15DAPS','GLUCYS', 'GTHPe','GTHOr', 'GRXR',  'GTHPi','GGLUCT2'],
                        'Glycogen metabolism': ['GLBRAN2','GLCP','GLDBRAN2'],
                        'Galactose metabolism': ['DDGALK','T6PK','TAG1PK','UT6PT', 'CT6PT','MLTG1', 'MLTG2', 'MLTG3', 'MLTG4', 'MLTG5', 'GALCTND', 'DDPGALA', 'TGBPA', 'PFK_2', 'GLTPD'],
                        'Glyoxylate and dicarboxylate metabolism': ['GLUSfx','HPYRI','TRSARr','FORMCOAL',  'FCOAH2', 'PRDX', 'FDH5pp', 'LCARS','GLYCL_2','LCARSyi','ACOXT','FOMETRi','CAT', 'CATpp','XCDC', 'GLYCTO2', 'GLYCTO3', 'GLYCTO4', 'TARTD', 'DTARTD', 'FORCT', 'HPYRRy', 'LCARR', 'GOR1'], #ask if its dicarboxylate or dicarbonylate
                        'Glycolysis/Gluconeogenesis':['PGCM','GAPD_1','GalMr','PYK5','GNK'],
                        'Glycine, serine and threonine metabolism': ['DMGDH2','HSK_2', 'CYSGL', 'ASAD_1','BETALDHx',  'CHOLD',  'GBDM',  'GLYBCOAT',  'DMMGDH2',  'SARCOX',  'BETALDHy',  'CRNDH','AACTOOR','GLYDHDA_copy1', 'GLYDHDA_copy2', 'THRS_2', 'THRD_L_1', 'SER_AL'],
                        'Histidine metabolism':['IZPN_1', 'HISTP_1', 'FGLU_1','DM_ergoth_c', 'IG3PS_1', 'EGTD', 'HISTDa_1', 'HISTDb_1','PRATPP_1','LHISO', 'LHIST', 'IGPDH_1'],
                        'Iron binding by siderophores': ['FE3PYOVDL2','FETRANS', 'CMCBTFU','SALCHS2FEexs','SALCHS4FEexs', 'DHBSZ3FEexs'],
                        'Inositol phosphate metabolism': ['INS2D','INOSR','INSCR','MI1PP_1'], 
                        'Levulinate metabolism':['4POPCOAM','4OXPTCOADHy',	'4POPCOAM3','4OXPTCOADHx', 'LEVUACOAL','4HPTNCOAL','4HPTNCOAK','4POPCOAM2'],
                        'Lipopolysaccharide biosynthesis':['DM_gmhep1p_c','GMHEPPA', 'PPM1', 'S7PI', 'GMHEPK'],
                        'Lipoarabinomannan and lipomannan and lipoarabinomanan biosynthesis': ['GMT1','G16MTM8',  'MANAT2',  'G16MTM1',  'G16MTM9',  'MANAT3',  'G16MTM6',  'G16MTM7',  'G16MTM4',  'G16MTM5',  'MANAT1',  'G16MTM10'],
                        'Lipid transporters': ['2AGPE180tipp', '2AGPG161tipp', '2AGPG181tipp', '2AGPE160tipp', '2AGPG141tipp', '2AGPG180tipp', '2AGPG140tipp', '2AGPG120tipp'],
                        'Lysine degradation':['PTRCTA','BUTCT','BUTCT2', 'ACOAD1','AASAD3','AATA','APENTAMAH2','APTNAT','LYSMO','OXPTNDH'],
                        'Lysine biosynthesis':['THPAT','APATi', 'APTA1i', 'ADPDA','ADAPAT','DAPDA','DHDPRx_r'],	 
                        'Miscellaneous and Unassigned':['AALDH','GLXCBL','ALDD5','ALDD31_1','NADS1_1','HOXPRx','ALDD3','TARTDC','PPDOy', 'G2PPpp','LACZpp','NOX','CU1Opp','PEAOpp','CYO4pp','DNTAep','MDNTHAep','FNOR','ACTD2','LACZ'],
                        'Miscellaneous and spontaneous reactions':['CMCBTFexs','FE3DCITCH','ASP3H','MCBTFUexs','DATPHs','ALOX','FNOR','H2CO3D','H2CO3D2','METOX1s','METOX2s', 'sink_mobd_c', 'sink_hemeO_c', 'HIPEOX',  'NH3c', 'sink_2ohph_c', 'sink_sheme_c','sink_aacald_c'],
                        'Mycolic acid biosynthesis': ['ACCC','FAMPL5', 'FAMPL1', 'MYC1CYC5', 'MYCOE', 'MYC2CYC2', 'MYCSacp56', 'MYCON3', 'ACCA3', 'OTSB1', 'MYCON5', 'MYC1CYC2', 'MYCSacp58', 'MYC2CYC1', 'MYCON2', 'MYCON4', 'FAMPL3', 'MYC1M1', 'MYC2CYC3', 'MYC1M2', 'TMMYCO', 'MYCOPP', 'MYC1CYC4', 'MYCON1', 'FBPA', 'MYCSacp50', 'MYC1CYC1', 'MYC1CYC3'],
                        'Mycothiol biosynthesis': ['NAGINS',  'MYCTR',  'MYCTR2',  'CIGAMS'],

                        'Methane metabolism': ['F4RHi','DMHDRFS','DMSOR1pp', 'TMAOR1pp','ALCD1'],


                        'Nicotinate and nicotinamide metabolism': ['NAMNPP','NADS2','NAPRT', 'NNAM', 'NMNS', 'ASPO2', 'ASPO2y','NADN','NADDPp_1','NADTRHD','PYNP1','NMNAT_1', 'NADDP_1','PYRDOX', 'NACHY', 'PYRZAM', '6HNACMO',  'NFMLDF','NADHXD', 'NADPHXE', 'NADPHXD', 'NADPHHS', 'NADHPO', 'NADHHR', 'NADPHHR', 'NADPHQR2', 'NADHHS', 'NADHXE'],
                        'Nitrotoluene degradation': ['TNTR1', 'TNTR1x', 'TNTR2', 'TNTR2x', 'TNTR3', 'TNTR4', 'TNTHAC1', 'TNTHAC2'],
                        'Nitrogen metabolism': ['NTRIR2x','CBMD','NODOy', 'NODOx', 'NO3R2bpp', 'NO3R1bpp','NH4OHDs','NTRARz','HCO3E', 'NOR_syn1', 'CYNTAH_1','NMO', 'NO', 'NTRIRy', 'NO3R1pp'],                        
                        'One carbon pool by folate': ['MTHFC', 'MTHFD', 'FTHFD', 'MTHFRD', 'FTHFCL', 'THFAT','MTHFR','DHFR2i', 'MTHFC_1', 'MTHFR3'],
                        'Other Amino Acid Metabolism':['PRAMPC_1','METSR_S1','CHORS_1','HSTPTr','ORNTAC_1'],
                        'Oxidative phosphorylation': ['SUCD4','NADH10',  'QRr',  'ATPS4rpp',  'QMO2',  'AMMQT8_2',  'NADH9',  'NADH5',  'CYO1b',  'NADH2r',  'NADH18pp',  'TMAOR2pp',  'CYTBO3_4pp',  'NADPHQR3',  'CYTCAA3pp',  'CYTBD2pp',  'THD2pp'],
                        'Other carbon fixation pathways':['COCO2'],
                        'Other glycan degradation (aka murein recycling)': ['ALAALAD','MLTGY3pp',  'MLTGY1pp', 'AGMH',  'MDDEP1pp',  'UM4PL',  'MLTGY2pp',  'AGM4PApp',  'UM3PL',  'MDDEP2pp', 'ANHMK'],

                        'Porphyrin metabolism': ['HMBS_1','HEMEOS_1','FCLT_2','LHTRK','UPP3MT_3','FEROc','FEROpp','ALMPC', '5HBZIDS2',  '5HBZIDMT', 'CPPPGO2_1'],
                        'Phosphonate and phosphinate metabolism': ['PHHL', 'PHHL2', 'AEPPYRTA','PPTHpp'],
                        'Pyruvate metabolism': ['MALLAC','PYROX_1','POR2','ACYP','ACYP_3','FHL','L_LACD3', 'POX', 'ACS_1','PEPC','ALDD2x','ALDD2y', 'ALDD3y', 'PDC', 'ACLSb', 'L_LACDcm_1', 'ALCD2y'],
                        'Purine metabolism': ['GNNUC','RNTR3c', 'NTP10', 'ADD', 'XTSNH', 'NTP11', 'RNTR1c', 'NDP3', 'CBMK4','BUP','ADNUC','CBMKr2', 'NTP6', 'ADPRDP_1', 'NTP2' ,'RNTR2c', 'NDPK10','ADNK3',  'PUNP3',  'NTP3',  'ADNK4', 'RNDR2b',  'RNDR1b',  'PUNP4','PUNP5',  'PUNP6',  'INSK',  'PUNP2',  'PUNP7',  'GSNK',  'URIK2',  'UUPP',  'PYK6',  'NTPP10',  'NTPP1_1',  'ALLTAHr',  'PNP',  '3NTD7pp',  'LDGUNPD',  'ALLTN',  'NTPP9_1',  'NTPP2_1',  'NTPP9',  'ADK2_1',  'NTPP10_1',  'NTPP11',  'PYK3',  'RNTR2c',  '3NTD9pp',  'UGLYCH',  'ADKd',  'DGUNC','CPK1', 'CYTK2_1','R5PPpp', 'CDGUNPD','DM_psd5p_c', 'URIC', 'AIRC4', 'DADNK', 'DGNSK', 'RNTR2', 'RNTR1', 'AIRCra', 'FGFT', 'AIRCrb','NDPK8','NTD11','NTD10','ADPT','GMPS2','NTD6','NTD2','NTD2pp','NTD7','NTD8','NTD9','NTD9pp','NTD12','NTD7pp','ADK1','ADK2','ADK3','ADK4','ADNK1','GUAD','ADSS','GTPH1','NTPTP1','NTPP1','NTPP3','NTPP5','NTPP6','NTPP7','GDPDPK','PDE4','ADA','NDPK1','RNDR2','GUACYC','PRFGS','PPGPPDP','RNDR1','DGK1','NDPK9','PRAIS','GK1','PRASCSi','HXPRT','GLUPRT','PDE1','NTPTP2','ADSL2r','AIRCr','XPPT','IMPD','PRAGSr','INSH','DADA','DADK','AICART','GTPDPDP','ADSL1r','GARFT','AIRC1','IMPC','GTPDPK','GMPS','GUAPRT','NDPK5','ADNCYC','r0280','HXAND','AP4AS','AIRC3','AIRC2','GDTPR','GDPDPKE','ATAH_1','GTPHS','ADPRDP','GK2','PRFGCL','XAND','GART','AIRSK','AP4AH','ATPHs','PPGPPDPE','GMPR','PPRGL'],
	                    'Pyrimidine metabolism': ['3AMACHYD', 'POAACR', 'DURADx','PYROX','URACPAH','URIH','CYTDH','NTP9','DHORD6','NTP8','PYK6','RNDR3b',  'RNDR4b', 'NTP5', 'TMDK1', 'CYTKD2', 'DCYTD', '3NTD4pp', 'MCSNAH', 'PYNP2r_1', 'PYK2', 'UPPN', 'DHPM2', 'PYNP1_1', '3NTD2pp', 'DHORD7', 'PYK4', 'DHPM1', 'DURAD2', 'ASPCT_2', 'DURIPP_1', 'NTPP8_1', 'NTPP11_1', 'DURAD','DHORD3', 'RNTR3', 'CTPS1_2','CYTK1','CYTK2','NDPK2','NDPK3','URIDK2r','BARB','DHORTS','NTD1','NTD5','NTD3','NTD4','NTD4pp','DCTPD2','PYNP2r','DCTPD','NTPP2','NTPP8','CBPS','RNDR3','ASPCT','DHORDi','DUTPDP','UMPK','NTPP4','PUNP1','RNDR4','CTPS2','RNTR4','NDPK7','TMDPP','UPPRT','DHORD2','ORPT','DURIPP','NDPK4','CTPS1','TMDS','NDPK6','TMDS3','OMPDC','DTMPK','YUMPS','CYTD','PSUDS','DHPD','BUP2','CSND','DHR_1','DHORDful','BUPN','DHORD5'],
                        'Pentose and glucuronate interconversions':['5DKGR','ARABDI','DXYLK','5DGLCNR','XYLK', 'XYLI1', 'XYLI2', 'RBP4E', 'ARAI', 'RBK_L1', 'DXYLTD', 'ABFA', 'ALR3'],
                        'Phenylalanine metabolism':['DLYSAD','araphe1', 'araphe3', 'araphe2','PEAMNOpp','ALDD19xr', 'HADPCOADH3', 'COALCDH', 'PACCOAE',  'OXDHCOAT', 'REPHACCOAI', 'OXCOAHDH',  'PPYRDC','PACCOAL','DHCIND', 'DHCINDO', 'HPPPNDO', '3HPPPNH', '3HCINNMH', 'CINNDO', 'HKNDDH','AMID2_1', 'HPHACA', 'PHEOR'],
                        'Pentose phosphate pathway':['GNP','DRBK_1','RBK2','PPM2','GNK_1','2DGLCNRx', '2DGLCNRy', 'RU5PP','R5PP','PFK_4','2DGLCNRx','2DGLCNRy','RU5PP','R5PP','IDOND', '2DGULRGx', '2DGULRx', '2DGULRy', '2DGULRGy', 'DKGLCNR1', 'DKGLCNR2x', 'DKGLCNR2y'],
                        'Propanoate metabolism':['13PPDH2GLCRAL','POR3','13PPDH2','13PPDH2','ALDD4x','ACS2','PTA2','HACOADr','MMCD', 'PPCSCT', '13PPDH2_1', 'PCNO', 'ACCOAL'],
                        'Propanediol Catabolism':['PPDD','12PPDRDH'],
                        'Phenylalanine, tyrosine and tryptophan biosynthesis' : ['DHQS_2', 'SHKK_1', 'PTHRP','PPND2','PACCOAL3','PACCOAL2','TRPS2_1', 'ANSN', 'F4H2O','ALDD20x'],
                        'Peptidoglycan biosynthesis': ['MCTP1Bpp', 'MLDCP3App','4PCPpp', 'AGM4PCPpp','UAAGDS', 'UAMAS',  'UAMAGS',  'UGMDDS',  'UDCPDP',  'ACGAMT',  'UAGPT3',  'MDDCP5pp',  'AGM3PApp',  'UM4PCP',  'MCTP1App',  'MDDCP3pp',  'MLDCP1App',  'MPTG',  'MLDCP2Bpp',  'MLDCP1Bpp',  'MLDCP2App','UGLDDS2_2', 'UGLDDS1_1','PAPPT3',  'ALAALAr', 'UDCPDPex',  'MDDCP4pp',  'MDDCP1pp'],
                        'Phenylpropanoid biosynthesis': ['FERULCOAS','CINNMCOAH','COALDDH'],
                        'Phthiocerol biosynthesis': ['PHTHS', 'PREPTHS', 'PREPHACPH', 'PHTHDLS', 'PREPTTA', 'FACOALPREPH', 'PDIMAS'],
                        'Phosphotransferase system (PTS)': ['TREptspp', 'TAGptspp', 'FRUpts2pp', 'MALTptspp'],
                        'Pantothenate and CoA Metabolism': ['ACPS','PPNCL2_1', 'DPCOAK_1','PNCDC_1'],
                        'Porines and efflux systems': ['GLYMETtex', 'LEULEUtex', 'METDtex', 'NACtex', 'G3PCtex', 'ALATHRtex', '34dhpactex', '2DHGLCNtex', 'TOLtex', 'DOPAtex', 'GLYGLUtex', 'HPAtex', 'GLYPHEtex', 'THYMtex', 'ALAHIStex', 'HISGLYtex', 'DXYLUDtex', '4HPTNtex', 'CUtex', '5DGLCNtex', 'G3PStex', 'G3PItex', 'CGLYtex', 'PEAMNtex', 'GLYGLNtex', 'HAT40tex', 'INSTtex', 'TAGtex', 'IDONtex', 'ALATRPtex', 'GLYSERtex', 'HISHIStex', 'TOLt6', 'FALDtpp2', 'TOLt5', 'INDOLEt5','ASO3t8pp','DNTAep','MDNTHAep', 'NI2t3pp', 'CU2tex', 'PHEMEtiex', 'COBALTt4pp', 'ZN2t3pp', 'FEOXAMUtpp', 'FEOXAMUtex', 'COBALT2t3pp', 'MN2t3pp', 'CD2t4pp', 'FEENTERtpp', 'CD2t3pp', 'ZN2t4pp', 'FE2t3pp', 'PYOVD2tex', 'FE2tex','MLDEP2pp', 'MLDEP1pp'],
                        'Riboflavin metabolism': ['DHPPDA_2','FMNRx', 'FLVRx','FMNRy','FMNAT_1', 'RBFK_1','2MAHMP'],
                        'Ribose metabolism':['PPM','RBK', 'PRPPS'],
                        'Signal transduction / Two-component system': ['MALt5','ALAALAD_1'],   

                        'Styrene degradation': ['FUMAC', 'MACACI', 'HGNTOR'],

                        'Starch and Sucrose metabolism':['PGMT_B','AMALT1','TRE6PH','TRE6PP','TRE6PS','GALUi','GLGC','MTI','GLCGSD','MALT', 'MLTP1', 'MLTP2', 'MLTP3', 'MMCD', 'TREP', 'TREH', 'AMALT2', 'AMALT3', 'AMALT4', 'MOTS2', 'MOTS3', 'MOTH2', 'MOTH3','SUCR', 'MAL6PG', 'MOTS4', 'MOTH4', 'BG_MAGD','BG_MADG'],
                        'Steroid degradation': ['ECHA20_2', 'MBKTHIO', 'FAD6', 'FAD31', 'IPDF', 'IPDCF'],
                        'Selenocompound metabolism':['CYSLY3','THIORDXi','TRDR',	  'ACSERL',	 'SEAHCYSHYD', 'SELCYSLY3', 'SEAHCYSHYD_1'],
                        'Sugar Metabolism': ['MMSAD3'],               
                        'Sugar transport': ['LCTSt3ipp'],
         
                        'Sulfur metabolism': ['PAPSR2','CYANSTpp','ALDD4','HLR_6515','SIRA2','SULR_1'],
                        'Taurine and hypotaurine metabolism':['HPYRRx','GTMLT', 'DALAOX','TAUDO'],
                        'Transport via diffusion': ['PHEHPAtexi','6ATHAtexi', 'PHEPTtexi', 'PHEOCTAtexi', 'NONAtexi','3OXOADPt_ex','HG2tppi', 'HAT40tpp'],

                        'Transport':['CYTMQOR3pp','CYOO2pp','FDH5r2pp','FRD2rpp'],
                        'Thiamine metabolism': ['TMK','HMPK1_1','THZSN', 'GOX', 'AHMMPS_1', 'AHMMPDPT', 'AHMMPS2'],
                        'Toluene degradation': ['PCMEH3pp','VNDH_2', 'BZDH', 'BZALDH','BZCSD', 'HMPSCD', 'BSCT', 'PICH', 'BZSCT','sink_4crsol_c'],                        
                        'Tuberculosis': ['THIORDXi','SPODM','DM_tbsinol_c','TBSYNLDPS',  'DM_isotbsinol_13S_c',  'TBSNLS3',  'TBSNLS1',  'TBSNLS2',  'DM_isotbsinol_13R_c','THIORDXi', 'SPODM'],
                        'Tryptophan metabolism': ['TRPTA','TRYPTAOX', 'LTDCL'],

                        'Tyrosine metabolism': ['2HH24DDH1','aratyr2','aratry2',  'aratry1', 'aratyr1', 'aratyr4', 'aratyr3','42A12BOOXpp','OPETDC', '3HPAOX',  '4H2KPILY', 'OPTCCL', '2HH24DDH', 'CMHMI', '4HOXPACMOF', 'CMCMSAD', 'DHPDO','DHEDAA','HHDDI', 'OHEDH'],                        
                        #'Benzanoate degradation':['CMLDC'], #check with Benzoate
                        'Terpenoid backbone biosynthesis': ['CDPMEK_1','DPMVD', 'PMEVK', 'HMGCOAR','MECDPDH3_syn','PPTT', 'OPP', 'UDPDPS2', 'DPPS', 'HMEDS', 'DCPDPS',  'IPDPUDCT', 'MECDPDH', 'HMEDR', 'UDPDPS', 'FRTT'],
                        'Ubiquinone and other terpenoid-quinone biosynthesis': ['FADMQOR','NPHSa', 'NPHSb','DHNCOAT','SUCBZL_1','DHNCOAT', 'DHNCOAS','VCACT','TYRTA', '34HPPOR'],
                        'Vitamin B6 metabolism':['PDXPP','HPYRP','4HTHRS_1', 'PLPS_1', 'ATPX5P'],
                        'Vitamins and Cofactor Biosynthesis': ['CHRPL','THZPSN','FE3Ri'],
                        'Valine, leucine and isoleucine (degradation)': ['3AIBT2','VOR','AT_MBD2','AT_MBD2','ILEDHr','OIVD2E', 'MB2CFO','MMM', 'MMSAD1', 'HIBADH','IBTMr', 'KARI_1', 'OCOAT2r', 'MMM2r', 'HMGCOAS','3M2OPLOXRD','BKDA2','BKDC1','BKDC','ACOAD20'],
                        'Valine, leucine and isoleucine (biosynthesis)':['DHAD3','KARA4', '2ACLMM', 'ACLS_b', 'ACLSb_1', 'MME_1',  'KARI_23dhmp_1', 'ACHBS2'],
                        'Valine, Leucine, and Isoleucine Metabolism':['HIBD','FADE233','ACOADH1'],
                        'Xylene degradation': ['HKNTDH', '2HH24DDH1','HOPNTAL2_1'],

                        
                        }


    for idx, model in enumerate(model_list): 
        for rxn in model.reactions:
                
            if rxn in bigg_meth.reactions:
                if bigg_meth.reactions.get_by_id(rxn.id).subsystem and bigg_meth.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_meth.reactions.get_by_id(rxn.id).subsystem   
               
            if rxn in bigg_salmo.reactions:
                if bigg_salmo.reactions.get_by_id(rxn.id).subsystem and bigg_salmo.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_salmo.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_lacto.reactions:
                if bigg_lacto.reactions.get_by_id(rxn.id).subsystem and bigg_lacto.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_lacto.reactions.get_by_id(rxn.id).subsystem 
                
            if rxn in bigg_pest.reactions:
                if bigg_pest.reactions.get_by_id(rxn.id).subsystem and bigg_pest.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_pest.reactions.get_by_id(rxn.id).subsystem                   
            
            if rxn in bigg_bsub.reactions:
                if bigg_bsub.reactions.get_by_id(rxn.id).subsystem and bigg_bsub.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_bsub.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_aur.reactions:
                if bigg_aur.reactions.get_by_id(rxn.id).subsystem and bigg_aur.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_aur.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_aur2.reactions:
                if bigg_aur2.reactions.get_by_id(rxn.id).subsystem and bigg_aur2.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_aur2.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_scere.reactions:
                if bigg_scere.reactions.get_by_id(rxn.id).subsystem and bigg_scere.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_scere.reactions.get_by_id(rxn.id).subsystem       
                
            if rxn in bigg_pneu.reactions:
                if bigg_pneu.reactions.get_by_id(rxn.id).subsystem and bigg_pneu.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_pneu.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_sson.reactions:
                if bigg_sson.reactions.get_by_id(rxn.id).subsystem and bigg_sson.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_sson.reactions.get_by_id(rxn.id).subsystem   
            
            if rxn in bigg_cdif.reactions:
                if bigg_cdif.reactions.get_by_id(rxn.id).subsystem and bigg_cdif.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_cdif.reactions.get_by_id(rxn.id).subsystem 
            
            if rxn in bigg_hpil.reactions:
                if bigg_hpil.reactions.get_by_id(rxn.id).subsystem and bigg_hpil.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_hpil.reactions.get_by_id(rxn.id).subsystem   
            
            if rxn in bigg_pseud_2.reactions:
                if bigg_pseud_2.reactions.get_by_id(rxn.id).subsystem and bigg_pseud_2.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_pseud_2.reactions.get_by_id(rxn.id).subsystem       
                
            if rxn in bigg_gmet.reactions:
                if bigg_gmet.reactions.get_by_id(rxn.id).subsystem and bigg_gmet.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_gmet.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_clju.reactions:
                if bigg_clju.reactions.get_by_id(rxn.id).subsystem and bigg_clju.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_clju.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_syn_2.reactions:
                if bigg_syn_2.reactions.get_by_id(rxn.id).subsystem and bigg_syn_2.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_syn_2.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_selo.reactions:
                if bigg_selo.reactions.get_by_id(rxn.id).subsystem and bigg_selo.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_selo.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_syn.reactions:
                if bigg_syn.reactions.get_by_id(rxn.id).subsystem and bigg_syn.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_syn.reactions.get_by_id(rxn.id).subsystem   
                
            if rxn in bigg_ecoli_4.reactions:
                if bigg_ecoli_4.reactions.get_by_id(rxn.id).subsystem and bigg_ecoli_4.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_ecoli_4.reactions.get_by_id(rxn.id).subsystem  
            if rxn in bigg_pseud.reactions:
                if bigg_pseud.reactions.get_by_id(rxn.id).subsystem and bigg_pseud.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_pseud.reactions.get_by_id(rxn.id).subsystem   
            if rxn in bigg_abaum.reactions:
                if bigg_abaum.reactions.get_by_id(rxn.id).subsystem and bigg_abaum.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_abaum.reactions.get_by_id(rxn.id).subsystem   
            
            if rxn in bigg_ecoli_core.reactions:
                if bigg_ecoli_core.reactions.get_by_id(rxn.id).subsystem and bigg_ecoli_core.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_ecoli_core.reactions.get_by_id(rxn.id).subsystem   
                    
            if rxn in bigg_ecoli_2.reactions:
                if bigg_ecoli_2.reactions.get_by_id(rxn.id).subsystem and bigg_ecoli_2.reactions.get_by_id(rxn.id).subsystem != '':
                #print(bigg_ecoli_2.reactions.get_by_id(rxn.id).subsystem)
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_ecoli_2.reactions.get_by_id(rxn.id).subsystem
            
            if rxn in bigg_ecoli_5.reactions:
                if bigg_ecoli_5.reactions.get_by_id(rxn.id).subsystem and bigg_ecoli_5.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_ecoli_5.reactions.get_by_id(rxn.id).subsystem
            
            if rxn in Rv.reactions:
                if Rv.reactions.get_by_id(rxn.id).subsystem and Rv.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = Rv.reactions.get_by_id(rxn.id).subsystem
            
            if rxn in bigg_iEK1008.reactions:
                if bigg_iEK1008.reactions.get_by_id(rxn.id).subsystem and bigg_iEK1008.reactions.get_by_id(rxn.id).subsystem != '':
                    model.reactions.get_by_id(rxn.id).subsystem = bigg_iEK1008.reactions.get_by_id(rxn.id).subsystem

            
            if str(rxn.id).__contains__('tex') or str(rxn.id).__contains__('t2pp') or str(rxn.id).__contains__('t3pp') or str(rxn.id).__contains__('t4pp'):
                model.reactions.get_by_id(rxn.id).subsystem = 'Transport'
            if str(rxn.id).__contains__('abc'):
                model.reactions.get_by_id(rxn.id).subsystem = 'ABC transporters'

            if rxn.subsystem == '' or rxn.subsystem == None:
                rxn.subsystem = 'NaN'
       
    ### Adding taylored subsystem association
    
            for subsystem, reaction_list in taylored_subsystem.items():
                if rxn.id in reaction_list:
                    rxn.subsystem = subsystem
                    #print(rxn.id, rxn.subsystem, subsystem)
        count = 0        
        for rxn in model.reactions:
            if not rxn.subsystem:
                count +=1
        print(model,count)
            
        

if plotGenesRxns:
    def count_annotations_and_keywords(fasta_file):
        with open(fasta_file, 'r') as file:
            gene_count = 0  #  count '>'

            enzyme_count = 0  # count enzyme or ase in annotation

            for line in file:
                # Check if line starts with '>' (start of a new gene annotation)
                if line.startswith('>'):
                    gene_count += 1
                    
        print(f"Total gene annotations (>): {gene_count}")
        print(f"Total of 'enzyme' or 'ase':{enzyme_count}")
        
        return gene_count 

    for idx, genome in enumerate(genome_list):
        globals()[f'genes_{species_name[idx]}'] = count_annotations_and_keywords(genome_folder+genome)
        globals()[f'enzymes_{species_name[idx]}'] = np.array([str(x).count('.',0,2) for x in globals()[f'{species_name[idx]}_annotations']['EC']]).sum()
    

    gf_carve_data_myco= {'M. tuberculosis H37Rv': {'Reactions': len(Rv.reactions),'Genes': genes_Rv, 'Enzymes':enzymes_Rv}, 
            'M. marinum':{'Reactions': len(marinum.reactions),'Genes': genes_marinum, 'Enzymes':enzymes_marinum},
            'M. smegmatis':{'Reactions': len(smegmatis.reactions),'Genes': genes_smegmatis, 'Enzymes':enzymes_smegmatis},
            'M. abscessus':{'Reactions': len(abscessus.reactions),'Genes': genes_abscessus, 'Enzymes':enzymes_abscessus},
            'M. aromaticivorans':{'Reactions': len(aromaticivorans.reactions),'Genes': genes_aromaticivorans, 'Enzymes':enzymes_aromaticivorans}}
            
    

    reaction_list =[]
    gene_list = []
    enzyme_list = []
    for bact in gf_carve_data_myco:
        reaction_list.append(gf_carve_data_myco[bact]['Reactions'])
        gene_list.append(gf_carve_data_myco[bact]['Genes'])
        enzyme_list.append(gf_carve_data_myco[bact]['Enzymes'])
        
    vocabulary = [1 , 2 , 3, 4, 5, 6]
    my_colors = {1:'red',2:'green',3:'blue', 4:'orange', 5:'cyan', 6:'purple'}
    reg = LinearRegression(fit_intercept=False).fit(np.array(gene_list).reshape((-1, 1)), np.array(reaction_list))
    R_2 = reg.score(np.array(gene_list).reshape((-1, 1)), np.array(reaction_list).reshape((-1, 1)))
    x_regression= list(range(0,7000,100))
    y_regression =  pred_reaction_no = reg.coef_*(x_regression) + reg.intercept_

    plt.figure(figsize=(10,8))

    for i,j in enumerate(gene_list):
        plt.scatter(gene_list[i] , reaction_list[i], color ='red',s=50)
        plt.text(gene_list[i]-450,reaction_list[i]+40, list(gf_carve_data_myco)[i], fontsize =12)
        
    plt.plot(x_regression,y_regression,'b-')
    plt.xlabel('Number of genes',size= 30)
    plt.ylabel('Number of reactions',size =30)
    #plt.title('Number of reactions vs number of protein-coding genes', size= 36)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlim(left=0)
    plt.ylim(bottom =0)
    #plt.legend(['R2: '+str(R_2)[0:5]])
    plt.savefig(f'./gene_rxns_{date}_{model_type}.png')


    plt.show()
    
    reg = LinearRegression(fit_intercept=False).fit(np.array(enzyme_list).reshape((-1, 1)), np.array(reaction_list))
    x_regression= list(range(0,3000,100))
    y_regression =  pred_reaction_no = reg.coef_*(x_regression) + reg.intercept_

    plt.figure(figsize=(10,8))

    for i,j in enumerate(enzyme_list):
        plt.scatter(enzyme_list[i] , reaction_list[i], color ='red',s=200)
        plt.text(enzyme_list[i]-50,reaction_list[i]+40, list(gf_carve_data_myco)[i], fontsize =12)
        
    plt.plot(x_regression,y_regression,'b-')
    plt.xlabel('Number of enzymes',size= 30)
    plt.ylabel('Number of reactions',size =30)

    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlim(left=0)
    plt.ylim(bottom =0)
    plt.savefig(f'./enzyme_rxns_{date}_{model_type}.png')


    plt.show()

if plotSuperVenn:
    mari_rxns = {x.id for x in marinum.reactions}
    #Ra_rxns = {x.id for x in Ra.reactions}
    Rv_rxns = {x.id for x in Rv.reactions}
    absc_rxns = {x.id for x in abscessus.reactions}
    smeg_rxns = {x.id for x in smegmatis.reactions}
    arom_rxns = {x.id for x in aromaticivorans.reactions}
    
    
    slow_rxns_uni =  (set(mari_rxns) | set(Rv_rxns))
    fast_rxns_uni = (set(arom_rxns) | set(smeg_rxns) | set(absc_rxns))
    intersection = (set.intersection(set(mari_rxns),set(arom_rxns),set(smeg_rxns),set(absc_rxns),set(Rv_rxns)))
    slow_ven  = [mari_rxns,Rv_rxns]
    fast_ven = [absc_rxns,arom_rxns,smeg_rxns]
        
    sets =[Rv_rxns, absc_rxns, mari_rxns,  arom_rxns,smeg_rxns]
    slow_fast_names = ['Slow growers','Fast growers']

    all_set = [Rv_rxns, absc_rxns, mari_rxns ,  arom_rxns,smeg_rxns]
    venn_dict_gen_size = {'$M. tuberculosis$':set(Rv_rxns),'$M. abscessus$':set(absc_rxns),'$M. marinum$':set(mari_rxns), '$M. aromaticivorans$':set(arom_rxns), '$M. smegmatis$':set(smeg_rxns) }
    venn_dict_fast_slow = {'$M. tuberculosis$':set(Rv_rxns),'$M. marinum$':set(mari_rxns), '$M. abscessus$':set(absc_rxns),'$M. aromaticivorans$':set(arom_rxns), '$M. smegmatis$':set(smeg_rxns) }

    plt.figure(figsize=(20, 10))
    supervenn(sets, list(venn_dict_gen_size.keys()), 
            chunks_ordering='occurrence', 
            #widths_minmax_ratio= 0.1,
            min_width_for_annotation = 80, 
            rotate_col_annotations=False, 
            col_annotations_area_height=2,
            side_plot_width = 1.5,
            fontsize=18
            )
    plt.xlabel('Number of Reactions')
    plt.ylabel('Species')
    plt.savefig(f'./supervenn_{date}_{model_type}.png', bbox_inches='tight', dpi= 500)

    plt.show()

    plt.figure(figsize=(10, 10))
    venn(venn_dict_gen_size)
    plt.savefig(f'./venn_{date}_{model_type}.png')

    #plt.title("Venn Diagram")

if groupSubsystem:

    # Extract subsystems from all reactions in all models and
    
    # Returns susbsystems at level1 or high level 

    def key_by_value_level1(df, value):
        for key, values in df.items():
            for key2, values2 in values.items():
                if value in values2:
                    return key
    #   return None
    # Returns subsystem at level2 or detailed level
    def key_by_value_level2(df, value):
        for key, values in df.items():
            for key2, values2 in values.items():
                if value in values2:
                    return key2
    # Returns subsystem at level3 or raw level
    def key_by_value_level3(df, value):
        for key, values in df.items():
            for key2, values2 in values.items():
                if value in values2:
                    return value
                    
#        return None
    # Group subsystems
    subsys_group= {
                    'Aminoacid metabolism':{
                        'Arginine biosynthesis':['Arginine biosynthesis'],
                        'Arginine and Proline Metabolism':
                            ['S_Arginine_and_Proline_Metabolism',
                            'Amino acid metabolism: Proline and Arginine metabolism','Spermidine biosynthesis',
                            'S_Alternate_Carbon_and_Nitrogen_source__Amines_and_Polyamines_Metabolism',
                            'Arginine and Proline Metabolism','Arginine and proline metabolism','Amino acid metabolism: Proline and Arginine metabolism'], 
                                                                       
                        'Histidine Metabolism':
                            ['Histidine Metabolism','S_Histidine_Metabolism','Histidine metabolism','Ergothioneine biosynthesis'],                                                 
                    
                        'Glycine, serine and threonine metabolism':
                            ['Glycine, serine and threonine metabolism','S_Glycine__Serine_and_threonine_metabolism',
                            'Glycine and Serine Metabolism','S_Threonine_and_Lysine_Metabolism',
                            'Threonine and Lysine Metabolism',
                            'S_Lysine_Metabolism',
                            'S_Alternate_Carbon_and_Nitrogen_source__Ectoine_Metabolism'],

                        'Lysine metabolism':['Lysine biosynthesis','Lysine degradation','L-alpha-aminoadipic acid (L-AAA) biosynthesis'], #ask
                        'Tryptophan Metabolism':['Tryptophan metabolism'],
                        'Tyrosine Metabolism':
                            ['Tyrosine metabolism'],
                        'Phenylalanine, tyrosine and tryptophan biosynthesis': ['Tyrosine, Tryptophan, and Phenylalanine Metabolism', 
                            'S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism','S_Phenylalanine_Tyrosine_Tryptophan_Biosynthesis','Phenylalanine, tyrosine and tryptophan biosynthesis'],
                        'Phenylalanine metabolism': ['Phenylalanine metabolism','Phenylalanine Metabolism'],
                        'Valine, Leucine, and Isoleucine Metabolism':
                            ['Valine, Leucine, and Isoleucine Metabolism',
                            'Amino acid metabolism: Valine, Leucine and Isoleucine biosynthesis',
                            'S_Valine__Leucine__and_Isoleucine_Metabolism', 'Valine, leucine and isoleucine degradation','Valine, leucine and isoleucine (degradation)', 'Valine, leucine and isoleucine (biosynthesis)','Valine, leucine and isoleucine biosynthesis'], 
                        'Cysteine and methionine metabolism':
                            ['Cysteine Metabolism',
                            'Amino acid metabolism: Cysteine and Methionine metabolism',
                            'Cysteine metabolism',
                            'S_Cysteine_Metabolism',
                            'Cysteine and methionine metabolism',
                            'Methionine metabolism', 
                            'Methionine Metabolism',
                            'Methionine Salvage Pathway', 
                            'S_Methionine_Metabolism'],     
                        'Alanine, Aspartate, and Glutamate Metabolism':
                            ['Alanine and Aspartate Metabolism','S_Alanine_and_Aspartate_Metabolism',
                            'Alanine, Aspartate, and Glutamate Metabolism',
                            'Alanine, aspartate and glutamate metabolism', 
                            'S_Glutamate_Metabolism',
                            'Glutamate metabolism',
                            'Glutamate Metabolism','Amino acid metabolism: Alanine, Aspartate and Glutamate metabolism']},
               

                    
                    'Biomass and maintenance functions': {
                        'Biomass and maintenance functions':  [
                            'Biomass and maintenance functions']},
                    
                    'Carbohydrate metabolism':{
                        'Propanoate Metabolism':
                            ['Propanoate Metabolism', 
                            'Propanoate metabolism',
                            'S_Alternate_Carbon__Propanoate_Metabolism',
                            'Propanoate metabolism','Methylglyoxal'],
                        'Levulinate metabolism':['Levulinate metabolism','S_Alternate_Carbon__Levulinate_Metabolism'],
                        
                        'Pentose Phosphate Pathway':['Pentose Phosphate Pathway',
                            'S_Alternate_Carbon__Ribose_Metabolism',
                            'Ribose metabolism',
                            'Pentose phosphate pathway',
                            'S_Pentose_Phosphate_Pathway',
                            "[array(['Pentose Phosphate Pathway'], dtype='<U25')]",
                            'Calvin cycle/Pentose phosphate pathway'],
                        'Inositol phosphate metabolism':['Inositol phosphate metabolism','Inositol Metabolism'],
                        'Ascorbate and Aldarate Metabolism':
                            ['S_Alternate_Carbon__Ascorbate_and_Aldarate_Metabolism',
                            'Ascorbate and aldarate metabolism'],
                        'Glyoxylate and dicarboxylate metabolism': ['Glyoxylate and dicarboxylate metabolism', 'S_Glyoxylate_and_dicarboxylate_metabolism', "[array(['Glyoxylate Metabolism'], dtype='<U21')]",
                            'Glyoxylate Metabolism','Glyoxylate Biosynthesis','Photorespiration/Glyoxylate Degredation'],
                        'Pyruvate metabolism': ["[array(['Pyruvate Metabolism'], dtype='<U19')]",
                            'Pyruvate Metabolism',
                            'Methylglyoxal Metabolism',
                            'S_Pyruvate_Metabolism',
                            'Pyruvate metabolism','Anaplerotic Pathway','Anaplerotic Reactions'],
                        'Amino sugar and nucleotide sugar metabolism':['Amino sugar and nucleotide sugar metabolism','S_Alginate_biosynthesis','Biosynthesis of amino sugars and nucleotide sugars'],
                        'Carbohydrate metabolism':['S_Carbohydrates_and_related_molecules',
                            'Carbohydrate Metabolism',
                            'Carbohydrate Metabolism',
                            'Carbohydrate metabolism'],
                        'Galactose metabolism':['Galactose metabolism'],
                        'Starch and Sucrose metabolism':['Starch and Sucrose metabolism','S_Starch_and_Sucrose_Metabolism',
                            'Starch and sucrose metabolism','Glycogen and sucrose metabolism',
                            'Starch and Sucrose Metabolism','Glycogen metabolism','Sugar Metabolism'],
                        'Butanoate metabolism':['Butanoate metabolism','S_Butanoate_Metabolism'],
                        'Glycolysis/Gluconeogenesis': ['Glycolysis/ Gluconeogenesis',                            
                            'Glycolysis / Gluconeogenesis',
                            'S_Gluconeogenesis',
                            'S_Glycolysis',
                            'Glycolysis',
                            'Glycolysis/Gluconeogenesis',
                            'Glycolysis/Gluconeogensis/Butanoate metabolism'],
                        'Fructose and mannose metabolism':
                            ['Fructose and mannose metabolism','S_Alternate_Carbon__Fructose_Metabolism'],
                        'Citrate Cycle':  [
                            'Citric Acid Cycle',
                            'TCA cycle',
                            'Citrate cycle',
                            'S_TCA_Cycle',
                            "[array(['Citric Acid Cycle'], dtype='<U17')]"],
                        'Pentose and Glucuronate Interconversions':['Pentose and Glucuronate Interconversions','Pentose and glucuronate interconversions']
                        },
                    
                    'Cell envelop and membrane':{
                        'Cell envelop and membrane': [
                            'S_Membrane_bioenergetics',
                            'S_Cell_Envelope_Biosynthesis',
                            
                            'S_Cell_Envelope_Biosynthesis__Biosynthesis_of_L_glycero_D_manno_heptose__Hep_',
                            'S_Cell_Envelope_Biosynthesis__Cellulose_Metabolism',
                            'Cell Envelope Biosynthesis', 
                            'Cell Wall/Membrane/Envelope Metabolism',

                            'Lipid &amp; Cell Wall Metabolism',
                            'Membrane Lipid Metabolism',
                            'Membrane Metabolism',
                            "[array(['Membrane Metabolism'], dtype='<U19')]",
                            'S_Cell_Envelope_Biosynthesis__Cellulose_Metabolism']},
                    
                        

                    'Energy Metabolism': {  
                        'Methane metabolism':['Methane metabolism','Methane Metabolism', 'Methanogenesis','S_Formaldehyde_Metabolism'],
                        'Nitrogen metabolism':['Nitrogen metabolism','S_Nitrogen_Metabolism',
                            'Nitrogen Metabolism'],
                        'Oxidative Phosphorylation': ['S_Oxidative_Phosphorylation','Oxidative Phosphorylation',
                            'Oxidative phosphorylation','S_Inorganic_polyphosphates_metabolism','Hydrogen production'],
                        'Other carbon fixation pathways':['Other carbon fixation pathways'],
                        'Sulfur Metabolism':['Sulfur metabolism','S_Sulfur_Metabolism','S_Phosphate_and_sulfur','Sulfite metabolism','SOB reactions']},
                              
                    
                    'Exchange': {
                        'Exchange':['Extracellular exchange',
                            'Exchange']}, 
                    
                    'Glycan biosynthesis and metabolism': {
                        'O-Antigen biosynthesis':['S_Cell_Envelope_Biosynthesis__O_antigen_Biosynthesis','S_Cell_Envelope_Biosynthesis__O_antigen_Biosynthesis'],
                            
                        'Peptidoglycan biosynthesis': ['Peptidoglycan Metabolism', 
                            "[array(['Peptidoglycan Metabolism'], dtype='<U24')]",'Peptidoglycan biosynthesis','S_Cell_Envelope_Biosynthesis__Peptidoglycan_Biosynthesis'],  
                        
                        'Other glycan degradation': ['Other glycan degradation (aka murein recycling)','Murein Recycling','Murein Biosynthesis'],
                        'Lipopolysaccharide biosynthesis':['Lipopolysaccharide biosynthesis','Lipopolysaccharide Biosynthesis / Recycling'],
                        'Arabinogalactan biosynthesis':['Arabinogalactan biosynthesis'],
                        'Lipoarabinomannan and lipomannan and lipoarabinomanan biosynthesis':['Lipoarabinomannan and lipomannan and lipoarabinomanan biosynthesis']        
                        },
                    
                     
                    'Human disease':{
                        'Tuberculosis':['Tuberculosis'],
                        'Human disease and host interactions':
                        ['Human disease and host interactions']},
                    
                    'Genetic Information Processing':{
                        'Translation':
                        ['Aminoacyl-tRNA biosynthesis']
                    },
                    
                    
                    'Lipid Metabolism':{
                        'Phthiocerol biosynthesis':['Phthiocerol biosynthesis'],
                        'Glycolipid biosynthesis':['Glycolipid biosynthesis'],
                        'Fatty acid elongation':['Fatty acid elongation','S_PHAs_Metabolism'],
                        'Steroid biosynthesis':['Steroid metabolism','Biosynthesis of steroids'],
                        'Fatty acid degradation':['Fatty acid degradation','Fatty Acid Degradation','Beta oxidation of unsaturated fatty acids'],
                        'Fatty acid biosynthesis':['Fatty acid biosynthesis','S_Fatty_Acid__Biosynthesis'],
                        'Fatty Acid Metabolism':["[array(['Fatty Acid Metabolism'], dtype='<U21')]",'Fatty Acid Metabolism','S_Fatty_Acid_Metabolism'],
                        'Biosynthesis of unsaturated fatty acids':['Biosynthesis of unsaturated fatty acids'],
                        'Mycolic acid biosynthesis':['Mycolic acid biosynthesis','Mycolic Acid Biosynthesis','Mycolic acid pathway'],
                        'Glycerophospholipid metabolism':['Glycerophospholipid metabolism','Glycerophospholipid Metabolism',"[array(['Glycerophospholipid metabolism'], dtype='<U30')]",
                            'S_Glycerophospholipid_Metabolism'],
                        'Glycerolipid metabolism':['Glycerolipid metabolism','Glycerolipid Synthesis','Glycerolipid Metabolism','1 2 Propanediol Catabolism','Propanediol Catabolism']
                    },                           
                    
                    'Metabolism of other amino acids':{
                            'Mycothiol biosynthesis':['Mycothiol biosynthesis'],
                            'D-Amino acids metabolism':['D-Amino acids metabolism','S_Alternate_Carbon_and_Nitrogen_source__D_Amino_acids_Metabolism'],
                            'Taurine and hypotaurine metabolism':['Taurine and hypotaurine metabolism'],
                            'Other Amino Acid Metabolism':['Other Amino Acid Metabolism','S_Alternate_Carbon_and_Nitrogen_source__D_Amino_acids_Metabolism'],
                            'beta-Alanine metabolism':['beta-Alanine metabolism'],
                            'Selenocompound metabolism': ['Selenocompound metabolism'],
                            'Glutathione metabolism':['Glutathione metabolism'],
                            'Phosphonate and phosphinate metabolism':['Phosphonate and phosphinate metabolism'],
                            'Other dipeptide Metabolism':['S_Alternate_Carbon_and_Nitrogen_source__Dipeptide_Metabolism']},


                    'Miscellaneous, unassigned, sink, and spontaneous reactions': {   
                        'Miscellaneous, unassigned, sink, and spontaneous reactions': [
                            'Miscellaneous', 
                            'Miscellaneous and spontaneous reactions',
                            'Unassigned', 
                            'Others', 
                            'Other',
                            'S_Formaldehyde_Metabolism',
                            'Energy Production and Conversion',
                            'Redox Metabolism', 
                            'Energy metabolism',
                            'Energy Metabolism',
                            'Oxidoreduction of electron transfer metabolites',
                            'NaN',
                            'Lipid metabolism',
                            'Lipid Metabolism',
                            'S_Lipids',
                            'Amino Acid Metabolism',
                            'S_Amino_acids_and_related_molecules',
                            'Central Metabolism',
                            'Carbon metabolism',
                            'Alternate Carbon Metabolism', 
                            'Alternate carbon metabolism',
                            'S_Alternate_Carbon',
                            '',
                            'Urea Cycle','Urea cycle and metabolism of amino groups',
                            'S_Urea_cycleamino_group_metabolism',
                            'Urea cycle',
                            'Miscellaneous and Unassigned',
                            'S_Plant_growth_promoting',
                            'Plant Growth promoting', # check indole-3-acetaldehyde
                            'Metabolite Repair',
                            'Intracellular demand']},
                    
                        #'Urea cycle':
                        #    ['Urea Cycle','Urea cycle and metabolism of amino groups',
                        #    'S_Urea_cycleamino_group_metabolism',
                        #    'Urea cycle']},
                    
                    'Nucleotide Metabolism': {
                        'Nucleotide Metabolism': ['Purine and Pyrimidine Biosynthesis',
                            'S_Nucleotides_and_nucleic_acids',
                            'Nucleotide Metabolism','S_Alternate_Carbon_and_Nitrogen_source__Nucleotide_Metabolism'],
                        'Purine metabolism': ['S_Nucleotide_Salvage_Pathway','Nucleotide Salvage Pathway','Nucleotide salvage pathway','Purine metabolism (De novo)','Purine metabolism',     
                            'S_Purine_Metabolism','Purine metabolism A G I Orn'],
                        'Pyrimidine metabolism':['Pyrimidine metabolism','Pyrimidine metabolism (De novo)',
                            'S_Pyrimidine_Metabolism']},
                
                    
                    #'Poliamine metabolism': {
                    #    'Poliamine metabolism':['Spermidine biosynthesis', 
                    #        'S_Alternate_Carbon_and_Nitrogen_source__Amines_and_Polyamines_Metabolism']},
                    
                    'Biosynthesis of other secondary metabolites':{
                        'Phenylpropanoid biosynthesis':['Phenylpropanoid biosynthesis'],
                        'Biosynthesis of various alkaloids':['Biosynthesis of various alkaloids']

                    },

                    
                    'Metabolism of terpenoids and polyketides':{  
                        'Biosynthesis of siderophore group nonribosomal peptides':['Biosynthesis of siderophore group nonribosomal peptides','Biosynthesis of siderophore','Mycobactin biosynthesis','siderophore biosynthesis'],  
                        'Terpenoid backbone biosynthesis':['Terpenoid backbone biosynthesis', 'S_Cofactor_and_Prosthetic_Group_Biosynthesis__Terpenoid_backbone_biosynthesis']},    
                   
                    
                    'Metabolism of cofactor and vitamins': {
                        'Vitamin and Cofactor Biosynthesis':['Vitamins &amp; Cofactor Biosynthesis',
                            'Vitamins and Cofactor Biosynthesis',                            
                            'Cofactor biosynthesis: Vitamin B12',
                            'Cofactor biosynthesis: Lipoate biosynthesis'], 
                        'Cofactor and Prosthetic Group Biosynthesis':['Cofactor and Prosthetic Group Biosynthesis','S_Cofactor_and_Prosthetic_Group_Biosynthesis','Cofactor and Prosthetic Group Metabolism',"[array(['Cofactor and Prosthetic Group Biosynthesis'], dtype='<U42')]",'S_Coenzymes_and_prosthetic_groups'],
                        'Polyprenyl Metabolism':['Polyprenyl Metabolism'],
                        'Vitamin B6 metabolism':['Vitamin B6 metabolism','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Vitamin_B6_Metabolism'],
                        'Biotin metabolism':['Biotin metabolism','Biotin Biosynthesis'],
                        'Riboflavin metabolism':['Riboflavin Metabolism', 'Riboflavin metabolism','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Riboflavin_Metabolism'],
                        'Pantothenate and CoA Biosynthesis':['Pantothenate and CoA Metabolism', 'Cofactor biosynthesis: Pantothenate and CoA biosynthesis','Pantothenate and CoA biosynthesis','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Pantothenate_and_CoA_Biosynthesis'],
                        'Ubiquinone and other terpenoid-quinone biosynthesis': ['Ubiquinone and other terpenoid-quinone biosynthesis','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Ubiquinone_biosynthesis','Ubiquinone biosynthesis'],    
                        'Porphyrin metabolism':['Porphyrin metabolism','Porphyrin and chlorophyll metabolism','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Porphyrin_and_Chlorophyll_Metabolism'],
                        'Thiamine Metabolism':['Thiamine Metabolism','Thiamine metabolism','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Thiamine_Metabolism'],
                        'Nicotinate and nicotinamide metabolism':['S_Cofactor_and_Prosthetic_Group_Biosynthesis__Nicotinamide_Biosynthesis', 'Nicotinate and nicotinamide metabolism',"[array(['Nicotinate and nicotinamide metabolism'], dtype='<U38')]"],
                        'One carbon pool by folate':['One carbon pool by folate',
                            'Folate Metabolism','Folate Biosynthesis',
                            'Folate biosynthesis','S_Cofactor_and_Prosthetic_Group_Biosynthesis__Folate_Biosynthesis','S_Cofactor_and_Prosthetic_Group_Biosynthesis__One_Carbon_pool_by_folate']
                        },

                    'Environmental Information Processing':{'Arsenic resistance':['Arsenic resistance'],
                        'Signal transduction / Two-component system':['Signal transduction / Two-component system'],
                        'Iron binding by siderophores':['Iron binding by siderophores'],
                        'Bacterial secretion system':['Bacterial secretion system'],
                        'Heavy Metal Tolerance':['S_Heavy_Metal_Tolerance'],
                        'Cellular processes':['Cellular community - prokaryotes'],
                        'Transport': ['Transport of unknown mechanism',
                            'S_Iron_uptake_and_metabolism',
                            'Transport',
                            'Inorganic Ion Transport and Metabolism', 
                            'S_Transport__solvent_extrusion',
                            'GLYCOLATE_transport'],
                        'Membrane transport':['Porines and efflux systems',
                            'ABC transporters',
                            'S_Transport__ABC_system',
                            'Lipid transporters',
                            'Phosphotransferase system (PTS)',
                            'Sugar transport',
                            'Electrochemical potential-driven transport',
                            'Transport via diffusion','S_Transport__ABC_system',
                            'Transport Outer Membrane Porin',
                            'Transport Outer Membrane',
                            'Transport: Inner Membrane',
                            'Transport, Inner Membrane', 
                            'Tranpsort, Inner Membrane',
                            'Transport, Outer Membrane', 
                            'Transport; Inner Membrane',
                            'Transport, Outer Membrane Porin','S_Transport__Outer_Membrane',
                            'S_Transport__Inner_Membrane',
                            'Tranpsort, Inner Membrane',
                            'Transport; Inner Membrane']

                        
                              },

                    'Xenobiotics degradation and metabolism':{
                        'Steroid degradation':['Cholesterol degradation', # Steroid metabolism
                            "[array(['Cholesterol degradation'], dtype='<U23')]",'Steroid degradation'],
                        'Styrene degradation':['Styrene degradation'],
                        'Degradation of polycyclic aromatic hydrocarbon':['Degradation of polycyclic aromatic hydrocarbon',
                            'S_Aromatic_Compounds_Degradation__Homogentisate_pathway', 
                            'Aromatic Compound degradation',
                            'S_Aromatic_Compounds_Degradation__Phenylacetyl_CoA_Catabolom',
                            'Phenylacetyl-CoA pathway',
                            'S_Aromatic_Compounds_Degradation__B_Ketoadipate_pathway'],
                        'Benzoate degradation':['Benzoate degradation','Benzanoate degradation'],
                        'Aminobenzoate degradation':['Aminobenzoate degradation'],
                        'Toluene degradation':['Toluene degradation','S_Aromatic_Compounds_Degradation__Toluene_Pathway'],
                        'Xylene degradation':['Xylene degradation'],
                        'Nitrotoluene degradation':['Nitrotoluene degradation']
                        },
                         
                    }

    subsys_group_df = pd.DataFrame(subsys_group)
    subsys_group_df.to_csv(f'subsys_classification_{date}.csv')
    
    subsystems = []

    for model in model_list:
        for reaction in model.reactions:
            if reaction.subsystem:
                subsystems.append(reaction.subsystem)
            else:
                subsystems.append('NaN')
            if reaction.id == 'PPDOy' or reaction.id == 'LCARSyi':
                model.remove_reactions([reaction])
        
    
    subsystems_df = pd.DataFrame(subsystems, columns=["subsystem"])

    # save unique subsystem dataframe
    unique_subsystems = subsystems_df["subsystem"].unique()
    unique_subsystems_df = pd.DataFrame(subsystems_df["subsystem"].unique())
    unique_subsystems_df.to_csv(f'all_subsystems_{model_type}_{date}.csv')
    
    # update the subsystems in the models
    
    if level_subsystem == 'detailed':
        for model in model_list:
            for reaction in model.reactions:
                if reaction.subsystem in unique_subsystems:
                    reaction.subsystem = key_by_value_level2(subsys_group, reaction.subsystem) #more detailed subsystems
    
    if level_subsystem == 'high':
        for model in model_list:
            for reaction in model.reactions:
                if reaction.subsystem in unique_subsystems:
                    reaction.subsystem = key_by_value_level1(subsys_group, reaction.subsystem) #higher hierarchy subsystems
    
    if level_subsystem == 'raw':
        for model in model_list:
            for reaction in model.reactions:
                reaction.subsystem = key_by_value_level3(subsys_group, reaction.subsystem) #higher hierarchy subsystems
    
    
    for idx, (model, file) in enumerate(zip(model_list,list(selected_models.values()))):
        write_sbml_model(model,f'{file}_{level_subsystem}_subsystem.xml') 
        save_json_model(model,f'{file}_{level_subsystem}_subsystem.json') 
        save_matlab_model(model,f'{file}_{level_subsystem}_subsystem.mat') 



    # start subsystem analysis    
    p_value_cutoff = 0.05
    

    Reaction = namedtuple('Reaction', ['name', 'subsystem','gene'])
    if arom == 'arom':
        all_models_rxns = list(marinum.reactions) + list(smegmatis.reactions) + list(abscessus.reactions) + list(Rv.reactions) + list(aromaticivorans.reactions)
    else:
        all_models_rxns = list(marinum.reactions) + list(smegmatis.reactions) + list(abscessus.reactions) + list(Rv.reactions) #+ list(aromaticivorans.reactions)

    mari_rxns = [Reaction(x.id, x.subsystem, x.gene_reaction_rule) for x in marinum.reactions ]
    absc_rxns =  [Reaction(x.id, x.subsystem, x.gene_reaction_rule) for x in abscessus.reactions ]
    smeg_rxns =  [Reaction(x.id, x.subsystem, x.gene_reaction_rule) for x in smegmatis.reactions ]
    if arom == 'arom': 
        arom_rxns = [Reaction(x.id, x.subsystem, x.gene_reaction_rule) for x in aromaticivorans.reactions]
    
    rv_rxns = [Reaction(x.id, x.subsystem, x.gene_reaction_rule) for x in Rv.reactions ]
    
    if arom == 'arom' : 
        all_reactions = [reaction for reaction in mari_rxns + arom_rxns + smeg_rxns + absc_rxns + rv_rxns if isinstance(reaction, Reaction)] #+ arom_rxns 
    else:
        all_reactions = [reaction for reaction in mari_rxns + smeg_rxns + absc_rxns + rv_rxns if isinstance(reaction, Reaction)] #+ arom_rxns 

    reaction_dict = {reaction.name: reaction for reaction in all_reactions}

    # extract reaction names from each list
    reaction_names_mari = {reaction.name for reaction in mari_rxns}
    reaction_names_absc = {reaction.name for reaction in absc_rxns}
    reaction_names_smeg = {reaction.name for reaction in smeg_rxns}
    if arom == 'arom':
        reaction_names_arom = {reaction.name for reaction in arom_rxns}

    reaction_names_mtb = {reaction.name for reaction in rv_rxns}

    subsystem_mari = [reaction.subsystem for reaction in mari_rxns]
    subsystem_absc = [reaction.subsystem for reaction in absc_rxns]
    subsystem_smeg = [reaction.subsystem for reaction in smeg_rxns]
    if arom == 'arom':
        subsystem_arom = [reaction.subsystem for reaction in arom_rxns]
        subsystem_arom_uni = [reaction.subsystem for reaction in arom_rxns if reaction not in [mari_rxns,absc_rxns,smeg_rxns,rv_rxns]]
    subsystem_mtb = [reaction.subsystem for reaction in rv_rxns]

    
    # Find similar reactions (those that appear in more than one list)
    if arom == 'arom':
        all_reaction_names = [reaction_names_mari, reaction_names_absc,reaction_names_mtb, reaction_names_smeg, reaction_names_arom]
    else:
        all_reaction_names = [reaction_names_mari, reaction_names_absc,reaction_names_mtb, reaction_names_smeg]
    
    similar_reactions_set = set.intersection(*all_reaction_names)

    # Find different reactions (those that appear in only one list)
    if arom == 'arom':
        all_reaction_names_flat = reaction_names_mari | reaction_names_absc | reaction_names_mtb | reaction_names_smeg | reaction_names_arom
    else:
        all_reaction_names_flat = reaction_names_mari | reaction_names_absc | reaction_names_mtb | reaction_names_smeg

    different_reactions_set = all_reaction_names_flat - similar_reactions_set
    different_reactions_set_mari = reaction_names_mari - similar_reactions_set
    different_reactions_set_absc = reaction_names_absc - similar_reactions_set
    different_reactions_set_mtb = reaction_names_mtb - similar_reactions_set
    different_reactions_set_smeg = reaction_names_smeg - similar_reactions_set
    if arom == 'arom':
        different_reactions_set_arom = reaction_names_arom - similar_reactions_set
    
    unique_reactions_set_arom = reaction_names_arom - reaction_names_mari - reaction_names_absc- reaction_names_mtb - reaction_names_smeg
    # Convert sets to lists of Reaction objects
    similar_reactions = [reaction_dict[name] for name in similar_reactions_set]
    different_reactions = [reaction_dict[name] for name in different_reactions_set]
    different_reactions_arom = [reaction_dict[name] for name in unique_reactions_set_arom]

    # Extract subsystems for similar and different reactions
    similar_subsystems = [reaction.subsystem for reaction in similar_reactions]
    different_subsystems = [reaction.subsystem for reaction in different_reactions]

    # Extract genes for similar and different reactions

    similar_genes = [1 if reaction.gene != '' else 0 for reaction in similar_reactions ]
    different_genes = [1 if reaction.gene != '' else 0 for reaction in different_reactions ]

    subsystem_int_mari = [reaction.subsystem for reaction in mari_rxns if reaction.name in similar_reactions_set]
    subsystem_int_absc = [reaction.subsystem for reaction in absc_rxns if reaction.name in similar_reactions_set]
    subsystem_int_smeg = [reaction.subsystem for reaction in smeg_rxns if reaction.name in similar_reactions_set]
    if arom == 'arom':
        subsystem_int_arom = [reaction.subsystem for reaction in arom_rxns if reaction.name in similar_reactions_set]
    subsystem_int_mtb = [reaction.subsystem for reaction in rv_rxns if reaction.name in similar_reactions_set]

    subsystem_dif_mari = [reaction.subsystem for reaction in mari_rxns if reaction.name in different_reactions_set_mari]
    subsystem_dif_absc = [reaction.subsystem for reaction in absc_rxns if reaction.name in different_reactions_set_absc]
    subsystem_dif_smeg = [reaction.subsystem for reaction in smeg_rxns if reaction.name in different_reactions_set_smeg]
    if arom == 'arom':
        subsystem_dif_arom = [reaction.subsystem for reaction in arom_rxns if reaction.name in different_reactions_set_arom]
    subsystem_dif_mtb = [reaction.subsystem for reaction in rv_rxns if reaction.name in different_reactions_set_mtb]
    rxn_subsystem_sim = [{reaction.subsystem: reaction.name} for reaction in mari_rxns if reaction.name in similar_reactions_set]

    rxn_subsystem_dif_mari = [{reaction.subsystem: reaction.name} for reaction in mari_rxns if reaction.name in different_reactions_set_mari]
    rxn_subsystem_dif_absc = [{reaction.subsystem: reaction.name} for reaction in absc_rxns if reaction.name in different_reactions_set_absc]
    rxn_subsystem_dif_smeg = [{reaction.subsystem: reaction.name} for reaction in smeg_rxns if reaction.name in different_reactions_set_smeg]
    if arom == 'arom':
        rxn_subsystem_dif_arom = [{reaction.subsystem: reaction.name} for reaction in arom_rxns if reaction.name in different_reactions_set_arom]
    rxn_subsystem_dif_mtb = [{reaction.subsystem: reaction.name} for reaction in rv_rxns if reaction.name in different_reactions_set_mtb]

    def detailRxns(df_rxn, name):
        processed_data = {}
        for item in df_rxn:
            for key, value in item.items():
                if key not in processed_data:
                    processed_data[key] = []
                processed_data[key].append(value)
        df = pd.DataFrame.from_dict(processed_data, orient='index').transpose()
        return df
    
    
    df1 = detailRxns(rxn_subsystem_dif_mari, 'marinum')
    df2 = detailRxns(rxn_subsystem_dif_absc, 'abscessus')
    df3 = detailRxns(rxn_subsystem_dif_arom, 'aromaticivorans')
    df4 = detailRxns(rxn_subsystem_dif_smeg, 'smegmatis')
    df5 = detailRxns(rxn_subsystem_dif_mtb, 'tuberculosis')
    df6 = detailRxns(rxn_subsystem_sim, 'shared')
    #df1.to_excel(pd.ExcelWriter(f'./non_shared_and_shared_rxn_subsystems_detail_species_{model_type}_{date}.xls', engine='openpyxl'), sheet_name='marinum', index=False)
    #df2.to_excel(pd.ExcelWriter(f'./non_shared_and_shared_rxn_subsystems_detail_species_{model_type}_{date}.xls', engine='openpyxl'), sheet_name='abscessus', index=False)
    #df3.to_excel(pd.ExcelWriter(f'./non_shared_and_shared_rxn_subsystems_detail_species_{model_type}_{date}.xls', engine='openpyxl'), sheet_name='aromaticivorans', index=False)
    #df4.to_excel(pd.ExcelWriter(f'./non_shared_and_shared_rxn_subsystems_detail_species_{model_type}_{date}.xls', engine='openpyxl'), sheet_name='smegmatis', index=False)
    #df5.to_excel(pd.ExcelWriter(f'./non_shared_and_shared_rxn_subsystems_detail_species_{model_type}_{date}.xls', engine='openpyxl'), sheet_name='tuberculosis', index=False)
    #df6.to_excel(pd.ExcelWriter(f'./non_shared_and_shared_rxn_subsystems_detail_species_{model_type}_{date}.xls', engine='openpyxl'), sheet_name='shared', index=False)
    df1.to_csv(f'./non_shared_and_shared_rxn_subsystems_detail_mari_{model_type}_{date}.csv')
    df2.to_csv(f'./non_shared_and_shared_rxn_subsystems_detail_absc_{model_type}_{date}.csv')
    df3.to_csv(f'./non_shared_and_shared_rxn_subsystems_detail_arom_{model_type}_{date}.csv')
    df4.to_csv(f'./non_shared_and_shared_rxn_subsystems_detail_smeg_{model_type}_{date}.csv')
    df5.to_csv(f'./non_shared_and_shared_rxn_subsystems_detail_mtb_{model_type}_{date}.csv')
    df6.to_csv(f'./non_shared_and_shared_rxn_subsystems_detail_shared_{model_type}_{date}.csv')


    
# Dis
    all_data = {
    "mari": rxn_subsystem_dif_mari,
    "absc": rxn_subsystem_dif_absc,
    "smeg": rxn_subsystem_dif_smeg,
    "arom": rxn_subsystem_dif_arom,
    "mtb": rxn_subsystem_dif_mtb
    }

    rows = []

    for list_name, items in all_data.items():
        for item in items:
            for subsystem, reaction in item.items():
                rows.append({
                    "Subsystem": subsystem,
                    list_name: reaction
                })

    # Create the DataFrame
    dif_rxn_subsystem_df = pd.DataFrame(rows)

    #dif_rxn_subsystem_df = dif_rxn_subsystem_df.groupby("Subsystem").agg(lambda x: ', '.join(filter(None, x)))



    # save to a CSV file
    dif_rxn_subsystem_df.to_csv(f"dif_rxn_subsystem_df{arom}_{date}_{model_type}.csv", index=False)

    data_intersection= {
        'reaction': [reaction.name for reaction in similar_reactions],
        'subsystem': similar_subsystems
    }
    data_aromaticivorans_dif = {
        'reaction': [reaction.name for reaction in different_reactions_arom],
        'subsystem':  [reaction.subsystem for reaction in different_reactions_arom]

    }

    data_different = {
        'reaction': [reaction.name for reaction in different_reactions],
        'subsystem': different_subsystems
    }

    data_gene_intersection= {
        'reaction': [reaction.name for reaction in similar_reactions],
        'gene': similar_genes
    }

    data_gene_different = {
        'reaction': [reaction.name for reaction in different_reactions],
        'gene': different_genes
    }

    # combine subsystems from both sets
    all_subsystems = set(data_intersection['subsystem']).union(set(data_different['subsystem']))
    all_genes = set(data_gene_intersection['gene']).union(set(data_gene_different['gene']))
    total_subsystems = data_intersection['subsystem'] + data_different['subsystem']
    data_all = {
        'reaction': data_intersection['reaction'] + data_different['reaction'],
        'subsystem': total_subsystems
    }
    data_aromaticivorans = {
        'reaction': data_aromaticivorans_dif['reaction'],
        'subsystem': data_aromaticivorans_dif['subsystem']
    }
    
    # count occurrences of subsystems in each set
    counts_inter = dict((sub,data_intersection['subsystem'].count(sub)) for sub in all_subsystems)
    counts_dif = dict((sub,data_different['subsystem'].count(sub)) for sub in all_subsystems)
    counts_all = dict((sub,total_subsystems.count(sub)) for sub in all_subsystems)


    counts_int_sub_mari = dict((sub,subsystem_int_mari.count(sub)) for sub in all_subsystems) # before all_subsystems
    counts_int_sub_absc = dict((sub,subsystem_int_absc.count(sub)) for sub in all_subsystems)
    if arom == 'arom':
        counts_int_sub_arom = dict((sub,subsystem_int_arom.count(sub)) for sub in all_subsystems)
    counts_int_sub_smeg = dict((sub,subsystem_int_smeg.count(sub)) for sub in all_subsystems)
    counts_int_sub_mtb = dict((sub,subsystem_int_mtb.count(sub)) for sub in all_subsystems)
    
    counts_dif_sub_mari = dict((sub,subsystem_dif_mari.count(sub)) for sub in all_subsystems) # before all_subsystems
    counts_dif_sub_absc = dict((sub,subsystem_dif_absc.count(sub)) for sub in all_subsystems)
    if arom == 'arom':
        counts_dif_sub_arom = dict((sub,subsystem_dif_arom.count(sub)) for sub in all_subsystems)
    counts_dif_sub_smeg = dict((sub,subsystem_dif_smeg.count(sub)) for sub in all_subsystems)
    counts_dif_sub_mtb = dict((sub,subsystem_dif_mtb.count(sub)) for sub in all_subsystems)

    counts_sub_mari = dict((sub,subsystem_mari.count(sub)) for sub in all_subsystems) # before all_subsystems
    counts_sub_absc = dict((sub,subsystem_absc.count(sub)) for sub in all_subsystems)
    if arom == 'arom':
        counts_sub_arom = dict((sub,subsystem_arom.count(sub)) for sub in all_subsystems)
        counts_sub_arom_uni = dict((sub,subsystem_arom_uni.count(sub)) for sub in all_subsystems)
    counts_sub_smeg = dict((sub,subsystem_smeg.count(sub)) for sub in all_subsystems)
    counts_sub_mtb = dict((sub,subsystem_mtb.count(sub)) for sub in all_subsystems)

    mari_ratios = {sub: counts_sub_mari[sub] / counts_all[sub] if sub in counts_all and counts_all[sub] != 0 else 0 
          for sub in counts_sub_mari}
    absc_ratios = {sub: counts_sub_absc[sub] / counts_all[sub] if sub in counts_all and counts_all[sub] != 0 else 0 
          for sub in counts_sub_absc}
    if arom == 'arom':
        arom_ratios = {sub: counts_sub_arom[sub] / counts_all[sub] if sub in counts_all and counts_all[sub] != 0 else 0 
              for sub in counts_sub_arom}
    smeg_ratios = {sub: counts_sub_smeg[sub] / counts_all[sub] if sub in counts_all and counts_all[sub] != 0 else 0 
          for sub in counts_sub_smeg}
    mtb_ratios = {sub: counts_sub_mtb[sub] / counts_all[sub] if sub in counts_all and counts_all[sub] != 0 else 0 
          for sub in counts_sub_mtb}
    
    # ratios for different reactions with respect to total reactions in subsystem in each species
    mari_dif_ratios = {sub: counts_dif_sub_mari[sub] / counts_sub_mari[sub] if sub in counts_sub_mari and counts_sub_mari[sub] != 0 else 0 
          for sub in counts_sub_mari}
    absc_dif_ratios = {sub: counts_dif_sub_absc[sub] / counts_sub_absc[sub] if sub in counts_sub_absc and counts_sub_absc[sub] != 0 else 0 
          for sub in counts_sub_absc}
    if arom == 'arom':
        arom_dif_ratios = {sub: counts_dif_sub_arom[sub] / counts_sub_arom[sub] if sub in counts_sub_arom and counts_sub_arom[sub] != 0 else 0 
              for sub in counts_sub_arom}
    smeg_dif_ratios = {sub: counts_dif_sub_smeg[sub] / counts_sub_smeg[sub] if sub in counts_sub_smeg and counts_sub_smeg[sub] != 0 else 0 
          for sub in counts_sub_smeg}
    mtb_dif_ratios = {sub: counts_dif_sub_mtb[sub] / counts_sub_mtb[sub] if sub in counts_sub_mtb and counts_sub_mtb[sub] != 0 else 0 
          for sub in counts_sub_mtb}

    for subsystem in counts_all:
        counts_all[subsystem] = {"count": counts_all[subsystem], "reactions": []}

    for subsystem in counts_sub_arom_uni:
        counts_sub_arom_uni[subsystem] = {"count": counts_sub_arom_uni[subsystem], "reactions": []}

    for reaction, subsystem in zip(data_all['reaction'], data_all['subsystem']):
        if subsystem in counts_all:
            counts_all[subsystem]["reactions"].append(reaction)
    
    for reaction, subsystem in zip(data_aromaticivorans['reaction'], data_aromaticivorans['subsystem']):
        if subsystem in counts_sub_arom_uni:
            counts_sub_arom_uni[subsystem]["reactions"].append(reaction)


    data_to_csv = {
        "Subsystem": [],
        "Reaction Count": [],
        "Reactions": []
    }
    for subsystem, details in counts_all.items(): #  for all counts_sub_arom_uni
        data_to_csv["Subsystem"].append(subsystem)
        data_to_csv["Reaction Count"].append(details["count"])
        data_to_csv["Reactions"].append("; ".join(details["reactions"]))  # Join reactions list into a single string

    # convert to DataFrame
    
    df = pd.DataFrame(data_to_csv)
    csv_file = f'subsystems_reactions_details_{model_type}_{arom}_{date}_{level_subsystem}.csv'
    df.to_csv(csv_file, index=False)
    # count occurrences of genes in each set
    counts_gene_inter = data_gene_intersection['gene'].count(1)
    counts_gene_dif = data_gene_different['gene'].count(1)
    not_in_gene_inter = len(data_gene_intersection['gene']) - counts_gene_inter
    not_in_gene_dif = len(data_gene_different['gene']) - counts_gene_dif

    #plt.figure(figsize=(8,10))
    mari_sub = pd.DataFrame(list(mari_ratios.values()), index=mari_ratios.keys(), columns=['$M. marinum$'])
    absc_sub = pd.DataFrame(list(absc_ratios.values()), index=absc_ratios.keys(), columns=['$M. abscessus$'])
    if arom == 'arom':
        arom_sub = pd.DataFrame(list(arom_ratios.values()), index=arom_ratios.keys(), columns=['$M. aromaticivorans$'])
    smeg_sub = pd.DataFrame(list(smeg_ratios.values()), index=smeg_ratios.keys(), columns=['$M. smegmatis$'])
    mtb_sub = pd.DataFrame(list(mtb_ratios.values()), index=mtb_ratios.keys(), columns=['$M. tuberculosis$'])

    mari_dif_sub = pd.DataFrame(list(mari_dif_ratios.values()), index=mari_dif_ratios.keys(), columns=['$M. marinum$'])
    absc_dif_sub = pd.DataFrame(list(absc_dif_ratios.values()), index=absc_dif_ratios.keys(), columns=['$M. abscessus$'])
    if arom == 'arom':
        arom_dif_sub = pd.DataFrame(list(arom_dif_ratios.values()), index=arom_dif_ratios.keys(), columns=['$M. aromaticivorans$'])
    smeg_dif_sub = pd.DataFrame(list(smeg_dif_ratios.values()), index=smeg_dif_ratios.keys(), columns=['$M. smegmatis$'])
    mtb_dif_sub = pd.DataFrame(list(mtb_dif_ratios.values()), index=mtb_dif_ratios.keys(), columns=['$M. tuberculosis$'])
    
    #    combined_df = pd.concat([mari_sub, absc_sub, arom_sub,smeg_sub, mtb_sub], axis=1)
    if arom == 'arom':
        combined_df = pd.concat([mari_sub, absc_sub, smeg_sub, mtb_sub,arom_sub], axis=1) 
        combined_dif_sub_df = pd.concat([mari_dif_sub, absc_dif_sub, smeg_dif_sub, mtb_dif_sub,arom_dif_sub], axis=1) 
    else:
        combined_df = pd.concat([mari_sub, absc_sub, smeg_sub, mtb_sub], axis=1) 
        combined_dif_sub_df = pd.concat([mari_dif_sub, absc_dif_sub, smeg_dif_sub, mtb_dif_sub], axis=1) 


    clustermap = sns.clustermap(combined_dif_sub_df, cmap='YlOrRd', annot=False, cbar=True, linewidths=0.5, col_cluster=False, row_cluster=True, figsize= (10,35)) #combined_df
    plt.setp(clustermap.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='right')
    plt.savefig(f'subsystem_heatmap_dif_{model_type}_{arom}.png', bbox_inches = 'tight', dpi = 300)
    
    
    ## Fisher test for genes
            
    if fisherTest:
        contingency_table_gene = [[counts_gene_inter, counts_gene_dif],
                                [not_in_gene_inter, not_in_gene_dif]]

        # Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(contingency_table_gene) #alternative = 'two-sided')
        print(odds_ratio, p_value)
        # p-value 0.093
        
        # Perform Fisher's exact test for each subsystem
        if arom == 'arom':
            counts_int_sub_myco = [counts_int_sub_absc,counts_int_sub_mari,counts_int_sub_smeg, counts_int_sub_arom,counts_int_sub_mtb]
            counts_dif_sub_myco = [counts_dif_sub_absc,counts_dif_sub_mari,counts_dif_sub_smeg, counts_dif_sub_arom,counts_dif_sub_mtb]
            subsystem_dif_all = [subsystem_dif_absc,subsystem_dif_mari,subsystem_dif_smeg,subsystem_dif_arom, subsystem_dif_mtb]
            #print(counts_sub_myco)
            sub_myco = [subsystem_absc, subsystem_mari,subsystem_smeg,subsystem_arom,subsystem_mtb]
            #print(sub_myco)
        else:
            counts_sub_myco = [counts_sub_absc,counts_sub_mari,counts_sub_smeg,counts_sub_mtb] #counts_sub_arom
            sub_myco = [subsystem_absc, subsystem_mari,subsystem_smeg,subsystem_mtb] #subsystem_arom

        for idx, species in enumerate(model_list):  #or for all model_list 

            globals()[f'results_{species_name[idx]}_less'] = []
            globals()[f'results_{species_name[idx]}_greater'] = []

            for subsystem in all_subsystems:         # remove exchange
                if subsystem != 'Exchange':
                    
                    # Construct contingency table
                    count_inter = counts_int_sub_myco[idx][subsystem]
                    count_dif = counts_dif_sub_myco[idx][subsystem]
                    not_in_inter = len(data_intersection['subsystem']) - count_inter
                    not_in_dif = len(subsystem_dif_all[idx]) - count_dif

                    contingency_table = [[count_inter, count_dif],
                                        [not_in_inter, not_in_dif]]
                    print(subsystem,contingency_table)

                    # Perform Fisher's exact test "less" for odds ratio <1 (more represented in dif reactions) and "greater" for odds ratio >1 
                    odds_ratio_less, p_value_less = fisher_exact(contingency_table, alternative = 'less')
                    odds_ratio_greater, p_value_greater = fisher_exact(contingency_table, alternative = 'greater')
                    globals()[f'results_{species_name[idx]}_less'].append({'subsystem': subsystem, 'odds_ratio': odds_ratio_less, 'p_value': p_value_less, 'contingency_table':contingency_table,'reactions_in_shared':count_inter, 'reactions_in_non_shared':count_dif})

                    globals()[f'results_{species_name[idx]}_greater'].append({'subsystem': subsystem, 'odds_ratio': odds_ratio_greater, 'p_value': p_value_greater, 'contingency_table':contingency_table,'reactions_in_shared':count_inter, 'reactions_in_non_shared':count_dif})

        # Convert results to a DataFrame for better visualization
        results_ab_df_less = pd.DataFrame(results_abscessus_less)
        results_ma_df_less = pd.DataFrame(results_marinum_less)
        results_sm_df_less = pd.DataFrame(results_smegmatis_less)
        results_rv_df_less = pd.DataFrame(results_Rv_less)
        results_ab_df_greater = pd.DataFrame(results_abscessus_greater)
        results_ma_df_greater = pd.DataFrame(results_marinum_greater)
        results_sm_df_greater = pd.DataFrame(results_smegmatis_greater)
        results_rv_df_greater = pd.DataFrame(results_Rv_greater)

        if arom == 'arom':
            results_ar_df_less = pd.DataFrame(results_aromaticivorans_less)
            results_ar_df_greater = pd.DataFrame(results_aromaticivorans_greater)

        for specie, result in zip(['ab_less','ma_less','sm_less','rv_less','ab_great','ma_great','sm_great','rv_great'],[results_ab_df_less,results_ma_df_less,results_sm_df_less,results_rv_df_less,results_ab_df_greater,results_ma_df_greater,results_sm_df_greater,results_rv_df_greater]):
            result.to_csv(f'Fisher_test_{specie}_{date}_{model_type}_{arom}_{level_subsystem}.csv')
        if arom == 'arom':    
            results_ar_df_less.to_csv(f'Fisher_test_ar_less_{date}_{model_type}_{arom}_{level_subsystem}.csv')
            results_ar_df_greater.to_csv(f'Fisher_test_ar_great_{date}_{model_type}_{arom}_{level_subsystem}.csv')
     
        # merging two levels for less 
        results_ab_df_less_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ab_less_{date}_{model_type}_{arom}_high.csv')
        results_ab_df_less_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ab_less_{date}_{model_type}_{arom}_detailed.csv')
        results_ar_df_less_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ar_less_{date}_{model_type}_{arom}_high.csv')
        results_ar_df_less_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ar_less_{date}_{model_type}_{arom}_detailed.csv')
        results_sm_df_less_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_sm_less_{date}_{model_type}_{arom}_high.csv')
        results_sm_df_less_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_sm_less_{date}_{model_type}_{arom}_detailed.csv')
        results_ma_df_less_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ma_less_{date}_{model_type}_{arom}_high.csv')
        results_ma_df_less_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ma_less_{date}_{model_type}_{arom}_detailed.csv')
        results_rv_df_less_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_rv_less_{date}_{model_type}_{arom}_high.csv')
        results_rv_df_less_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_rv_less_{date}_{model_type}_{arom}_detailed.csv')
        
        # merging two levels for greater 
        results_ab_df_great_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ab_great_{date}_{model_type}_{arom}_high.csv')
        results_ab_df_great_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ab_great_{date}_{model_type}_{arom}_detailed.csv')
        results_ar_df_great_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ar_great_{date}_{model_type}_{arom}_high.csv')
        results_ar_df_great_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ar_great_{date}_{model_type}_{arom}_detailed.csv')
        results_sm_df_great_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_sm_great_{date}_{model_type}_{arom}_high.csv')
        results_sm_df_great_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_sm_great_{date}_{model_type}_{arom}_detailed.csv')
        results_ma_df_great_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ma_great_{date}_{model_type}_{arom}_high.csv')
        results_ma_df_great_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_ma_great_{date}_{model_type}_{arom}_detailed.csv')
        results_rv_df_great_high = pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_rv_great_{date}_{model_type}_{arom}_high.csv')
        results_rv_df_great_detailed =  pd.read_csv(f'/home/icancino/Documents/Carveme/Genomes/ModelsRef1241024Egg/Fisher_test_rv_great_{date}_{model_type}_{arom}_detailed.csv')
        

        merge_ab_detail_high_pvalue_less  = pd.concat([results_ab_df_less_high,results_ab_df_less_detailed], ignore_index=True)        
        merge_ar_detail_high_pvalue_less  = pd.concat([results_ar_df_less_high,results_ar_df_less_detailed], ignore_index=True)      
        merge_sm_detail_high_pvalue_less  = pd.concat([results_sm_df_less_high,results_sm_df_less_detailed], ignore_index=True)
        merge_ma_detail_high_pvalue_less  = pd.concat([results_ma_df_less_high,results_ma_df_less_detailed], ignore_index=True)
        merge_rv_detail_high_pvalue_less  = pd.concat([results_rv_df_less_high,results_rv_df_less_detailed], ignore_index=True)        

        merge_ab_detail_high_pvalue_great  = pd.concat([results_ab_df_great_high,results_ab_df_great_detailed], ignore_index=True)        
        merge_ar_detail_high_pvalue_great  = pd.concat([results_ar_df_great_high,results_ar_df_great_detailed], ignore_index=True)      
        merge_sm_detail_high_pvalue_great  = pd.concat([results_sm_df_great_high,results_sm_df_great_detailed], ignore_index=True)
        merge_ma_detail_high_pvalue_great  = pd.concat([results_ma_df_great_high,results_ma_df_great_detailed], ignore_index=True)
        merge_rv_detail_high_pvalue_great  = pd.concat([results_rv_df_great_high,results_rv_df_great_detailed], ignore_index=True)        


        ## Multitest with p-value for less and greater 
       
        merge_ab_detail_high_pvalue_less['adjusted_p_value'] = false_discovery_control(np.array(merge_ab_detail_high_pvalue_less['p_value']),axis=None, method='bh')
        merge_ma_detail_high_pvalue_less['adjusted_p_value'] = false_discovery_control(np.array(merge_ma_detail_high_pvalue_less['p_value']),axis=None, method='bh')
        merge_sm_detail_high_pvalue_less['adjusted_p_value'] = false_discovery_control(np.array(merge_sm_detail_high_pvalue_less['p_value']),axis=None, method='bh')
        merge_rv_detail_high_pvalue_less['adjusted_p_value'] = false_discovery_control(np.array(merge_rv_detail_high_pvalue_less['p_value']),axis=None, method='bh')
        
        merge_ab_detail_high_pvalue_great['adjusted_p_value'] = false_discovery_control(np.array(merge_ab_detail_high_pvalue_great['p_value']),axis=None, method='bh')
        merge_ma_detail_high_pvalue_great['adjusted_p_value'] = false_discovery_control(np.array(merge_ma_detail_high_pvalue_great['p_value']),axis=None, method='bh')
        merge_sm_detail_high_pvalue_great['adjusted_p_value'] = false_discovery_control(np.array(merge_sm_detail_high_pvalue_great['p_value']),axis=None, method='bh')
        merge_rv_detail_high_pvalue_great['adjusted_p_value'] = false_discovery_control(np.array(merge_rv_detail_high_pvalue_great['p_value']),axis=None, method='bh')
        
        if arom == 'arom':
            merge_ar_detail_high_pvalue_less['adjusted_p_value'] = false_discovery_control(np.array(merge_ar_detail_high_pvalue_less['p_value']),axis=None, method='bh')
            merge_ar_detail_high_pvalue_great['adjusted_p_value'] = false_discovery_control(np.array(merge_ar_detail_high_pvalue_great['p_value']),axis=None, method='bh')

        
        significant_results_ab_df = merge_ab_detail_high_pvalue_less.sort_values(by=['adjusted_p_value'])
        significant_results_ab_df.to_csv(f'Significant_results_abscessus_less_{date}_{model_type}_{arom}_high_low.csv')
        significant_results_ab_df = merge_ab_detail_high_pvalue_great.sort_values(by=['adjusted_p_value'])
        significant_results_ab_df.to_csv(f'Significant_results_abscessus_great_{date}_{model_type}_{arom}_high_low.csv')


        significant_results_ma_df = merge_ma_detail_high_pvalue_less.sort_values(by=['adjusted_p_value'])
        significant_results_ma_df.to_csv(f'Significant_results_marinum_less_{date}_{model_type}_{arom}_high_low.csv')
        significant_results_ma_df = merge_ma_detail_high_pvalue_great.sort_values(by=['adjusted_p_value'])
        significant_results_ma_df.to_csv(f'Significant_results_marinum_great_{date}_{model_type}_{arom}_high_low.csv')
        

        significant_results_sm_df = merge_sm_detail_high_pvalue_less.sort_values(by=['adjusted_p_value'])
        significant_results_sm_df.to_csv(f'Significant_results_smegmatis_less_{date}_{model_type}_{arom}_high_low.csv')
        significant_results_sm_df = merge_sm_detail_high_pvalue_great.sort_values(by=['adjusted_p_value'])
        significant_results_sm_df.to_csv(f'Significant_results_smegmatis_great_{date}_{model_type}_{arom}_high_low.csv')
        

        significant_results_rv_df = merge_rv_detail_high_pvalue_less.sort_values(by=['adjusted_p_value'])
        significant_results_rv_df.to_csv(f'Significant_results_mtb_less_{date}_{model_type}_{arom}_high_low.csv')
        significant_results_rv_df = merge_rv_detail_high_pvalue_great.sort_values(by=['adjusted_p_value'])
        significant_results_rv_df.to_csv(f'Significant_results_mtb_great_{date}_{model_type}_{arom}_high_low.csv')
        
        if arom == 'arom':
            significant_results_ar_df = merge_ar_detail_high_pvalue_less.sort_values(by=['adjusted_p_value'])
            significant_results_ar_df.to_csv(f'Significant_results_aromaticivorans_less_{date}_{model_type}_{arom}_high_low.csv')
            significant_results_ar_df = merge_ar_detail_high_pvalue_great.sort_values(by=['adjusted_p_value'])
            significant_results_ar_df.to_csv(f'Significant_results_aromaticivorans_great_{date}_{model_type}_{arom}_high_low.csv')
        
        
        
        dict_list_sim = []
        row_dict_sim = {}
        dict_list_dif = []
        row_dict_dif = {}
        
        for rxn in Rv.reactions:
            if rxn.subsystem in list(significant_results_rv_df['subsystem']) and rxn.id in [reaction.name for reaction in similar_reactions]:
                row_dict_sim = {'Subsystem':rxn.subsystem, 'Similar Reaction':'Yes','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_sim.append(row_dict_sim)

            elif rxn.subsystem in list(significant_results_rv_df['subsystem']) and rxn.id not in [reaction.name for reaction in similar_reactions]:
                #if model == model_list_noarom[0]:
                row_dict_dif = {'Subsystem':rxn.subsystem, 'Similar Reaction':'No','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_dif.append(row_dict_dif)
        for rxn in abscessus.reactions:
            if rxn.subsystem in list(significant_results_ab_df['subsystem']) and rxn.id in [reaction.name for reaction in similar_reactions]:
                row_dict_sim = {'Subsystem':rxn.subsystem, 'Similar Reaction':'Yes','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_sim.append(row_dict_sim)

            elif rxn.subsystem in list(significant_results_ab_df['subsystem']) and rxn.id not in [reaction.name for reaction in similar_reactions]:
                #if model == model_list_noarom[0]:
                row_dict_dif = {'Subsystem':rxn.subsystem, 'Similar Reaction':'No','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_dif.append(row_dict_dif)
    
        for rxn in smegmatis.reactions:
            if rxn.subsystem in list(significant_results_sm_df['subsystem']) and rxn.id in [reaction.name for reaction in similar_reactions]:
                row_dict_sim = {'Subsystem':rxn.subsystem, 'Similar Reaction':'Yes','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_sim.append(row_dict_sim)

            elif rxn.subsystem in list(significant_results_sm_df['subsystem']) and rxn.id not in [reaction.name for reaction in similar_reactions]:
                row_dict_dif = {'Subsystem':rxn.subsystem, 'Similar Reaction':'No','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_dif.append(row_dict_dif)
    
        for rxn in marinum.reactions:
            if rxn.subsystem in list(significant_results_ma_df['subsystem']) and rxn.id in [reaction.name for reaction in similar_reactions]:
                row_dict_sim = {'Subsystem':rxn.subsystem, 'Similar Reaction':'Yes','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_sim.append(row_dict_sim)

            elif rxn.subsystem in list(significant_results_ma_df['subsystem']) and rxn.id not in [reaction.name for reaction in similar_reactions]:
                #if model == model_list_noarom[0]:
                row_dict_dif = {'Subsystem':rxn.subsystem, 'Similar Reaction':'No','Reaction': rxn.id, 'Name':rxn.name}
                dict_list_dif.append(row_dict_dif)
        
        if arom == 'arom':
            for rxn in aromaticivorans.reactions:
                if rxn.subsystem in list(significant_results_ar_df['subsystem']) and rxn.id in [reaction.name for reaction in similar_reactions]:
                    row_dict_sim = {'Subsystem':rxn.subsystem, 'Similar Reaction':'Yes','Reaction': rxn.id, 'Name':rxn.name}
                    dict_list_sim.append(row_dict_sim)
                elif rxn.subsystem in list(significant_results_ar_df['subsystem']) and rxn.id not in [reaction.name for reaction in similar_reactions]:
                    #if model == model_list_noarom[0]:
                    row_dict_dif = {'Subsystem':rxn.subsystem, 'Similar Reaction':'No','Reaction': rxn.id, 'Name':rxn.name}
                    dict_list_dif.append(row_dict_dif)
    
        if arom == 'arom':
            data = {
                'Reaction': [reaction.name for reaction in different_reactions],
                'Marinum': [1 if reaction.name in [x.name for x in mari_rxns] else 0 for reaction in different_reactions],
                'Smegmatis': [1 if reaction.name in [x.name for x in smeg_rxns] else 0 for reaction in different_reactions],
                'Aromaticivorans': [1 if reaction.name in [x.name for x in arom_rxns] else 0 for reaction in different_reactions],
                'Abscessus': [1 if reaction.name in [x.name for x in absc_rxns] else 0 for reaction in different_reactions],
                'MTB': [1 if reaction.name in [x.name for x in rv_rxns] else 0 for reaction in different_reactions]
                    }
            data_subsys = {
                'Reaction': [reaction.subsystem for reaction in different_reactions],
                'Marinum': [1 if reaction.subsystem in [x.subsystem for x in mari_rxns] else 0 for reaction in different_reactions],
                'Smegmatis': [1 if reaction.subsystem in [x.subsystem for x in smeg_rxns] else 0 for reaction in different_reactions],
                'Aromaticivorans': [1 if reaction.subsystem in [x.subsystem for x in arom_rxns] else 0 for reaction in different_reactions],
                'Abscessus': [1 if reaction.subsystem in [x.subsystem for x in absc_rxns] else 0 for reaction in different_reactions],
                'MTB': [1 if reaction.subsystem in [x.subsystem for x in rv_rxns] else 0 for reaction in different_reactions]
            }
        else:
            data = {
                'Reaction': [reaction.name for reaction in different_reactions],
                'Marinum': [1 if reaction.name in [x.name for x in mari_rxns] else 0 for reaction in different_reactions],
                'Smegmatis': [1 if reaction.name in [x.name for x in smeg_rxns] else 0 for reaction in different_reactions],
                #'Aromaticivorans': [1 if reaction.name in [x.name for x in arom_rxns] else 0 for reaction in different_reactions],
                'Abscessus': [1 if reaction.name in [x.name for x in absc_rxns] else 0 for reaction in different_reactions],
                'MTB': [1 if reaction.name in [x.name for x in rv_rxns] else 0 for reaction in different_reactions]
                    }
            data_subsys = {
                'Reaction': [reaction.subsystem for reaction in different_reactions],
                'Marinum': [1 if reaction.subsystem in [x.subsystem for x in mari_rxns] else 0 for reaction in different_reactions],
                'Smegmatis': [1 if reaction.subsystem in [x.subsystem for x in smeg_rxns] else 0 for reaction in different_reactions],
                #'Aromaticivorans': [1 if reaction.subsystem in [x.subsystem for x in arom_rxns] else 0 for reaction in different_reactions],
                'Abscessus': [1 if reaction.subsystem in [x.subsystem for x in absc_rxns] else 0 for reaction in different_reactions],
                'MTB': [1 if reaction.subsystem in [x.subsystem for x in rv_rxns] else 0 for reaction in different_reactions]

                    }
        data_rxn_dif = {
                'Reaction': [reaction.name for reaction in different_reactions],
                'Marinum': [reaction.name if reaction.name in [x.name for x in mari_rxns] else 0 for reaction in different_reactions],
                'Smegmatis': [reaction.name if reaction.name in [x.name for x in smeg_rxns] else 0 for reaction in different_reactions],
                #'Aromaticivorans': [reaction.name if reaction.name in [x.name for x in arom_rxns] else 0 for reaction in different_reactions],
                'Abscessus': [reaction.name if reaction.name in [x.name for x in absc_rxns] else 0 for reaction in different_reactions],
                'MTB': [reaction.name if reaction.name in [x.name for x in rv_rxns] else 0 for reaction in different_reactions]
                    }
        
        Dif_rxns_location_filtered = pd.DataFrame(data)
        Dif_rxns_location_filtered = Dif_rxns_location_filtered[~Dif_rxns_location_filtered['Reaction'].str.contains('EX_')]
        Dif_rxns_location_filtered.to_csv(f'all_dif_rxn_{model_type}_{arom}_noex.csv')
        Dif_rxns_location_myco = pd.DataFrame(data_rxn_dif)
        path_rxns_location = pd.DataFrame({})
        non_path_rxns_location = pd.DataFrame({})
        slow_rxns_location = pd.DataFrame({})
        fast_rxns_location = pd.DataFrame({})
        
        Dif_ss_location = pd.DataFrame(data_subsys)
        if model_type == 'max_gf':
            filter_patho = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']))})
            filter_slow = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
            filter_env = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))})
            filter_fast = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']))})
            if arom == 'arom':
                Dif_rxns_location_filtered['Count'] = Dif_rxns_location_filtered[['Marinum', 'Smegmatis', 'Aromaticivorans' ,'Abscessus','MTB']].sum(axis=1) 
            else: 
                Dif_rxns_location_filtered['Count'] = Dif_rxns_location_filtered[['Marinum', 'Smegmatis', 'Abscessus','MTB']].sum(axis=1) #, 'Aromaticivorans' (before abscessus)
            
            if arom == 'arom':

                filter_patho['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction'])))
                filter_patho['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)& (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))
                filter_patho['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)& (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                filter_slow['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_slow['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)& (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))
                filter_slow['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)& (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                
                filter_fast['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction'])))
                filter_fast['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1)& (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))
                filter_fast['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1)& (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))

                
                filter_env['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) + [np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) &  (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))
                filter_env['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0)& (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))
                filter_env['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0)& (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))

                filter_patho_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0)] #& (Dif_ss_location['Aromaticivorans'] == 0)
                filter_slow_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0) ]
                filter_fast_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1) ]
                filter_env_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1) ]
                filter_tb_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0)]
                filter_ntm_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)& (Dif_ss_location['Aromaticivorans'] == 1)]


                #filter_patho['Subsystem'] =  filter_patho['Reaction'].apply(lambda x: marinum.reactions.get_by_id(str(x)).subsystem , axis=1)
                #filter_patho['Name'] = marinum.reactions.get_by_id(str(filter_patho['Reaction'])).name
                #filter_patho['Reaction'] = list(filter_patho['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum']==1]) - len(filter_patho['Reaction']))

                filter_patho['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                
                filter_patho['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                
                
                filter_patho['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                
                print('PATHO',filter_patho)
                #filter_env['Reaction'] = list(filter_env['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_env['Reaction']))
                filter_env['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                filter_env['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))

                print('NONPATHO',filter_env)

                #filter_slow['Subsystem'] = marinum.reactions.get_by_id(filter_slow['Reaction']).subsystem
                #filter_slow['Name'] = marinum.reactions.get_by_id(filter_slow['Reaction']).name
                filter_slow['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                filter_slow['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                print('SLOW',filter_slow)
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_fast['Reaction']))
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + (([np.nan] * (len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])- len(filter_fast))))
                filter_fast['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']]  + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
    
                filter_fast['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                
                filter_fast['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                print('FAST',filter_fast )
                
            else:
                
                filter_patho['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction'])))
                filter_patho['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))
                filter_patho['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                filter_slow['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_slow['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))
                filter_slow['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                
                filter_fast['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction'])))
                filter_fast['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))
                filter_fast['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))

                
                filter_env['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) + [np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))
                filter_env['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))
                filter_env['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']])))
                
                #filter_tb['Reaction'] = Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                #    (Dif_rxns_location_filtered['Smegmatis'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)]['Reaction']
                #filter_ntm['Reaction'] = Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                #    (Dif_rxns_location_filtered['Smegmatis'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction']


                filter_patho_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)] #& (Dif_ss_location['Aromaticivorans'] == 0)
                filter_slow_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)]
                filter_fast_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)]
                filter_env_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)]
                filter_tb_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)]
                filter_ntm_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)]


                #filter_patho['Subsystem'] =  filter_patho['Reaction'].apply(lambda x: marinum.reactions.get_by_id(str(x)).subsystem , axis=1)
                #filter_patho['Name'] = marinum.reactions.get_by_id(str(filter_patho['Reaction'])).name
                #filter_patho['Reaction'] = list(filter_patho['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum']==1]) - len(filter_patho['Reaction']))

                filter_patho['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                
                filter_patho['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                
                
                filter_patho['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                
                print('PATHO',filter_patho)
                #filter_env['Reaction'] = list(filter_env['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_env['Reaction']))
                filter_env['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
        
                print('NONPATHO',filter_env)

                #filter_slow['Subsystem'] = marinum.reactions.get_by_id(filter_slow['Reaction']).subsystem
                #filter_slow['Name'] = marinum.reactions.get_by_id(filter_slow['Reaction']).name
                filter_slow['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                filter_slow['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                print('SLOW',filter_slow)
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_fast['Reaction']))
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + (([np.nan] * (len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])- len(filter_fast))))
                filter_fast['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']]  + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                filter_fast['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                print('FAST',filter_fast )
                
        elif model_type == 'min_gf':  
            if arom == 'arom':  
                filter_patho = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
                filter_slow = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
                filter_env = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))})
                filter_fast = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))})

                Dif_rxns_location_filtered['Count'] = Dif_rxns_location_filtered[['Marinum', 'Smegmatis', 'Aromaticivorans', 'Abscessus','MTB']].sum(axis=1) 
                
                filter_patho['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_patho['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))
                filter_patho['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1)& (Dif_rxns_location_filtered['Aromaticivorans'] == 0)  & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1)& (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                filter_slow['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_slow['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)  & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))
                filter_slow['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                
                filter_fast['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 1) &(Dif_rxns_location_filtered['Aromaticivorans'] == 1)  ]['Reaction'])))
                filter_fast['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_fast['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))

                
                filter_env['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + [np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0)  & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))
                filter_env['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) &(Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_env['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) &(Dif_rxns_location_filtered['Aromaticivorans'] == 1)& (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))


                filter_patho_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0)]
                filter_slow_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0)]
                filter_fast_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1)]
                filter_env_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1)]
                filter_tb_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0)]
                filter_ntm_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1)]

                filter_patho['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                
                filter_patho['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                
                
                filter_patho['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                
                print('PATHO',filter_patho)
                #filter_env['Reaction'] = list(filter_env['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_env['Reaction']))
                filter_env['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                filter_env['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))


                print('NONPATHO',filter_env)

                #filter_slow['Subsystem'] = marinum.reactions.get_by_id(filter_slow['Reaction']).subsystem
                #filter_slow['Name'] = marinum.reactions.get_by_id(filter_slow['Reaction']).name
                filter_slow['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                filter_slow['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                print('SLOW',filter_slow)
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_fast['Reaction']))
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + (([np.nan] * (len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])- len(filter_fast))))
                filter_fast['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']]  + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                
                filter_fast['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                
                filter_fast['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                print('FAST',filter_fast )


            else:
                filter_patho = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
                filter_slow = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
                filter_env = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))})
                filter_fast = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))})

                Dif_rxns_location_filtered['Count'] = Dif_rxns_location_filtered[['Marinum', 'Smegmatis', 'Abscessus','MTB']].sum(axis=1) 
                
                filter_patho['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_patho['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))
                filter_patho['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                filter_slow['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_slow['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))
                filter_slow['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']])))

                
                filter_fast['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction'])))
                filter_fast['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_fast['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))

                
                filter_env['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + [np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']))
                filter_env['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_env['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))


                filter_patho_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)]
                filter_slow_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)]
                filter_fast_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)]
                filter_env_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)]
                filter_tb_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)]
                filter_ntm_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)]


                #filter_patho['Subsystem'] =  filter_patho['Reaction'].apply(lambda x: marinum.reactions.get_by_id(str(x)).subsystem , axis=1)
                #filter_patho['Name'] = marinum.reactions.get_by_id(str(filter_patho['Reaction'])).name
                #filter_patho['Reaction'] = list(filter_patho['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum']==1]) - len(filter_patho['Reaction']))

                filter_patho['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                
                filter_patho['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                
                
                filter_patho['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                
                print('PATHO',filter_patho)
                #filter_env['Reaction'] = list(filter_env['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_env['Reaction']))
                filter_env['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                #filter_env['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                #filter_env['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                #filter_env['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                #filter_env['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))


                print('NONPATHO',filter_env)

                #filter_slow['Subsystem'] = marinum.reactions.get_by_id(filter_slow['Reaction']).subsystem
                #filter_slow['Name'] = marinum.reactions.get_by_id(filter_slow['Reaction']).name
                filter_slow['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                filter_slow['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                print('SLOW',filter_slow)
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_fast['Reaction']))
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + (([np.nan] * (len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])- len(filter_fast))))
                filter_fast['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']]  + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                
                #filter_fast['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                #filter_fast['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                #filter_fast['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                #filter_fast['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                
                filter_fast['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                print('FAST',filter_fast )
                

        elif model_type == 'max_genes_perc_gf':
            if arom == 'arom':
                filter_patho = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']))})
                filter_slow = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
                filter_env = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction']))})
                filter_fast = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction']))})

                Dif_rxns_location_filtered['Count'] = Dif_rxns_location_filtered[['Marinum', 'Smegmatis', 'Aromaticivorans','Abscessus','MTB']].sum(axis=1) #'Aromaticivorans'
                
                filter_patho['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)& (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_patho['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1)& (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0)& (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))
                filter_patho['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))

                filter_slow['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction'])))
                filter_slow['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))
                filter_slow['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1)& (Dif_rxns_location_filtered['Aromaticivorans'] == 0)  & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))

                
                filter_fast['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 1)  & (Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction'])))
                filter_fast['Common Name'] = [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0)& (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_fast['Common Subsystem'] = [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1)& (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))

                
                filter_env['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + [np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Aromaticivorans'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']))
                filter_env['Common Name'] = [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Aromaticivorans'] == 1)& (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_env['Common Subsystem'] = [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0)& (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0)& (Dif_rxns_location_filtered['Aromaticivorans'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))

                filter_patho_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0) ]
                filter_slow_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0) ]
                filter_fast_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1) ]
                filter_env_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1)]
                filter_tb_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0) & (Dif_ss_location['Aromaticivorans'] == 0)]
                filter_ntm_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1) & (Dif_ss_location['Aromaticivorans'] == 1)]

                filter_patho['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                
                filter_patho['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                
                
                filter_patho['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                
                print('PATHO',filter_patho)
                #filter_env['Reaction'] = list(filter_env['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_env['Reaction']))
                
                filter_env['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_env['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_env) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))

                filter_env['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                print('NONPATHO',filter_env)

                #filter_slow['Subsystem'] = marinum.reactions.get_by_id(filter_slow['Reaction']).subsystem
                #filter_slow['Name'] = marinum.reactions.get_by_id(filter_slow['Reaction']).name
                filter_slow['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                filter_slow['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                print('SLOW',filter_slow)
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_fast['Reaction']))
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + (([np.nan] * (len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])- len(filter_fast))))
                filter_fast['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']]  + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                
                filter_fast['Aromaticivorans']= [(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Name']= [(aromaticivorans.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Subsystem']= [(aromaticivorans.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                filter_fast['Ar. Gene']= [(aromaticivorans.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_fast) - len([(aromaticivorans.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Aromaticivorans'] ==1 ]['Reaction']])))
                
                filter_fast['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                print('FAST',filter_fast )
            else:
                filter_patho = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']))})
                filter_slow = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']))})
                filter_env = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']))})
                filter_fast = pd.DataFrame({'Common Reaction' : range(len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']))})

                Dif_rxns_location_filtered['Count'] = Dif_rxns_location_filtered[['Marinum', 'Smegmatis', 'Abscessus','MTB']].sum(axis=1) #'Aromaticivorans'
                
                filter_patho['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction'])))
                filter_patho['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0)]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))
                filter_patho['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))

                filter_slow['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction'])))
                filter_slow['Common Name'] = [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))
                filter_slow['Common Subsystem'] = [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 1) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 1) & (Dif_rxns_location_filtered['Smegmatis'] == 0) ]['Reaction']])))

                
                filter_fast['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + ([np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Abscessus'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                        (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction'])))
                filter_fast['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_fast['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 1) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))

                
                filter_env['Common Reaction'] = list(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']) + [np.nan] * (len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Smegmatis'] == 1)]['Reaction']) - len(Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & 
                    (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']))
                filter_env['Common Name'] = [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))
                filter_env['Common Subsystem'] = [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[(Dif_rxns_location_filtered['Marinum'] == 0) & (Dif_rxns_location_filtered['Abscessus'] == 0) & (Dif_rxns_location_filtered['MTB'] == 0) & (Dif_rxns_location_filtered['Smegmatis'] == 1) ]['Reaction']])))

                filter_patho_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)  ]
                filter_slow_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)  ]
                filter_fast_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)  ]
                filter_env_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)  ]
                filter_tb_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 0) & (Dif_ss_location['Abscessus'] == 0) & (Dif_ss_location['MTB'] == 1) & 
                    (Dif_ss_location['Smegmatis'] == 0)  ]
                filter_ntm_ss = Dif_ss_location[(Dif_ss_location['Marinum'] == 1) & (Dif_ss_location['Abscessus'] == 1) & (Dif_ss_location['MTB'] == 0) & 
                    (Dif_ss_location['Smegmatis'] == 1)  ]

                filter_patho['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_patho['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                
                filter_patho['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_patho['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                
                
                filter_patho['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_patho['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_patho) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                
                print('PATHO',filter_patho)
                #filter_env['Reaction'] = list(filter_env['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_env['Reaction']))
                filter_env['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_env['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_env) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))

                print('NONPATHO',filter_env)

                #filter_slow['Subsystem'] = marinum.reactions.get_by_id(filter_slow['Reaction']).subsystem
                #filter_slow['Name'] = marinum.reactions.get_by_id(filter_slow['Reaction']).name
                filter_slow['Marinum']= [(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Name']= [(marinum.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Subsystem']= [(marinum.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                filter_slow['M. Gene']= [(marinum.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(marinum.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Marinum'] ==1 ]['Reaction']])))
                
                filter_slow['MTB']= [(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Name']= [(Rv.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Subsystem']= [(Rv.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))
                filter_slow['TB. Gene']= [(Rv.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']] + ([np.nan]* (len(filter_slow) - len([(Rv.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['MTB'] ==1 ]['Reaction']])))

                print('SLOW',filter_slow)
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + [0] * (len(Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis']==1]) - len(filter_fast['Reaction']))
                #filter_fast['Reaction'] = list(filter_fast['Reaction']) + (([np.nan] * (len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])- len(filter_fast))))
                filter_fast['Smegmatis']= [(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']]  + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Name']= [(smegmatis.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Subsystem']= [(smegmatis.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                filter_fast['S. Gene']= [(smegmatis.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(smegmatis.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Smegmatis'] ==1 ]['Reaction']])))
                
                
                filter_fast['Abscessus']= [(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Name']= [(abscessus.reactions.get_by_id(str(x)).name) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Subsystem']= [(abscessus.reactions.get_by_id(str(x)).subsystem) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                filter_fast['A. Gene']= [(abscessus.reactions.get_by_id(str(x)).genes) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']] + ([np.nan] * (len(filter_fast) - len([(abscessus.reactions.get_by_id(str(x)).id) for x in Dif_rxns_location_filtered[Dif_rxns_location_filtered['Abscessus'] ==1 ]['Reaction']])))
                print('FAST',filter_fast )
            
        Dif_rxns_location_filtered.to_csv(f'dif_rxns_{model_type}_{arom}.csv')
        filter_patho.to_csv(f'common_patho_{model_type}_{arom}.csv')
        filter_env.to_csv(f'common_env_{model_type}_{arom}.csv')
        filter_slow.to_csv(f'common_slow_{model_type}_{arom}.csv')
        filter_fast.to_csv(f'common_fast_{model_type}_{arom}.csv')
               
        detailed_sig_subsyst_sim = pd.DataFrame.from_dict(dict_list_sim)
        detailed_sig_subsyst_dif = pd.DataFrame.from_dict(dict_list_dif)
        #detailed_sig_subsyst_dif['Reaction'] = detailed_sig_subsyst_dif
        #df_combined = pd.merge(detailed_sig_subsyst_dif, Dif_rxns_location_filtered, on='Reaction', how='inner')
        
        
        detailed_sig_subsyst_dif.to_csv(f'Detail_dif_reactions_in_subsystem_fisher_{model_type}_{arom}.csv')
        detailed_sig_subsyst_sim.to_csv(f'Detail_sim_reactions_in_subsystem_fisher_{model_type}_{arom}.csv')
if heatmapDifMin:
    bigg_df = pd.read_csv('./bigg_reactions.txt', sep='\t', usecols=['bigg_id', 'database_links'])
    dif_min = pd.read_csv(f'./dif_rxns_{model_type}_{arom}.csv', usecols = [1,2,3,4,5,6])
    dif_min.loc[:,'ReactionName',] = [universe.reactions.get_by_id(x).name for x in dif_min['Reaction']]
    #dif_min_trim = dif_min[dif_min['ReactionName'] != '']

    # extract EC number from the bigg database file (database_links column)
    def extract_ec_number(links):
        if pd.notnull(links):  # Check if the link entry is not NaN
            # Search for the 'EC Number:' substring
            start = links.find('EC Number: http://identifiers.org/ec-code/')
            if start != -1:
                # Extract the EC number after the specified substring
                start += len('EC Number: http://identifiers.org/ec-code/')
                end = links.find(';', start)
                if end == -1:
                    end = len(links)  # if ';' is not found, take the rest of the string
                return links[start:end]
        return '-'  # Return '-' if no EC number is found

    # Apply the function to create a new 'EC_Number' column
    bigg_df['EC_Number'] = bigg_df['database_links'].apply(extract_ec_number)
    merged_df = dif_min.merge(bigg_df[['bigg_id', 'EC_Number']], how='left', left_on='Reaction', right_on='bigg_id')

    # Fill NaN in EC_Number with '-' to indicate missing EC numbers
    merged_df['EC_Number'] = merged_df['EC_Number'].fillna('-')
    merged_df = merged_df.drop(columns=['bigg_id'])

    merged_df.loc[merged_df['Reaction'] == 'DHBSZ3FEexs', 'ReactionName'] = '2-3-dihydroxybenzoylserine trimer Fe III sequestration'
    merged_df.loc[merged_df['Reaction'] == '4CMLCL_kt', 'ReactionName'] = '4-carboxymuconolactone-decarboxylase'
    merged_df.loc[merged_df['Reaction'] == '13PPDH2', 'ReactionName'] = '1,3-propanediol dehydrogenase (3-hydroxypropanal)'
    merged_df.loc[merged_df['Reaction'] == '13PPDH2_1', 'ReactionName'] = '1,3-propanediol dehydrogenase (L-lactaldehyde)'
    merged_df.loc[merged_df['Reaction'] == 'MCTC_1', 'ReactionName'] = 'Methylcrotonoyl-coa carboxylase'
    merged_df.loc[merged_df['Reaction'] == 'ACTD_1', 'ReactionName'] = 'Acetoin dehydrogenase (NADH)'
    merged_df.loc[merged_df['Reaction'] == 'ODH2E', 'ReactionName'] = '2-oxoadipate dehydrogenase'
    merged_df.loc[merged_df['Reaction'] == 'ATAH_1', 'ReactionName'] = 'Allantoate amidinohydrolase'
    merged_df.loc[merged_df['Reaction'] == 'ACTDa', 'ReactionName'] = 'Acetoin dehydrogenase (NADP)'
    merged_df.loc[merged_df['Reaction'] == 'ACTD', 'ReactionName'] = 'Acetoin dehydrogenase (NAD)'
    merged_df.loc[merged_df['Reaction'] == 'GLUDy', 'ReactionName'] = 'Glutamate dehydrogenase (NADP)'
    merged_df.drop(merged_df[merged_df['Reaction'] == "GUAt2pp"].index, inplace = True)
    merged_df.drop(merged_df[merged_df['Reaction'] == "GLCRAL_1"].index, inplace = True)
    merged_df.drop(merged_df[merged_df['Reaction'] == "MCTC_1"].index, inplace = True)
    merged_df.drop(merged_df[merged_df['Reaction'] == "GLCRD_1"].index, inplace = True)
    merged_df.drop(merged_df[merged_df['Reaction'] == "4CMLCL_kt"].index, inplace = True)



    

    merged_df.to_csv(f'bigg_EC_{model_type}_{arom}.csv')
    merged_df_latex = merged_df.drop(columns = ['Marinum','Smegmatis','Aromaticivorans','Abscessus','MTB'])
    merged_df_latex.loc[:,'Formula'] = [universe.reactions.get_by_id(x).reaction for x in merged_df['Reaction']]
    merged_df_latex.rename(columns= {'Reaction':'BiGG ID','ReactionName': 'Reaction Name','EC_number':'EC number'}, inplace =True)

    if model_type == 'min_gf':
        merged_df_latex.to_csv(f'latex_csv_test.csv')
        merged_df_latex.to_latex(f'table_{model_type}_{arom}_{date}_latex.txt',index= False)
    # open txt bigg_id ;'EC Number:'; 
    #dif_min.loc[:,'EC',] = [ for x in dif_min['Reaction']]
    # remove clustering band
    if arom == 'arom':
        clmap = sns.clustermap(merged_df.iloc[:, 1:6], cmap='Greens', annot=False, cbar=None, linewidths=0.5, yticklabels = merged_df['ReactionName']+' ('+ merged_df['EC_Number']+ ')', row_cluster = True, col_cluster = False, figsize=(10, 17))
        clmap.cax.set_visible(False)
        custom_xticks = ['$M. marinun$','$M. smegmatis$','$M. aromaticivorans$', '$M. abscessus$','$M. tuberculosis$']  
    else:
        clmap = sns.clustermap(merged_df.iloc[:, 1:5], cmap='Greens', annot=False, cbar=None, linewidths=0.5, yticklabels = merged_df['ReactionName']+' ('+ merged_df['EC_Number']+ ')', row_cluster = True, col_cluster = False, figsize=(10, 17))
        clmap.cax.set_visible(False)
        custom_xticks = ['$M. marinun$','$M. smegmatis$','$M. abscessus$','$M. tuberculosis$']  

    clmap.ax_heatmap.set_xticklabels(custom_xticks, rotation=45, ha='center')
    clmap.ax_row_dendrogram.set_visible(False)
    plt.savefig(f'./heatmap_{model_type}_reactions_{arom}.png', bbox_inches='tight', dpi=400)
    plt.show()

    # search similar reactions in other models 
    if model_type == 'min_gf':
        find_min_rxn_other = pd.DataFrame({})
        for rxn in merged_df['Reaction']:
            for model in model_list: 
                if rxn in [rxn.id for rxn in model.reactions]:
                    for model_others in model_list:
                        if model.reactions.get_by_id(rxn).reactants[0].id in model_others.metabolites:
                            find_min_rxn_other.loc[rxn,model_others] = str([item.id for item in iter(model_others.metabolites.get_by_id(model.reactions.get_by_id(rxn).reactants[0].id).reactions)])

        find_min_rxn_other.to_csv(f'similar_reactions_in_other_models_{model_type}_{date}.csv')
    
    if model_type == 'min_gf':
        find_min_rxn_other = pd.DataFrame({})
        for rxn in merged_df['Reaction']:
            for model in model_list: 
                if rxn in [rxn.id for rxn in model.reactions]:
                    for model_others in model_list:
                        if model.reactions.get_by_id(rxn).reactants[0].id in model_others.metabolites:
                            find_min_rxn_other.loc[rxn,model_others] = str([item.id for item in iter(model_others.metabolites.get_by_id(model.reactions.get_by_id(rxn).reactants[0].id).reactions)])

        find_min_rxn_other.to_csv(f'similar_reactions_in_other_models_{model_type}_{date}.csv')    
    orphan_met = pd.DataFrame({})
    for model in model_list:
        for met in model.metabolites:
            if met.id.endswith('_p'):
                if len([x for x in iter(model.metabolites.get_by_id(met.id).reactions)]) < 2:
                    orphan_met.loc[met.id,model] = str(met.id) + ', in reactions : '+str([x for x in iter(model.metabolites.get_by_id(met.id).reactions)])
            #if model.reactions.get_by_id(rxn).reactants[0].id in model_others.metabolites:
        #    find_min_rxn_other.loc[rxn,model_others] = str([item.id for item in iter(model_others.metabolites.get_by_id(model.reactions.get_by_id(rxn).reactants[0].id).reactions)])
    orphan_met.to_csv(f'orphan_met_periplasm_{model_type}_{date}.csv')
if plotSubsystemMaxgen:

    # Plot fisher test and multitest result for all species and both sides for max gene models
    significant_results_ab_df_less = pd.read_csv(f'Significant_results_abscessus_less_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_ab_df_great = pd.read_csv(f'Significant_results_abscessus_great_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_ma_df_less = pd.read_csv(f'Significant_results_marinum_less_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_ma_df_great = pd.read_csv(f'Significant_results_marinum_great_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_sm_df_less = pd.read_csv(f'Significant_results_smegmatis_less_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_sm_df_great = pd.read_csv(f'Significant_results_smegmatis_great_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_rv_df_less = pd.read_csv(f'Significant_results_mtb_less_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_rv_df_great = pd.read_csv(f'Significant_results_mtb_great_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_ar_df_less = pd.read_csv(f'Significant_results_aromaticivorans_less_{date}_max_genes_perc_gf_{arom}_high_low.csv')
    significant_results_ar_df_great = pd.read_csv(f'Significant_results_aromaticivorans_great_{date}_max_genes_perc_gf_{arom}_high_low.csv')



    print(significant_results_rv_df_less)
    dataframes_less = [significant_results_rv_df_less, significant_results_ab_df_less, significant_results_ma_df_less,significant_results_ar_df_less,significant_results_sm_df_less]
    dataframes_great = [significant_results_rv_df_great, significant_results_ab_df_great, significant_results_ma_df_great,significant_results_ar_df_great,significant_results_sm_df_great]

    colors =['blue','red','orange', 'violet','green']
    labels =  ['$M. tuberculosis$','$M. abscessus$','$M. marinum$','$M. aromaticivorans$', '$M. smegmatis$']
    pvalue_cutoff = 0.05

    # plot
    #fig, ax = plt.subplots()
    #fig.canvas.draw()

    plt.figure(figsize=(12, 8))

    #fig, axes = plt.subplots(1, 2, figsize=(14, 10), sharey=True)

    plt.title('Overrepresented in non-shared reactions')
    tick_array = ['Ascorbate and Aldarate Metabolism','Butanoate metabolism','Galactose metabolism','Levulinate metabolism','Phenylalanine metabolism','Tyrosine Metabolism','Other glycan degradation','Nitrotoluene degradation','Toluene degradation' ,'Lipid Metabolism','Fatty acid elongation','Fatty acid degradation','Glycerophospholipid metabolism','Environmental Information Processing','Membrane transport']
    #significant_results_sm_df_less['subsystem'] = pd.Categorical(significant_results_sm_df_less['subsystem'], categories=tick_array, ordered=True)
    #print(significant_results_sm_df_less['subsystem'])

    for df, color, label in zip(dataframes_less, colors, labels):
        #df = df[df['subsystem'].isin(tick_array)]
        #df['subsystem'] = df['subsystem'].fillna('Unknown').astype(str)
        #df.sort_values(['subsystem'], inplace=True)
        print(df['subsystem'])

        filtered_df = df[df['adjusted_p_value'] < pvalue_cutoff]
        filtered_df['subsystem'] =  pd.Categorical(filtered_df['subsystem'], categories=tick_array, ordered=True)
        #filtered_df['subsystem'].cat.set_categories(tick_array)

        #print(df['subsystem'].cat.categories)  # Should match 
        #print(filtered_df)
        if subsystem_plot_ratio:
            size_ratio = filtered_df['reactions_in_non_shared']*100/(filtered_df['reactions_in_non_shared'] +filtered_df['reactions_in_shared'])
        else: 
            size_ratio = filtered_df['reactions_in_non_shared']

        plt.scatter(filtered_df['adjusted_p_value'], filtered_df['subsystem'], color=color,s = size_ratio)
    
    
    if subsystem_plot_ratio:
        title_fig = 'Reactions in subsystem (%)'
        sizes_legend = [1,10,50,100]
        for size in sizes_legend:
            plt.scatter([], [], c='k', alpha=0.3, s=size, label=str(size))
    
    else: 
        title_fig = 'Number of Reactions'
        sizes_legend = [10,100,300]
        for size in sizes_legend:
            plt.scatter([], [], c='k', alpha=0.3, s=size,label=str(size))
    
    #ax.set_yticks(ticks = list(range(len(tick_array))))
    #ax.set_yticklabels(tick_array)

    plt.legend(scatterpoints=1, frameon=True, labelspacing=1, title= title_fig)
    plt.xlim(1e-14, 0.1)

    plt.xlabel('p-value')
    plt.ylabel('Subsystem')
    plt.xscale('log')

    #lgnd = plt.legend()
    #for handle in lgnd.legend_handles:
    #    handle.set_sizes([12.0])f
    
    plt.grid(True)
    plt.tight_layout()
    if subsystem_plot_ratio:
        plt.savefig(f'./subsystem_ratio_pvalue_{model_type}_ratio_non_shared_reactions_{date}_{arom}.png', bbox_inches='tight', dpi=400)
    else: 
        plt.savefig(f'./subsystem_absolute_pvalue_{model_type}_ratio_non_shared_reactions_{date}_{arom}.png', bbox_inches='tight', dpi=400)

    plt.show()

    plt.figure(figsize=(12, 8))

    plt.title('Overrepresented in shared reactions')
    for df, color, label in zip(dataframes_great, colors, labels):
        df.sort_values(['subsystem'], inplace=True)
        filtered_df = df[df['adjusted_p_value'] < pvalue_cutoff]
        if subsystem_plot_ratio:
            size_ratio = filtered_df['reactions_in_shared']*100/(filtered_df['reactions_in_non_shared'] +filtered_df['reactions_in_shared'])
        else: 
            size_ratio = filtered_df['reactions_in_shared']

        #size_ratio = filtered_df['reactions_in_shared']*100/(filtered_df['reactions_in_non_shared'] +filtered_df['reactions_in_shared'])
        plt.scatter(filtered_df['adjusted_p_value'], filtered_df['subsystem'], color=color,s = size_ratio)

    plt.xlabel('p-value')
    plt.ylabel('Subsystem')
    if subsystem_plot_ratio:
        title_fig = 'Reactions in subsystem (%)'
        sizes_legend = [1,10,50,100]
        for size in sizes_legend:
            plt.scatter([], [], c='k', alpha=0.3, s=size,label=str(size))
    
    else: 
        title_fig = 'Number of Reactions'
        sizes_legend = [10,100,300]
        for size in sizes_legend:
            plt.scatter([], [], c='k', alpha=0.3, s=size,label=str(size))
    

    plt.legend(scatterpoints=1, frameon=True, labelspacing=1, title=title_fig)
    plt.grid(True)
    plt.xlim(1e-14, 0.1)

    #plt.tick_params(axis='y', which='both', left=False, labelleft=False)
    plt.xscale('log')

    plt.tight_layout()
    if subsystem_plot_ratio:
        plt.savefig(f'./subsystem_ratio_pvalue_{model_type}_ratio_shared_reactions_{date}_{arom}.png', bbox_inches='tight', dpi=400)
    else: 
        plt.savefig(f'./subsystem_absolute_pvalue_{model_type}_ratio_shared_reactions_{date}_{arom}.png', bbox_inches='tight', dpi=400)

    plt.show()
if plotSubsystemMin:
    # Plot fisher test and multitest result for all species and both sides for min models

    significant_results_ab_df_less = pd.read_csv(f'Significant_results_abscessus_less_{date}_min_gf_{arom}_high_low.csv')
    significant_results_ab_df_great = pd.read_csv(f'Significant_results_abscessus_great_{date}_min_gf_{arom}_high_low.csv')
    significant_results_ma_df_less = pd.read_csv(f'Significant_results_marinum_less_{date}_min_gf_{arom}_high_low.csv')
    significant_results_ma_df_great = pd.read_csv(f'Significant_results_marinum_great_{date}_min_gf_{arom}_high_low.csv')
    significant_results_sm_df_less = pd.read_csv(f'Significant_results_smegmatis_less_{date}_min_gf_{arom}_high_low.csv')
    significant_results_sm_df_great = pd.read_csv(f'Significant_results_smegmatis_great_{date}_min_gf_{arom}_high_low.csv')
    significant_results_rv_df_less = pd.read_csv(f'Significant_results_mtb_less_{date}_min_gf_{arom}_high_low.csv')
    significant_results_rv_df_great = pd.read_csv(f'Significant_results_mtb_great_{date}_min_gf_{arom}_high_low.csv')
    significant_results_ar_df_less = pd.read_csv(f'Significant_results_aromaticivorans_less_{date}_min_gf_{arom}_high_low.csv')
    significant_results_ar_df_great = pd.read_csv(f'Significant_results_aromaticivorans_great_{date}_min_gf_{arom}_high_low.csv')

    dataframes_less = [significant_results_rv_df_less, significant_results_ab_df_less, significant_results_ma_df_less,significant_results_ar_df_less,significant_results_sm_df_less]
    dataframes_great = [significant_results_rv_df_great, significant_results_ab_df_great, significant_results_ma_df_great,significant_results_ar_df_great,significant_results_sm_df_great]

    colors =['blue','red','orange', 'violet','green']
    labels =  ['$M. tuberculosis$','$M. abscessus$','$M. marinum$','$M. aromaticivorans$', '$M. smegmatis$']
    pvalue_cutoff = 0.05
    

    # plot
    
    plt.figure(figsize=(7, 4))
    

    plt.title('Overrepresented in non-shared reactions')
    for df, color, label in zip(dataframes_less, colors, labels):
        filtered_df = df[df['adjusted_p_value'] < pvalue_cutoff]
        #flitered_df order 
        if subsystem_plot_ratio:
            size_ratio = filtered_df['reactions_in_non_shared']*100/(filtered_df['reactions_in_non_shared'] +filtered_df['reactions_in_shared'])
        else: 
            size_ratio = filtered_df['reactions_in_non_shared']

        plt.scatter(filtered_df['adjusted_p_value'], filtered_df['subsystem'], color=color, s = size_ratio)

    plt.xlabel('p-value')
    plt.ylabel('Subsystem')
    plt.xscale('log')

    if subsystem_plot_ratio:
        title_fig = 'Reactions in subsystem (%)'
        sizes_legend = [1,10,100]
        for size in sizes_legend:
            plt.scatter([], [], c='k', alpha=0.3, s=size, label=str(size))
    
    else: 
        title_fig = 'Number of Reactions'
        sizes_legend = [1,15,30]
        for size in sizes_legend:
            plt.scatter([], [], c='k', alpha=0.3, s=size, label=str(size))
    

    plt.legend(scatterpoints=1, frameon=True, labelspacing=1, title=title_fig)

#    lgnd = plt.legend(loc="best", scatterpoints=1, fontsize=10)
#    for handle in lgnd.legend_handles:
#        handle.set_sizes([12.0])
    #lgnd.legend_handles[1].set_sizes([6.0])
    plt.grid(True)

    #axes[1].set_title('Overrepresented in shared reactions')
    #for df, color, label in zip(dataframes_great, colors, labels):
    #    df.sort_values(['subsystem'], inplace=True)
    #    filtered_df = df[df['adjusted_p_value'] < pvalue_cutoff]
        
    #    axes[1].scatter(filtered_df['adjusted_p_value'], filtered_df['subsystem'], color=color)

    #axes[1].set_xlabel('p-value')
    #axes[1].legend()
    #axes[1].grid(True)
    #axes[1].tick_params(axis='y', which='both', left=False, labelleft=False)
    #axes[1].set_xscale('log')
    plt.tight_layout()
    if subsystem_plot_ratio:
        plt.savefig(f'./subsystem_ratio_adjusted_pvalue_{model_type}_reactions_{date}_{arom}.svg', bbox_inches='tight', dpi=400)
    else: 
        plt.savefig(f'./subsystem_absolute_adjusted_pvalue_{model_type}_reactions_{date}_{arom}.svg', bbox_inches='tight', dpi=400)

    plt.show()
if latexHierarchySubsystem:
     
    # export to latex table the hierarchy 
    rows = []
    for key, values in subsys_group.items():
        for value in values:
            rows.append({'Parent group': key, 'Subsystem': value})

    subsys_group_edit = pd.DataFrame(rows)

    
    subsys_group_edit.to_latex('subsys_hierarchy.txt')
    
    