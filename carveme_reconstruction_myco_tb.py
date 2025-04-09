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

## Specify working directory
dir_cmd = "/home/icancino/Documents/Models"
new_dir_cmd = "."
## Reference model for carveme in xml format
ref = '/home/icancino/Documents/Models/Ref/ref_iEK1008_30723.xml'

## Universe database selection
config = configparser.ConfigParser()
config.read('/home/icancino/anaconda3/envs/reconstruction/lib/python3.7/site-packages/carveme/config.cfg')
config['generated']['default_universe'] = 'data/generated/universe_bact_arch_3.xml.gz'
universe= 'bact_arch_3'
uni_file = '/home/icancino/Documents/Models/universe/universe_bact_arch_3.xml'
## Specifying genomes folder (must be fasta)
genome_folder ='/eggnog_fasta_5'


## Calling carveme
def call_carve(filename,out,diamond,piden,evalue,ref, universe_file):
        
        cmd = 'carve '+str(filename)+' --output '+str(out)+' --diamond-args='+str(diamond) + ' --piden ' + str(piden) + ' --evalue ' + str(evalue) +' --reference ' + str(ref) + ' --universe-file '+ universe_file + ' --reference-score 1' 
        print(cmd)
        process = subprocess.Popen(cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True
                    )
        stdout, stderr = process.communicate()
        return stdout, stderr         

## Calling MEMOTE
def call_memote(name,model):
    cmd  = 'memote report snapshot --exclusive test_consistency --filename '+ str(name)+'.html '+ str(model) 
    print(cmd)
    process = subprocess.Popen(cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True
                    ) 
    stdout, stderr = process.communicate()
    return stdout, stderr         



## Genomes to use for the models. Taking models from folder
genome_list = []
for f in os.listdir(os.fsencode(dir_cmd + genome_folder)):
    filename = os.fsdecode(f)
    if filename.__contains__('H37Rv'):
        genome_list.append(filename)

## Diamond alignment parameters
args = ["--fast","--sensitive","--more-sensitive"]
## 
evalues = [1,1E-10,1E-50]
piden = [30,40,45,50,55,60,65]

## Dictionary to save carveme parameters
param = {}

## models to use

bigg_mtb = cobra.io.read_sbml_model('/home/icancino/Documents/Models/GEMS/iEK1008(1).xml')
bigg_ecoli =cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iJO1366.json') 
bigg_ecoli_2 = cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iML1515.json')
bigg_ecoli_3 = cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iY75_1357.json')
bigg_meth = cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iAF692.json')
bigg_pseud = cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iJN1463.json')
bigg_pseud_2 = cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iJN746.json')
bigg_syn_2 = cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iJN678.json')
bigg_salmo =cobra.io.load_json_model('/home/icancino/Documents/Models/GEMS/iYS1720.json')


metabolites_to_change_dict = {'fdxo_2_2_c':'FeS',
                        'fmnRD_c':'C17H21N4O9P',
                        'mobd_p':'H2MoO4',
                        'mobd_e':'H2MoO4', #
                         'mobd_c':'H2MoO4', #
                        'dmso2_p':'C2H6O2S',
                        'dmso2_e':'C2H6O2S',
                        'm3hdp_c':'C5H11O8P2',
                        '23dhbzs_c':'C10H10NO6',
                        'cholc5coa_c':'C45H70N7O18P3S',
                        '14glucan_c':'C6H10O5',
                        'mphen_c':'C37N2O1H50',
                        'mphenh2_c':'C37N2O1H52',
                        '14glucan_p':'C6H10O5',
                        '14glucan_e':'C6H10O5',
                        'pdima_e' :'C96H190O5',
                        'pdima_c' :'C96H190O5',
                        'hexccoa_c' :'C47H83N7O17P3S', #as iEK1011_2 model
                        'fdxrd_c':'FeS',
                        'fdxox_c':'FeS',
                        '23dhbzs3_c':'C30H29N3O16',
                        'fe3dhbzs3_p':'C30FeH29N3O16',
                        'fe3dhbzs3_e':'C30FeH29N3O16',
                        'fe3dhbzs3_c':'C30FeH29N3O16',
                        'pdima_p':'C96H190O5',
                        'hexccoa_c':'C47H82N7O17P3S',
                        'fe3dhbzs_p':'C10H10NO6Fe',
                        'fe3dhbzs_e':'C10H10NO6Fe',
                        'fe3dhbzs_c':'C10H10NO6Fe',
                        'ocholc5coa_c':'C45H68N7O19P3S',
                        'hcholc5coa_c':'C45H70N7O19P3S',
                        'cholenec5coa_c':'C45H68N7O18P3S',
                        'cholc5coa_c':'C45H70N7O18P3S',
                        'ocholc8coa_c':'C48H74N7O19P3S',
                        'hcholc8coa_c':'C48H76N7O19P3S',
                        'ochol_c':'C27H42O2',
                        'chol4en3one_c':'C27H44O',
                        'cholenec8coa_c':'C48H74N7O18P3S',
                        'hchol_c':'C27H44O2'} 

## Mass balance
def is_mass_balanced(reaction):
    balance = defaultdict(int)
    for metabolite, coefficient in iteritems(reaction.metabolites):
        if metabolite.elements is None or len(metabolite.elements) == 0:
            return False
        for element, amount in iteritems(metabolite.elements):
            balance[element] += coefficient * amount
    return all(amount == 0 for amount in itervalues(balance))

## Function to change met formula
def met_change_formula(model, metabolites_to_change):
    model_var = cobra.io.read_sbml_model(model) 
    met_keys =  model_var.metabolites.__dict__['_dict'].keys()

    for met_key in metabolites_to_change.keys():
        if met_key in met_keys:
            try:
                model_var.metabolites.get_by_id(met_key).formula = metabolites_to_change[met_key]    
            except:
                continue
    cobra.io.write_sbml_model(model_var,str(model))
    return 

## Function to change formula of metabolites from other models
def met_change_other_models(model, bigg_mtb, bigg_ecoli, bigg_meth, bigg_ecoli_2, bigg_pseud, bigg_salmo, bigg_ecoli_3,bigg_syn_2):
    model_var = cobra.io.read_sbml_model(model) 
    length =[]
    cond = True

    while cond: 
        rxn_mass_unbal =[]
        rxn_charge_unbal = []
        met_mass_change=[]
        for rxn in model_var.reactions:
            if rxn.check_mass_balance() and list(rxn.check_mass_balance().keys())[0] == 'charge' and len(list(rxn.check_mass_balance().keys())) == 1:
                if str(rxn.id)[:2] != 'EX' and str(rxn.id)[:4] != 'sink' and str(rxn.id)[:2] != 'DM':
                    rxn_charge_unbal.append(rxn)
            if is_mass_balanced(rxn) == False and str(rxn.id)[:2] != 'EX' and str(rxn.id)[:4] != 'sink' and str(rxn.id)[:2] != 'DM':
                rxn_mass_unbal.append(rxn)                        

        for x in rxn_mass_unbal:
            for met in model_var.reactions.get_by_id(x.id).metabolites:
                met_mass_change.append(met) 
        met_mass_change_nodup =list(set(met_mass_change))
        length.append(len(rxn_mass_unbal))
        remaining_mass_met =[]
        for met in met_mass_change_nodup:
            if met in bigg_mtb.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_mtb.metabolites.get_by_id(str(met.id)).formula            #print(met)
            elif met in bigg_meth.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_meth.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_pseud.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_pseud.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_salmo.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_salmo.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_ecoli.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_ecoli.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_ecoli_2.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_ecoli_2.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_ecoli_3.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_ecoli_3.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_pseud_2.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_pseud_2.metabolites.get_by_id(str(met.id)).formula
            elif met in bigg_syn_2.metabolites:
                model_var.metabolites.get_by_id(str(met.id)).formula = bigg_syn_2.metabolites.get_by_id(str(met.id)).formula
            else:
                remaining_mass_met.append(met)
        
        if len(length) >2 and length[-1] == length[-2]:
            cond = False
    cobra.io.write_sbml_model(model_var,str(model))
    return
##

#directory = os.fsencode(dir_cmd + '/') 
x=0 #Counter for failed models (carveme optimization problem fail)
m =0 #counter for model creation
for genome in genome_list:
    if Path(str(new_dir_cmd)+'/'+str(genome).split('.')[0]).exists():
        os.chdir(str(new_dir_cmd)+'/'+str(genome).split('.')[0])
        print('Repository exists')
    else:
        print('Creating rep. '+str(genome).split('.')[0])
        os.mkdir(str(new_dir_cmd)+'/'+str(genome).split('.')[0])
        os.chdir(str(new_dir_cmd)+'/'+str(genome).split('.')[0])
    with open("CARVEME_PARAMETERS.json", "w") as readme:
        for arg in args:
            for evalue in evalues:
    	        for pidentity in piden:
                    gen = dir_cmd + genome_folder + '/'+ genome
                    outname = genome.split('.')[0] +str(arg)[2:]+ str(pidentity)+str(evalue)+'.xml'
                    print(gen)

                    if Path('./'+outname).exists() == False:
                        print('Reconstructing '+outname)
                        try:
                            stdout, stderr = call_carve(gen,outname,'"'+arg+'"',pidentity, evalue,ref, uni_file)
                            param['Date'] = str(datetime.datetime.now())
                            param['Genome'] = genome
                            param['Model name'] = outname
                            param['Diamond alignment sensitivity'] = str(arg)
                            param['e-value of alignment cutoff'] = str(evalue)
                            param['Percentage of identity of alignment cutoff (%)'] = str(pidentity)
                            param['Universe database'] = str(config['generated']['default_universe'].split('/')[2][:-4])
                            param['Reference model'] = str(ref)
                            param['Genome folder'] = str(genome_folder)
                            param['Reference score'] = '1'

                            param['Error'] = 'No error'
                            json.dump(param, readme)

                        except:
                            param['Date'] = str(datetime.datetime.now())
                            param['Genome'] = genome
                            param['Model name'] = outname
                            param['Diamond alignment sensitivity'] = str(arg)
                            param['e-value of alignment cutoff'] = str(evalue)
                            param['Percentage of identity of alignment cutoff (%)'] = str(pidentity)
                            param['Universe database'] = str(config['generated']['default_universe'].split('/')[2][:-4])
                            param['Reference model'] = str(ref)
                            param['Genome folder'] = str(genome_folder)
                            param['Reference score'] = '1'

                            param['Error'] = 'Error in creation. No model'
                            json.dump(param,readme)
                            continue
                    
                        model_file = Path('./'+ outname)
                    
                        met_change_other_models(outname, bigg_mtb, bigg_ecoli, bigg_meth, bigg_ecoli_2, bigg_pseud, bigg_salmo, bigg_ecoli_3,bigg_syn_2)
                        met_change_formula(outname,metabolites_to_change_dict)
                        try:
                            stdout_m, stderr_m = call_memote(outname[:-4],outname)
                        except:
                            print('Memote did not run')
                        m += 1
                    

                        print('models ',m)
    os.chdir('../')
 
    
# directory_xml = dir_cmd + '/' + fold_mod
# directory_tsv = dir_cmd + '/' + fold_tsv
# directory_csv = dir_cmd + '/' + fold_csv
# for f in os.listdir(directory):
#         filename = os.fsdecode(f)
#         if filename.endswith('.xml'):
#                 model = dir_cmd + '/' + filename
#                 shutil.move(model,directory_xml + '/' + filename)
#         elif filename.endswith('.tsv'):
#                 align = dir_cmd + '/' + filename
#                 shutil.move(align,directory_tsv + '/' + filename)
#         elif filename.endswith('.csv'):
#                 csv = dir_cmd + '/' + filenameFwd:
#                 shutil.move(csv,directory_csv + '/' + filename)
#         else:
#                 continue
            
