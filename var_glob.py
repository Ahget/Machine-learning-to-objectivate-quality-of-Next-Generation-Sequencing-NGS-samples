import pandas as pd
import numpy as np
# des fonctions pour extraire les variables globales

def Extract_glob_feat(df, feat, longueur):
    nume = []
    nume.append(np.array([1 for i in df[feat].values if i>0.3]).sum()/longueur)
    nume.append(np.array([1 for i in df[feat].values if i>0.5]).sum()/longueur)
    nume.append(np.array([1 for i in df[feat].values if i>0.7]).sum()/longueur)
    nume.append(np.array([1 for i in df[feat].values if i>0.9]).sum()/longueur)
    df_res = pd.DataFrame(data = [[df[feat].mean(), df[feat].median(), nume[0], nume[1], nume[2], nume[3]]] , columns = [feat+"_mean", feat+"_median", feat+"_"+str(0.3), feat+"_"+str(0.5), feat+"_"+str(0.7), feat+"_"+str(0.9)])
    return df_res
    
def extract_all_var(df, variables, longueur):
    list_df = []
    for feat in variables:
        df_t = Extract_glob_feat(df, feat, longueur)
        list_df.append(df_t)
    df_tot = pd.concat(list_df, axis = 1)
    return df_tot

def create_var_glob_amplicon(df_amp):
    var_percent = ['gc_percent', 'cov20_percent', 'cov100_percent', 'cov500_percent', 'fwd_e2e_percent', 'rev_e2e_percent']

    longueur = df_amp["longueur"].values[0]
    # information des variables en pourcent
    df_var = extract_all_var(df_amp, var_percent, longueur)
    df_var.loc[:, "longueur"] = longueur

    # varible globale % de OPT
    OPT = np.array([1 for i in range(len(df_amp)) if df_amp["Pred2class"].values[i]=="OPT"]).sum()
    OPT_percent = OPT/longueur

    df_var.loc[:, "OPT_percent"] = OPT_percent

    return df_var

def create_var_glob_tsv(df_tsv):
    var_percent = ['Frequency_norm', 'Quality_norm', 'Coverage_norm', 'Allele Cov_norm', 'longueur_norm']

    longueur = df_tsv["longueur"].values[0]
# information des variables en pourcent
    df_var = extract_all_var(df_tsv, var_percent, longueur)

    return df_var

def create_var_glob_tot(df_amp, df_tsv):
    df_amp_glob = create_var_glob_amplicon(df_amp)
    df_tsv_glob = create_var_glob_tsv(df_tsv)
    df_global = pd.concat([df_amp_glob, df_tsv_glob], axis = 1)
    return df_global