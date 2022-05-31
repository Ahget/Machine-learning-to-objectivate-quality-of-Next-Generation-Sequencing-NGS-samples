import pandas as pd
from joblib import load
from var_glob import create_var_glob_tot

DT = load("DT_final_erasme.joblib")
RF = load("RF_final_erasme.joblib")
scaler = load("Scaler_final_erasme.joblib")


features = ['gc_percent_mean', 'gc_percent_median',
       'gc_percent_0.3', 'gc_percent_0.5', 'gc_percent_0.7', 'gc_percent_0.9',
       'cov20_percent_mean', 'cov20_percent_median', 'cov20_percent_0.3',
       'cov20_percent_0.5', 'cov20_percent_0.7', 'cov20_percent_0.9',
       'cov100_percent_mean', 'cov100_percent_median', 'cov100_percent_0.3',
       'cov100_percent_0.5', 'cov100_percent_0.7', 'cov100_percent_0.9',
       'cov500_percent_mean', 'cov500_percent_median', 'cov500_percent_0.3',
       'cov500_percent_0.5', 'cov500_percent_0.7', 'cov500_percent_0.9',
       'fwd_e2e_percent_mean', 'fwd_e2e_percent_median', 'fwd_e2e_percent_0.3',
       'fwd_e2e_percent_0.5', 'fwd_e2e_percent_0.7', 'fwd_e2e_percent_0.9',
       'rev_e2e_percent_mean', 'rev_e2e_percent_median', 'rev_e2e_percent_0.3',
       'rev_e2e_percent_0.5', 'rev_e2e_percent_0.7', 'rev_e2e_percent_0.9', 'longueur', 'OPT_percent',
       'Frequency_norm_mean', 'Frequency_norm_median', 'Frequency_norm_0.3',
       'Frequency_norm_0.5', 'Frequency_norm_0.7', 'Frequency_norm_0.9',
       'Quality_norm_mean', 'Quality_norm_median', 'Quality_norm_0.3',
       'Quality_norm_0.5', 'Quality_norm_0.7', 'Quality_norm_0.9',
       'Coverage_norm_mean', 'Coverage_norm_median', 'Coverage_norm_0.3',
       'Coverage_norm_0.5', 'Coverage_norm_0.7', 'Coverage_norm_0.9',
       'Allele Cov_norm_mean', 'Allele Cov_norm_median', 'Allele Cov_norm_0.3',
       'Allele Cov_norm_0.5', 'Allele Cov_norm_0.7', 'Allele Cov_norm_0.9',
       'longueur_norm_mean', 'longueur_norm_median', 'longueur_norm_0.3',
       'longueur_norm_0.5', 'longueur_norm_0.7', 'longueur_norm_0.9']


# préprocessing ######################################################

# traitement tsv

def prep_tsv(df_tsv):
    DF_tsv_norm = []
    if str(df_tsv["Barcode"][0])!="nan":
        df_tsv.fillna(0)
        df1 = df_tsv[df_tsv["Barcode"]==df_tsv["Barcode"][0]]
        df1.index = range(len(df1))
        df1.loc[:,"longueur"] = len(df1)
        
        DF_tsv_norm = scaler.transform(df1[["Frequency", "Quality", "Coverage", "Allele Cov", "longueur"]])
        DF_tsv_norm = pd.concat([df1, pd.DataFrame(DF_tsv_norm, columns = ["Frequency_norm", "Quality_norm", "Coverage_norm", "Allele Cov_norm", "longueur_norm"])], axis = 1)
        is_NC = False
    else:
        is_NC = True
        
    return DF_tsv_norm, is_NC




# traitement amplicon

def prep_amp(df_amplicon):
    df_amplicon.loc[:,"longueur"] = len(df_amplicon)
    df_amplicon.loc[:, "contig_length"] = df_amplicon["contig_end"] - df_amplicon["contig_srt"]+1
    df_amplicon.loc[:, "gc_percent"] = df_amplicon["gc_count"]/df_amplicon["contig_length"]
    df_amplicon.loc[:, "cov20_percent"] = df_amplicon["cov20x"]/df_amplicon["contig_length"]
    df_amplicon.loc[:, "cov100_percent"] = df_amplicon["cov100x"]/df_amplicon["contig_length"]
    df_amplicon.loc[:, "cov500_percent"] = df_amplicon["cov500x"]/df_amplicon["contig_length"]
    df_amplicon.loc[:, "fwd_e2e_percent"] = df_amplicon["fwd_e2e"]/df_amplicon["fwd_reads"]
    df_amplicon.loc[:, "rev_e2e_percent"] = df_amplicon["rev_e2e"]/df_amplicon["rev_reads"]
    return df_amplicon

def pred2class(DF_amp):
    New_var1 = ['gc_percent', 'cov500_percent', 'fwd_e2e_percent', 'rev_e2e_percent', 'cov20_percent', 'cov100_percent', 'contig_length']
    pred = DT.predict(DF_amp[New_var1])
    DF_amp.loc[:,"Pred2class"]=pred
    return DF_amp


# classification ######################################################

def classif(df):
    pred = RF.predict(df[features])
    pourcentages = pd.DataFrame(RF.predict_proba(df[features]), columns= RF.classes_)
    OPT_percent = df["OPT_percent"]

    return [pred, pourcentages, OPT_percent]

# fonction principale ######################################################

def function(liste):
    results = []
    df_amp = pd.read_csv([i for i in liste if i.find("amplicon")>-1][0], sep='\t').fillna(0)
    df_tsv = pd.read_csv([i for i in liste if i.find("tsv")>-1][0], sep='\t')
    
    df_tsv_norm, is_NC = prep_tsv(df_tsv)
    if is_NC == False:
        df_ampli = prep_amp(df_amp)
        df_ampli = pred2class(df_ampli)
        df = create_var_glob_tot(df_ampli, df_tsv_norm)
        results = classif(df)
        results[0] = results[0][0]
        results[1] = {"NC":results[1]["NC"].values[0], "NVA":results[1]["NVA"].values[0], "OPT":results[1]["OPT"].values[0], "SOPT":results[1]["SOPT"].values[0]}
        results[2] = str(results[2].values[0])[:4]
        
    return results
    