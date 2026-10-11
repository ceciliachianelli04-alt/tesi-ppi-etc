# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 11:57:21 2026

@author: cecil
"""

import pandas as pd

file_name = "C:/Users/cecil/OneDrive - alcampus.it/UNI/tirocinio/BIOGRID-ALL-5.0.250.tab.txt"

prova = pd.read_csv(file_name,sep="\t") #aprire il file in pandas colonnato (sto mettendo il file come input)

prova = prova[["OFFICIAL_SYMBOL_A","OFFICIAL_SYMBOL_B"]] # sto estraendo le colonne che ci interessano
prova=prova.set_axis(["proteina1","proteina2"],axis=1) #con axis=1 vuol dire che sto rinominando le colonne


file_name = "C:/Users/cecil/OneDrive - alcampus.it/UNI/tirocinio/listapro_neurodegeneration.txt"
lista_pro= pd.read_csv(file_name, names=["proteine"])

#devo selezionare le righe dellla matrice [28190102x2] che hanno entrambe le proteine che fanno parte della lista proteine
proteine = set(lista_pro["proteine"])#ho trasformato il dataframe in una lista 
mat_def = prova[prova["proteina1"].isin(proteine) & prova["proteina2"].isin(proteine)].reset_index(drop=True) #ho selezionato le protein dalla matrice prova degli indici che ho trovato

#dobbiamo fare la lista definitiva con le interazioni bidirezionali, ad esempio se c'è un'interazione proteina x e proteina y, ci deve essere anche un'interazione proteina y proteina x, quindi dovremo avere due righe con le stesse proteine, ma in ordine invertito
#dobbiamo quindi creare un'altra mat def con gli stessi elementi della mat def originale ma con le colonne invertite, per poi fare un confronto tra la nuova lista def e quella originale e elimare le righe uguali (perché nella lista def originale possono comunque esserci delle proteine con un'interazione bidirezionale)
mat_def2=mat_def[mat_def.columns[::-1]]#sto scambiando le colonne della matrice mat_def
mat_def2=mat_def2.set_axis(["proteina1","proteina2"],axis=1)

mat_conc=pd.concat([mat_def,mat_def2], ignore_index=True)#ho ottenuto una matrice con il doppio delle righe di mat_def

modello_biogrid_malato=mat_conc.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_biogrid_malato=modello_biogrid_malato[modello_biogrid_malato["proteina1"] != modello_biogrid_malato["proteina2"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)
# Rimuovi le interazioni duplicate mantenendo solo una direzione
modello_biogrid_malato_def = modello_biogrid_malato[modello_biogrid_malato.apply(lambda row: row['proteina1'] < row['proteina2'], axis=1)].reset_index(drop=True)

modello_biogrid_malato_def.to_csv('BIOGRID_NEURODEGENERATION_DEDUP.csv', sep=',', index=False)