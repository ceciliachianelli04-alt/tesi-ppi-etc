# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 15:18:10 2025

@author: hp
"""

import pandas as pd
import re

file_name = "C:/Users/cecil/OneDrive - alcampus.it/UNI/tirocinio/MINT-species_human.txt"

prova_min= pd.read_csv(file_name,sep="\t", header=None, dtype=str) #aprire il file in pandas colonnato (sto mettendo il file come input)
# tengo solo le colonne utili
prova_min= prova_min[[0,1,4,5]].copy()
prova_min.columns = ["idA", "idB", "aliasA", "aliasB"]

#Funzione per estrarre gene name
def extract_gene_name(alias_field):
    if pd.isna(alias_field):
        return None
    parts = alias_field.split("|")
    for p in parts:
        m = re.search(r"uniprotkb:([A-Za-z0-9\-]+)\(gene name\)", p)
        if m:
            return m.group(1)
        m2 = re.search(r"uniprotkb:([A-Za-z0-9\-]+)\(gene name synonym\)", p)
        if m2:
            return m2.group(1)
    return None

# Estrazione gene name
prova_min["geneA"] = prova_min["aliasA"].apply(extract_gene_name)
prova_min["geneB"] = prova_min["aliasB"].apply(extract_gene_name)
prova_min= prova_min[["geneA", "geneB"]].dropna()
prova_min=prova_min.set_axis(["proteina1","proteina2"],axis=1) #con axis=1 vuol dire che sto rinominando le colonne

file_name = "C:/Users/cecil/OneDrive - alcampus.it/UNI/tirocinio/mitochondrial_ETC-lista-proteine.txt"
lista_pro= pd.read_csv(file_name, names=["proteine"])

#devo selezionare le righe della matrice che hanno entrambe le proteine che fanno parte della lista proteine
proteine = set(lista_pro["proteine"])#ho trasformato il dataframe in una lista 
mat_def = prova_min[prova_min["proteina1"].isin(proteine) & prova_min["proteina2"].isin(proteine)].reset_index(drop=True) #ho selezionato le protein dalla matrice prova degli indici che ho trovato

#dobbiamo fare la lista definitiva con le interazioni bidirezionali, ad esempio se c'è un'interazione proteina x e proteina y, ci deve essere anche un'interazione proteina y proteina x, quindi dovremo avere due righe con le stesse proteine, ma in ordine invertito
#dobbiamo quindi creare un'altra mat def con gli stessi elementi della mat def originale ma con le colonne invertite, per poi fare un confronto tra la nuova lista def e quella originale e elimare le righe uguali (perché nella lista def originale possono comunque esserci delle proteine con un'interazione bidirezionale)
mat_def2=mat_def[mat_def.columns[::-1]]#sto scambiando le colonne della matrice mat_def
mat_def2=mat_def2.set_axis(["proteina1","proteina2"],axis=1)

mat_conc=pd.concat([mat_def,mat_def2], ignore_index=True)#ho ottenuto una matrice con il doppio delle righe di mat_def

modello_min=mat_conc.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_min=modello_min[modello_min["proteina1"] != modello_min["proteina2"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)
modello_min.to_csv('MINT_ETC.csv', sep=',', index=False)


import networkx as nx
import matplotlib.pyplot as plt

G1_min=nx.from_pandas_edgelist(modello_min,source="proteina1", target="proteina2") #ho creato il grafo dalla matrice modello 
G1_min.remove_nodes_from(["CYCS", "MT-CO1"])

nodes_min = list(G1_min.nodes);#lista delle proteine che fanno parte del grafo
deg_min= [list(G1_min.degree())[i][1] for i in range(0, len(nodes_min))]; #grado di ogni nodo
comm_net_min= nx.community.louvain_communities(G1_min, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_min= {n: i for i, com in enumerate(comm_net_min) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_min);
node_color = [community_index_min[n] for n in G1_min];
node_size = [v * 500 for v in deg_min];
nx.draw_networkx(
     G1_min,
     pos=pos,
     with_labels=False,
     node_color=node_color,
     node_size=node_size,
     edge_color="grey",
     alpha=0.5,
     width=5,
     cmap=plt.cm.RdYlBu,
)
plt.show()

num_archi_min= G1_min.number_of_edges()
diam_min= nx.diameter(G1_min)
mod_min= nx.community.modularity(G1_min, comm_net_min)
dens_min= nx.density(G1_min)
grado_medio_min=sum(deg_min)/len(deg_min)
lista_clustering_min=[nx.clustering(G1_min, node) for node in G1_min.nodes()]
clustering_medio_min=sum(lista_clustering_min)/len(lista_clustering_min)
clustering_medio_min=sum(lista_clustering_min)/len(lista_clustering_min)
betweenness_min = nx.betweenness_centrality(G1_min, normalized=True) 
valori_bet_min = list(betweenness_min.values()) 
betweenness_medio_min=sum(valori_bet_min)/len(valori_bet_min)
ecc_min= nx.eccentricity(G1_min)         # lista delle eccentricità dei nodi
r_min= nx.radius(G1_min)                 # raggio (min eccentricity)
center_nodes_min= nx.center(G1_min)      # lista dei nodi baricentro
closeness_min= nx.closeness_centrality(G1_min)
valori_clos_min= list(closeness_min.values())
closeness_media_min= sum(valori_clos_min)/len(valori_clos_min)
