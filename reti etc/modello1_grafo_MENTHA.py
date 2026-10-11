# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:45:10 2025

@author: hp
"""

import pandas as pd


file_name = "C:/Users/hp/CampusOneDrive/OneDrive - alcampus.it/UNI/tirocinio/MENTHA-relations_result_22_10_25.tsv"

with open(file_name, "r", encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

righe_pulite = []

for line in lines:
    parts = line.strip().split(";")
    # prendere solo le righe con esattamente 6 colonne
    if len(parts) == 6:
        righe_pulite.append(parts)

tabella_men = pd.DataFrame(righe_pulite, columns=["Protein A", "Gene A", "Protein B", "Gene B", "Score", "PMID"])


# Estrarre le due colonne
tabella_men = tabella_men[["Gene A", "Gene B"]]

# Rinominare le colonne
tabella_men.columns = ["proteina1_men", "proteina2_men"]

file_name = "C:/Users/hp/CampusOneDrive/OneDrive - alcampus.it/UNI/tirocinio/mitochondrial_ETC-lista-proteine.txt"
lista_pro= pd.read_csv(file_name, names=["proteine"])

# Convertire tutto in stringhe e rimuovere spazi
lista_pro["proteine"] = lista_pro["proteine"].astype(str).str.strip().str.upper()

tabella_men["proteina1_men"] = tabella_men["proteina1_men"].astype(str).str.strip().str.upper()
tabella_men["proteina2_men"] = tabella_men["proteina2_men"].astype(str).str.strip().str.upper()


proteine = set(lista_pro["proteine"])

mask1 = tabella_men["proteina1_men"].isin(proteine)
mask2 = tabella_men["proteina2_men"].isin(proteine)

mat_def_men = tabella_men[mask1 & mask2]

#dobbiamo fare la lista definitiva con le interazioni bidirezionali, ad esempio se c'è un'interazione proteina x e proteina y, ci deve essere anche un'interazione proteina y proteina x, quindi dovremo avere due righe con le stesse proteine, ma in ordine invertito
#dobbiamo quindi creare un'altra mat def con gli stessi elementi della mat def originale ma con le colonne invertite, per poi fare un confronto tra la nuova lista def e quella originale e elimare le righe uguali (perché nella lista def originale possono comunque esserci delle proteine con un'interazione bidirezionale)
mat_def2_men=mat_def_men[mat_def_men.columns[::-1]]#sto scambiando le colonne della matrice mat_def
mat_def2_men=mat_def2_men.set_axis(["proteina1_men","proteina2_men"],axis=1)

mat_conc_men=pd.concat([mat_def_men,mat_def2_men], ignore_index=True)#ho ottenuto una matrice con il doppio delle righe di mat_def

modello_men=mat_conc_men.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_men=modello_men[modello_men["proteina1_men"] != modello_men["proteina2_men"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)
modello_men.to_csv('MENTHA_ETC.csv', sep=',', index=False)

import networkx as nx
import matplotlib.pyplot as plt

G1_men=nx.from_pandas_edgelist(modello_men,source="proteina1_men", target="proteina2_men") #ho creato il grafo dalla matrice modello 

nodes_men = list(G1_men.nodes);#lista delle proteine che fanno parte del grafo
deg_men = [list(G1_men.degree())[i][1] for i in range(0, len(nodes_men))]; #grado di ogni nodo
comm_net_men = nx.community.louvain_communities(G1_men, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_men = {n: i for i, com in enumerate(comm_net_men) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_men);
node_color = [community_index_men[n] for n in G1_men];
node_size = [v * 500 for v in deg_men];
nx.draw_networkx(
     G1_men,
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

num_archi_men = G1_men.number_of_edges()
diam_men= nx.diameter(G1_men)
mod_men= nx.community.modularity(G1_men, comm_net_men)
dens_men= nx.density(G1_men)
grado_medio_men=sum(deg_men)/len(deg_men)
lista_clustering_men=[nx.clustering(G1_men, node) for node in G1_men.nodes()]
clustering_medio_men=sum(lista_clustering_men)/len(lista_clustering_men)
betweenness_men = nx.betweenness_centrality(G1_men, normalized=True)
valori_bet_men= list(betweenness_men.values())
betweenness_medio_men=sum(valori_bet_men)/len(valori_bet_men)
ecc_men= nx.eccentricity(G1_men)         # lista delle eccentricità dei nodi
r_men= nx.radius(G1_men)                 # raggio (min eccentricity)
center_nodes_men= nx.center(G1_men)      # lista dei nodi baricentro
closeness_men= nx.closeness_centrality(G1_men)
valori_clos_men= list(closeness_men.values())
closeness_media_men = sum(valori_clos_men)/len(valori_clos_men)



