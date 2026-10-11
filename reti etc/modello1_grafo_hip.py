# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 14:46:48 2025

@author: hp
"""

import pandas as pd

file_name = "C:/Users/hp/CampusOneDrive/OneDrive - alcampus.it/UNI/tirocinio/HIPPIE-subnetwork.tsv"

prova = pd.read_csv(file_name,sep="\t") #aprire il file in pandas colonnato (sto mettendo il file come input)
print(prova.columns)

prova = prova[["SDHB","SDHA"]] # sto estraendo le colonne che ci interessano
prova=prova.set_axis(["proteina1","proteina2"],axis=1) #con axis=1 vuol dire che sto rinominando le colonne

file_name = "C:/Users/hp/CampusOneDrive/OneDrive - alcampus.it/UNI/tirocinio/mitochondrial_ETC-lista-proteine.txt"
lista_pro= pd.read_csv(file_name, names=["proteine"])

#devo selezionare le righe della matrice che hanno entrambe le proteine che fanno parte della lista proteine
proteine = set(lista_pro["proteine"])#ho trasformato il dataframe in una lista 
mat_def = prova[prova["proteina1"].isin(proteine) & prova["proteina2"].isin(proteine)].reset_index(drop=True) #ho selezionato le protein dalla matrice prova degli indici che ho trovato

#dobbiamo fare la lista definitiva con le interazioni bidirezionali, ad esempio se c'è un'interazione proteina x e proteina y, ci deve essere anche un'interazione proteina y proteina x, quindi dovremo avere due righe con le stesse proteine, ma in ordine invertito
#dobbiamo quindi creare un'altra mat def con gli stessi elementi della mat def originale ma con le colonne invertite, per poi fare un confronto tra la nuova lista def e quella originale e elimare le righe uguali (perché nella lista def originale possono comunque esserci delle proteine con un'interazione bidirezionale)
mat_def2=mat_def[mat_def.columns[::-1]]#sto scambiando le colonne della matrice mat_def
mat_def2=mat_def2.set_axis(["proteina1","proteina2"],axis=1)

mat_conc=pd.concat([mat_def,mat_def2], ignore_index=True)#ho ottenuto una matrice con il doppio delle righe di mat_def

modello_hip=mat_conc.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_hip=modello_hip[modello_hip["proteina1"] != modello_hip["proteina2"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)


import networkx as nx
import matplotlib.pyplot as plt

G1_hip=nx.from_pandas_edgelist(modello_hip,source="proteina1", target="proteina2") #ho creato il grafo dalla matrice modello 

nodes_hip = list(G1_hip.nodes);#lista delle proteine che fanno parte del grafo
deg_hip= [list(G1_hip.degree())[i][1] for i in range(0, len(nodes_hip))]; #grado di ogni nodo
comm_net_hip= nx.community.louvain_communities(G1_hip, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_hip= {n: i for i, com in enumerate(comm_net_hip) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_hip);
node_color = [community_index_hip[n] for n in G1_hip];
node_size = [v * 500 for v in deg_hip];
nx.draw_networkx(
     G1_hip,
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

num_archi_hip= G1_hip.number_of_edges()
diam_hip= nx.diameter(G1_hip)
mod_hip= nx.community.modularity(G1_hip, comm_net_hip)
dens_hip= nx.density(G1_hip)
grado_medio_hip=sum(deg_hip)/len(deg_hip)
lista_clustering_hip=[nx.clustering(G1_hip, node) for node in G1_hip.nodes()]
clustering_medio_hip=sum(lista_clustering_hip)/len(lista_clustering_hip)
clustering_medio_hip=sum(lista_clustering_hip)/len(lista_clustering_hip)
betweenness_hip = nx.betweenness_centrality(G1_hip, normalized=True) 
valori_bet_hip = list(betweenness_hip.values()) 
betweenness_medio_hip=sum(valori_bet_hip)/len(valori_bet_hip)
ecc_hip= nx.eccentricity(G1_hip)         # lista delle eccentricità dei nodi
r_hip= nx.radius(G1_hip)                 # raggio (min eccentricity)
center_nodes_hip= nx.center(G1_hip)      # lista dei nodi baricentro
closeness_hip= nx.closeness_centrality(G1_hip)
valori_clos_hip= list(closeness_hip.values())
closeness_media_hip= sum(valori_clos_hip)/len(valori_clos_hip)
