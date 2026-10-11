# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 15:01:22 2025

@author: hp
"""

# voglio costruire sempre il grafo modello 1 ma questa volta basandomi sul database IID, quello fatto precedentemente era rispetto al database BIOGRID

import pandas as pd

file_name = r"C:\Users\hp\CampusOneDrive\OneDrive - alcampus.it\UNI\tirocinio\IID-PPIs.txt"
tabella_iid= pd.read_csv(file_name,sep="\t") #aprire il file in pandas colonnato (sto mettendo il file come input)

tabella_iid = tabella_iid[["symbol1","symbol2"]] # sto estraendo le colonne che ci interessano
tabella_iid=tabella_iid.set_axis(["proteina1_iid","proteina2_iid"],axis=1) #con axis=1 vuol dire che sto rinominando le colonne

file_name = "C:/Users/hp/CampusOneDrive/OneDrive - alcampus.it/UNI/tirocinio/mitochondrial_ETC-lista-proteine.txt"
lista_pro= pd.read_csv(file_name, names=["proteine"])

#devo selezionare le righe della matrice che hanno entrambe le proteine che fanno parte della lista proteine
proteine = set(lista_pro["proteine"])#ho trasformato il dataframe in una lista 
mat_def_iid = tabella_iid[tabella_iid["proteina1_iid"].isin(proteine) & tabella_iid["proteina2_iid"].isin(proteine)].reset_index(drop=True) #ho selezionato le protein dalla matrice prova degli indici che ho trovato

#dobbiamo fare la lista definitiva con le interazioni bidirezionali, ad esempio se c'è un'interazione proteina x e proteina y, ci deve essere anche un'interazione proteina y proteina x, quindi dovremo avere due righe con le stesse proteine, ma in ordine invertito
#dobbiamo quindi creare un'altra mat def con gli stessi elementi della mat def originale ma con le colonne invertite, per poi fare un confronto tra la nuova lista def e quella originale e elimare le righe uguali (perché nella lista def originale possono comunque esserci delle proteine con un'interazione bidirezionale)
mat_def2_iid=mat_def_iid[mat_def_iid.columns[::-1]]#sto scambiando le colonne della matrice mat_def
mat_def2_iid=mat_def2_iid.set_axis(["proteina1_iid","proteina2_iid"],axis=1)

mat_conc_iid=pd.concat([mat_def_iid,mat_def2_iid], ignore_index=True)#ho ottenuto una matrice con il doppio delle righe di mat_def

modello_iid=mat_conc_iid.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_iid=modello_iid[modello_iid["proteina1_iid"] != modello_iid["proteina2_iid"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)


import networkx as nx
import matplotlib.pyplot as plt

G1_iid=nx.from_pandas_edgelist(modello_iid,source="proteina1_iid", target="proteina2_iid") #ho creato il grafo dalla matrice modello 

nodes_iid= list(G1_iid.nodes);#lista delle proteine che fanno parte del grafo
deg_iid = [list(G1_iid.degree())[i][1] for i in range(0, len(nodes_iid))]; #grado di ogni nodo
comm_net_iid= nx.community.louvain_communities(G1_iid, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_iid = {n: i for i, com in enumerate(comm_net_iid) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_iid);
node_color = [community_index_iid[n] for n in G1_iid];
node_size = [v * 500 for v in deg_iid];
nx.draw_networkx(
     G1_iid,
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

num_archi_iid= G1_iid.number_of_edges()
diam_iid= nx.diameter(G1_iid)
mod_iid= nx.community.modularity(G1_iid, comm_net_iid)
dens_iid= nx.density(G1_iid)
grado_medio_iid=sum(deg_iid)/len(deg_iid)
betweenness_iid = nx.betweenness_centrality(G1_iid, normalized=True)
valori_bet_iid= list(betweenness_iid.values()) 
betweenness_medio_iid=sum(valori_bet_iid)/len(valori_bet_iid)
lista_clustering_iid=[nx.clustering(G1_iid, node) for node in G1_iid.nodes()]
clustering_medio_iid=sum(lista_clustering_iid)/len(lista_clustering_iid)
ecc_iid= nx.eccentricity(G1_iid)         # lista delle eccentricità dei nodi
r_iid= nx.radius(G1_iid)                 # raggio (min eccentricity)
center_nodes_iid= nx.center(G1_iid)      # lista dei nodi baricentro
closeness_iid= nx.closeness_centrality(G1_iid)
valori_clos_iid= list(closeness_iid.values())
closeness_media_iid = sum(valori_clos_iid)/len(valori_clos_iid)
