# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 14:55:51 2025

@author: hp
"""

import pandas as pd

file_name = "C:/Users/hp/CampusOneDrive/OneDrive - alcampus.it/UNI/tirocinio/STRING_interactions.tsv"

prova = pd.read_csv(file_name,sep="\t") #aprire il file in pandas colonnato (sto mettendo il file come input)

prova = prova[["#node1","node2"]] # sto estraendo le colonne che ci interessano
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

modello_string=mat_conc.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_string=modello_string[modello_string["proteina1"] != modello_string["proteina2"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)
modello_string.to_csv('STRING_ETC.csv', sep=',', index=False)

import networkx as nx
import matplotlib.pyplot as plt

G1_string=nx.from_pandas_edgelist(modello_string,source="proteina1", target="proteina2") #ho creato il grafo dalla matrice modello 

nodes_string = list(G1_string.nodes);#lista delle proteine che fanno parte del grafo
deg_string= [list(G1_string.degree())[i][1] for i in range(0, len(nodes_string))]; #grado di ogni nodo
comm_net_string = nx.community.louvain_communities(G1_string, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_string = {n: i for i, com in enumerate(comm_net_string) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_string);
node_color = [community_index_string[n] for n in G1_string];
node_size = [v * 500 for v in deg_string];
nx.draw_networkx(
     G1_string,
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

num_archi_string = G1_string.number_of_edges()
diam_string= nx.diameter(G1_string)
mod_string= nx.community.modularity(G1_string, comm_net_string)
dens_string= nx.density(G1_string)
grado_medio_string=sum(deg_string)/len(deg_string)
lista_clustering_string=[nx.clustering(G1_string, node) for node in G1_string.nodes()]
clustering_medio_string=sum(lista_clustering_string)/len(lista_clustering_string)
betweenness_string = nx.betweenness_centrality(G1_string, normalized=True)
valori_bet_string= list(betweenness_string.values())
betweenness_medio_string=sum(valori_bet_string)/len(valori_bet_string)
ecc_string= nx.eccentricity(G1_string)         # lista delle eccentricità dei nodi
r_string= nx.radius(G1_string)                 # raggio (min eccentricity)
center_nodes_string= nx.center(G1_string)      # lista dei nodi baricentro
closeness_string= nx.closeness_centrality(G1_string)
valori_clos_string= list(closeness_string.values())
closeness_media_string = sum(valori_clos_string)/len(valori_clos_string)