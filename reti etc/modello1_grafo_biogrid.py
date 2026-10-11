# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""
import pandas as pd

file_name = "C:/Users/cecil/OneDrive - alcampus.it/UNI/tirocinio/BIOGRID-ALL-5.0.250.tab.txt"

prova = pd.read_csv(file_name,sep="\t") #aprire il file in pandas colonnato (sto mettendo il file come input)

prova = prova[["OFFICIAL_SYMBOL_A","OFFICIAL_SYMBOL_B"]] # sto estraendo le colonne che ci interessano
prova=prova.set_axis(["proteina1","proteina2"],axis=1) #con axis=1 vuol dire che sto rinominando le colonne


file_name = "C:/Users/cecil/OneDrive - alcampus.it/UNI/tirocinio/mitochondrial_ETC-lista-proteine.txt"
lista_pro= pd.read_csv(file_name, names=["proteine"])

#devo selezionare le righe della matrice che hanno entrambe le proteine che fanno parte della lista proteine
proteine = set(lista_pro["proteine"])#ho trasformato il dataframe in una lista 
mat_def = prova[prova["proteina1"].isin(proteine) & prova["proteina2"].isin(proteine)].reset_index(drop=True) #ho selezionato le protein dalla matrice prova degli indici che ho trovato

#dobbiamo fare la lista definitiva con le interazioni bidirezionali, ad esempio se c'è un'interazione proteina x e proteina y, ci deve essere anche un'interazione proteina y proteina x, quindi dovremo avere due righe con le stesse proteine, ma in ordine invertito
#dobbiamo quindi creare un'altra mat def con gli stessi elementi della mat def originale ma con le colonne invertite, per poi fare un confronto tra la nuova lista def e quella originale e elimare le righe uguali (perché nella lista def originale possono comunque esserci delle proteine con un'interazione bidirezionale)
mat_def2=mat_def[mat_def.columns[::-1]]#sto scambiando le colonne della matrice mat_def
mat_def2=mat_def2.set_axis(["proteina1","proteina2"],axis=1)

mat_conc=pd.concat([mat_def,mat_def2], ignore_index=True)#ho ottenuto una matrice con il doppio delle righe di mat_def

modello_biogrid=mat_conc.drop_duplicates(ignore_index=True)#sto eliminando le righe uguali
modello_biogrid=modello_biogrid[modello_biogrid["proteina1"] != modello_biogrid["proteina2"]].reset_index(drop=True)#sto elimando le righe con gli stessi elemnti (self interactions)
modello_biogrid.to_csv('BIOGRID_ETC.csv', sep=',', index=False)


import networkx as nx
import matplotlib.pyplot as plt

G1_biogrid=nx.from_pandas_edgelist(modello_biogrid,source="proteina1", target="proteina2") #ho creato il grafo dalla matrice modello 

nodes_biogrid = list(G1_biogrid.nodes);#lista delle proteine che fanno parte del grafo
deg_biogrid= [list(G1_biogrid.degree())[i][1] for i in range(0, len(nodes_biogrid))]; #grado di ogni nodo
comm_net_biogrid= nx.community.louvain_communities(G1_biogrid, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_biogrid = {n: i for i, com in enumerate(comm_net_biogrid) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_biogrid);
node_color = [community_index_biogrid[n] for n in G1_biogrid];
node_size = [v * 500 for v in deg_biogrid];
nx.draw_networkx(
     G1_biogrid,
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

num_archi_biogrid= G1_biogrid.number_of_edges()
diam_biogrid= nx.diameter(G1_biogrid)
mod_biogrid= nx.community.modularity(G1_biogrid, comm_net_biogrid)
dens_biogrid= nx.density(G1_biogrid)
grado_medio_biogrid=sum(deg_biogrid)/len(deg_biogrid)
lista_clustering_biogrid=[nx.clustering(G1_biogrid, node) for node in G1_biogrid.nodes()]
clustering_medio_biogrid=sum(lista_clustering_biogrid)/len(lista_clustering_biogrid)
clustering_medio_biogrid=sum(lista_clustering_biogrid)/len(lista_clustering_biogrid)
betweenness_biogrid = nx.betweenness_centrality(G1_biogrid, normalized=True) 
valori_bet_biogrid= list(betweenness_biogrid.values()) 
betweenness_medio_biogrid=sum(valori_bet_biogrid)/len(valori_bet_biogrid)
ecc_biogrid= nx.eccentricity(G1_biogrid)         # lista delle eccentricità dei nodi
r_biogrid= nx.radius(G1_biogrid)                 # raggio (min eccentricity)
center_nodes_biogrid= nx.center(G1_biogrid)      # lista dei nodi baricentro
closeness_biogrid= nx.closeness_centrality(G1_biogrid)
valori_clos_biogrid= list(closeness_biogrid.values())
closeness_media_biogrid = sum(valori_clos_biogrid)/len(valori_clos_biogrid)






