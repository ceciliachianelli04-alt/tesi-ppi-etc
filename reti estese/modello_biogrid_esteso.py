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

#devo selezionare le righe dellla matrice [28190102x2] che hanno entrambe le proteine che fanno parte della lista proteine
proteine = set(lista_pro["proteine"])#ho trasformato il dataframe in una lista 
modello_bio_etc_esteso= prova[prova["proteina1"].isin(proteine) | prova["proteina2"].isin(proteine)].reset_index(drop=True) #ho selezionato le protein dalla matrice prova degli indici che ho trovato
modello_bio_etc_esteso.to_csv('BIOGRID_ETC_ALL.csv', sep=',', index=False)

import networkx as nx
import matplotlib.pyplot as plt

G1_bio_etc_esteso=nx.from_pandas_edgelist(modello_bio_etc_esteso,source="proteina1", target="proteina2") #ho creato il grafo dalla matrice modello 

nodes_bio_etc_esteso = list(G1_bio_etc_esteso.nodes);#lista delle proteine che fanno parte del grafo
deg_bio_etc_esteso= [list(G1_bio_etc_esteso.degree())[i][1] for i in range(0, len(nodes_bio_etc_esteso))]; #grado di ogni nodo
comm_net_bio_etc_esteso= nx.community.louvain_communities(G1_bio_etc_esteso, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_bio_etc_esteso = {n: i for i, com in enumerate(comm_net_bio_etc_esteso) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_bio_etc_esteso, k=0.03);
node_color = [community_index_bio_etc_esteso[n] for n in G1_bio_etc_esteso];
node_size = [v * 500 for v in deg_bio_etc_esteso];
nx.draw_networkx(
     G1_bio_etc_esteso,
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

num_archi_bio_etc_esteso= G1_bio_etc_esteso.number_of_edges()
diam_bio_etc_esteso= nx.diameter(G1_bio_etc_esteso)
mod_bio_etc_esteso= nx.community.modularity(G1_bio_etc_esteso, comm_net_bio_etc_esteso)
dens_bio_etc_esteso= nx.density(G1_bio_etc_esteso)
grado_medio_bio_etc_esteso=sum(deg_bio_etc_esteso)/len(deg_bio_etc_esteso)
lista_clustering_bio_etc_esteso=[nx.clustering(G1_bio_etc_esteso, node) for node in G1_bio_etc_esteso.nodes()]
clustering_medio_bio_etc_esteso=sum(lista_clustering_bio_etc_esteso)/len(lista_clustering_bio_etc_esteso)
clustering_medio_bio_etc_esteso=sum(lista_clustering_bio_etc_esteso)/len(lista_clustering_bio_etc_esteso)
betweenness_bio_etc_esteso = nx.betweenness_centrality(G1_bio_etc_esteso, normalized=True) 
valori_bet_bio_etc_esteso= list(betweenness_bio_etc_esteso.values()) 
betweenness_medio_bio_etc_esteso=sum(valori_bet_bio_etc_esteso)/len(valori_bet_bio_etc_esteso)
ecc_bio_etc_esteso= nx.eccentricity(G1_bio_etc_esteso)         # lista delle eccentricità dei nodi
r_bio_etc_esteso= nx.radius(G1_bio_etc_esteso)                 # raggio (min eccentricity)
center_nodes_bio_etc_esteso= nx.center(G1_bio_etc_esteso)      # lista dei nodi baricentro
closeness_bio_etc_esteso= nx.closeness_centrality(G1_bio_etc_esteso)
valori_clos_bio_etc_esteso= list(closeness_bio_etc_esteso.values())
closeness_media_bio_etc_esteso = sum(valori_clos_bio_etc_esteso)/len(valori_clos_bio_etc_esteso) 

with open("lista_centri_etc_esteso.txt", "w") as f:
    for elemento in center_nodes_bio_etc_esteso:
        f.write(str(elemento) + "\n")

lista_pro_clean = lista_pro.iloc[:, 0].str.strip().tolist()
center_nodes_clean = [x.strip() for x in center_nodes_bio_etc_esteso]

# Trova gli elementi in comune
centri_etc_esteso = list(set(lista_pro_clean) & set(center_nodes_clean))

with open("lista_centri_etc_esteso_etc.txt", "w") as f:
    f.write(", ".join(str(elemento) for elemento in centri_etc_esteso))