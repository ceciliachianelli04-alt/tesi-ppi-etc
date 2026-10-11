# -*- coding: utf-8 -*-
"""
Created on Sun May 10 09:40:25 2026

@author: cecil
"""

import pandas as pd
# Carico le due matrici
matrice_neuro = pd.read_csv("C:/PPI prediction/BIOGRID_NEURODEGENERATION_dedup.csv")
matrice_etc_estesa = pd.read_csv("C:/PPI prediction/BIOGRID_ETC_ALL_dedup.csv")

# Prendo le proteine della rete neuro
proteine_neuro = set(matrice_neuro["proteina1"]) | set(matrice_neuro["proteina2"])

# Prendo le proteine ETC
proteine_etc_estesa = set(matrice_etc_estesa["source"]) | set(matrice_etc_estesa["target"])

# Seleziono dalla rete ETC estesa le righe dove una proteina è neuro e l'altra è ETC
nuove_interazioni = matrice_etc_estesa[
    (matrice_etc_estesa["source"].isin(proteine_neuro) & matrice_etc_estesa["target"].isin(proteine_etc_estesa)) |
    (matrice_etc_estesa["source"].isin(proteine_etc_estesa) & matrice_etc_estesa["target"].isin(proteine_neuro))
].rename(columns={"source": "proteina1", "target": "proteina2"})

# Unisco le due matrici
matrice_neuro_estesa = pd.concat([matrice_neuro, nuove_interazioni]).drop_duplicates().reset_index(drop=True)
matrice_neuro_estesa[["proteina1", "proteina2"]] = pd.DataFrame(
    matrice_neuro_estesa.apply(lambda row: sorted([row["proteina1"], row["proteina2"]]), axis=1).tolist(),
    index=matrice_neuro_estesa.index
)

# elimino le righe uguali
matrice_neuro_estesa = matrice_neuro_estesa.drop_duplicates().reset_index(drop=True)
matrice_neuro_estesa.to_csv("C:/PPI prediction/BIOGRID_NEURO_ESTESO.csv", index=False)

import networkx as nx
import matplotlib.pyplot as plt

G1_bio_neuro_esteso1=nx.from_pandas_edgelist(matrice_neuro_estesa,source="proteina1", target="proteina2") #ho creato il grafo dalla matrice modello 

nodes_bio_neuro_esteso1 = list(G1_bio_neuro_esteso1.nodes);#lista delle proteine che fanno parte del grafo
deg_bio_neuro_esteso1= [list(G1_bio_neuro_esteso1.degree())[i][1] for i in range(0, len(nodes_bio_neuro_esteso1))]; #grado di ogni nodo
comm_net_bio_neuro_esteso1= nx.community.louvain_communities(G1_bio_neuro_esteso1, weight='weight', resolution=1, threshold=1e-07, seed=None); #identifica le comunità del grafo
community_index_bio_neuro_esteso1 = {n: i for i, com in enumerate(comm_net_bio_neuro_esteso1) for n in com}; 

 
fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G1_bio_neuro_esteso1);
node_color = [community_index_bio_neuro_esteso1[n] for n in G1_bio_neuro_esteso1];
node_size = [v * 500 for v in deg_bio_neuro_esteso1];
nx.draw_networkx(
     G1_bio_neuro_esteso1,
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

num_archi_bio_neuro_esteso1= G1_bio_neuro_esteso1.number_of_edges()
diam_bio_neuro_esteso1= nx.diameter(G1_bio_neuro_esteso1)
mod_bio_neuro_esteso1= nx.community.modularity(G1_bio_neuro_esteso1, comm_net_bio_neuro_esteso1)
dens_bio_neuro_esteso1= nx.density(G1_bio_neuro_esteso1)
grado_medio_bio_neuro_esteso1=sum(deg_bio_neuro_esteso1)/len(deg_bio_neuro_esteso1)
lista_clustering_bio_neuro_esteso1=[nx.clustering(G1_bio_neuro_esteso1, node) for node in G1_bio_neuro_esteso1.nodes()]
clustering_medio_bio_neuro_esteso1=sum(lista_clustering_bio_neuro_esteso1)/len(lista_clustering_bio_neuro_esteso1)
clustering_medio_bio_neuro_esteso1=sum(lista_clustering_bio_neuro_esteso1)/len(lista_clustering_bio_neuro_esteso1)
betweenness_bio_neuro_esteso1 = nx.betweenness_centrality(G1_bio_neuro_esteso1, normalized=True) 
valori_bet_bio_neuro_esteso1= list(betweenness_bio_neuro_esteso1.values()) 
betweenness_medio_bio_neuro_esteso1=sum(valori_bet_bio_neuro_esteso1)/len(valori_bet_bio_neuro_esteso1)
ecc_bio_neuro_esteso1= nx.eccentricity(G1_bio_neuro_esteso1)         # lista delle eccentricità dei nodi
r_bio_neuro_esteso1= nx.radius(G1_bio_neuro_esteso1)                 # raggio (min eccentricity)
center_nodes_bio_neuro_esteso1= nx.center(G1_bio_neuro_esteso1)      # lista dei nodi baricentro
closeness_bio_neuro_esteso1= nx.closeness_centrality(G1_bio_neuro_esteso1)
valori_clos_bio_neuro_esteso1= list(closeness_bio_neuro_esteso1.values())
closeness_media_bio_neuro_esteso1 = sum(valori_clos_bio_neuro_esteso1)/len(valori_clos_bio_neuro_esteso1) 

with open("lista_centri_neuro_esteso_etcest.txt", "w") as f:
    for elemento in center_nodes_bio_neuro_esteso1:
        f.write(str(elemento) + "\n")


