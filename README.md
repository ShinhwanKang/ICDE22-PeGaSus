# Personalized Graph Summarization: Formulation, Scalable Algorithms, and Applications

## Overview
- Problem Formulation: We introduce a new problem, personalized graph summarization, and demonstrate the usefulness of personalizing summary graphs.

- Algorithm Design: We propose PeGaSus (**Pe**rsonalized **G**r**a**ph **Su**mmarization with **S**calability), a linear-time algorithm for the problem. We show empirically that it scales to a graph with one billion edges.

- Extensive Experiments: We exhibit the effectiveness of PeGaSus, and its applicability to distributed multi-query answering using six real-world graphs.

## Code
The algorithms used in the paper is available at ```./PeGaSus/```.


## Datasets

|Name|#Nodes|#Edges|Summary|Source|Download|
|:---:|:---:|:---:|:---:|:---:|:---:|
|LastFM-Asia (LA)|7,624|27,806|Social|[SNAP](https://snap.stanford.edu/data/feather-lastfm-social.html)|[LINK](https://snap.stanford.edu/data/lastfm_asia.zip)|
|Caida (CA)|26,475|53,381|Internet|[SNAP](https://snap.stanford.edu/data/as-Caida.html)|[LINK](https://snap.stanford.edu/data/as-caida.tar.gz)|
|DBLP (DB)|317,080|1,049,866|Collaboration|[SNAP](https://snap.stanford.edu/data/com-DBLP.html)|[LINK](https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz)|
|Amazon0601 (A6)|403,364|2,443,311|Co-purchase|[SNAP](https://snap.stanford.edu/data/amazon0601.html)|[LINK](https://snap.stanford.edu/data/amazon0601.txt.gz)|
|Skitter (SK)|1,694,616|11,094,209|Internet|[SNAP](https://snap.stanford.edu/data/as-Skitter.html)|[LINK](https://snap.stanford.edu/data/as-skitter.txt.gz)|
|Wikipedia (WK)|3,174,745|103,310,688|Hyperlinks|[KONNECT](http://konect.cc/networks/wikipedia_link_sr/)|[LINK](http://konect.cc/files/download.tsv.wikipedia_link_sr.tar.bz2)|

