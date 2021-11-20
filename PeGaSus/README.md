# PeGaSus

## Requirements

* \>= OpenJDK 12

## Input Format

PeGaSus assumes that the input graph *G = (**V**, **E**)* is undirected without self-loops.
The format of the example file is given below.
Each line represents an edge.
Each edge *{u, v} ∈ **E*** joins two distinct nodes *u != v ∈ **V***, separated by a tab.
Each node *v ∈ **V*** is assigned to a unique integer id.

### Example Input Format
```
    0   1
    2   3
    3   4
```
- The example consists of 5 nodes with 3 edges.

## Target Node Set Format

PeGaSus personalizes a summary graph for a target node set *T $\subset$ **V***.
The format of the example file is given below.
Each line represents a node *v ∈ **V*** with a unique integer id.
### Example Input Format
```
    4721
    274
    2195
```
- The example target node set consists of 3 nodes.

## Execution

```
java -jar PeGaSus.jar [data path] [target compression ratio] [target node set file path] [checking personalized error] [saving summary graph] [alpha]
```

### argument
- data path: path to the input text file
- target compression ratio [0,1]: desired size of a summary graph compared relative to the input graph size in bits
- target node set file path: path to the target node set file
- checking personalized error: if *true* compute the personalized error. Otherwise, not.
- saving summary graph: if *true* save a summary graph. Otherwise, not.
- alpha: the degree of personalization

```
java -jar PeGaSus.jar ./data/lastfm_asia.txt 0.5 ./data/lastfm_asia_TNS.txt true true 1.25   
```
## Output Format

The output file contains information about nodes (nodes in *G*) belonging to each
supernode *s ∈ **S*** of the output summary graph *$\bar{G}$ = (**S**, **P**)* and information about each superedge *p ∈ **P***. 
The first integer on each line following the line "\<Node of Each Supernode\>" represents the id of the supernode, and the following integers separated by tabs represent the ids of the subnodes belonging to that supernode.
Each line following  the line "\<Superedge Info\>" represents a single superedge. 
The two integers separated by tabs represent the id of the source supernode and the id of the destination supernode.

* Output
```
    dataPath: ./lastfm_asia.txt
    Target Ratio: 0.1
    Target Node Set File Name: ./lastfm_asia_TNS.txt
    is Checking Personalized Error: true
    is Saving Summary Graph: true
    Alpha: 1.25
    |V|:    7624
    |E|:    27806
    TopNDrop
    Compression Ratio       0.1     Execution time   0.788s
    Personalized Error      0.001089218172509844
```


* Output file
```
    <Node of Each Supernode>
    2569	3145	2225	3975	2306    498 428 1464
    2775	5346	2162	4012	3774
    .
    .
    <Superedge Info>
    2569	2569
    2775	2775
    1481	1357
    .
    .
```




