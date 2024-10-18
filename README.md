# Density_Decomposition
This project encompasses an efficient algorithm for computing precise density decomposition, the implementation of various efficient iterative algorithms for hypergraph density decomposition, data collection for assessing approximation qualities, measurement of time and maximum memory usage, and visualization of figures based on the collected data. 

Iterative algorithms implemnted: Elist++, FISTA*, Greedy++, PR, PR_EXP, PR_LIN. 

## Table of Contents

- [Data Preprocessing](#Data-Preprocessing)
- [Exact Decomposition](#Exact-Decomposition)
- [Approximation Algorithms](#Approximation)
- [Data Visualization](#Visualization)
- [License](#license)

## Data-Preprocessing
The code for data preprocessing can be found in the 'data-preprocessing' folder. The input data should be in a bipartite representation of a (hyper)graph, following this format: each line should consist of 'u v', representing the edge (u, v) in the bipartite graph.

Additionally, for a normal graph provided in the format where each line is 'u v', denoting the edge (u, v) in the normal graph, we offer scripts to convert the normal graph into its bipartite representation and its bipartite double cover.

To convert a normal graph file named 'norma.txt' into its bipartite representation, run the 'norm_to_bi.py' script. This will produce a file named 'bipartite_mark.txt'. To generate the bipartite double cover, run the 'norm_to_hyper.py' script, resulting in a file named 'bipartite_hyper.txt'.


To install the project, follow these steps:
1. Clone the repository.
2. Run `npm install` to install dependencies.

## Usage

Here is an example code snippet:
```javascript
console.log('Hello, World!');
