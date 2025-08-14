# Density_Decomposition
This project encompasses an efficient algorithm for computing precise density decomposition, the implementation of various efficient iterative algorithms for hypergraph density decomposition, data collection for assessing approximation qualities, measurement of time and maximum memory usage, and visualization of figures based on the collected data. 

Iterative algorithms implemnted: Elist++, FISTA*, Greedy++, PR, PR_EXP, PR_LIN. 

The codes were tested on Linux system. 

## Table of Contents

- [Data Preprocessing](#Data-Preprocessing)
- [Exact Decomposition](#Exact-Decomposition)
- [Approximation Algorithms](#Approximation)
- [Data Visualization](#Visualization)
- [Time and Memory](#Time_Mem)
- [License](#license)

## Data Preprocessing
The code for data preprocessing can be found in the 'data-preprocessing' folder. Our data is sourced from https://snap.stanford.edu/. While certain graphs are immediately usable, others require the removal of irrelevant text from the initial section of the downloaded file. To achieve this, you can utilize a Python script located at data-preprocessing/process.py. 

The original data downloaded from https://snap.stanford.edu/ is a normal graph. After processing, it should be in a bipartite representation of a (hyper)graph, following this format: each line should consist of 'u v', representing the edge (u, v) in the bipartite graph.

For a normal graph provided in the format where each line is 'u v', denoting the edge (u, v) in the normal graph, we offer scripts to convert the normal graph into its bipartite representation and its bipartite double cover. 

To convert a normal graph file named 'graph.txt'(we have provided a test sample) into its bipartite representation, run the 'norm_to_bi.py' script. This will produce a file named 'bipartite_mark.txt'. To generate the bipartite double cover, run the 'norm_to_hyper.py' script, resulting in a file named 'bipartite_hyper.txt'.

## Exact Decomposition
Bipartite representation: To calculate the exact decomposition for bipartite representation of normal graphs, refer to the code in /normal/decomposition.cpp. This file is written in C++. Prior to compiling the code, ensure that you have installed the Boost package (available at https://www.boost.org/). Place the input file named "bipartite_mark.txt" in the same directory as the compiled code, then proceed to execute the code. Upon execution, two files will be generated: "maxflow_res_large.txt," which contains the decomposition, and "Exact_normal.txt," detailing the information about the density vector of each hypernode. "Exact_normal.txt" will be utilized in evaluating the approximation quality.

Bipartite double cover: For calculating the exact decomposition for bipartite double cover of normal graphs, locate the code in /hyper/decomposition.cpp. This file is also written in C++. Before compiling the code, ensure that you have the Boost package installed (accessible at https://www.boost.org/). Place the input file named "bipartite_hyper.txt" in the same directory as the compiled code, then run the code. Upon completion, two files will be created: "maxflow_res_large.txt," which includes the decomposition, and "Exact_hyper.txt," which provides information about the density vector of each hypernode. "Exact_hyper.txt" will be essential for assessing the quality of the approximation.


## Approximation Algorithms
Bipartite Representation:
To execute the iterative approximation algorithms and assess the quality of the approximation, compile the code located at /normal/iterative_normal.cpp. Place the original input file named "bipartite_mark.txt" and the density vector file named "Exact_normal.txt" in the same directory as the compiled file. Subsequently, create a folder named "data_normal". You can then run the compiled file, and the resulting data will be stored inside the "data_normal" folder. The data for all the algorithms will be created. 

Bipartite Double Cover:
For running the iterative approximation algorithms and evaluating the approximation quality, compile the code found at /hyper/iterative_hyper.cpp. Ensure that the original input file named "bipartite_hyper.txt" and the density vector file named "Exact_hyper.txt" are in the same directory as the compiled file. Create a folder named "data_hyper", and after that, execute the compiled file. The data generated will be stored within the "data_hyper" folder. The data for all the algorithms will be created. 

Within the data folder, the following conventions apply: "abs" denotes the global error, "mul" signifies the multiplicative error, and "inv" indicates the number of inversions.

## Data Visualization
Bipartite Representation:
For generating figures, place the density vector file "Exact_normal.txt" in the same directory as "draw_all_normal.py" located at /normal/draw_all_normal.py. Create a folder named "figures" in this directory, then execute draw_all_normal.py. The resulting figures will be displayed in the "figures" folder.

Bipartite Double Cover:
To create figures, position the density vector file "Exact_hyper.txt" in the same directory as "draw_all_hyper.py" located at /hyper/draw_all_hyper.py. Establish a folder named "figures" in this directory, then run draw_all_hyper.py. The figures will be presented in the "figures" folder.

## Time and Memory
For any given graph, place its bipartite representation in the folder /time_mem/normal/ and its bipartite double cover in the folder /time_mem/hyper/ with the same file name.

To evaluate the average time, compile the file located at /time_mem/normal/test_time.cpp as test_time without changing its location. Create a directory named "data_normal" at /time_mem/normal/. Similarly, compile the file /time_mem/hyper/test_time_hyper.cpp as test_time_hyper without relocating it. Create a folder named "data_hyper" at /time_mem/hyper/. Next, open the file time_mem/run_program_time.sh and modify the variable params to match the dataset file names inside /time_mem/normal/ or /time_mem/hyper/. Run it, subsequently, the data will appear in "data_normal" and "data_hyper", respectively. Create a "figures" folder within "data_normal" or "data_hyper," execute plot.py, and the figures will be displayed inside the "figures" directory.

To assess the maximum memory usage, compile the file found at /time_mem/normal/test_memory.cpp as test_mem without changing its location. Establish a folder named "data_normal_mem" at /time_mem/normal/. Execute the file located at /time/mem/normal/run_memory.sh. In this file, the variable params stores the list of file names of the datasets, while algos holds the list of algorithms to test; these can be adjusted as necessary. Upon completion, the data will be visible within the "data_normal_mem" folder. Create a "figures" directory within it, run the script plot.py, and the resulting figures will be displayed. The process for the bipartite double cover version follows a similar sequence.

