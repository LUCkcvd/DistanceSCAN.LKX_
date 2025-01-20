# DistanceSCAN

DistanceSCAN is an efficient algorithm for distance-based structural graph clustering, implementing novel techniques for processing large-scale graph data. This implementation is based on the research paper published in SIGMOD 2023.

## Features

- **Distance-based Structural Graph Clustering**: Implements advanced clustering techniques that consider both structural similarities and distances between nodes
- **Graph Processing Capabilities**:
  - Support for weighted/unweighted and directed/undirected graphs
  - Dynamic graph operations (edge insertion, deletion, updates)
  - Multiple similarity measures:
    * Jaccard similarity
    * Cosine similarity
    * Set containment measures
  - Advanced graph reweighting schemes:
    * Logarithmic scaling
    * Sum-weight normalization
    * Max-weight normalization
  - Memory-efficient graph processing:
    * Batched edge processing
    * Pre-allocated data structures
    * Automatic memory monitoring
    * Dynamic cleanup for large graphs
- **Multiple Algorithm Implementations**:
  - `distancescan`: Our primary method using optimized sketches with max-heaps
  - `distancescan_pst`: Variant using persist search trees for all-distances sketches
  - `scan`: Basic SCAN implementation for the distance-SCAN problem
  - `pscan`: PSCAN implementation adapted for distance-SCAN
  - `exact`: Exact algorithm implementation
  - `wscan`: Weighted SCAN variant

- **Operations**:
  1. **Sketch Construction**: Build all-distances bottom-k sketches for efficient similarity computation
     - Optimized bottom-k sketch construction and serialization
     - Histogram-based approximation with configurable bin sizes
     - Double hashing for efficient storage
     - Batch processing with automatic memory management
     - Chunked histogram storage for large graphs
     - Real-time progress tracking with ETA estimation
  2. **Clustering Queries**: Execute clustering operations with various parameters
  3. **Cluster Validation**: Evaluate clustering quality using modularity metrics
  4. **Quality Validation**: Compare results between different algorithms
  5. **Graph Conversion**: Transform graphs to undirected similarity-weighted format
  6. **Performance Analysis**: Built-in performance monitoring
     - Memory usage tracking
     - Operation timing measurements
     - Efficiency metrics collection

## Requirements

- Ubuntu (or Linux-based OS)
- C++ 14 or higher
- GCC 4.8 or higher
- Boost libraries
- CMake 3.9 or higher
- OpenMP support

## Installation

```sh
$ cmake .
$ make
```

## Usage

### 1. Constructing Sketches

```sh
./Distance_SCAN_SIGMOD --operation construct_sketches --algo <algorithm> [options]
```

Algorithm options:
- `distancescan_pst`: Uses persist search trees for storing all-distances sketches
- `distancescan`: Uses max-heaps for storing all-distances sketches

Required parameters:
- `--prefix <prefix>`: Output prefix for generated files
- `--dataset <dataset>`: Input dataset name
- `-d <value>`: Maximum distance threshold
- `-k <value>`: Number of samples in bottom-k sketches

Example:
```sh
$ ./Distance_SCAN_SIGMOD --operation construct_sketches --dataset ego-facebook --algo distancescan -k 16 -d 0.4 
```

Note: The sketch construction process includes several optimizations:
- Batch processing of nodes to control memory usage
- Automatic memory cleanup after each batch
- Progress tracking with edge processing rate
- Peak memory usage monitoring
- OpenMP parallel processing support
- Efficient disk I/O with chunked storage

### 2. Querying Clusters

```sh
./Distance_SCAN_SIGMOD --operation query --algo <algorithm> [options]
```

Algorithm options:
- `scan`: Basic SCAN implementation
- `pscan`: PSCAN implementation
- `exact`: Exact algorithm
- `distancescan`: Our optimized method
- `wscan`: Weighted SCAN variant

Required parameters:
- `--prefix <prefix>`: Input/output prefix
- `--dataset <dataset>`: Dataset name
- `-d <value>`: Distance threshold
- `-m <value>`: Maximum distance threshold for sketches
- `-k <value>`: Number of samples in bottom-k sketches
- `-u <value>`: Threshold for structurally similar neighbors (μ)
- `-e <value>`: Similarity threshold (ε)

Example:
```sh
$ ./Distance_SCAN_SIGMOD --operation query --dataset ego-facebook --algo distancescan -k 16 -u 5 -e 0.2 -d 0.3 -m 0.4 
```

Performance optimizations in the query process:
- Early termination optimizations for similarity computations
- Efficient bottom-k sketch comparison algorithms
- Optimized Jaccard similarity calculations
- Memory-efficient histogram-based filtering

### 3. Validation Operations

The implementation supports two types of validation:
- `cluster_validation`: Evaluates clustering quality using modularity metrics
- `quality_validation`: Compares results between different algorithms

Advanced validation features include:
- Adjusted Rand Index (ARI) computation for cluster comparison
- Modularity-based quality assessment
- Comprehensive performance metrics
- Core/non-core node classification analysis
- Cluster quality comparison across different algorithms

## Input Graph Format

The implementation accepts graph input files in the following formats:
- Edge list format (source_node target_node [weight])
- Undirected graph files (automatically handles symmetry)
- Weighted graph files (supports both similarity and distance weights)
- Attribute files for node metadata (optional)

Implementation optimizations:
- Buffered I/O for efficient graph loading
- OpenMP parallel processing for edge processing
- Chunked memory allocation to control memory usage
- Efficient data structures for sparse graphs
- Dynamic memory cleanup for processing large graphs

## Supported Datasets

The implementation works with graph datasets from:
- [SNAP](https://snap.stanford.edu/data/)
- [LAW](http://law.di.unimi.it/datasets.php)
- [AMiner](https://www.aminer.cn/data/?nav=openData#Topic-coauthor)

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@inproceedings{DBLP:***,
  author    = {Kaixin Liu and
               Sibo Wang and
               Yong Zhang and
               Chunxiao Xing},
  title     = {An Efficient Algorithm for Distance-based Structural Graph Clustering},
  journal   = {PACMMOD},
  volume    = {1},
  number    = {45},
  pages     = {**--**},
  year      = {2023},
  doi       = {10.1145/3588725},
}
```

## Contact

For questions or issues, please contact Kaixin Liu at lkx17@mails.tsinghua.edu.cn.
