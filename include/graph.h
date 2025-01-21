#ifndef SCAN_GRAPH_H
#define SCAN_GRAPH_H

#include "mylib.h"
#include "config.h"
#include <memory>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

// Memory pool for efficient small allocations
class MemoryPool {
private:
    static constexpr size_t CHUNK_SIZE = 1024 * 1024; // 1MB chunks
    vector<unique_ptr<char[]>> chunks;
    size_t current_offset = 0;
    size_t current_chunk = 0;

public:
    void* allocate(size_t size) {
        if (current_chunk >= chunks.size() || current_offset + size > CHUNK_SIZE) {
            chunks.push_back(make_unique<char[]>(CHUNK_SIZE));
            current_chunk = chunks.size() - 1;
            current_offset = 0;
        }
        void* ptr = &chunks[current_chunk][current_offset];
        current_offset += size;
        return ptr;
    }
    
    void reset() {
        chunks.clear();
        current_offset = 0;
        current_chunk = 0;
    }
};

// Compressed Sparse Row format for efficient graph storage
struct CSRGraph {
    vector<size_t> offsets;  // Size n+1, stores start positions of each vertex's edges
    vector<int> edges;       // Size m, stores destination vertices
    vector<double> weights;  // Size m, stores edge weights (if weighted)
    
    void reserve(size_t num_vertices, size_t num_edges) {
        offsets.reserve(num_vertices + 1);
        edges.reserve(num_edges);
        weights.reserve(num_edges);
    }
    
    void clear() {
        vector<size_t>().swap(offsets);
        vector<int>().swap(edges);
        vector<double>().swap(weights);
    }
};

class Graph {
private:
    static constexpr size_t MAX_MEMORY_BYTES = 45ULL * 1024 * 1024 * 1024; // 45GB
    size_t current_memory_usage = 0;
    
    // Core graph structure in CSR format
    CSRGraph csr_graph;
    
    // Memory mapped file support
    int mmap_fd = -1;
    void* mmap_data = nullptr;
    size_t mmap_size = 0;
    
    // Memory pool for temporary allocations
    MemoryPool mem_pool;
    
    // Lazy loaded data structures
    class LazyMatrix {
        vector<unordered_map<int, double>> data;
        bool is_loaded = false;
        Graph* graph = nullptr;
        
    public:
        LazyMatrix() = default;
        
        LazyMatrix& operator=(const vector<unordered_map<int, double>>& other) {
            data = other;
            is_loaded = true;
            return *this;
        }
        
        void init(size_t size) {
            if (!is_loaded) {
                data.clear();
                data.resize(size);
                is_loaded = true;
            }
        }
        
        void clear() {
            vector<unordered_map<int, double>>().swap(data);
            is_loaded = false;
        }
        
        operator const vector<unordered_map<int, double>>&() const {
            if (!is_loaded) {
                throw runtime_error("Accessing uninitialized LazyMatrix");
            }
            return data;
        }
        
        unordered_map<int, double>& operator[](size_t idx) {
            if (!is_loaded) {
                init(graph->n);
            }
            return data[idx];
        }
        
        const unordered_map<int, double>& operator[](size_t idx) const {
            if (!is_loaded) {
                throw runtime_error("Accessing uninitialized LazyMatrix");
            }
            return data[idx];
        }
        
        bool empty() const {
            return !is_loaded || data.empty();
        }
        
        void set_graph(Graph* g) { graph = g; }

        friend void swap(LazyMatrix& first, LazyMatrix& second) noexcept {
            using std::swap;
            swap(first.data, second.data);
            swap(first.is_loaded, second.is_loaded);
            swap(first.graph, second.graph);
        }
    };
    
    LazyMatrix path_weight;
    LazyMatrix similarity_matrix;
    LazyMatrix jac_res;
    
    // Similarity matrix for boolean values
    vector<unordered_map<int, bool>> similarity;
    
    // Accessors for CSR graph
    vector<int> get_neighbors(int v) const {
        size_t start = csr_graph.offsets[v];
        size_t end = csr_graph.offsets[v + 1];
        vector<int> neighbors;
        neighbors.reserve(end - start);
        for (size_t i = start; i < end; i++) {
            neighbors.push_back(csr_graph.edges[i]);
        }
        return neighbors;
    }
    
    double get_edge_weight(int u, int v) const {
        size_t start = csr_graph.offsets[u];
        size_t end = csr_graph.offsets[u + 1];
        for (size_t i = start; i < end; i++) {
            if (csr_graph.edges[i] == v) {
                return csr_graph.weights[i];
            }
        }
        return 0.0;
    }
    
    void set_edge_weight(int u, int v, double w) {
        size_t start = csr_graph.offsets[u];
        size_t end = csr_graph.offsets[u + 1];
        for (size_t i = start; i < end; i++) {
            if (csr_graph.edges[i] == v) {
                csr_graph.weights[i] = w;
                break;
            }
        }
    }
    
    bool check_and_update_memory(size_t additional_bytes);
    void track_memory_usage();
    void release_memory();
    
    // Memory mapping helpers
    bool map_file(const string& filename);
    void unmap_file();

public:
    int n;
    long long m;
    
    string data_folder;
    
    bool directed = false;
    bool weighted = true;
    
    // Core vertex data
    vector<int> clusterID;  // -1 not classified, -2 hub, -3 outlier
    vector<int> is_core;    // -1 init, 0 non-member, 1 core, 2 hub, 3 outlier
    vector<double> weight_degree;
    double whole_weight;
    
    // Temporary working sets
    vector<vector<int>> d_neighbors;  // Only allocated during specific algorithms
    
    // Compatibility layer for old code
    class AdjacencyList {
        Graph& graph;
        mutable vector<vector<int>> temp_lists;
    public:
        AdjacencyList() = delete;
        AdjacencyList(Graph& g) : graph(g) {
            temp_lists.resize(g.n);
        }
        
        AdjacencyList& operator=(const vector<vector<int>>& other) {
            temp_lists = other;
            return *this;
        }
        
        void resize(size_t size) {
            temp_lists.resize(size);
        }
        
        vector<int>& operator[](int v) {
            temp_lists[v] = graph.get_neighbors(v);
            return temp_lists[v];
        }
        
        const vector<int>& operator[](int v) const {
            temp_lists[v] = graph.get_neighbors(v);
            return temp_lists[v];
        }
        
        size_t size() const { return graph.n; }

        friend void swap(AdjacencyList& first, vector<vector<int>>& second) noexcept {
            using std::swap;
            swap(first.temp_lists, second);
        }
        
        friend void swap(vector<vector<int>>& first, AdjacencyList& second) noexcept {
            using std::swap;
            swap(first, second.temp_lists);
        }
        
        friend void swap(AdjacencyList& first, AdjacencyList& second) noexcept {
            using std::swap;
            // Note: Don't swap graph references as they should stay with their original objects
            swap(first.temp_lists, second.temp_lists);
        }
    };
    
    class EdgeWeights {
        Graph& graph;
        mutable vector<unordered_map<int, double>> temp_weights;
    public:
        EdgeWeights() = delete;
        EdgeWeights(Graph& g) : graph(g) {
            temp_weights.resize(g.n);
        }
        
        EdgeWeights& operator=(const vector<unordered_map<int, double>>& other) {
            temp_weights = other;
            return *this;
        }
        
        void resize(size_t size) {
            temp_weights.resize(size);
        }
        
        void init() {
            for (int u = 0; u < graph.n; u++) {
                auto neighbors = graph.get_neighbors(u);
                for (int v : neighbors) {
                    temp_weights[u][v] = graph.get_edge_weight(u, v);
                }
            }
        }
        
        unordered_map<int, double>& operator[](int u) {
            return temp_weights[u];
        }
        
        const unordered_map<int, double>& operator[](int u) const {
            return temp_weights[u];
        }
        
        size_t size() const { return graph.n; }
        
        void sync() {
            for (int u = 0; u < graph.n; u++) {
                for (const auto& [v, w] : temp_weights[u]) {
                    graph.set_edge_weight(u, v, w);
                }
            }
        }

        friend void swap(EdgeWeights& first, EdgeWeights& second) noexcept {
            using std::swap;
            // Note: Don't swap graph references as they should stay with their original objects
            swap(first.temp_weights, second.temp_weights);
        }
    };
    
    // Make these non-const references to support mutation
    AdjacencyList adj_list;
    EdgeWeights edge_weight;
    
    Graph() : adj_list(*this), edge_weight(*this) {
        path_weight.set_graph(this);
        similarity_matrix.set_graph(this);
        jac_res.set_graph(this);
    }
    
    explicit Graph(const string &graph_path);
    
    void init(const string &graph_path);
    void init_nm();
    void init_d_neighbors();
    void init_similarity();
    
    double get_path_weight(int u, int v);
    void set_path_weight(int u, int v, double w);
    
    int get_similairty(int u, int v);
    void set_similarity(int u, int v, bool s);
    
    double get_jac_res(int u, int v);
    void set_jac_res(int u, int v, double s);
    
    unordered_map<int, double> compute_path_weight(int u, double dis_thre);
    
    void edge_del(int u, int v);
    void edge_ins(int u, int v, double w);
    void edge_update(int u, int v, double w);
    
    const vector<unordered_map<int, double>> &getPathWeight() const;
    
    void get_undirected_weighted_graph_snap_uniform();
    void get_undirected_weighted_graph_snap_exponent(double exp = -1);
    void get_undirected_similarity_weighted();
    void save_graph(const string &new_graph_name);
    void handle_LFR_graph(const string &graph_path, int nodes = -1);
    
    double jaccard_raw(int u, int v);
    double jaccard_raw_wscan(int u, int v);
    
    void reweighted(int type);
    
    ipair random_choose_edge();
    void convert_to_undirected_graph(string graph_path);

    // Enable swapping of Graph objects
    friend void swap(Graph& first, Graph& second) noexcept {
        using std::swap;
        
        // Swap basic members
        swap(first.n, second.n);
        swap(first.m, second.m);
        swap(first.data_folder, second.data_folder);
        swap(first.directed, second.directed);
        swap(first.weighted, second.weighted);
        
        // Swap CSR structure
        swap(first.csr_graph.offsets, second.csr_graph.offsets);
        swap(first.csr_graph.edges, second.csr_graph.edges);
        swap(first.csr_graph.weights, second.csr_graph.weights);
        
        // Swap memory management
        swap(first.current_memory_usage, second.current_memory_usage);
        swap(first.mmap_fd, second.mmap_fd);
        swap(first.mmap_data, second.mmap_data);
        swap(first.mmap_size, second.mmap_size);
        
        // Swap core data structures
        swap(first.clusterID, second.clusterID);
        swap(first.is_core, second.is_core);
        swap(first.weight_degree, second.weight_degree);
        swap(first.whole_weight, second.whole_weight);
        
        // Swap working sets
        swap(first.d_neighbors, second.d_neighbors);
        swap(first.similarity, second.similarity);
        
        // Swap lazy matrices
        swap(first.path_weight, second.path_weight);
        swap(first.similarity_matrix, second.similarity_matrix);
        swap(first.jac_res, second.jac_res);
        
        // Update graph pointers in lazy matrices
        first.path_weight.set_graph(&first);
        first.similarity_matrix.set_graph(&first);
        first.jac_res.set_graph(&first);
        second.path_weight.set_graph(&second);
        second.similarity_matrix.set_graph(&second);
        second.jac_res.set_graph(&second);
    }
};

extern Graph graph;
#endif //SCAN_GRAPH_H
