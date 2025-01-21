
#include "../include/graph.h"
#include <sys/resource.h>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <algorithm>

Graph graph;

bool Graph::map_file(const string& filename) {
    if (mmap_fd >= 0) {
        unmap_file();
    }
    
    mmap_fd = open(filename.c_str(), O_RDONLY);
    if (mmap_fd == -1) {
        return false;
    }
    
    struct stat sb;
    if (fstat(mmap_fd, &sb) == -1) {
        close(mmap_fd);
        mmap_fd = -1;
        return false;
    }
    
    mmap_size = sb.st_size;
    mmap_data = mmap(nullptr, mmap_size, PROT_READ, MAP_PRIVATE, mmap_fd, 0);
    
    if (mmap_data == MAP_FAILED) {
        close(mmap_fd);
        mmap_fd = -1;
        mmap_data = nullptr;
        return false;
    }
    
    return true;
}

void Graph::unmap_file() {
    if (mmap_data) {
        munmap(mmap_data, mmap_size);
        mmap_data = nullptr;
    }
    if (mmap_fd >= 0) {
        close(mmap_fd);
        mmap_fd = -1;
    }
    mmap_size = 0;
}

bool Graph::check_and_update_memory(size_t additional_bytes) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    size_t current_mem = static_cast<size_t>(usage.ru_maxrss) * 1024ULL;
    
    if (current_mem + additional_bytes > MAX_MEMORY_BYTES) {
        stringstream ss;
        ss << "Memory limit reached: " << (current_mem / (1024.0 * 1024 * 1024)) 
           << "GB used, cannot allocate additional " 
           << (additional_bytes / (1024.0 * 1024 * 1024)) << "GB";
        INFO(ss.str());
        
        release_memory();
        
        getrusage(RUSAGE_SELF, &usage);
        current_mem = static_cast<size_t>(usage.ru_maxrss) * 1024ULL;
        
        if (current_mem + additional_bytes > MAX_MEMORY_BYTES) {
            return false;
        }
    }
    
    current_memory_usage = current_mem;
    return true;
}

void Graph::track_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    current_memory_usage = static_cast<size_t>(usage.ru_maxrss) * 1024ULL;
    
    if (current_memory_usage > MAX_MEMORY_BYTES) {
        stringstream ss;
        ss << "Memory usage exceeds limit: " 
           << (current_memory_usage / (1024.0 * 1024 * 1024)) << "GB";
        INFO(ss.str());
        release_memory();
    }
}

void Graph::release_memory() {
    // Clear lazy-loaded matrices
    path_weight.clear();
    similarity_matrix.clear();
    jac_res.clear();
    
    // Reset memory pool
    mem_pool.reset();
    
    // Clear temporary working sets
    d_neighbors.clear();
    d_neighbors.shrink_to_fit();
    
    // Clear CSR buffers if possible
    if (m == 0) {
        csr_graph.clear();
    }
    
    // Force garbage collection
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void Graph::init_d_neighbors() {
    Timer timer(D_NEIGHBOR_TIME);
    
    // Initialize data structures
    path_weight.init(n);
    d_neighbors.clear();
    d_neighbors.resize(n);
    
    // Process each vertex
    for (int n_id = 0; n_id < n; ++n_id) {
        // Initialize with self-loop
        d_neighbors[n_id].push_back(n_id);
        
        // Setup priority queue for Dijkstra's algorithm
        priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
        pq.push(make_pair(n_id, 0));
        set_path_weight(n_id, n_id, 0);
        while (!pq.empty()) {
            auto cur_node = pq.top();
            pq.pop();
            if (cmp_double(cur_node.second, get_path_weight(n_id, cur_node.first)) == 1)continue;
            if (cur_node.first > n_id) {
                d_neighbors[n_id].emplace_back(cur_node.first);
                d_neighbors[cur_node.first].emplace_back(n_id);
            }

            for (int nei_id: adj_list[cur_node.first]) {
                if (nei_id < n_id) continue;
                double nei_dis = cur_node.second + edge_weight[cur_node.first][nei_id];
                if (cmp_double(nei_dis, config.distance) < 1) {
                    double old_dis = get_path_weight(n_id, nei_id);
                    if ((cmp_double(old_dis, -1) == 0 || cmp_double(old_dis, nei_dis) > 0)) {
                        set_path_weight(n_id, nei_id, nei_dis);
                        pq.push(make_pair(nei_id, nei_dis));
                    }
                }
            }
        }
    }
}


Graph::Graph(const string &graph_path) : adj_list(*this), edge_weight(*this) {
    path_weight.set_graph(this);
    similarity_matrix.set_graph(this);
    jac_res.set_graph(this);
    init(graph_path);
}

void Graph::init(const string &graph_path) {
    Timer timer(READ_GRAPH_TIME);
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    init_nm();
    
    // Calculate estimated memory for CSR format and auxiliary data
    size_t estimated_mem = 
        (n + 1) * sizeof(size_t) + // offsets array
        m * sizeof(int) +          // edges array
        (weighted ? m * sizeof(double) : 0) + // weights array
        n * sizeof(int) * 2 +     // clusterID and is_core
        n * sizeof(unordered_map<int, bool>); // similarity matrix
        
    if (!check_and_update_memory(estimated_mem)) {
        throw runtime_error("Insufficient memory to load graph of size n=" + 
                          to_string(n) + ", m=" + to_string(m));
    }
    
    // Initialize CSR structure with minimal memory footprint
    csr_graph.offsets.reserve(n + 1);
    csr_graph.offsets.push_back(0);
    csr_graph.edges.reserve(m);
    if (weighted) {
        csr_graph.weights.reserve(m);
    }
    
    // Initialize core data structures with exact sizes
    clusterID.resize(n, -1);
    is_core.resize(n, -1);
    
    // Initialize compatibility layers with minimal initial capacity
    adj_list.resize(n);
    edge_weight.resize(n);
    
    // Initialize similarity matrix only when needed
    if (!similarity.empty()) {
        similarity.clear();
    }
    similarity.resize(n);
    
    // Initialize lazy matrices only when needed
    if (config.operation == CLUSTER_VALIDATION) {
        jac_res.init(n);
    }
    
    // First pass: Count edges per node
    vector<size_t> edge_counts(n, 0);
    string graph_file = data_folder;
    if (directed) {
        graph_file += "graph.txt";
    } else if (config.operation == CONVERT_GRAPH) {
        graph_file += "undirect_graph.txt";
        weighted = false;
    } else {
        if (data_folder.find("graphs_coauthors") != data_folder.npos) {
            graph_file += "undirect_graph2.txt";
        } else if (data_folder.find("LFR") != data_folder.npos || config.operation == EXPONLFR) {
            graph_file += "undirect_graph.txt";
        } else {
            graph_file += "jac_graph.txt";
        }
    }

    FILE *fin = fopen(graph_file.c_str(), "r");
    if (!fin) {
        throw runtime_error("Failed to open graph file: " + graph_file);
    }

    // First pass: Count edges and build CSR structure
    if (weighted) {
        int t1, t2;
        double w;
        while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
            if (t1 == t2 || t1 < 0 || t2 < 0 || t1 >= n || t2 >= n) continue;
            edge_counts[t1]++;
            if (!directed) {
                edge_counts[t2]++; // Count both directions for undirected graphs
            }
        }
    } else {
        int t1, t2;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            if (t1 == t2 || t1 < 0 || t2 < 0 || t1 >= n || t2 >= n) continue;
            edge_counts[t1]++;
            if (!directed) {
                edge_counts[t2]++;
            }
        }
    }

    // Build CSR offsets
    for (int i = 0; i < n; i++) {
        csr_graph.offsets.push_back(csr_graph.offsets.back() + edge_counts[i]);
    }

    // Pre-allocate edges and weights arrays
    csr_graph.edges.resize(csr_graph.offsets.back());
    if (weighted) {
        csr_graph.weights.resize(csr_graph.offsets.back(), 1.0); // Default weight 1.0
    }

    // Reset edge counts for second pass
    vector<size_t> current_idx(n, 0);

    // Second pass: Fill CSR arrays
    fseek(fin, 0, SEEK_SET);
    if (weighted) {
        int t1, t2;
        double w;
        while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
            if (t1 == t2 || t1 < 0 || t2 < 0 || t1 >= n || t2 >= n) continue;
            
            // Add forward edge
            size_t idx1 = csr_graph.offsets[t1] + current_idx[t1]++;
            csr_graph.edges[idx1] = t2;
            csr_graph.weights[idx1] = w;
            
            // Add reverse edge for undirected graphs
            if (!directed) {
                size_t idx2 = csr_graph.offsets[t2] + current_idx[t2]++;
                csr_graph.edges[idx2] = t1;
                csr_graph.weights[idx2] = w;
            }
        }
    } else {
        int t1, t2;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            if (t1 == t2 || t1 < 0 || t2 < 0 || t1 >= n || t2 >= n) continue;
            
            // Add forward edge
            size_t idx1 = csr_graph.offsets[t1] + current_idx[t1]++;
            csr_graph.edges[idx1] = t2;
            
            // Add reverse edge for undirected graphs
            if (!directed) {
                size_t idx2 = csr_graph.offsets[t2] + current_idx[t2]++;
                csr_graph.edges[idx2] = t1;
            }
        }
    }
    fclose(fin);
    result.n = this->n;
    result.m = this->m;
    cout << "init graph graphn: " << this->n << " m: " << this->m << endl;

    clusterID = vector<int>(n, -1);
    is_core = vector<int>(n, -1);
    similarity = vector<unordered_map<int, bool >>(n, unordered_map<int, bool>{});
    if (config.algo == W_SCAN) {
        weight_degree = vector<double>(n, 0);
        whole_weight = 0;
        if (config.similarityType == Config::cos) {
            for (int i = 0; i < n; ++i) {
                for (int nei: adj_list[i]) {
                    weight_degree[i] += edge_weight[i][nei] * edge_weight[i][nei];
                }
                weight_degree[i] = sqrt(weight_degree[i] + 1);
                whole_weight += weight_degree[i] + 1;
            }
        }
    }
    if (data_folder.find("graphs_coauthors") != data_folder.npos || data_folder.find("LFR") != data_folder.npos) {
        reweighted(config.type);
    }
    if (config.operation == CLUSTER_VALIDATION) {
        jac_res = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    }
}


void Graph::init_nm() {
    string attribute_file = data_folder + "attribute.txt";
    assert_file_exist("attribute file", attribute_file);
    ifstream attr(attribute_file);
    string line;
    while (getline(attr, line)) {
        vector<string> output;
        split_string(line, output, "=");
        if (output.size() != 2)continue;
        if (output[0] == "n") {
            n = stoi(output[1]);
        } else if (output[0] == "m") {
            m = stoi(output[1]);
        } else if (output[0] == "weighted") {
            if (output[1] != "0") {
                weighted = true;
            } else {
                weighted = false;
            }
        }
    }
    attr.close();
}

double Graph::get_path_weight(int u, int v) {
    if (u > v) swap(u, v);
    if (path_weight[u].find(v) == path_weight[u].end()) {
        return -1;
    }
    return path_weight[u][v];
}

void Graph::set_path_weight(int u, int v, double w) {
    if (u > v) swap(u, v);
    path_weight[u][v] = w;
}


int Graph::get_similairty(int u, int v) {
    if (u > v) swap(u, v);
    if (similarity[u].find(v) == similarity[u].end()) {
        return -1;
    }
    return similarity[u][v];
}

void Graph::set_similarity(int u, int v, bool s) {
    if (u > v) swap(u, v);
    similarity[u][v] = s;
}

void Graph::init_similarity() {
    vector<unordered_map<int, bool >>(n, unordered_map<int, bool>{}).swap(similarity);
}


unordered_map<int, double> Graph::compute_path_weight(int u, double dis_thre) {
    unordered_map<int, double> nei_dis_table;
    nei_dis_table[u] = 0;
    priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
    pq.push(make_pair(u, 0));
    while (!pq.empty()) {
        auto cur_node = pq.top();
        pq.pop();
        if (cmp_double(cur_node.second, nei_dis_table[cur_node.first]) == 1)continue;
#ifdef _DEBUG_
        assert(cmp_double(cur_node.second, nei_dis_table[cur_node.first]) == 0);
#endif
        for (int nei_id: adj_list[cur_node.first]) {
            double nei_dis = cur_node.second + edge_weight[cur_node.first][nei_id];
            if (cmp_double(nei_dis, dis_thre) < 1) {
                if (nei_dis_table.find(nei_id) == nei_dis_table.end() ||
                    cmp_double(nei_dis_table[nei_id], nei_dis) > 0) {
                    nei_dis_table[nei_id] = nei_dis;
                    pq.push(make_pair(nei_id, nei_dis));
                }
            }
        }
    }
    return nei_dis_table;
}

void Graph::edge_del(int u, int v) {
    // Check memory usage before operation
    if (!check_and_update_memory(sizeof(size_t) * 2)) {
        throw runtime_error("Insufficient memory for edge deletion");
    }
    
    m--;
    
    // Find and mark edge in CSR structure
    {
        size_t start_u = csr_graph.offsets[u];
        size_t end_u = csr_graph.offsets[u + 1];
        bool found = false;
        for (size_t i = start_u; i < end_u && !found; i++) {
            if (csr_graph.edges[i] == v) {
                // Mark edge as deleted by setting weight to 0
                if (weighted) {
                    csr_graph.weights[i] = 0.0;
                }
                found = true;
            }
        }
        
        if (!directed && found) {
            size_t start_v = csr_graph.offsets[v];
            size_t end_v = csr_graph.offsets[v + 1];
            for (size_t i = start_v; i < end_v; i++) {
                if (csr_graph.edges[i] == u) {
                    if (weighted) {
                        csr_graph.weights[i] = 0.0;
                    }
                    break;
                }
            }
        }
    }
    
    // Update compatibility layers efficiently
    {
        // Remove edge from adjacency lists
        auto& adj_u = adj_list[u];
        auto it_u = find(adj_u.begin(), adj_u.end(), v);
        if (it_u != adj_u.end()) {
            *it_u = adj_u.back();
            adj_u.pop_back();
        }
        
        if (!directed) {
            auto& adj_v = adj_list[v];
            auto it_v = find(adj_v.begin(), adj_v.end(), u);
            if (it_v != adj_v.end()) {
                *it_v = adj_v.back();
                adj_v.pop_back();
            }
        }
        
        // Remove edge weights
        edge_weight[u].erase(v);
        if (!directed) {
            edge_weight[v].erase(u);
        }
    }
    
    // Update similarity if needed
    if (!similarity.empty()) {
        if (u < v) {
            similarity[u].erase(v);
        } else {
            similarity[v].erase(u);
        }
    }
    
    // Track memory usage after operation
    track_memory_usage();
}

void Graph::edge_ins(int u, int v, double w) {
    // Calculate memory needed for new edge
    size_t additional_mem = sizeof(int) * 2;  // For CSR edges
    if (weighted) {
        additional_mem += sizeof(double) * 2;  // For weights
    }
    if (!directed) {
        additional_mem *= 2;  // Double for undirected edges
    }
    
    // Check memory availability
    if (!check_and_update_memory(additional_mem)) {
        throw runtime_error("Insufficient memory for edge insertion");
    }
    
    m++;
    
    // Update CSR structure efficiently
    {
        // Pre-allocate space to avoid multiple reallocations
        size_t new_size = csr_graph.edges.size() + (directed ? 1 : 2);
        csr_graph.edges.reserve(new_size);
        if (weighted) {
            csr_graph.weights.reserve(new_size);
        }
        
        // Insert forward edge
        size_t idx_u = csr_graph.offsets[u + 1] - 1;
        csr_graph.edges.insert(csr_graph.edges.begin() + idx_u, v);
        if (weighted) {
            csr_graph.weights.insert(csr_graph.weights.begin() + idx_u, w);
        }
        
        // Insert reverse edge for undirected graphs
        if (!directed) {
            size_t idx_v = csr_graph.offsets[v + 1] - 1;
            csr_graph.edges.insert(csr_graph.edges.begin() + idx_v, u);
            if (weighted) {
                csr_graph.weights.insert(csr_graph.weights.begin() + idx_v, w);
            }
        }
        
        // Update offsets efficiently
        for (size_t i = u + 1; i < csr_graph.offsets.size(); i++) {
            csr_graph.offsets[i]++;
        }
        if (!directed) {
            for (size_t i = v + 1; i < csr_graph.offsets.size(); i++) {
                csr_graph.offsets[i]++;
            }
        }
    }
    
    // Update compatibility layers efficiently
    {
        // Pre-reserve space in adjacency lists
        if (adj_list[u].capacity() == adj_list[u].size()) {
            adj_list[u].reserve(adj_list[u].size() * 2);
        }
        adj_list[u].push_back(v);
        edge_weight[u][v] = w;
        
        if (!directed) {
            if (adj_list[v].capacity() == adj_list[v].size()) {
                adj_list[v].reserve(adj_list[v].size() * 2);
            }
            adj_list[v].push_back(u);
            edge_weight[v][u] = w;
        }
    }
    
    // Track memory usage after operation
    track_memory_usage();
}

void Graph::edge_update(int u, int v, double w) {
    // Check memory usage before operation
    if (!check_and_update_memory(sizeof(double) * (directed ? 1 : 2))) {
        throw runtime_error("Insufficient memory for edge update");
    }
    
    // Update CSR structure efficiently
    {
        size_t start_u = csr_graph.offsets[u];
        size_t end_u = csr_graph.offsets[u + 1];
        bool found = false;
        
        // Use binary search if edge list is large enough
        if (end_u - start_u > 32) {
            auto it = lower_bound(csr_graph.edges.begin() + start_u, 
                                csr_graph.edges.begin() + end_u, v);
            if (it != csr_graph.edges.begin() + end_u && *it == v) {
                csr_graph.weights[it - csr_graph.edges.begin()] = w;
                found = true;
            }
        } else {
            // Linear search for small lists
            for (size_t i = start_u; i < end_u && !found; i++) {
                if (csr_graph.edges[i] == v) {
                    csr_graph.weights[i] = w;
                    found = true;
                }
            }
        }
        
        if (!directed && found) {
            size_t start_v = csr_graph.offsets[v];
            size_t end_v = csr_graph.offsets[v + 1];
            
            // Use same search strategy for reverse edge
            if (end_v - start_v > 32) {
                auto it = lower_bound(csr_graph.edges.begin() + start_v,
                                    csr_graph.edges.begin() + end_v, u);
                if (it != csr_graph.edges.begin() + end_v && *it == u) {
                    csr_graph.weights[it - csr_graph.edges.begin()] = w;
                }
            } else {
                for (size_t i = start_v; i < end_v; i++) {
                    if (csr_graph.edges[i] == u) {
                        csr_graph.weights[i] = w;
                        break;
                    }
                }
            }
        }
    }
    
    // Update compatibility layer efficiently
    {
        edge_weight[u][v] = w;
        if (!directed) {
            edge_weight[v][u] = w;
        }
    }
    
    // Track memory usage after operation
    track_memory_usage();
}

const vector<unordered_map<int, double>> &Graph::getPathWeight() const {
    return path_weight;
}

void Graph::save_graph(const string &new_graph_name) {
    string w_file = data_folder + new_graph_name;
    ofstream out(w_file, ios::out);
    if (out.is_open()) {
        for (int i = 0; i < n; ++i) {
            for (int nei: adj_list[i]) {
                out << i << "\t" << nei << "\t";
                if (edge_weight[i].find(nei) != edge_weight[i].end()) {
                    out << edge_weight[i][nei];
                }
                out << "\n";
            }
        }
    }
}

void Graph::get_undirected_similarity_weighted() {
    string w_file_jac = data_folder + "jac_graph.txt";
    string w_file_cos = data_folder + "cos_graph.txt";
    ofstream out_jac(w_file_jac, ios::out);
    ofstream out_cos(w_file_cos, ios::out);
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            int inter = 2, i_index = 0, nei_index = 0;
            while (i_index < adj_list[i].size() && nei_index < adj_list[nei].size()) {
                if (adj_list[i][i_index] == adj_list[nei][nei_index]) {
                    inter++;
                    i_index++;
                    nei_index++;
                } else if (adj_list[i][i_index] < adj_list[nei][nei_index]) {
                    i_index++;
                } else {
                    nei_index++;
                }
            }
            double jac = double(inter) / double(adj_list[i].size() + adj_list[nei].size() + 2 - inter);
            double cos = double(inter) / sqrt(double((adj_list[i].size() + 1) * (adj_list[nei].size() + 1)));
            jac = 1 - (0.9 * jac);
            cos = 1 - (0.9 * cos);
            if (out_jac.is_open()) {
                out_jac << i << "\t" << nei << "\t" << jac << "\n";
            }
            if (out_cos.is_open()) {
                out_cos << i << "\t" << nei << "\t" << cos << "\n";
            }
        }
    }
}

void Graph::get_undirected_weighted_graph_snap_uniform() {
    static default_random_engine generator(time(NULL));
    static uniform_real_distribution<double> dis(0, 0.9);
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            if (edge_weight[i].find(nei) == edge_weight[i].end()) {
                edge_weight[i][nei] = 1 - dis(generator);
                edge_weight[nei][i] = edge_weight[i][nei];
            }
        }
    }
    save_graph("uniform_weighted_graph.txt");
}

void Graph::get_undirected_weighted_graph_snap_exponent(double log_base) {
    std::default_random_engine generator(time(NULL));
    std::uniform_real_distribution<double> dis(0.1, 1);
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            double w = -1 * log(1 - dis(generator));
            if (log_base > 1) {
                w /= log(log_base);
            }
            if (edge_weight[i].find(nei) == edge_weight[i].end()) {
                edge_weight[i][nei] = w;
                edge_weight[nei][i] = edge_weight[i][nei];
            }
        }
    }
    if (log_base > 1) {
        save_graph("exponent_" + to_string(log_base) + "_weighted_graph.txt");
    } else {
        save_graph("exponent_e_weighted_graph.txt");
    }
}

void Graph::handle_LFR_graph(const string &graph_path, int nodes) {
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    INFO(nodes);
    if (nodes == -1) {
        string attribute_file = data_folder + "flags.dat";
        INFO(attribute_file);
        ifstream attr(attribute_file);
        string line;
        while (getline(attr, line)) {
            vector<string> output;
            split_string(line, output, " ");
            INFO(line);
            INFO(output);
            if (output.size() != 2)continue;
            if (output[0] == "-N") {
                n = stoi(output[1]) + 1;
                break;
            }
        }
        attr.close();
    } else {
        n = nodes;
    }
    INFO(n);
    adj_list = vector<vector<int >>(n, vector<int>());
    edge_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>());
    string graph_file = data_folder + "network.dat";

    m = 0;
    FILE *fin = fopen(graph_file.c_str(), "r");
    int t1, t2;
    double w;
    while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
        if (t1 == t2 || cmp_double(w, 0) == 0)continue;
        adj_list[t1].push_back(t2);
        if (edge_weight[t1].find(t2) == edge_weight[t1].end()) {
            edge_weight[t1][t2] = w;
            edge_weight[t2][t1] = w;
        } else {
            edge_weight[t1][t2] += w;
            edge_weight[t2][t1] += w;
            //edge_weight[t1][t2] /= 2;
            //edge_weight[t2][t1] /= 2;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            if (cmp_double(edge_weight[i][nei], 0) < 1) {
                edge_del(i, nei);
            }
        }
    }
    m = 0;
    for (int i = 0; i < n; ++i) {
        m += adj_list[i].size();
    }
    fclose(fin);
    cout << "init graph graph n: " << this->n << " m: " << this->m << endl;
    save_graph("raw_undirect_graph.txt");
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            edge_weight[i][nei] = 1 / edge_weight[i][nei];
        }
    }
    m = 0;
    for (int i = 0; i < n; ++i) {
        m += adj_list[i].size();
    }
    save_graph("undirect_graph.txt");
    string w_file = data_folder + "attribute.txt";
    ofstream out(w_file, ios::out);
    if (out.is_open()) {
        out << "n=" << n << "n" << endl;
        out << "m=" << m << "n" << endl;
    }
    out.close();
}

double Graph::get_jac_res(int u, int v) {
    if (u > v) swap(u, v);
    if (jac_res[u].find(v) == jac_res[u].end()) {
        return -1;
    }
    return jac_res[u][v];
}

void Graph::set_jac_res(int u, int v, double s) {
    if (u > v) swap(u, v);
    jac_res[u][v] = s;
}

double Graph::jaccard_raw_wscan(int u, int v) {
    double inter = 0;
    unordered_map<int, bool> union_table;
    for (int nei: adj_list[u]) {
        union_table[nei] = true;
    }
    for (int nei: adj_list[v]) {
        if (union_table.find(nei) != union_table.end()) {
            if (u == nei || v == nei) {
                inter += exp(1 / edge_weight[u][v] - 1);
            } else {
                double real_w1 = exp(1 / edge_weight[u][nei] - 1);
                double real_w2 = exp(1 / edge_weight[v][nei] - 1);
                inter += real_w1 * real_w2;
            }
        }
    }
    assert(config.similarityType == Config::cos);
    return inter / weight_degree[u] / weight_degree[v];
}

double Graph::jaccard_raw(int u, int v) {
    unordered_map<int, bool> union_table;
    for (int nei: adj_list[u]) {
        union_table[nei] = true;
    }
    for (int nei: adj_list[v]) {
        union_table[nei] = true;
    }
    int inter_size = adj_list[u].size() + adj_list[v].size() - union_table.size();
    double similarity_result;
    if(config.similarityType == config.jac){
        similarity_result = double(inter_size) / double(union_table.size());
    } else if (config.similarityType == config.cos) {
        similarity_result = inter_size / sqrt(adj_list[u].size()) / sqrt(adj_list[v].size());
    } else if(config.similarityType == config.set_containment1){
        similarity_result = double(inter_size) / double (adj_list[u].size());
    } else if(config.similarityType == config.set_containment2){
        similarity_result = double(inter_size) / double (adj_list[v].size());
    }
    if (config.operation == CLUSTER_VALIDATION) {
        set_jac_res(u, v, similarity_result);
    }
    return similarity_result;
}

void Graph::reweighted(int type) {
    int max_weight = 0, sum_weight = 0;
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            if (max_weight < edge_weight[i][nei]) {
                max_weight = edge_weight[i][nei];
            }
            sum_weight += edge_weight[i][nei];
        }
    }
    sum_weight /= 2;
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            int x = edge_weight[i][nei];
            switch (type) {
                case 0:
                    // 1/(log(x) + 1)
                    edge_weight[i][nei] = 1 / (log(x) + 1);
                    //INFO(edge_weight[i][nei]);
                    break;
                case 1:
                    // log(sum_weight/x)/log(sum_weight) 归一化后的 inverse document frequency
                    edge_weight[i][nei] = log(sum_weight / x) / log(sum_weight);
                    break;
                case 2:
                    // (log(sum_weight/(1+x))+1)/(log(sum_weight/2)+1) 归一化的 inverse document frequency smooth
                    edge_weight[i][nei] = (log(sum_weight / (1 + x)) + 1) / (log(sum_weight / 2) + 1);
                    break;
                case 3:
                    // log((max_weight+1)/x)/log(max_weight+1) 归一化的 inverse document frequency max
                    edge_weight[i][nei] = log((max_weight + 1) / x) / log(max_weight + 1);
                    break;
                case 4:
                    // log((sum_weight-x)/x)/log(max_weight-1) 归一化的 probabilistic inverse document frequency
                    edge_weight[i][nei] = log((sum_weight - x) / x) / log(sum_weight - 1);
                    break;

            }
        }
    }
}

ipair Graph::random_choose_edge() {
    static std::default_random_engine generator(time(NULL));
    std::uniform_int_distribution<int> dis(0, n - 1);
    int source = -1, target = -1;
    do {
        source = dis(generator);
    } while (adj_list[source].empty());
    int max_idx = adj_list[source].size() - 1;
    std::uniform_int_distribution<int> dis2(0, max_idx);
    target = dis2(generator);
    return make_pair(source, adj_list[source][target]);
}

void Graph::convert_to_undirected_graph(string graph_path) {
    data_folder = config.graph_location + graph_path;
    init_nm();
    vector<set<int>> undirect_adj_list(n);
    string graph_file = data_folder + FILESEP + "graph.txt";
    assert_file_exist("graph file", graph_file);
    FILE *fin = fopen(graph_file.c_str(), "r");
    char line[1000];
    long long t1, t2;
    int min_id = INT_MAX, max_id = -1, final_n = 0, final_m = 0, tmp = 0;
    unordered_map<long, int> str_index;
    fin = fopen(graph_file.c_str(), "r");
    while (fgets(line, 1000, fin) != NULL) {
        if (line[0] > '9' || line[0] < '0') { continue; }
        sscanf(line, "%lld\t%lld\n", &t1, &t2);
        if (t1 == t2)continue;
        min_id = min(min_id, (int) t1);
        min_id = min(min_id, (int) t2);
        max_id = max(max_id, (int) t1);
        max_id = max(max_id, (int) t2);
        undirect_adj_list[t1].insert(t2);
        undirect_adj_list[t2].insert(t1);
    }
    string undirect_graph_file = data_folder + FILESEP + "undirect_graph.txt";
    FILE *fout = fopen(undirect_graph_file.c_str(), "w");
    for (int j = 0; j < undirect_adj_list.size(); ++j) {
        if (!undirect_adj_list[j].empty()) {
            final_n++;
            for (int des: undirect_adj_list[j]) {
                final_m++;
                fprintf(fout, "%d\t%d\n", j, des);
            }
        }
    }
    fclose(fout);
    string undirect_graph_attr_file = data_folder + FILESEP + "undirect_graph_attribute.txt";
    fout = fopen(undirect_graph_attr_file.c_str(), "w");
    fprintf(fout, "n=%d\nm=%lld\n", n, m);
    fprintf(fout, "final_n=\t%d\nfinal_m=\t%d\n", final_n, final_m);
    fprintf(fout, "min_node=\t%d\nmax_node=\t%d\n", min_id, max_id);
    for (auto item: str_index) {
        fprintf(fout, "%ld\t%d\n", item.first, item.second);
    }
    fclose(fin);
    fclose(fout);
}
