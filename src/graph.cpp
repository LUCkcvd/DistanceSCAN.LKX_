#include "../include/graph.h"

Graph graph;

void Graph::init_d_neighbors() {
    Timer timer(D_NEIGHBOR_TIME);
    path_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    d_neighbors = vector<vector<int >>(n, vector<int>());
    
    // Use average degree to estimate connectivity
    double avg_degree = (2.0 * m) / n;
    int reserve_size = min(static_cast<int>(avg_degree * 3), n/100); // Reserve triple avg degree, max 1% of n
    
    INFO("Initializing d_neighbors with reserve size: " + to_string(reserve_size));
    
    // Pre-allocate memory in chunks
    const int CHUNK_SIZE = 10000;
    for (int i = 0; i < n; i += CHUNK_SIZE) {
        int chunk_end = min(i + CHUNK_SIZE, n);
        for (int j = i; j < chunk_end; j++) {
            path_weight[j].reserve(reserve_size);
            d_neighbors[j].reserve(reserve_size);
        }
        
        // Monitor memory less frequently
        if ((i/CHUNK_SIZE) % 100 == 0) {
            long long curr_mem = get_proc_memory();
            if (curr_mem > 2000000) { // ~2GB
                INFO("Memory: " + format_memory_usage(curr_mem * 1024) + " - cleaning up");
                for (int k = 0; k < i; k++) {
                    d_neighbors[k].shrink_to_fit();
                    path_weight[k].rehash(0);
                }
                std::system("sync");
            }
        }
    }
    
    // Process nodes
    for (int n_id = 0; n_id < n; ++n_id) {
        d_neighbors[n_id].emplace_back(n_id);
        priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
        pq.push(make_pair(n_id, 0));
        set_path_weight(n_id, n_id, 0);
        
        for (int nei_id: d_neighbors[n_id]) {
            pq.push(make_pair(nei_id, get_path_weight(n_id, nei_id)));
        }
        
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
        
        // Clear queue memory
        priority_queue<idpair, vector<idpair>, cmp_idpair>().swap(pq);
    }
}

void Graph::init(const string &graph_path) {
    Timer timer(READ_GRAPH_TIME);
    INFO("Reading graph...");
    this->data_folder = graph_path;
    init_nm();
    
    // Monitor initial memory usage
    long long initial_mem = get_proc_memory();
    INFO("Initial memory: " + format_memory_usage(initial_mem * 1024));
    
    // Calculate average degree from total edges
    double avg_degree = (2.0 * m) / n;  // For undirected graph
    int reserve_size = min(static_cast<int>(avg_degree * 2), n/100);  // Reserve double avg degree, max 1% of n
    
    // Pre-allocate vectors with more conservative sizing
    const int CHUNK_SIZE = 10000;  // Smaller fixed chunk size
    adj_list.clear();
    edge_weight.clear();
    adj_list.resize(n);
    edge_weight.resize(n);
    
    // Allocate in chunks
    for (int i = 0; i < n; i += CHUNK_SIZE) {
        int chunk_end = min(i + CHUNK_SIZE, n);
        for (int j = i; j < chunk_end; j++) {
            adj_list[j].reserve(reserve_size);
            edge_weight[j].reserve(reserve_size);
        }
        
        // Monitor memory less frequently
        if ((i/CHUNK_SIZE) % 100 == 0) {
            long long curr_mem = get_proc_memory();
            if (curr_mem > 2000000) {  // ~2GB
                INFO("Memory: " + format_memory_usage(curr_mem * 1024) + " - cleaning up");
                for (int k = 0; k < i; k++) {
                    adj_list[k].shrink_to_fit();
                    edge_weight[k].rehash(0);
                }
                std::system("sync");
            }
        }
    }
    
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

    // Open file with error checking
    FILE *fin = fopen(graph_file.c_str(), "r");
    if (!fin) {
        throw std::runtime_error("Failed to open graph file: " + graph_file);
    }
    
    // Use buffered I/O for better performance
    char buffer[65536];
    if (setvbuf(fin, buffer, _IOFBF, sizeof(buffer)) != 0) {
        fclose(fin);
        throw std::runtime_error("Failed to set file buffer");
    }

    try {
        if (weighted) {
            int t1, t2;
            double w;
            int batch_size = 0;
            const int MAX_BATCH = 100000;  // Larger batch size for better performance
            vector<tuple<int,int,double>> edge_batch;
            edge_batch.reserve(MAX_BATCH);
            
            long long edges_read = 0;
            
            while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) == 3) {
                if (t1 == t2 || t1 >= n || t2 >= n) continue;
                
                edge_batch.emplace_back(t1, t2, w);
                batch_size++;
                edges_read++;
                
                if (batch_size >= MAX_BATCH) {
                    // Process batch
                    for (const auto& [src, dst, weight] : edge_batch) {
                        adj_list[src].push_back(dst);
                        adj_list[dst].push_back(src);
                        edge_weight[src][dst] = weight;
                        edge_weight[dst][src] = weight;
                    }
                    
                    // Clear batch and check memory
                    edge_batch.clear();
                    edge_batch.shrink_to_fit();
                    batch_size = 0;
                    
                    // Monitor memory less frequently
                    if (edges_read % 1000000 == 0) {
                        long long curr_mem = get_proc_memory();
                        if (curr_mem > 2000000) {  // ~2GB
                            INFO("Memory: " + format_memory_usage(curr_mem * 1024) + " - cleaning up");
                            for (auto& vec : adj_list) vec.shrink_to_fit();
                            for (auto& map : edge_weight) map.rehash(0);
                            std::system("sync");
                        }
                    }
                }
            }
            
            // Process remaining edges
            if (!edge_batch.empty()) {
                for (const auto& [src, dst, weight] : edge_batch) {
                    adj_list[src].push_back(dst);
                    adj_list[dst].push_back(src);
                    edge_weight[src][dst] = weight;
                    edge_weight[dst][src] = weight;
                }
                edge_batch.clear();
                edge_batch.shrink_to_fit();
            }
            
        } else {
            int t1, t2;
            int batch_size = 0;
            const int MAX_BATCH = 10000;
            vector<pair<int,int>> edge_batch;
            edge_batch.reserve(MAX_BATCH);
            
            while (fscanf(fin, "%d%d", &t1, &t2) == 2) {
                if (t1 == t2 || t1 >= n || t2 >= n) continue;
                
                edge_batch.emplace_back(t1, t2);
                batch_size++;
                
                if (batch_size >= MAX_BATCH) {
                    // Process batch
                    for (const auto& [src, dst] : edge_batch) {
                        adj_list[src].push_back(dst);
                        adj_list[dst].push_back(src);
                    }
                    edge_batch.clear();
                    batch_size = 0;
                }
            }
            
            // Process remaining edges
            for (const auto& [src, dst] : edge_batch) {
                adj_list[src].push_back(dst);
                adj_list[dst].push_back(src);
            }
        }
    } catch (...) {
        fclose(fin);
        throw;
    }
    
    fclose(fin);

    result.n = this->n;
    result.m = this->m;
    INFO("Graph loaded: " + to_string(n) + " nodes, " + to_string(m) + " edges");

    // Initialize other data structures with pre-allocation
    clusterID = vector<int>(n, -1);
    is_core = vector<int>(n, -1);
    similarity = vector<unordered_map<int, bool>>(n);
    
    // Use average degree for similarity map reservation
    int sim_reserve_size = min(static_cast<int>(avg_degree * 3), n/100);
    
    // Allocate similarity maps in chunks
    for (int i = 0; i < n; i += CHUNK_SIZE) {
        int chunk_end = min(i + CHUNK_SIZE, n);
        for (int j = i; j < chunk_end; j++) {
            similarity[j].reserve(sim_reserve_size);
        }
        
        // Monitor memory less frequently
        if ((i/CHUNK_SIZE) % 100 == 0) {
            long long curr_mem = get_proc_memory();
            if (curr_mem > 2000000) {
                INFO("Memory: " + format_memory_usage(curr_mem * 1024) + " - cleaning up");
                for (int k = 0; k < i; k++) {
                    similarity[k].rehash(0);
                }
                std::system("sync");
            }
        }
    }

    if (config.algo == W_SCAN) {
        weight_degree = vector<double>(n, 0);
        whole_weight = 0;
        if (config.similarityType == Config::cos) {
            #pragma omp parallel for reduction(+:whole_weight)
            for (int i = 0; i < n; ++i) {
                for (int nei: adj_list[i]) {
                    weight_degree[i] += edge_weight[i][nei] * edge_weight[i][nei];
                }
                weight_degree[i] = sqrt(weight_degree[i] + 1);
                whole_weight += weight_degree[i] + 1;
            }
        }
    }

    if (data_folder.find("graphs_coauthors") != data_folder.npos || 
        data_folder.find("LFR") != data_folder.npos) {
        reweighted(config.type);
    }

    if (config.operation == CLUSTER_VALIDATION) {
        jac_res = vector<unordered_map<int, double>>(n);
        
        // Allocate jac_res maps in chunks
        for (int i = 0; i < n; i += CHUNK_SIZE) {
            int chunk_end = min(i + CHUNK_SIZE, n);
            for (int j = i; j < chunk_end; j++) {
                jac_res[j].reserve(sim_reserve_size);
            }
            
            // Monitor memory less frequently
            if ((i/CHUNK_SIZE) % 100 == 0) {
                long long curr_mem = get_proc_memory();
                if (curr_mem > 2000000) {
                    INFO("Memory: " + format_memory_usage(curr_mem * 1024) + " - cleaning up");
                    for (int k = 0; k < i; k++) {
                        jac_res[k].rehash(0);
                    }
                    std::system("sync");
                }
            }
        }
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
        if (output.size() != 2) continue;
        if (output[0] == "n") {
            n = stoi(output[1]);
        } else if (output[0] == "m") {
            m = stoi(output[1]);
        } else if (output[0] == "weighted") {
            weighted = (output[1] != "0");
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
    union_table.reserve(adj_list[u].size() + adj_list[v].size());
    
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

void Graph::get_undirected_similarity_weighted() {
    string w_file_jac = data_folder + "jac_graph.txt";
    string w_file_cos = data_folder + "cos_graph.txt";
    ofstream out_jac(w_file_jac, ios::out);
    ofstream out_cos(w_file_cos, ios::out);
    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        stringstream jac_stream, cos_stream;
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
            
            jac_stream << i << "\t" << nei << "\t" << jac << "\n";
            cos_stream << i << "\t" << nei << "\t" << cos << "\n";
        }
        
        #pragma omp critical
        {
            if (out_jac.is_open()) out_jac << jac_stream.str();
            if (out_cos.is_open()) out_cos << cos_stream.str();
        }
    }
}

double Graph::jaccard_raw(int u, int v) {
    unordered_map<int, bool> union_table;
    union_table.reserve(adj_list[u].size() + adj_list[v].size());
    
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
    int max_weight = 0;
    long long sum_weight = 0;
    
    #pragma omp parallel for reduction(max:max_weight) reduction(+:sum_weight)
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            double w = edge_weight[i][nei];
            max_weight = max(max_weight, static_cast<int>(w));
            sum_weight += static_cast<long long>(w);
        }
    }
    
    sum_weight /= 2; // Adjust for undirected graph
    
    // Pre-compute common values
    double log_sum_weight = log(sum_weight);
    double log_max_weight = log(max_weight + 1);
    
    // Update edge weights in parallel
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int nei: adj_list[i]) {
            int x = edge_weight[i][nei];
            switch (type) {
                case 0: // 1/(log(x) + 1)
                    edge_weight[i][nei] = 1.0 / (log(x) + 1.0);
                    break;
                case 1: // log(sum_weight/x)/log(sum_weight)
                    edge_weight[i][nei] = log(sum_weight/x) / log_sum_weight;
                    break;
                case 2: // (log(sum_weight/(1+x))+1)/(log(sum_weight/2)+1)
                    edge_weight[i][nei] = (log(sum_weight/(1.0+x)) + 1.0) / 
                                        (log(sum_weight/2.0) + 1.0);
                    break;
                case 3: // log((max_weight+1)/x)/log(max_weight+1)
                    edge_weight[i][nei] = log((max_weight + 1.0)/x) / log_max_weight;
                    break;
                case 4: // log((sum_weight-x)/x)/log(sum_weight-1)
                    edge_weight[i][nei] = log((sum_weight - x)/x) / log(sum_weight - 1);
                    break;
            }
        }
    }
}

void Graph::edge_del(int u, int v) {
    if (u == v) return;
    
    m--;
    
    // Remove edges from adjacency lists using swap-and-pop
    auto remove_edge = [](vector<int>& adj, int target) {
        auto it = find(adj.begin(), adj.end(), target);
        if (it != adj.end()) {
            *it = adj.back();
            adj.pop_back();
        }
    };
    
    remove_edge(adj_list[u], v);
    remove_edge(adj_list[v], u);
    
    // Remove weights
    edge_weight[u].erase(v);
    edge_weight[v].erase(u);
    
    // Update similarity if needed
    if (!similarity.empty()) {
        if (u < v) {
            similarity[u].erase(v);
        } else {
            similarity[v].erase(u);
        }
    }
}

void Graph::edge_ins(int u, int v, double w) {
    if (u == v) return;
    
    m++;
    
    // Add edges to adjacency lists
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
    
    // Add weights
    edge_weight[u][v] = w;
    edge_weight[v][u] = w;
}

void Graph::edge_update(int u, int v, double w) {
    if (u == v) return;
    
    edge_weight[u][v] = w;
    edge_weight[v][u] = w;
}
