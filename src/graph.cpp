#include "../include/graph.h"

Graph graph;

// Previous optimized init_d_neighbors() implementation remains unchanged...
void Graph::init_d_neighbors() {
    Timer timer(D_NEIGHBOR_TIME);
    path_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    d_neighbors = vector<vector<int >>(n, vector<int>());
    
    // Pre-allocate memory for path weights
    for (auto& map : path_weight) {
        map.reserve(n/10); // Estimate average connectivity
    }
    
    for (int n_id = 0; n_id < n; ++n_id) {
        d_neighbors[n_id].reserve(n/10); // Pre-allocate estimated size
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
    }
}

// Previous optimized init() implementation remains unchanged...
void Graph::init(const string &graph_path) {
    Timer timer(READ_GRAPH_TIME);
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    init_nm();
    
    // Pre-allocate vectors with estimated sizes
    adj_list = vector<vector<int>>(n);
    edge_weight = vector<unordered_map<int, double>>(n);
    for (auto& vec : adj_list) {
        vec.reserve(n/10); // Estimate average degree
    }
    for (auto& map : edge_weight) {
        map.reserve(n/10);
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
            while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) == 3) {
                if (t1 == t2 || t1 >= n || t2 >= n) continue;
                
                // Add edges in both directions for undirected graph
                adj_list[t1].push_back(t2);
                adj_list[t2].push_back(t1);
                
                edge_weight[t1][t2] = w;
                edge_weight[t2][t1] = w;
            }
        } else {
            int t1, t2;
            while (fscanf(fin, "%d%d", &t1, &t2) == 2) {
                if (t1 == t2 || t1 >= n || t2 >= n) continue;
                adj_list[t1].push_back(t2);
                adj_list[t2].push_back(t1); // Add reverse edge for undirected graph
            }
        }
    } catch (...) {
        fclose(fin);
        throw;
    }
    
    fclose(fin);

    result.n = this->n;
    result.m = this->m;
    cout << "init graph graphn: " << this->n << " m: " << this->m << endl;

    // Initialize other data structures
    clusterID = vector<int>(n, -1);
    is_core = vector<int>(n, -1);
    similarity = vector<unordered_map<int, bool>>(n);
    for (auto& map : similarity) {
        map.reserve(n/10);
    }

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

    if (data_folder.find("graphs_coauthors") != data_folder.npos || 
        data_folder.find("LFR") != data_folder.npos) {
        reweighted(config.type);
    }

    if (config.operation == CLUSTER_VALIDATION) {
        jac_res = vector<unordered_map<int, double>>(n);
        for (auto& map : jac_res) {
            map.reserve(n/10);
        }
    }
}

// Optimized jaccard_raw implementation using two-pointer technique
double Graph::jaccard_raw(int u, int v) {
    vector<int>& u_adj = adj_list[u];
    vector<int>& v_adj = adj_list[v];
    
    // Sort adjacency lists if needed
    if (!is_sorted(u_adj.begin(), u_adj.end())) {
        sort(u_adj.begin(), u_adj.end());
    }
    if (!is_sorted(v_adj.begin(), v_adj.end())) {
        sort(v_adj.begin(), v_adj.end());
    }
    
    // Calculate intersection size using two-pointer technique
    int i = 0, j = 0;
    int intersection_size = 0;
    while (i < u_adj.size() && j < v_adj.size()) {
        if (u_adj[i] == v_adj[j]) {
            intersection_size++;
            i++;
            j++;
        } else if (u_adj[i] < v_adj[j]) {
            i++;
        } else {
            j++;
        }
    }
    
    double similarity_result;
    if (config.similarityType == config.jac) {
        int union_size = u_adj.size() + v_adj.size() - intersection_size;
        similarity_result = static_cast<double>(intersection_size) / union_size;
    } else if (config.similarityType == config.cos) {
        similarity_result = intersection_size / sqrt(u_adj.size() * v_adj.size());
    } else if (config.similarityType == config.set_containment1) {
        similarity_result = static_cast<double>(intersection_size) / u_adj.size();
    } else {
        similarity_result = static_cast<double>(intersection_size) / v_adj.size();
    }
    
    if (config.operation == CLUSTER_VALIDATION) {
        set_jac_res(u, v, similarity_result);
    }
    
    return similarity_result;
}

// Optimized reweighted implementation with parallel processing
void Graph::reweighted(int type) {
    // Pre-compute max and sum weights in one pass
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
            double x = edge_weight[i][nei];
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

// Optimized edge operations
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

// Rest of the implementation remains unchanged...
