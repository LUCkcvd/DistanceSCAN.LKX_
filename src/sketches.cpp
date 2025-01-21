
#include "sketches.h"
#include <omp.h>
#include <algorithm>
#include <queue>
#include <mutex>
#include <thread>
#include <chrono>


vector<unsigned long long> SKETCHES::pick_random_coffs(int k, const unsigned long long &P) {
    vector<unsigned long long> coffs;
    while (k--) {
        coffs.emplace_back(rand_ulong(0, P - 1));
    }
    return coffs;
}

int
SKETCHES::double_hashing(const vector<unsigned long long int> &coff_a, const vector<unsigned long long int> &coff_b,
                         int key, int i) {
    // 2-universal hash h(key) = (a * key + b) % P % n.
    unsigned long long int h_0 = (coff_a[0] * key + coff_b[0]) % P % (2 * graph.n);
    unsigned long long int h_1 = (coff_a[1] * key + coff_b[1]) % P % (2 * graph.n);
    //double hashing  h(key, i) = (h1(key) + i * h2(key)) mod n
    return int((h_0 + i * h_1) % (2 * graph.n));
}

void SKETCHES::init_hash_table() {
    vector<unsigned long long> coff_a, coff_b;
    coff_a = pick_random_coffs(2, P);
    coff_b = pick_random_coffs(2, P);

#ifdef _DEBUG_
    coff_a = {409230109, 2524636865};
    coff_b = {92027918, 2638940465};
#endif
    key2value = vector<int>(graph.n, -1);
    value2key = vector<int>(2 * graph.n, -1);

    // Use OpenMP for parallel hash table initialization
    #pragma omp parallel
    {
        vector<int> local_key2value(graph.n, -1);
        vector<int> local_value2key(2 * graph.n, -1);

        #pragma omp for schedule(dynamic, BATCH_SIZE)
        for (int key = 0; key < key2value.size(); ++key) {
            for (int j = 0; j < value2key.size(); ++j) {
                int value = double_hashing(coff_a, coff_b, key, j);
                if (local_value2key[value] == -1) {
                    local_value2key[value] = key;
                    local_key2value[key] = value;
                    break;
                }
            }
        }

        // Merge local results into global arrays
        #pragma omp critical
        {
            for (int i = 0; i < graph.n; ++i) {
                if (local_key2value[i] != -1) {
                    key2value[i] = local_key2value[i];
                    value2key[local_key2value[i]] = i;
                }
            }
        }
    }
}

// Thread-local storage for intermediate results
struct ThreadLocalData {
    vector<int> current_histogram;  // Histogram for current node
    vector<idpair> local_dis_source_vec;
    priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
    Treap treap;
    ADS ads;
    
    ThreadLocalData(int n, int num_bins) : 
        current_histogram(num_bins, 0) {
        try {
            // Pre-allocate with proper size and default values
            local_dis_source_vec.resize(n, make_pair(-1, std::numeric_limits<double>::max()));
        } catch (const std::bad_alloc& e) {
            throw std::runtime_error("Failed to allocate memory for thread local data: " + std::string(e.what()));
        }
    }
    
    void reset() {
        // Reset histogram
        std::fill(current_histogram.begin(), current_histogram.end(), 0);
        
        // Reset distance vector for next node
        std::fill(local_dis_source_vec.begin(), local_dis_source_vec.end(), 
                 make_pair(-1, std::numeric_limits<double>::max()));
        
        // Clear priority queue
        while (!pq.empty()) {
            pq.pop();
        }
        
        // Reset data structures
        treap = Treap();
        ads = ADS();
    }
    
    // Prevent copying to avoid memory issues
    ThreadLocalData(const ThreadLocalData&) = delete;
    ThreadLocalData& operator=(const ThreadLocalData&) = delete;
    
    // Allow moving
    ThreadLocalData(ThreadLocalData&&) = default;
    ThreadLocalData& operator=(ThreadLocalData&&) = default;
};

void SKETCHES::construct_sketches() {
    Timer timer(CONSTRUCT_SKETCH_TIME);
    INFO("constructing sketches...");
    
    const size_t total_nodes = static_cast<size_t>(graph.n);
    std::atomic<size_t> processed_nodes{0};
    auto start_time = std::chrono::high_resolution_clock::now();
    
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string idx = config.graph_location + prefix + 
                 (config.algo == BOTK_SCAN ? "_multi_ads_bot.idx" : "_my_ads.idx");

    // Initialize file output and synchronization
    std::ofstream ofs(idx, std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    std::mutex oa_mutex;
    
    // Initialize shared data structures
    int num_bins = get_bin_id(config.max_distance, false) + 1;
    histogram = vector<vector<int>>(graph.n, vector<int>(num_bins));
    
    // Configure OpenMP and pre-allocate thread local storage
    const int num_threads = omp_get_max_threads();
    vector<unique_ptr<ThreadLocalData>> thread_data;
    thread_data.reserve(num_threads);
    
    try {
        for (int i = 0; i < num_threads; ++i) {
            thread_data.emplace_back(make_unique<ThreadLocalData>(graph.n, num_bins));
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to allocate thread local storage: " + std::string(e.what()));
    }
    
    // Create synchronization primitives for progress monitoring
    std::mutex progress_mutex;
    std::condition_variable progress_cv;
    std::atomic<bool> processing_complete{false};
    std::atomic<bool> progress_reported_100{false};
    
    // Create progress monitoring thread
    std::thread progress_thread([&]() {
        std::unique_lock<std::mutex> lock(progress_mutex);
        
        while (!processing_complete) {
            // Wait for notification or timeout
            if (progress_cv.wait_for(lock, std::chrono::seconds(1)) == std::cv_status::timeout) {
                auto current_time = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
                
                if (elapsed > 0) {
                    double progress = static_cast<double>(processed_nodes) / total_nodes * 100;
                    
                    // Only report if not at 100% or haven't reported 100% yet
                    if (progress < 100.0 || !progress_reported_100) {
                        double nodes_per_second = static_cast<double>(processed_nodes) / elapsed;
                        int eta_seconds = static_cast<int>((total_nodes - processed_nodes) / nodes_per_second);
                        
                        std::ostringstream oss;
                        oss << "Progress: " << std::fixed << std::setprecision(2) << progress << "% - ";
                        
                        // Only show ETA if not at 100%
                        if (progress < 100.0) {
                            oss << "ETA: " << (eta_seconds / 3600) << "h "
                                << ((eta_seconds % 3600) / 60) << "m "
                                << (eta_seconds % 60) << "s";
                        } else {
                            oss << "Complete!";
                            progress_reported_100 = true;
                        }
                        
                        INFO(oss.str());
                    }
                }
            }
        }
    });
    
    // Process nodes in parallel with dynamic scheduling for better load balancing
    #pragma omp parallel
    {
        const int thread_id = omp_get_thread_num();
        auto& local_data = *thread_data[thread_id];
        
        #pragma omp for schedule(dynamic, 64)
        for (size_t source_nid = 0; source_nid < static_cast<size_t>(graph.n); ++source_nid) {
            // Update progress counter
            if (static_cast<int>((source_nid) % 100) == 99) {
                // Update in batches of 100
                processed_nodes.fetch_add(100, std::memory_order_relaxed);
            } else if (source_nid == total_nodes - 1) {
                // Handle final batch
                size_t remaining = total_nodes - processed_nodes.load(std::memory_order_relaxed);
                if (remaining > 0) {
                    processed_nodes.fetch_add(remaining, std::memory_order_relaxed);
                }
            }
            
            // Initialize data structures for this node
            if (config.algo == BOTK_SCAN) {
                local_data.treap = Treap();
            } else if (config.algo == MY_ADS) {
                local_data.ads = ADS();
            }
    
            // Process source node
            idpair source(source_nid, 0);
            local_data.local_dis_source_vec[source_nid] = source;
            
            // Clear and initialize priority queue
            while (!local_data.pq.empty()) {
                local_data.pq.pop();
            }
            local_data.pq.push(source);
            
            // Process node's neighbors
            while (!local_data.pq.empty()) {
                auto cur_node = local_data.pq.top();
                local_data.pq.pop();
                
                if (cmp_double(cur_node.second, local_data.local_dis_source_vec[cur_node.first].second) == 1) 
                    continue;
                
                // Update histogram and data structures
                int bin_id = get_bin_id(cur_node.second, false);
                if (bin_id >= 0 && bin_id < num_bins) {
                    local_data.current_histogram[bin_id]++;
                }
                
                if (config.algo == BOTK_SCAN) {
                    local_data.treap.insert_botk(key2value[cur_node.first], cur_node.second);
                } else if (config.algo == MY_ADS) {
                    local_data.ads.insert_botk(key2value[cur_node.first], cur_node.second);
                }
                
                // Process neighbors
                const auto& neighbors = graph.adj_list[cur_node.first];
                
                for (int nei_id : neighbors) {
                    // Get edge weight using the safe EdgeWeights accessor
                    const auto& edge_weights = graph.edge_weight[cur_node.first];
                    auto weight_it = edge_weights.find(nei_id);
                    if (weight_it == edge_weights.end()) {
                        continue;  // Skip if no edge weight exists
                    }
                    
                    // Calculate path weight from source to neighbor through current node
                    double nei_dis = cur_node.second + weight_it->second;
                    
                    if (cmp_double(nei_dis, config.max_distance) < 1) {
                        if (local_data.local_dis_source_vec[nei_id].first != source_nid ||
                            cmp_double(local_data.local_dis_source_vec[nei_id].second, nei_dis) > 0) {
                            local_data.local_dis_source_vec[nei_id] = {source_nid, nei_dis};
                            local_data.pq.push({nei_id, nei_dis});
                        }
                    }
                }
            }
            
            // Prepare histogram data locally first
            vector<int> cumulative_histogram = local_data.current_histogram;
            for (int j = 1; j < cumulative_histogram.size(); ++j) {
                cumulative_histogram[j] += cumulative_histogram[j - 1];
            }
            
            // Batch update global data with minimal critical section
            #pragma omp critical
            {
                // Update global histogram
                histogram[source_nid] = std::move(cumulative_histogram);
                
                // Serialize data structures
                Timer timer(SAVE_BOTK_TIME);
                if (config.algo == BOTK_SCAN) {
                    oa << local_data.treap;
                } else if (config.algo == MY_ADS) {
                    oa << local_data.ads;
                }
            }
            
            // Reset data for next node
            local_data.reset();
        }
    }

    
    // Process nodes in parallel...
    // [existing parallel processing code]
    
    // Signal completion and wait for progress thread
    {
        std::lock_guard<std::mutex> lock(progress_mutex);
        processing_complete = true;
    }
    progress_cv.notify_one();
    progress_thread.join();
}

int SKETCHES::update_histogram(int source_nid, double path_weight) {
    int bin_id = get_bin_id(path_weight, false);
    histogram[source_nid][bin_id]++;
    return bin_id;
}

int SKETCHES::get_bin_id(double dis, bool is_lb) {
    int bin_id = floor(dis / config.bin);
    if (is_lb) {
        bin_id--;
    }
    return bin_id;
}

int SKETCHES::query_histogram(int source_nid, double dis, bool is_lb) {
    int bin_id = get_bin_id(dis, is_lb);
    if (bin_id < 0) {
        return 0;
    } else {
        return histogram[source_nid][bin_id];
    }
}

void SKETCHES::deserialize_sketches() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string file_name = config.graph_location + prefix + "_multi_ads_bot.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_multi_ads_bot.idx";
    }
    std::ostringstream oss;
    oss << file_name;
    INFO(oss.str());
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    Treap tmp_treap;
    for (int i = 0; i < graph.n; ++i) {
        ia >> tmp_treap;
        ads_bot.emplace_back(tmp_treap);
    }

    file_name = config.graph_location + prefix + "_histogram.idx";
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(file_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> histogram >> key2value >> value2key;
}

void SKETCHES::deserialize_sketches2() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string file_name;
    if (config.algo == BOTK_SCAN) {
        file_name = config.graph_location + prefix + "_multi_ads_bot.idx";
        if (!exists_test(file_name)) {
            file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_multi_ads_bot.idx";
        }
    } else if (config.algo == MY_ADS) {
        file_name = config.graph_location + prefix + "_my_ads.idx";
        if (!exists_test(file_name)) {
            file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_my_ads.idx";
        }
    }

    std::ostringstream oss;
    oss << file_name;
    INFO(oss.str());
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    Treap tmp_treap;
    ADS ads;
    bot_k = vector<vector<int>>(graph.n, vector<int>{});
    if (config.algo == MY_ADS) {
        bot_k_dis = vector<vector<double>>(graph.n, vector<double>{});
    }
    for (int i = 0; i < graph.n; ++i) {
        if (config.algo == BOTK_SCAN) {
            ia >> tmp_treap;
            bot_k[i] = tmp_treap.get_bot_k();
        } else if (config.algo == MY_ADS) {
            ia >> ads;
            bot_k[i] = ads.get_bot_k(bot_k_dis[i]);
        }
    }
    file_name = config.graph_location + prefix + "_hashmap.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_hashmap.idx";
    }
    std::ostringstream oss3;
    oss3 << file_name;
    INFO(oss3.str());
    std::ifstream info_ifs2(file_name);
    boost::archive::binary_iarchive info_ia2(info_ifs2);
    info_ia2 >> key2value >> value2key;

    prefix += "_bin_" + to_string(config.bin);
    file_name = config.graph_location + prefix + "_histogram.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_bin_" +
                    to_string(config.bin) + "_histogram.idx";
    }
    std::ostringstream oss4;
    oss4 << file_name;
    INFO(oss4.str());
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(file_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> histogram;
    std::ostringstream oss2;
    oss2 << histogram.size() << " " << key2value.size() << " " << value2key.size();
    INFO(oss2.str());
}

void SKETCHES::serialize_sketches() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string idx = config.graph_location + prefix + "_hashmap.idx";
    std::ofstream info_ofs2(idx);
    boost::archive::binary_oarchive info_oa2(info_ofs2);
    info_oa2 << key2value << value2key;
    prefix += +"_bin_" + to_string(config.bin);
    idx = config.graph_location + prefix + "_histogram.idx";
    std::ofstream info_ofs(idx);
    boost::archive::binary_oarchive info_oa(info_ofs);
    info_oa << histogram;
}


void SKETCHES::get_approx_neis() {
    //neis_in_dis = vector<double>(graph.n);
    neis_in_dis_lb = vector<int>(graph.n);
    neis_in_dis_ub = vector<int>(graph.n);
    int bin_id = get_bin_id(config.distance, false);
    for (int i = 0; i < graph.n; ++i) {
        neis_in_dis_lb[i] = query_histogram(i, config.distance, true);
        if (bot_k[i].size() < config.hash_k) {
            neis_in_dis_ub[i] = bot_k[i].size();
        } else {
            neis_in_dis_ub[i] = query_histogram(i, config.distance, false);
        }
    }
}

double SKETCHES::jaccard_with_botk(int u, int v) {
    int u_cur = 0, v_cur = 0, union_size = 0, inter_size = 0;
    while (u_cur != bot_k[u].size() && v_cur != bot_k[v].size() && union_size < config.hash_k) {
        union_size++;
        if (bot_k[u][u_cur] < bot_k[v][v_cur]) {
            u_cur++;
        } else if (bot_k[u][u_cur] > bot_k[v][v_cur]) {
            v_cur++;
        } else {
            inter_size++;
            u_cur++;
            v_cur++;
        }
    }
    if (union_size < config.hash_k - 0.1) {
        union_size += bot_k[u].size() - u_cur + bot_k[v].size() - v_cur;
    }
    union_size = min(config.hash_k, union_size);
    return double(inter_size) / union_size;
}

double SKETCHES::jaccard_with_sketches(int u, int v, double dis) {
#ifdef  _DEBUG_
    assert(dis > -0.5);
#endif
    if (config.sketches_optimize > 1) {
        int u_ub_size = neis_in_dis_ub[u], u_lb_size = neis_in_dis_lb[u];
        int v_ub_size = neis_in_dis_ub[v], v_lb_size = neis_in_dis_lb[v];
        //int u_ub_size = query_histogram(u, config.distance, false), u_lb_size = query_histogram(u, config.distance, true);
        //int v_ub_size = query_histogram(v, config.distance, false), v_lb_size = query_histogram(v, config.distance, true);
        double sim_ub = 0;
        if (u_ub_size < v_lb_size) {
            sim_ub = max(sim_ub, double(u_ub_size) / v_lb_size);
        }
        if (v_ub_size < u_lb_size) {
            sim_ub = max(sim_ub, double(v_ub_size) / u_lb_size);
        }
        if (sim_ub != 0 && cmp_double(sim_ub, config.epsilon) == -1) {
            return false;
        }
        int u_intersection = query_histogram(u, config.distance - dis, true);
        int v_intersection = query_histogram(v, config.distance - dis, true);
        if (u_intersection < v_intersection) {
            swap(u, v);
            swap(u_intersection, v_intersection);
        }
        double sim_lb = double(u_intersection) / (u_ub_size + v_ub_size - u_intersection);
        if (cmp_double(sim_lb, config.epsilon) > -1) {
            return true;
        }
    }
    if (config.sketches_optimize > 0) {
        int u_cur = 0, v_cur = 0, union_size = 0, inter_size = 0;
        while (u_cur != bot_k[u].size() && v_cur != bot_k[v].size() && union_size < config.hash_k) {
            union_size++;
            if (bot_k[u][u_cur] < bot_k[v][v_cur]) {
                u_cur++;
            } else if (bot_k[u][u_cur] > bot_k[v][v_cur]) {
                v_cur++;
            } else {
                inter_size++;
                u_cur++;
                v_cur++;
            }
            if (inter_size > intersection_ub) {
                return true;
            } else if (inter_size < intersection_lb[union_size]) {
                return false;
            }
        }
        if (union_size < config.hash_k) {
            union_size += bot_k[u].size() - u_cur + bot_k[v].size() - v_cur;
        }
        union_size = min(config.hash_k, union_size);
        return cmp_double(double(inter_size) / union_size, config.epsilon) > -1;
    } else {
        return cmp_double(jaccard_with_botk(u, v), config.epsilon) > -1;
    }
}
