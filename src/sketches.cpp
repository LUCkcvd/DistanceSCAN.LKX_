#include <sstream>
#include <iomanip>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <mutex>
#include <chrono>
#include "sketches.h"

// Helper functions
string format_memory_size(long long bytes) {
    if (bytes < 1024) return to_string(bytes) + "B";
    if (bytes < 1024*1024) return to_string(bytes/1024) + "KB";
    if (bytes < 1024*1024*1024) return to_string(bytes/(1024*1024)) + "MB";
    return to_string(bytes/(1024*1024*1024)) + "GB";
}

template<typename T>
void clear_vector(vector<T>& vec) {
    vector<T>().swap(vec);
}

string format_time(double seconds) {
    int hours = seconds / 3600;
    int minutes = (seconds - hours * 3600) / 60;
    int secs = seconds - hours * 3600 - minutes * 60;
    stringstream ss;
    if (hours > 0) ss << hours << "h ";
    if (minutes > 0) ss << minutes << "m ";
    ss << secs << "s";
    return ss.str();
}

string format_eta(double progress, double elapsed) {
    if (progress <= 0) return "N/A";
    double total = elapsed / progress;
    double remaining = total - elapsed;
    return format_time(remaining);
}

vector<unsigned long long> SKETCHES::pick_random_coffs(int k, const unsigned long long &P) {
    vector<unsigned long long> coffs;
    while (k--) {
        coffs.emplace_back(rand_ulong(0, P - 1));
    }
    return coffs;
}

int SKETCHES::double_hashing(const vector<unsigned long long int> &coff_a, const vector<unsigned long long int> &coff_b,
                         int key, int i) {
    unsigned long long int h_0 = (coff_a[0] * key + coff_b[0]) % P % (2 * graph.n);
    unsigned long long int h_1 = (coff_a[1] * key + coff_b[1]) % P % (2 * graph.n);
    return int((h_0 + i * h_1) % (2 * graph.n));
}

void SKETCHES::init_hash_table() {
    auto start_time = std::chrono::steady_clock::now();
    
    vector<unsigned long long> coff_a = pick_random_coffs(2, P);
    vector<unsigned long long> coff_b = pick_random_coffs(2, P);

    key2value = vector<int>(graph.n, -1);
    value2key = vector<int>(2 * graph.n, -1);

    int collisions = 0;
    for (int key = 0; key < key2value.size(); ++key) {
        for (int j = 0; j < value2key.size(); ++j) {
            int value = double_hashing(coff_a, coff_b, key, j);
            if (value2key[value] == -1) {
                value2key[value] = key;
                key2value[key] = value;
                if (j > 0) collisions++;
                break;
            }
        }
        
        // Progress update every 25%
        if (key % (graph.n/4) == 0) {
            auto current_time = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
            double progress = (double)key / graph.n;
            INFO("Hash table: " + to_string((int)(progress * 100)) + "% - ETA: " + 
                 format_eta(progress, elapsed));
        }
    }
    
    clear_vector(coff_a);
    clear_vector(coff_b);

    auto end_time = std::chrono::steady_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    INFO("Hash table initialization complete in " + format_time(total_duration) + " with " + 
         to_string(collisions) + " collisions");
}

void SKETCHES::update_memory_stats() {
    double current_mem = get_proc_memory();
    progress.last_memory_usage = current_mem;
    if (current_mem > progress.peak_memory_usage) {
        progress.peak_memory_usage = current_mem;
    }
    peak_memory_mb = std::max(peak_memory_mb, current_mem);
}

void SKETCHES::report_progress(bool force) {
    static std::mutex progress_mutex;
    std::lock_guard<std::mutex> lock(progress_mutex);

    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - progress.start_time).count();
    
    if (!force && progress.nodes_processed % PROGRESS_UPDATE_NODES != 0) {
        return;
    }

    double percent = (double)progress.nodes_processed / graph.n * 100;
    double edges_per_sec = progress.edges_processed / (elapsed > 0 ? elapsed : 1);
    double eta = (graph.n - progress.nodes_processed) / (progress.nodes_processed / (double)elapsed);
    
    update_memory_stats();
    
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2)
       << "Progress: " << percent << "% "
       << "(" << progress.nodes_processed << "/" << graph.n << " nodes) - "
       << format_memory_size(progress.last_memory_usage * 1024) << " RAM "
       << "(" << format_memory_size(progress.peak_memory_usage * 1024) << " peak) - "
       << (int)edges_per_sec << " edges/s - "
       << "ETA: " << format_time(eta);
    
    INFO(ss.str());
}

void SKETCHES::construct_sketches() {
    Timer timer(CONSTRUCT_SKETCH_TIME);
    progress.reset();
    
    const string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    const string idx = config.graph_location + prefix + (config.algo == BOTK_SCAN ? "_multi_ads_bot.idx" : "_my_ads.idx");
    
    std::ofstream ofs(idx, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Failed to open output file: " + idx);
    }
    boost::archive::binary_oarchive oa(ofs);

    // Initialize construction parameters
    const int num_bins = get_bin_id(config.max_distance, false) + 1;
    const int total_nodes = graph.n;
    
    // Pre-allocate histogram with reserve to prevent reallocation
    histogram.resize(total_nodes);
    
    // Configure OpenMP
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    INFO("Using " + to_string(num_threads) + " threads");
    
    // Process nodes in chunks
    for (int chunk_start = 0; chunk_start < total_nodes; chunk_start += HISTOGRAM_CHUNK_SIZE) {
        const int chunk_end = min(chunk_start + HISTOGRAM_CHUNK_SIZE, total_nodes);
            
        // Allocate histogram bins for current chunk
        for (int i = chunk_start; i < chunk_end; i++) {
            histogram[i].resize(num_bins);
            fill(histogram[i].begin(), histogram[i].end(), 0);
        }

        // Process nodes in parallel batches
        #pragma omp parallel
        {
            // Thread-local storage
            vector<idpair> dis_source_vec(graph.n, make_pair(-1, -1));
            priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
            size_t local_edges_processed = 0;

            #pragma omp for schedule(dynamic) nowait
            for (int batch_start = chunk_start; batch_start < chunk_end; batch_start += BATCH_SIZE) {
                const int batch_end = min(batch_start + BATCH_SIZE, chunk_end);
            
                for (int source_nid = batch_start; source_nid < batch_end; ++source_nid) {
                    // Reset distance vector for new node
                    fill(dis_source_vec.begin(), dis_source_vec.end(), make_pair(-1, -1));
                
                    // Initialize data structures based on algorithm
                    Treap treap;
                    ADS ads;
                    if (config.algo == BOTK_SCAN) {
                        treap = Treap();
                    } else if (config.algo == MY_ADS) {
                        ads = ADS();
                    }
                
                    // Initialize source node and priority queue
                    idpair source = make_pair(source_nid, 0);
                    dis_source_vec[source_nid] = source;
                    while (!pq.empty()) pq.pop();
                    pq.push(source);
                    
                    size_t edges_processed = 0;
                    size_t nodes_visited = 0;
                
                    // Process nodes in priority queue
                    while (!pq.empty()) {
                        auto cur_node = pq.top();
                        pq.pop();
                        nodes_visited++;
                        
                        if (cmp_double(cur_node.second, dis_source_vec[cur_node.first].second) == 1) continue;

                        // Update histogram and data structures
                        double path_distance = cur_node.second;
                        update_histogram(source_nid, path_distance);
                        
                        // Update sketch data structures
                        if (config.algo == BOTK_SCAN) {
                            treap.insert_botk(key2value[cur_node.first], path_distance);
                        } else if (config.algo == MY_ADS) {
                            ads.insert_botk(key2value[cur_node.first], path_distance);
                        }

                        // Process neighbors
                        for (int nei_id : graph.adj_list[cur_node.first]) {
                            double nei_dis = cur_node.second + graph.edge_weight[cur_node.first][nei_id];
                            if (cmp_double(nei_dis, config.max_distance) < 1) {
                                if (dis_source_vec[nei_id].first != source_nid ||
                                    cmp_double(dis_source_vec[nei_id].second, nei_dis) > 0) {
                                    dis_source_vec[nei_id].first = source_nid;
                                    dis_source_vec[nei_id].second = nei_dis;
                                    pq.push(make_pair(nei_id, nei_dis));
                                    edges_processed++;
                                }
                            }
                        }
                    }
                    
                    {
                        Timer timer(SAVE_BOTK_TIME);
                        if (config.algo == BOTK_SCAN) {
                            oa << treap;
                        } else if (config.algo == MY_ADS) {
                            oa << ads;
                        }
                    }
                    
                    {
                        Timer timer(UPDATE_HISTO_TIME);
                        for (int j = 1; j < histogram[source_nid].size(); ++j) {
                            histogram[source_nid][j] += histogram[source_nid][j - 1];
                        }
                    }

                    #pragma omp critical
                    {
                        progress.nodes_processed++;
                        progress.edges_processed += edges_processed;
                        
                        if (progress.nodes_processed % PROGRESS_UPDATE_NODES == 0) {
                            report_progress();
                        }
                    }
                    
                    // Clear queue memory
                    priority_queue<idpair, vector<idpair>, cmp_idpair>().swap(pq);
                }
            }
            
            // Force memory cleanup after each batch
            clear_vector(dis_source_vec);
        }
        
        // Save and clear histogram data for this chunk
        {
            Timer timer(SAVE_BOTK_TIME);
            string chunk_file = config.graph_location + prefix + "_histogram_chunk_" + 
                              to_string(chunk_start/HISTOGRAM_CHUNK_SIZE) + ".tmp";
            std::ofstream chunk_ofs(chunk_file, std::ios::binary);
            boost::archive::binary_oarchive chunk_oa(chunk_ofs);
            
            for (int i = chunk_start; i < chunk_end; i++) {
                chunk_oa << histogram[i];
                vector<int>().swap(histogram[i]);
            }
        }
        
        // Use fsync instead of system("sync") to avoid vfork
        int fd = open(config.graph_location.c_str(), O_RDONLY);
        if (fd >= 0) {
            fsync(fd);
            close(fd);
        }
    }
    
    // Merge histogram chunks
    string final_histogram_file = config.graph_location + prefix + "_histogram.idx";
    std::ofstream final_ofs(final_histogram_file, std::ios::binary);
    boost::archive::binary_oarchive final_oa(final_ofs);
    
    int num_chunks = (total_nodes + HISTOGRAM_CHUNK_SIZE - 1) / HISTOGRAM_CHUNK_SIZE;
    for (int i = 0; i < num_chunks; i++) {
        string chunk_file = config.graph_location + prefix + "_histogram_chunk_" + to_string(i) + ".tmp";
        std::ifstream chunk_ifs(chunk_file, std::ios::binary);
        boost::archive::binary_iarchive chunk_ia(chunk_ifs);
        
        for (int j = 0; j < HISTOGRAM_CHUNK_SIZE && (i * HISTOGRAM_CHUNK_SIZE + j) < total_nodes; j++) {
            vector<int> hist_entry;
            chunk_ia >> hist_entry;
            final_oa << hist_entry;
        }
        
        chunk_ifs.close();
        std::remove(chunk_file.c_str());
    }
    
    auto end_time = std::chrono::steady_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(
        end_time - progress.start_time).count();
    
    INFO("Construction complete in " + format_time(total_duration) + " - " + 
         format_memory_size(get_proc_memory() * 1024) + " peak RAM");
}

void SKETCHES::deserialize_sketches() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string file_name = config.graph_location + prefix + "_multi_ads_bot.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_multi_ads_bot.idx";
    }
    
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    Treap tmp_treap;
    for (int i = 0; i < graph.n; ++i) {
        ia >> tmp_treap;
        ads_bot.emplace_back(tmp_treap);
    }

    file_name = config.graph_location + prefix + "_histogram.idx";
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
    std::ifstream info_ifs2(file_name);
    boost::archive::binary_iarchive info_ia2(info_ifs2);
    info_ia2 >> key2value >> value2key;

    prefix += "_bin_" + to_string(config.bin);
    file_name = config.graph_location + prefix + "_histogram.idx";
    if (!exists_test(file_name)) {
        file_name = config.graph_location + "dmax_" + to_str(config.max_distance) + "_k_65536_bin_" +
                    to_string(config.bin) + "_histogram.idx";
    }
    std::ifstream info_ifs(file_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> histogram;
}

void SKETCHES::serialize_sketches() {
    string prefix = "dmax_" + to_str(config.max_distance) + "_k_" + to_string(config.hash_k);
    string idx = config.graph_location + prefix + "_hashmap.idx";
    std::ofstream info_ofs2(idx);
    boost::archive::binary_oarchive info_oa2(info_ofs2);
    info_oa2 << key2value << value2key;
    
    prefix += "_bin_" + to_string(config.bin);
    idx = config.graph_location + prefix + "_histogram.idx";
    std::ofstream info_ofs(idx);
    boost::archive::binary_oarchive info_oa(info_ofs);
    info_oa << histogram;
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

void SKETCHES::get_approx_neis() {
    neis_in_dis_lb = vector<int>(graph.n);
    neis_in_dis_ub = vector<int>(graph.n);
    
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
    if (config.sketches_optimize > 1) {
        int u_ub_size = neis_in_dis_ub[u], u_lb_size = neis_in_dis_lb[u];
        int v_ub_size = neis_in_dis_ub[v], v_lb_size = neis_in_dis_lb[v];
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
