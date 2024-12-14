#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <fstream>
#include "libs/json.hpp"

#include <sys/resource.h> // Para medir consumo de memória e tempos de CPU
#include <sys/time.h>     // Para tempos adicionais
#include <mpi.h>
#include <cstring>
#include <stack>

using json = nlohmann::json;

// Estruturas para passar dados para as threads
struct WCCData {
    const std::unordered_map<int, std::vector<int>>* adj_list;
    const std::set<int>* nodes;
    int edge_count;
    int largest_wcc_size = 0;
    int largest_wcc_edges = 0;
    double fraction_nodes = 0;
    double fraction_edges = 0;
};

struct SCCData {
    const std::unordered_map<int, std::vector<int>>* adj_list;
    const std::set<int>* nodes;
    int edge_count;
    int largest_scc_size = 0;
    int largest_scc_edges = 0;
    double fraction_nodes = 0;
    double fraction_edges = 0;
};

struct ClusteringData {
    const std::unordered_map<int, std::vector<int>>* adj_list;
    const std::set<int>* nodes;
    double average_clustering_coefficient = 0.0;
};

struct TriangleData {
    const std::unordered_map<int, std::vector<int>>* adj_list;
    const std::set<int>* nodes;
    int total_triangles = 0;
};

// Função para obter métricas de execução
void collectExecutionMetrics(json &output) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // Tempo de CPU ocupado (em microssegundos)
    long user_time_ms = usage.ru_utime.tv_sec * 1000 + usage.ru_utime.tv_usec / 1000;
    long system_time_ms = usage.ru_stime.tv_sec * 1000 + usage.ru_stime.tv_usec / 1000;

    // Consumo de memória máxima (em KB)
    long memory_kb = usage.ru_maxrss;

    // Adiciona ao JSON
    output["execution"] = {
        {"cpu_time_user_ms", user_time_ms},
        {"cpu_time_system_ms", system_time_ms},
        {"memory_kb", memory_kb}
    };

    // Mostra no terminal
    std::cout << "\nExecution Metrics:\n";
    std::cout << "CPU Time (User): " << user_time_ms << " ms\n";
    std::cout << "CPU Time (System): " << system_time_ms << " ms\n";
    std::cout << "Memory Usage: " << memory_kb << " KB\n";
}

// Função para realizar BFS e encontrar o componente conectado
void bfs(const std::unordered_map<int, std::vector<int>> &adj_list, int start_node, std::set<int> &visited, std::set<int> &component) {
    std::queue<int> to_visit;
    to_visit.push(start_node);
    visited.insert(start_node);
    component.insert(start_node);

    while (!to_visit.empty()) {
        int current = to_visit.front();
        to_visit.pop();

        for (int neighbor : adj_list.at(current)) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                component.insert(neighbor);
                to_visit.push(neighbor);
            }
        }
    }
}

// Função para carregar o grafo ordenado e sem duplicação de arestas
bool loadGraphDirected(const std::string& filename, std::unordered_map<int, std::vector<int>>& adj_list, std::set<int>& nodes, int& edge_count) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << "\n";
        return false;
    }

    std::unordered_map<int, int> degree_map;
    std::string line;
    edge_count = 0;

    while (std::getline(file, line)) {
        // Ignorar linhas que começam com '#'
        if (!line.empty() && line[0] == '#') {
            std::cout << "Linha ignorada: " << line << "\n";
            continue;
        }

        std::istringstream iss(line);
        int node1, node2;
        if (iss >> node1 >> node2) {
            if (std::find(adj_list[node1].begin(), adj_list[node1].end(), node2) == adj_list[node1].end()) {
                adj_list[node1].push_back(node2);
                degree_map[node1]++;
                degree_map[node2]++;
                nodes.insert(node1);
                nodes.insert(node2);
                edge_count++;
            }
        } else {
            std::cout << "Linha inválida: " << line << "\n";
        }
    }

    // Ordenar adj_list em ordem decrescente de grau
    std::vector<std::pair<int, int>> sorted_nodes(degree_map.begin(), degree_map.end());
    std::sort(sorted_nodes.begin(), sorted_nodes.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    std::unordered_map<int, std::vector<int>> sorted_adj_list;
    for (const auto& [node, _] : sorted_nodes) {
        sorted_adj_list[node] = adj_list[node];
    }
    adj_list = std::move(sorted_adj_list);

    file.close();
    return true;
}

// Função para criar subgrafo com os top_k nós mais conectados
void createSubgraph(const std::unordered_map<int, std::vector<int>>& adj_list,
                    std::unordered_map<int, std::vector<int>>& subgraph,
                    std::set<int>& subgraph_nodes,
                    int& subgraph_edge_count,
                    int top_k) {
    // Calcular grau de cada nó
    std::vector<std::pair<int, int>> node_degrees;
    for (const auto& pair : adj_list) {
        int degree = pair.second.size();
        node_degrees.emplace_back(pair.first, degree);
    }

    // Ordenar nós por grau (decrescente)
    std::sort(node_degrees.begin(), node_degrees.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    // Selecionar os top_k nós
    subgraph_nodes.clear();
    for (int i = 0; i < std::min(top_k, static_cast<int>(node_degrees.size())); ++i) {
        subgraph_nodes.insert(node_degrees[i].first);
    }

    // Construir subgrafo
    subgraph_edge_count = 0;
    for (int node : subgraph_nodes) {
        if (adj_list.find(node) != adj_list.end()) {
            for (int neighbor : adj_list.at(node)) {
                if (subgraph_nodes.find(neighbor) != subgraph_nodes.end()) {
                    subgraph[node].push_back(neighbor);
                    ++subgraph_edge_count;
                }
            }
        }
    }
}


// Função para calcular o maior componente fracamente conectado (WCC)
void computeLargestWCCOpenMPI(WCCData& data, int rank, int size) {
    int num_nodes = data.nodes->size();
    std::vector<int> node_indices(data.nodes->begin(), data.nodes->end());

    int local_start = rank * (num_nodes / size);
    int local_end = (rank == size - 1) ? num_nodes : (rank + 1) * (num_nodes / size);

    std::unordered_set<int> visited;

    // Identificar componentes conexos locais
    for (int i = local_start; i < local_end; ++i) {
        int node = node_indices[i];
        if (visited.find(node) == visited.end()) {
            std::queue<int> to_visit;
            to_visit.push(node);
            visited.insert(node);
            int component_size = 0;

            while (!to_visit.empty()) {
                int current = to_visit.front();
                to_visit.pop();
                component_size++;

                if (data.adj_list->find(current) != data.adj_list->end()) {
                    for (int neighbor : data.adj_list->at(current)) {
                        if (visited.find(neighbor) == visited.end()) {
                            visited.insert(neighbor);
                            to_visit.push(neighbor);
                        }
                    }
                }
            }

            data.largest_wcc_size = std::max(data.largest_wcc_size, component_size);
        }
    }

    // Coletar componentes globais
    MPI_Allreduce(MPI_IN_PLACE, &data.largest_wcc_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    data.fraction_nodes = static_cast<double>(data.largest_wcc_size) / num_nodes;
}



//Função para calcular o maior SCC
void computeLargestSCCOpenMPI(SCCData& data, int rank, int size) {
    int num_nodes = data.nodes->size();
    std::vector<int> node_indices(data.nodes->begin(), data.nodes->end());

    int local_start = rank * (num_nodes / size);
    int local_end = (rank == size - 1) ? num_nodes : (rank + 1) * (num_nodes / size);

    std::unordered_map<int, std::vector<int>> reversed_graph;

    // Criar grafo reverso local
    for (int i = local_start; i < local_end; ++i) {
        int node = node_indices[i];
        if (data.adj_list->find(node) != data.adj_list->end()) {
            for (int neighbor : data.adj_list->at(node)) {
                reversed_graph[neighbor].push_back(node);
            }
        }
    }

    // Coletar o grafo reverso global
    std::unordered_map<int, std::vector<int>> global_reversed_graph;
    for (const auto& [key, neighbors] : reversed_graph) {
        MPI_Allreduce(MPI_IN_PLACE, global_reversed_graph[key].data(), neighbors.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    std::vector<bool> visited(num_nodes, false);
    std::stack<int> finish_stack;

    // Primeira DFS para ordem de término
    for (int i = local_start; i < local_end; ++i) {
        if (!visited[i]) {
            std::stack<int> dfs_stack;
            dfs_stack.push(i);

            while (!dfs_stack.empty()) {
                int current = dfs_stack.top();
                dfs_stack.pop();
                visited[current] = true;

                for (int neighbor : data.adj_list->at(node_indices[current])) {
                    int neighbor_index = std::distance(node_indices.begin(), std::find(node_indices.begin(), node_indices.end(), neighbor));
                    if (!visited[neighbor_index]) {
                        dfs_stack.push(neighbor_index);
                    }
                }
                finish_stack.push(current);
            }
        }
    }

    std::fill(visited.begin(), visited.end(), false);

    // Segunda DFS para identificar SCCs
    while (!finish_stack.empty()) {
        int node = finish_stack.top();
        finish_stack.pop();

        if (!visited[node]) {
            std::stack<int> dfs_stack;
            dfs_stack.push(node);
            int scc_size = 0;

            while (!dfs_stack.empty()) {
                int current = dfs_stack.top();
                dfs_stack.pop();
                visited[current] = true;
                scc_size++;

                for (int neighbor : global_reversed_graph[node_indices[current]]) {
                    int neighbor_index = std::distance(node_indices.begin(), std::find(node_indices.begin(), node_indices.end(), neighbor));
                    if (!visited[neighbor_index]) {
                        dfs_stack.push(neighbor_index);
                    }
                }
            }

            data.largest_scc_size = std::max(data.largest_scc_size, scc_size);
        }
    }

    data.fraction_nodes = static_cast<double>(data.largest_scc_size) / num_nodes;
}

// Função para calcular as frações de nós e arestas
void computeFractions(int largest_wcc_size, int largest_wcc_edges, int total_nodes, int total_edges, double &fraction_nodes, double &fraction_edges) {
    fraction_nodes = static_cast<double>(largest_wcc_size) / total_nodes;
    fraction_edges = static_cast<double>(largest_wcc_edges) / total_edges;
}

// Adaptar funções para operar com o subgrafo
void computeFractionsForSubgraph(int largest_component_size, int largest_component_edges, int total_nodes, int total_edges, double &fraction_nodes, double &fraction_edges) {
    fraction_nodes = static_cast<double>(largest_component_size) / total_nodes;
    fraction_edges = static_cast<double>(largest_component_edges) / total_edges;
}

// Função para calcular o coeficiente de agrupamento de um nó
double computeAverageClusteringCoefficientOpenMPI(const ClusteringData& data, int rank, int size) {
    double total_clustering = 0.0;
    int node_count = 0;

    std::vector<int> nodes_vector(data.nodes->begin(), data.nodes->end());

    int num_nodes = nodes_vector.size();
    int local_start = rank * (num_nodes / size);
    int local_end = (rank == size - 1) ? num_nodes : (rank + 1) * (num_nodes / size);

    for (int i = local_start; i < local_end; ++i) {
        int node = nodes_vector[i];
        if (data.adj_list->find(node) == data.adj_list->end()) continue;
        const auto& neighbors = data.adj_list->at(node);
        int degree = neighbors.size();
        if (degree < 2) continue;

        int triangle_count = 0;
        std::set<int> neighbor_set(neighbors.begin(), neighbors.end());
        for (size_t j = 0; j < neighbors.size(); ++j) {
            for (size_t k = j + 1; k < neighbors.size(); ++k) {
                if (neighbor_set.find(neighbors[k]) != neighbor_set.end() &&
                    data.adj_list->find(neighbors[j]) != data.adj_list->end() &&
                    std::find(data.adj_list->at(neighbors[j]).begin(), data.adj_list->at(neighbors[j]).end(), neighbors[k]) != data.adj_list->at(neighbors[j]).end()) {
                    ++triangle_count;
                }
            }
        }

        total_clustering += static_cast<double>(2 * triangle_count) / (degree * (degree - 1));
        ++node_count;
    }

    double global_clustering;
    MPI_Reduce(&total_clustering, &global_clustering, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_clustering / node_count;
}


// Função para contar triângulos no grafo
int countTrianglesOpenMPI(const TriangleData& data, int rank, int size) {
    int total_triangles = 0;

    std::vector<int> nodes_vector(data.nodes->begin(), data.nodes->end());

    int num_nodes = nodes_vector.size();
    int local_start = rank * (num_nodes / size);
    int local_end = (rank == size - 1) ? num_nodes : (rank + 1) * (num_nodes / size);

    for (int i = local_start; i < local_end; ++i) {
        int node = nodes_vector[i];
        if (data.adj_list->find(node) == data.adj_list->end()) continue;
        const auto& neighbors = data.adj_list->at(node);
        std::set<int> neighbor_set(neighbors.begin(), neighbors.end());

        for (size_t j = 0; j < neighbors.size(); ++j) {
            for (size_t k = j + 1; k < neighbors.size(); ++k) {
                if (neighbor_set.find(neighbors[k]) != neighbor_set.end() &&
                    data.adj_list->find(neighbors[j]) != data.adj_list->end() &&
                    std::find(data.adj_list->at(neighbors[j]).begin(), data.adj_list->at(neighbors[j]).end(), neighbors[k]) != data.adj_list->at(neighbors[j]).end()) {
                    ++total_triangles;
                }
            }
        }
    }

    int global_triangles;
    MPI_Reduce(&total_triangles, &global_triangles, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_triangles / 3;
}



// Função para calcular a Fração de Triângulos Fechados
double calculateFractionOfClosedTriangles(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes) {
    int closed_triangles = 0; // Número de triângulos fechados
    int open_triangles = 0;  // Número de triângulos abertos

    for (int node : nodes) {
        const auto &neighbors = adj_list.at(node);

        // Usar conjunto para verificar pares de vizinhos
        std::set<int> neighbor_set(neighbors.begin(), neighbors.end());

        // Itera sobre todos os pares de vizinhos
        for (size_t i = 0; i < neighbors.size(); ++i) {
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                if (neighbor_set.find(neighbors[j]) != neighbor_set.end() &&
                    std::find(adj_list.at(neighbors[i]).begin(), adj_list.at(neighbors[i]).end(), neighbors[j]) != adj_list.at(neighbors[i]).end()) {
                    closed_triangles++; // Triângulo fechado
                } else {
                    open_triangles++; // Triângulo aberto
                }
            }
        }
    }

    // Cada triângulo fechado é contado três vezes, então divide por 3
    closed_triangles /= 3;

    // Total de triângulos possíveis (fechados + abertos)
    int total_triangles = closed_triangles + open_triangles;

    // Retorna a fração de triângulos fechados
    if (total_triangles == 0) return 0.0;
    return static_cast<double>(closed_triangles) / total_triangles;
}

//triangulo fechado com subgrafo (não consegui fazer de outra maneira / não sei como fazer ;( sorry alysson)
double closedTrianglesSubgraph(const std::unordered_map<int, std::vector<int>>& subgraph, const std::set<int>& subgraph_nodes) {
    int closed_triangles = 0;
    int open_triangles = 0;

    for (int node : subgraph_nodes) {
        if (subgraph.find(node) == subgraph.end()) continue; // Verifique se o nó existe
        const auto& neighbors = subgraph.at(node);
        std::set<int> neighbor_set(neighbors.begin(), neighbors.end());

        for (size_t i = 0; i < neighbors.size(); ++i) {
            if (subgraph.find(neighbors[i]) == subgraph.end()) continue; // Verifique se o nó existe
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                if (subgraph.find(neighbors[j]) == subgraph.end()) continue; // Verifique se o nó existe
                if (neighbor_set.find(neighbors[j]) != neighbor_set.end() &&
                    std::find(subgraph.at(neighbors[i]).begin(), subgraph.at(neighbors[i]).end(), neighbors[j]) != subgraph.at(neighbors[i]).end()) {
                    closed_triangles++;
                } else {
                    open_triangles++;
                }
            }
        }
    }

    closed_triangles /= 3; // Cada triângulo é contado 3 vezes.
    int total_triangles = closed_triangles + open_triangles;
    if (total_triangles == 0) return 0.0;
    return static_cast<double>(closed_triangles) / total_triangles;
}

// Função para calcular o diâmetro do grafo
int calculateGraphDiameter(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes) {
    int diameter = 0;

    for (int node : nodes) {
        std::queue<int> to_visit;
        std::unordered_map<int, int> distances;
        to_visit.push(node);
        distances[node] = 0;

        while (!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();

            for (int neighbor : adj_list.at(current)) {
                if (distances.find(neighbor) == distances.end()) {
                    distances[neighbor] = distances[current] + 1;
                    diameter = std::max(diameter, distances[neighbor]);
                    to_visit.push(neighbor);
                }
            }
        }
    }

    return diameter;
}

//diametro do subgrafo 
int graphDiameterSubgraph(const std::unordered_map<int, std::vector<int>>& subgraph, const std::set<int>& subgraph_nodes) {
    int diameter = 0;

    for (int node : subgraph_nodes) {
        if (subgraph.find(node) == subgraph.end()) continue; // Verifique se o nó existe
        std::queue<int> to_visit;
        std::unordered_map<int, int> distances;
        to_visit.push(node);
        distances[node] = 0;

        while (!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();

            if (subgraph.find(current) == subgraph.end()) continue; // Verifique se o nó existe
            for (int neighbor : subgraph.at(current)) {
                if (distances.find(neighbor) == distances.end()) {
                    distances[neighbor] = distances[current] + 1;
                    diameter = std::max(diameter, distances[neighbor]);
                    to_visit.push(neighbor);
                }
            }
        }
    }

    return diameter;
}

// Função para calcular o diâmetro efetivo a 90%
double calculateEffectiveDiameter(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes, double percentile = 0.9) {
    std::vector<int> distances;

    // Calcular distâncias entre todos os pares de nós
    for (int node : nodes) {
        std::queue<int> to_visit;
        std::unordered_map<int, int> dist;
        to_visit.push(node);
        dist[node] = 0;

        while (!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();

            for (int neighbor : adj_list.at(current)) {
                if (dist.find(neighbor) == dist.end()) {
                    dist[neighbor] = dist[current] + 1;
                    to_visit.push(neighbor);
                    distances.push_back(dist[neighbor]);
                }
            }
        }
    }

    // Remover distâncias inválidas
    distances.erase(std::remove_if(distances.begin(), distances.end(), [](int d) { return d <= 0; }), distances.end());

    if (distances.empty()) {
        return 0.0; // Retornar 0 se não houver pares válidos
    }

    // Ordenar as distâncias
    std::sort(distances.begin(), distances.end());

    // Calcular distribuição cumulativa
    size_t total_pairs = distances.size();
    double exact_index = percentile * total_pairs;

    // Encontrar o percentil exato
    size_t lower_index = static_cast<size_t>(std::floor(exact_index));
    size_t upper_index = static_cast<size_t>(std::ceil(exact_index));

    if (lower_index == upper_index || upper_index >= total_pairs) {
        return distances[lower_index];
    }

    // Interpolação linear entre índices
    double fraction = exact_index - lower_index;
    return distances[lower_index] + fraction * (distances[upper_index] - distances[lower_index]);
}

//diametro efetivo do subgrafo
double effectiveDiameterSubgraph(const std::unordered_map<int, std::vector<int>>& subgraph, const std::set<int>& subgraph_nodes, double percentile = 0.9) {
    std::vector<int> distances;

    for (int node : subgraph_nodes) {
        if (subgraph.find(node) == subgraph.end()) continue; // Verifique se o nó existe
        std::queue<int> to_visit;
        std::unordered_map<int, int> dist;
        to_visit.push(node);
        dist[node] = 0;

        while (!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();

            if (subgraph.find(current) == subgraph.end()) continue; // Verifique se o nó existe
            for (int neighbor : subgraph.at(current)) {
                if (dist.find(neighbor) == dist.end()) {
                    dist[neighbor] = dist[current] + 1;
                    to_visit.push(neighbor);
                    distances.push_back(dist[neighbor]);
                }
            }
        }
    }

    distances.erase(std::remove_if(distances.begin(), distances.end(), [](int d) { return d <= 0; }), distances.end());
    if (distances.empty()) return 0.0;

    std::sort(distances.begin(), distances.end());
    size_t total_pairs = distances.size();
    double exact_index = percentile * total_pairs;

    size_t lower_index = static_cast<size_t>(std::floor(exact_index));
    size_t upper_index = static_cast<size_t>(std::ceil(exact_index));

    if (lower_index == upper_index || upper_index >= total_pairs) {
        return distances[lower_index];
    }

    double fraction = exact_index - lower_index;
    return distances[lower_index] + fraction * (distances[upper_index] - distances[lower_index]);
}

int main(int argc, char** argv) {
    std::unordered_map<int, std::vector<int>> adj_list;
    std::set<int> nodes;
    int edge_count = 0;

    // Início do timer para carregar o grafo
    auto start_time = std::chrono::high_resolution_clock::now();

    // Carregar o grafo
    if (!loadGraphDirected("web-Google.txt", adj_list, nodes, edge_count)) {
        return 1;
    }


    auto loadGraph_end = std::chrono::high_resolution_clock::now();
    std::cout << "\nTime to load graph: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loadGraph_end - start_time).count() << " ms\n";

    // Criar subgrafo com os 1000 nós mais conectados
    /*
    */
    // Criar subgrafo
    std::unordered_map<int, std::vector<int>> subgraph;
    std::set<int> subgraph_nodes;
    int subgraph_edge_count = 0;
    int top_k = 1000;
    createSubgraph(adj_list, subgraph, subgraph_nodes, subgraph_edge_count, top_k);
    
    // Estruturas de dados para do grafo original
    /*
    WCCData wcc_data = {&adj_list, &nodes, edge_count};
    SCCData scc_data = {&adj_list, &nodes, edge_count};
    ClusteringData clustering_data = {&adj_list, &nodes};
    TriangleData triangle_data = {&adj_list, &nodes};
    */
    
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Estruturas de dados para o subgrafo
    WCCData wcc_data = {&subgraph, &subgraph_nodes, subgraph_edge_count};
    SCCData scc_data = {&subgraph, &subgraph_nodes, subgraph_edge_count};
    ClusteringData clustering_data = {&subgraph, &subgraph_nodes};
    TriangleData triangle_data = {&subgraph, &subgraph_nodes};

    double average_clustering_coefficient = 0;
    int total_triangles = 0;

    // NÓS
    auto nodes_start = std::chrono::high_resolution_clock::now();
    std::cout << "\nnodes: " << subgraph_nodes.size() << "\n";
    //std::cout << "\nnodes: " << nodes.size() << "\n";
    auto nodes_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for nodes: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(nodes_end - nodes_start).count() << " ms\n\n";

    // ARESTAS
    auto edges_start = std::chrono::high_resolution_clock::now();
    std::cout << "edges: " << subgraph_edge_count << "\n";
    //std::cout << "edges: " << edge_count << "\n";
    auto edges_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for edges: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(edges_end - edges_start).count() << " ms\n\n";

    auto task_start = std::chrono::high_resolution_clock::now();
    
    // Executar as funções simultaneamente com OpenMPI 
    // Execução paralela das funções
    /*
    if (rank == 0) {
        MPI_Request requests[4];

        MPI_Isend(&wcc_data, sizeof(WCCData), MPI_BYTE, 1, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(&scc_data, sizeof(SCCData), MPI_BYTE, 2, 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(&clustering_data, sizeof(ClusteringData), MPI_BYTE, 3, 0, MPI_COMM_WORLD, &requests[2]);
        MPI_Isend(&triangle_data, sizeof(TriangleData), MPI_BYTE, 4, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
    } else {
        if (rank == 1) computeLargestWCCOpenMPI(wcc_data, rank, size);
        if (rank == 2) computeLargestSCCOpenMPI(scc_data, rank, size);
        if (rank == 3) computeAverageClusteringCoefficientOpenMPI(clustering_data, rank, size);
        if (rank == 4) countTrianglesOpenMPI(triangle_data, rank, size);
    }
    */

    if (rank == 0) {
        computeLargestWCCOpenMPI(wcc_data, rank, size);
    } else if (rank == 1) {
        computeLargestSCCOpenMPI(scc_data, rank, size);
    }else if (rank == 2) {
        computeAverageClusteringCoefficientOpenMPI(clustering_data, rank, size);
    }
    else if (rank == 3) {
        countTrianglesOpenMPI(triangle_data, rank, size);
    } 
     else {
        printf("Rank %d não está atribuído a nenhuma tarefa.\n", rank);
    }
    
    MPI_Finalize();

    auto task_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for tasks: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(task_end - task_start).count() << " ms\n\n";

    /*
    */
    std::cout << "WCC Nodes: " << wcc_data.largest_wcc_size << "\n";
    std::cout << "WCC Fraction of Nodes: " << wcc_data.fraction_nodes << "\n";
    std::cout << "WCC Edges: " << wcc_data.largest_wcc_edges << "\n";
    std::cout << "WCC Fraction of Edges: " << wcc_data.fraction_edges << "\n";

    std::cout << "SCC Nodes: " << scc_data.largest_scc_size << "\n";
    std::cout << "SCC Fraction of Nodes: " << scc_data.fraction_nodes << "\n";
    std::cout << "SCC Edges: " << scc_data.largest_scc_edges << "\n";
    std::cout << "SCC Fraction of Edges: " << scc_data.fraction_edges << "\n";

    std::cout << "Average Clustering Coefficient: " << average_clustering_coefficient << "\n";

    std::cout << "Total Triangles: " << total_triangles << "\n";

    // FRAÇÃO DE TRIÂNGULOS FECHADOS
    auto closed_triangles_start = std::chrono::high_resolution_clock::now();
    //double fraction_closed_triangles = calculateFractionOfClosedTriangles(adj_list, nodes);
    double fraction_closed_triangles = closedTrianglesSubgraph(subgraph, subgraph_nodes);
    auto closed_triangles_end = std::chrono::high_resolution_clock::now();
    std::cout << "terminou " << "\n";
    std::cout << "fraction of closed triangles: " << std::fixed << std::setprecision(4) << fraction_closed_triangles << "\n";
    std::cout << "Time for fraction of closed triangles: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(closed_triangles_end - closed_triangles_start).count() << " ms\n\n";

    // DIÂMETRO DO GRAFO
    auto diameter_start = std::chrono::high_resolution_clock::now();
    int graph_diameter = 0;
    //graph_diameter = calculateGraphDiameter(adj_list, nodes);
    graph_diameter = graphDiameterSubgraph(subgraph, subgraph_nodes);
    auto diameter_end = std::chrono::high_resolution_clock::now();
    std::cout << "graph diameter: " << graph_diameter << "\n";
    std::cout << "Time for graph diameter: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(diameter_end - diameter_start).count() << " ms\n\n";

    // DIÂMETRO EFETIVO DO GRAFO
    auto effective_diameter_start = std::chrono::high_resolution_clock::now();
    double effective_diameter = 0;
    //effective_diameter = calculateEffectiveDiameter(adj_list, nodes);
    effective_diameter = effectiveDiameterSubgraph(subgraph, subgraph_nodes);
    auto effective_diameter_end = std::chrono::high_resolution_clock::now();
    std::cout << "effective diameter (90%): " << std::fixed << std::setprecision(4) << effective_diameter << "\n";
    std::cout << "Time for effective diameter: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(effective_diameter_end - effective_diameter_start).count() << " ms\n";

    // Fim do timer geral
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "\nTotal execution time: " << total_duration << " ms\n";

    // Criar o JSON com as métricas
    int node_count = nodes.size();
    json output;
    output["graph_metrics"]["nodes"] = nodes.size();
    output["graph_metrics"]["edges"] = edge_count;

    output["graph_metrics"]["largest_wcc"] = {
        {"nodes", wcc_data.largest_wcc_size},
        {"fraction_of_total_nodes", wcc_data.fraction_nodes},
        {"edges", wcc_data.largest_wcc_edges},
        {"fraction_of_total_edges", wcc_data.fraction_edges}
    };

    output["graph_metrics"]["largest_scc"] = {
        {"nodes", scc_data.largest_scc_size},
        {"fraction_of_total_nodes", scc_data.fraction_nodes},
        {"edges", (scc_data.largest_scc_edges/2)},
        {"fraction_of_total_edges", scc_data.fraction_edges}
    };

    output["graph_metrics"]["average_clustering_coefficient"] = clustering_data.average_clustering_coefficient;
    output["graph_metrics"]["triangles"] = triangle_data.total_triangles;
    output["graph_metrics"]["fraction_of_closed_triangles"] = fraction_closed_triangles;
    output["graph_metrics"]["diameter"] = graph_diameter;
    output["graph_metrics"]["effective_diameter_90_percentile"] = effective_diameter;

    output["execution_time_ms"] = total_duration;

    // Coletar métricas adicionais
    collectExecutionMetrics(output);

    // Salvar o JSON no arquivo
    std::ofstream file("dados_exportados.json");
    file << output.dump(4); // Indentação para facilitar a leitura
    file.close();

    std::cout << "Metrics exported to 'dados_exportados.json'.\n";
    return 0;
}
