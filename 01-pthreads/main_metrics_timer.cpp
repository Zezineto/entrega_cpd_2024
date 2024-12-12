#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <unordered_map>
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
#include <pthread.h>
#include <cstring>

using json = nlohmann::json;

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

// Função para carregar o grafo do arquivo
bool loadGraph(const std::string &filename, std::unordered_map<int, std::vector<int>> &adj_list, std::set<int> &nodes, int &edge_count) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erro ao abrir o arquivo " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int source, target;

        if (ss >> source >> target) {
            adj_list[source].push_back(target); // Aresta direta
            adj_list[target].push_back(source); // Aresta reversa (ignorar direção)
            nodes.insert(source);
            nodes.insert(target);
            edge_count++;
        } else {
            std::cerr << "Erro ao processar a linha: " << line << std::endl;
        }
    }
    file.close();
    return true;
}

// Função para calcular o maior componente fracamente conectado (WCC)
void computeLargestWCC(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes, int &largest_wcc_size, int &largest_wcc_edges) {
    std::set<int> visited;
    largest_wcc_size = 0;
    largest_wcc_edges = 0;

    for (int node : nodes) {
        if (visited.find(node) == visited.end()) {
            std::set<int> component;
            bfs(adj_list, node, visited, component);

            int component_edges = 0;
            for (int comp_node : component) {
                component_edges += adj_list.at(comp_node).size();
            }
            component_edges /= 2; // Divide por 2 porque as arestas foram contadas duas vezes

            if (component.size() > largest_wcc_size) {
                largest_wcc_size = component.size();
                largest_wcc_edges = component_edges;
            }
        }
    }
}

//new

// Função para calcular o maior SCC
void computeLargestSCC(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes, int &largest_scc_size, int &largest_scc_edges) {
    std::set<int> visited;
    std::vector<int> finish_order;

    // Passo 1: Preencher a ordem de término
    for (int node : nodes) {
        if (visited.find(node) == visited.end()) {
            std::set<int> dummy_component;
            bfs(adj_list, node, visited, dummy_component);
            finish_order.push_back(node);
        }
    }

    // Reverter a ordem de término
    std::reverse(finish_order.begin(), finish_order.end());

    // Passo 2: Criar o grafo reverso
    std::unordered_map<int, std::vector<int>> reversed_graph;
    for (const auto &pair : adj_list) {
        int source = pair.first;
        for (int target : pair.second) {
            reversed_graph[target].push_back(source);
        }
    }

    // Passo 3: Processar os nós no grafo reverso
    visited.clear();
    largest_scc_size = 0;
    largest_scc_edges = 0;

    for (int node : finish_order) {
        if (visited.find(node) == visited.end()) {
            std::set<int> component;
            bfs(reversed_graph, node, visited, component);

            // Contar as arestas no SCC atual
            int component_edges = 0;
            for (int u : component) {
                for (int v : adj_list.at(u)) {
                    if (component.find(v) != component.end()) {
                        component_edges++;
                    }
                }
            }

            // Atualizar se for o maior SCC
            if (component.size() > largest_scc_size) {
                largest_scc_size = component.size();
                largest_scc_edges = component_edges;
            }
        }
    }
}


//new

// Função para calcular as frações de nós e arestas
void computeFractions(int largest_wcc_size, int largest_wcc_edges, int total_nodes, int total_edges, double &fraction_nodes, double &fraction_edges) {
    fraction_nodes = static_cast<double>(largest_wcc_size) / total_nodes;
    fraction_edges = static_cast<double>(largest_wcc_edges) / total_edges;
}

// Função para calcular o coeficiente de agrupamento de um nó
double computeAverageClusteringCoefficient(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes) {
    double total_clustering = 0.0;

    for (int node : nodes) {
        const auto &neighbors = adj_list.at(node);
        int degree = neighbors.size();
        if (degree < 2) {
            total_clustering += 0.0; // Contribuição explícita para manter a média correta
            continue;
        }

        int triangle_count = 0;
        std::set<int> neighbor_set(neighbors.begin(), neighbors.end());

        for (size_t i = 0; i < neighbors.size(); ++i) {
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                if (neighbor_set.find(neighbors[j]) != neighbor_set.end() &&
                    std::find(adj_list.at(neighbors[i]).begin(), adj_list.at(neighbors[i]).end(), neighbors[j]) != adj_list.at(neighbors[i]).end()) {
                    triangle_count++;
                }
            }
        }

        double clustering = static_cast<double>(2 * triangle_count) / (degree * (degree - 1));
        total_clustering += clustering;
    }

    return total_clustering / nodes.size(); // Média sobre todos os nós
}


// Função para contar triângulos no grafo
int countTriangles(const std::unordered_map<int, std::vector<int>> &adj_list, const std::set<int> &nodes) {
    int total_triangles = 0;

    for (int node : nodes) {
        const auto &neighbors = adj_list.at(node);
        std::set<int> neighbor_set(neighbors.begin(), neighbors.end());

        for (size_t i = 0; i < neighbors.size(); ++i) {
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                if (neighbor_set.find(neighbors[j]) != neighbor_set.end() && 
                    std::find(adj_list.at(neighbors[i]).begin(), adj_list.at(neighbors[i]).end(), neighbors[j]) != adj_list.at(neighbors[i]).end()) {
                    total_triangles++;
                }
            }
        }
    }

    // Dividir por 3 porque cada triângulo é contado três vezes
    return total_triangles / 3;
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

int main() {
    std::unordered_map<int, std::vector<int>> adj_list;
    std::set<int> nodes;
    int edge_count = 0;

    // Início do timer para carregar o grafo
    auto start_time = std::chrono::high_resolution_clock::now();

    // Carregar o grafo
    if (!loadGraph("facebook_combined.txt", adj_list, nodes, edge_count)) {
        return 1;
    }


    auto loadGraph_end = std::chrono::high_resolution_clock::now();
    std::cout << "\nTime to load graph: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(loadGraph_end - start_time).count() << " ms\n";

    // NÓS
    auto nodes_start = std::chrono::high_resolution_clock::now();
    std::cout << "\nnodes: " << nodes.size() << "\n";
    auto nodes_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for nodes: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(nodes_end - nodes_start).count() << " ms\n\n";

    // ARESTAS
    auto edges_start = std::chrono::high_resolution_clock::now();
    std::cout << "edges: " << edge_count << "\n";
    auto edges_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for edges: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(edges_end - edges_start).count() << " ms\n\n";

    // WCC 
    auto wcc_start = std::chrono::high_resolution_clock::now();
    int largest_wcc_size = 0;
    int largest_wcc_edges = 0;
    computeLargestWCC(adj_list, nodes, largest_wcc_size, largest_wcc_edges);
    double fraction_nodes = 0.0, fraction_edges = 0.0;
    computeFractions(largest_wcc_size, largest_wcc_edges, nodes.size(), edge_count, fraction_nodes, fraction_edges);
    auto wcc_end = std::chrono::high_resolution_clock::now();
    std::cout << "nodes (WCC): " << largest_wcc_size << "\n";
    std::cout << "fraction of total nodes (WCC): " << std::fixed << std::setprecision(3) << fraction_nodes << "\n";
    std::cout << "edges (WCC): " << largest_wcc_edges << "\n";
    std::cout << "fraction of total edges (WCC): " << std::fixed << std::setprecision(3) << fraction_edges << "\n";
    std::cout << "Time for WCC: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(wcc_end - wcc_start).count() << " ms\n\n";

    // SCC
    auto scc_start = std::chrono::high_resolution_clock::now();
    int largest_scc_size = 0;
    int largest_scc_edges = 0;
    computeLargestSCC(adj_list, nodes, largest_scc_size, largest_scc_edges);
    double fraction_nodes_scc = 0.0, fraction_edges_scc = 0.0;
    computeFractions(largest_scc_size, largest_scc_edges / 2, nodes.size(), edge_count, fraction_nodes_scc, fraction_edges_scc);
    auto scc_end = std::chrono::high_resolution_clock::now();
    std::cout << "Nodes in largest SCC: " << largest_scc_size << "\n";
    std::cout << "fraction of total nodes (SCC): " << std::fixed << std::setprecision(3) << fraction_nodes_scc << "\n";
    std::cout << "Edges in largest SCC: " << (largest_scc_edges / 2) << "\n";
    std::cout << "fraction of total edges (SCC): " << std::fixed << std::setprecision(3) << fraction_edges_scc << "\n";
    std::cout << "Time for SCC: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(scc_end - scc_start).count() << " ms\n\n";

    

    // COEFICIENTE DE AGRUPAMENTO MÉDIO
    auto clustering_start = std::chrono::high_resolution_clock::now();
    double average_clustering_coefficient = computeAverageClusteringCoefficient(adj_list, nodes);
    auto clustering_end = std::chrono::high_resolution_clock::now();
    std::cout << "average clustering coefficient: " << std::fixed << std::setprecision(4) << average_clustering_coefficient << "\n";
    std::cout << "Time for clustering coefficient: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(clustering_end - clustering_start).count() << " ms\n\n";

    // TRIÂNGULOS
    auto triangles_start = std::chrono::high_resolution_clock::now();
    int total_triangles = countTriangles(adj_list, nodes);
    auto triangles_end = std::chrono::high_resolution_clock::now();
    std::cout << "triangles: " << total_triangles << "\n";
    std::cout << "Time for triangles: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(triangles_end - triangles_start).count() << " ms\n\n";

    // FRAÇÃO DE TRIÂNGULOS FECHADOS
    auto closed_triangles_start = std::chrono::high_resolution_clock::now();
    double fraction_closed_triangles = calculateFractionOfClosedTriangles(adj_list, nodes);
    auto closed_triangles_end = std::chrono::high_resolution_clock::now();
    std::cout << "fraction of closed triangles: " << std::fixed << std::setprecision(4) << fraction_closed_triangles << "\n";
    std::cout << "Time for fraction of closed triangles: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(closed_triangles_end - closed_triangles_start).count() << " ms\n\n";

    // DIÂMETRO DO GRAFO
    auto diameter_start = std::chrono::high_resolution_clock::now();
    int graph_diameter = calculateGraphDiameter(adj_list, nodes);
    auto diameter_end = std::chrono::high_resolution_clock::now();
    std::cout << "graph diameter: " << graph_diameter << "\n";
    std::cout << "Time for graph diameter: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(diameter_end - diameter_start).count() << " ms\n\n";

    // DIÂMETRO EFETIVO DO GRAFO
    auto effective_diameter_start = std::chrono::high_resolution_clock::now();
    double effective_diameter = calculateEffectiveDiameter(adj_list, nodes);
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
        {"nodes", largest_wcc_size},
        {"fraction_of_total_nodes", static_cast<double>(largest_wcc_size) / nodes.size()},
        {"edges", largest_wcc_edges},
        {"fraction_of_total_edges", static_cast<double>(largest_wcc_edges) / edge_count}
    };

    output["graph_metrics"]["largest_scc"] = {
        {"nodes", largest_scc_size},
        {"fraction_of_total_nodes", static_cast<double>(largest_scc_size) / nodes.size()},
        {"edges", (largest_scc_edges/2)},
        {"fraction_of_total_edges", (static_cast<double>(largest_scc_edges) / edge_count) / 2}
    };

    output["graph_metrics"]["average_clustering_coefficient"] = average_clustering_coefficient;
    output["graph_metrics"]["triangles"] = total_triangles;
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
