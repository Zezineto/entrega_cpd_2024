#include <pthread.h>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <atomic>
#include <mutex>

std::mutex adj_list_mutex;
std::mutex nodes_mutex;

// Estruturas globais
std::unordered_map<int, std::vector<int>> adj_list;
std::set<int> nodes;
std::atomic<int> edge_count(0);

// Função para carregar parte do grafo (usando pthreads)
void* loadGraphPart(void* arg) {
    int thread_id = *((int*)arg);
    int lines_per_thread = 10000;  // Número de linhas que cada thread vai processar (ajustar conforme necessário)
    int start_line = thread_id * lines_per_thread;
    int end_line = start_line + lines_per_thread;

    std::ifstream file("facebook_combined.txt");
    std::string line;
    int current_line = 0;

    // Encontrar a linha que a thread deve começar a processar
    while (std::getline(file, line)) {
        if (current_line >= start_line && current_line < end_line) {
            std::stringstream ss(line);
            int source, target;

            if (ss >> source >> target) {
                // Lock para acessar a lista de adjacência de forma segura
                std::lock_guard<std::mutex> lock(adj_list_mutex);
                adj_list[source].push_back(target);
                adj_list[target].push_back(source);

                // Lock para acessar o conjunto de nós de forma segura
                std::lock_guard<std::mutex> node_lock(nodes_mutex);
                nodes.insert(source);
                nodes.insert(target);

                // Contar as arestas
                edge_count.fetch_add(1, std::memory_order_relaxed);
            }
        }
        current_line++;
    }

    file.close();
    return nullptr;
}

int main() {
    int num_threads = 4; // Pode ser 4 ou 8
    pthread_t threads[num_threads];
    int thread_ids[num_threads];

    // Criar as threads
    for (int i = 0; i < num_threads; ++i) {
        thread_ids[i] = i;
        pthread_create(&threads[i], nullptr, loadGraphPart, (void*)&thread_ids[i]);
    }

    // Esperar todas as threads terminarem
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], nullptr);
    }

    // Exibir os resultados
    std::cout << "Nodos: " << nodes.size() << std::endl;
    std::cout << "Arestas: " << edge_count.load() << std::endl;

    return 0;
}
