/*
*
mpiexec - n 4 "C:\Users\mingf\Desktop\MPIBellmanFord\Debug\MPIBellmanFord.exe" "C:\Users\mingf\Desktop\MPIBellmanFord\MPIBellmanFord\input.txt"
*
*/

#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstring>
#include <vector>
#include <chrono>
#include <mpi.h>
#include <ctime> 

using std::string;
using std::cout;
using std::endl;
using std::cin;
#define INF 1000000

//define color
#define YELLOW "\033[93m"
#define AQUA "\033[96m"
#define CYAN "\033[36m"
#define MAGENTA "\033[35m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define GREY "\033[90m"
#define RESET "\033[0m"

string locationName[] = {
	"Main Entrance",
	"YumYum Cafeteria/Block L",
	"The Rimba",
	"Block M",
	"Bangungan KKB/Block A",
	"RedBricks Cafeteria",
	"Bangungan TSS",
	"CITC",
	"Block K",
	"Block D",
	"Library",
	"DTAR",
	"Sport Complex",
	"Hostel",
	"Casuarina Cafe",
	"Block DK",
	"Block AB"
};


/**
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and matrix dimension convert(2D->1D) function
 */


 /**
 *A structure to represent a path from the starting vertex to a destination vertex
 */
struct Path {
	std::vector<int> vertices; // Stores the vertices in the path
	int total_distance;        // Stores the total distance of the path

	Path() : total_distance(INF) {}
};

/**
 * A helper function to print a path
 */
int printPath(const Path& path, const string locationName[]) {
	int pathTime;
	for (size_t i = 0; i < path.vertices.size(); i++) {
		int vertex = path.vertices[i];
		std::string location = locationName[vertex];

		cout << location;

		if (i < path.vertices.size() - 1) {
			cout << " -> ";
		}
	}

	time_t currentTm = time(nullptr);
	struct tm localTime;

	// Use localtime_s to get the local time safely
	if (localtime_s(&localTime, &currentTm) != 0) {
		// Handle the error, e.g., print an error message or exit
		cout << RED;
		std::cerr << "Error getting local time" << std::endl;
		cout << RESET;
		return 1;
	}

	// Check for nullptr, although it's unlikely to be nullptr with this code

	char currentTimeStr[10];  // HH:MM format
	strftime(currentTimeStr, sizeof(currentTimeStr), "%I:%M %p", &localTime);

	// Output the formatted time
	//cout << GREEN << "\nCurrent time: " << currentTimeStr << RESET << endl;
	// Create a copy of the current time's tm structure
	struct tm estimatedTime = localTime;

	// Calculate estimated arrival time (current time + estimation time)
	estimatedTime.tm_min += path.total_distance;

	mktime(&estimatedTime); // Normalize the time structure after adding minutes


	// Format and output the estimated time as "hh:mm am/pm"
	strftime(currentTimeStr, sizeof(currentTimeStr), "%I:%M %p", &estimatedTime);
	cout << GREY << "Estimation Time (mins): " << MAGENTA << path.total_distance << RESET << endl;
	cout << GREY << "Estimated Arrival Time: " << MAGENTA << currentTimeStr << RESET << endl;

	pathTime = path.total_distance;

	return pathTime;
}



/**
 * A function to find and store all paths from the starting vertex to the destination vertex
 */
void findAllPaths(int current_vertex, int destination_vertex, std::vector<int>& current_path, std::vector<Path>& all_paths, int current_distance, int* mat, int n) {
	current_path.push_back(current_vertex);

	if (current_vertex == destination_vertex) {
		// Found a path
		Path path;
		path.vertices = current_path;
		path.total_distance = current_distance;
		all_paths.push_back(path);
	}
	else {
		for (int v = 0; v < n; v++) {
			int weight = mat[current_vertex * n + v];
			if (weight < INF && std::find(current_path.begin(), current_path.end(), v) == current_path.end()) {
				// Vertex v is a neighbor of the current_vertex and has not been visited yet
				findAllPaths(v, destination_vertex, current_path, all_paths, current_distance + weight, mat, n);
			}
		}
	}

	current_path.pop_back();
}

/**
 * A function to print all paths from the starting vertex to the destination vertex
 */
void printAllPaths(int starting_location, int destination_location, const std::vector<Path>& all_paths, const string locationName[], int bestTime) {
	int count = 1;
	Path bestpath;
	int pathTime;
	int tempPathTime = 10000000;
	int bestPathTime;
	if (all_paths.empty()) {
		cout << RED << "No paths found." << RESET << endl;
	}
	else {
		cout << AQUA << "\nAll Possible Paths from " << locationName[starting_location] << " to " << locationName[destination_location] << ": " << RESET << endl;
		for (const Path& path : all_paths) {
			cout << YELLOW << "Path " << count << ":" << RESET << endl;
			pathTime = printPath(path, locationName); // Pass locationName to printPath
			if (pathTime < tempPathTime) {
				tempPathTime = pathTime;
				bestpath = path;
			}

			cout << "\n";

			count++;
			bestPathTime = tempPathTime;
		}

		if (bestPathTime == bestTime) {
			cout << CYAN << "Shortest Path from " << locationName[starting_location] << " to " << locationName[destination_location] << "(Suggests):" << RESET << endl;
			printPath(bestpath, locationName);
		}

	}


}

void printUpdateMessage(int my_rank, int v, int new_dis) {
	cout << GREEN << "Process " << my_rank << RESET << ": Updating dist[" << v << "] = " << new_dis << endl;
}

void printCurrentTime() {
	// Get the current time
	time_t currentTime = time(nullptr);
	struct tm localTime;
	if (localtime_s(&localTime, &currentTime) != 0) {
		std::cerr << "Error getting local time." << endl;
		return;
	}

	char currentTimeStr[10];
	strftime(currentTimeStr, sizeof(currentTimeStr), "%I:%M %p", &localTime);

	cout << GREEN << "\nCurrent time: " << currentTimeStr << "\n" << RESET << endl;
}

/**
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and matrix dimension convert(2D->1D) function
 */
namespace utils {
	int N; //number of vertices
	int* mat; // the adjacency matrix

	void abort_with_error_message(string msg) {
		std::cerr << msg << endl;
		abort();
	}

	//translate 2-dimension coordinate to 1-dimension
	int convert_dimension_2D_1D(int x, int y, int n) {
		return x * n + y;
	}

	int read_file(string filename) {
		std::ifstream inputf(filename, std::ifstream::in);
		if (!inputf.good()) {
			abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
		}
		inputf >> N;
		//input matrix should be smaller than 20MB * 20MB (400MB, we don't have too much memory for multi-processors)
		assert(N < (1024 * 1024 * 20));
		mat = (int*)malloc(N * N * sizeof(int));
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) {
				inputf >> mat[convert_dimension_2D_1D(i, j, N)];
			}
		return 0;
	}

	int print_result(bool has_negative_cycle, int* dist) {
		std::ofstream outputf("output.txt", std::ofstream::out);
		if (!has_negative_cycle) {
			for (int i = 0; i < N; i++) {
				if (dist[i] > INF)
					dist[i] = INF;
				outputf << dist[i] << '\n';
			}
			outputf.flush();
		}
		else {
			outputf << "FOUND NEGATIVE CYCLE!" << endl;
		}
		outputf.close();
		return 0;
	}

	void print_result_to_console(bool has_negative_cycle, int* dist) {
		if (!has_negative_cycle) {
			for (int i = 0; i < N; i++) {
				if (dist[i] > INF)
					dist[i] = INF;
				//std::cout << "Vertex " << i << ": " << dist[i] << '\n';
			}
			std::cout << std::flush;
		}
		else {
			std::cout << "FOUND NEGATIVE CYCLE!" << std::endl;
		}
	}
}//namespace utils


void bellman_ford(int my_rank, int p, MPI_Comm comm, int n, int* mat, int* dist, bool* has_negative_cycle, int starting_vertex, int destination_vertex) {
	int loc_n; // need a local copy for N
	int loc_start, loc_end;
	int* loc_mat; //local matrix
	int* loc_dist; //local distance

	//step 1: broadcast N
	if (my_rank == 0) {
		loc_n = n;
	}
	MPI_Bcast(&loc_n, 1, MPI_INT, 0, comm);

	//step 2: find local task range
	int ave = loc_n / p;
	loc_start = ave * my_rank;
	loc_end = ave * (my_rank + 1);
	if (my_rank == p - 1) {
		loc_end = loc_n;
	}

	//step 3: allocate local memory
	loc_mat = (int*)malloc(loc_n * loc_n * sizeof(int));
	loc_dist = (int*)malloc(loc_n * sizeof(int));

	//step 4: broadcast matrix mat
	if (my_rank == 0)
		memcpy(loc_mat, mat, sizeof(int) * loc_n * loc_n);
	MPI_Bcast(loc_mat, loc_n * loc_n, MPI_INT, 0, comm);

	//step 5: bellman-ford algorithm
	for (int i = 0; i < loc_n; i++) {
		loc_dist[i] = INF;
	}

	loc_dist[starting_vertex] = 0;  // Set the starting vertex's distance to 0

	MPI_Barrier(comm);

	bool loc_has_change;
	int loc_iter_num = 0;
	for (int iter = 0; iter < loc_n - 1; iter++) {
		loc_has_change = false;
		loc_iter_num++;
		for (int u = loc_start; u < loc_end; u++) {
			for (int v = 0; v < loc_n; v++) {
				int weight = loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
				if (weight < INF) {
					if (loc_dist[u] + weight < loc_dist[v]) {
						loc_dist[v] = loc_dist[u] + weight;
						loc_has_change = true;
						printUpdateMessage(my_rank, v, loc_dist[v]);
					}
				}
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, &loc_has_change, 1, MPI_C_BOOL, MPI_LOR, comm);
		if (!loc_has_change)
			break;
		MPI_Allreduce(MPI_IN_PLACE, loc_dist, loc_n, MPI_INT, MPI_MIN, comm);
	}

	//check negative cycle
	if (loc_iter_num == loc_n - 1 || loc_iter_num == destination_vertex) {
		loc_has_change = false;
		for (int u = loc_start; u < loc_end; u++) {
			for (int v = 0; v < loc_n; v++) {
				int weight = loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
				if (weight < INF) {
					if (loc_dist[u] + weight < loc_dist[v]) {
						loc_dist[v] = loc_dist[u] + weight;
						loc_has_change = true;
						break;
					}
				}
			}
		}
		MPI_Allreduce(&loc_has_change, has_negative_cycle, 1, MPI_C_BOOL, MPI_LOR, comm);
	}

	//step 6: retrieve results back
	if (my_rank == 0)
		memcpy(dist, loc_dist, loc_n * sizeof(int));

	//step 7: remember to free memory
	free(loc_mat);
	free(loc_dist);

}
int main(int argc, char** argv) {
	if (argc <= 1) {
		utils::abort_with_error_message("INPUT FILE WAS NOT FOUND!");
	}
	string filename = argv[1];

	int* dist = nullptr; // Initialize dist pointer to nullptr
	bool has_negative_cycle = false;
	int initial_node = -1;
	int final_node = -1;

	// MPI initialization
	MPI_Init(&argc, &argv);
	MPI_Comm comm;

	int p; // number of processors
	int my_rank; // my global rank
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &my_rank);

	// Only rank 0 process does the I/O
	if (my_rank == 0) {
		assert(utils::read_file(filename) == 0);
		dist = (int*)malloc(sizeof(int) * utils::N);


		cout << AQUA;

		cout << "=======================================================\n";
		cout << "[ Bellman-Ford Algorithm for Campus Navigation System ]\n";
		cout << "=======================================================\n\n";

		cout << CYAN;

		cout << "Serial Method\n";
		cout << "............................\n\n";

		cout << YELLOW;

		cout << "No.\tLocation Name\n";
		cout << "---------------------------------\n";

		cout << RESET;

		for (int i = 0; i < sizeof(locationName) / sizeof(locationName[0]); i++) {
			cout << MAGENTA << i + 1 << RESET << "\t" << locationName[i] << endl;
		}

		cout << YELLOW;
		cout << "--------------------------------\n\n";

		cout << RESET;

		// Get user input only from rank 0
		cout << "Enter the starting point (1-17): ";
		cin >> initial_node;
		cout << "Enter the destination point (1-17): ";
		cin >> final_node;

		// Validate input and handle any errors
		while (initial_node < 1 || initial_node > 17 || final_node < 1 || final_node > 17) {
			cout << RED;
			cout << "Invalid start or end point. Please enter valid points.\n" << endl;
			cout << RESET;
			cout << "Enter the starting point (1-17): ";
			cin >> initial_node;
			cout << "Enter the destination point (1-17): ";
			cin >> final_node;
		}

		initial_node--;
		final_node--;
		if (initial_node < 0 || final_node < 0) {
			utils::abort_with_error_message("starting_vertex or destination_vertex WAS NOT FOUND!");
		}

		cout << CYAN;
		cout << "\n\n\nParallel Processing...\n\n";
		cout << RESET;
	}

	// Broadcast input values to all processes
	MPI_Bcast(&initial_node, 1, MPI_INT, 0, comm);
	MPI_Bcast(&final_node, 1, MPI_INT, 0, comm);


	// Time counter
	double t1, t2;
	MPI_Barrier(comm);
	t1 = MPI_Wtime();

	// Bellman-Ford algorithm
	bellman_ford(my_rank, p, comm, utils::N, utils::mat, dist, &has_negative_cycle, initial_node, final_node);
	MPI_Barrier(comm);

	// End timer
	t2 = MPI_Wtime();

	if (my_rank == 0) {

		cout << MAGENTA;
		cout << "\n" << std::setprecision(6) << "Processing Time(s): " << (t2 - t1) << "\n" << endl;
		cout << RESET;

		printCurrentTime();


		// Find all paths from initial_node to final_node
		std::vector<int> current_path;
		std::vector<Path> all_paths;
		findAllPaths(initial_node, final_node, current_path, all_paths, 0, utils::mat, utils::N);
		int bestTime = dist[final_node];
		// Print all paths
		printAllPaths(initial_node, final_node, all_paths, locationName, bestTime);
		//std::cerr.setf(std::ios::fixed);
		/*std::cerr << std::setprecision(6) << "Time(s): " << (t2 - t1) << endl;*/
		//utils::print_result(has_negative_cycle, dist);
		//utils::print_result_to_console(has_negative_cycle, dist);


		cout << AQUA << "\n========================== END ========================\n" << RESET << endl;

		free(dist);
		free(utils::mat);
	}

	MPI_Finalize();
	return 0;
}
