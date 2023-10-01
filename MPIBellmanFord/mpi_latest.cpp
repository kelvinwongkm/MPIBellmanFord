/*
*
mpiexec -n 4 "C:\Users\mingf\Desktop\MPIBellmanFord\Debug\MPIBellmanFord.exe" "C:\Users\mingf\Desktop\MPIBellmanFord\MPIBellmanFord\input.txt"
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
#include <vector>
#include <chrono>
#include <windows.h>
#include <psapi.h>

using namespace std;

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

struct Path {
	std::vector<int> vertices; // Stores the vertices in the path
	int total_distance;        // Stores the total distance of the path

	Path() : total_distance(INF) {}
};

// to print path and estimate arrival time
int printPath(const Path& path, const string locationName[]) {
	int pathTime;
	// iterate throught the vertices in the path
	for (size_t i = 0; i < path.vertices.size(); i++) {
		if (i == 0) {
			cout << "\t\t>   ";
		}
		else if (i == 3 || i == 6 || i == 9 || i == 12 || i == 15) {
			cout << "\n";
			cout << "\t\t>   ";
		}
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
		std::cerr << "Error getting local time" << endl;
		cout << RESET;
		return 1;
	}

	char currentTimeStr[10];  // HH:MM format
	strftime(currentTimeStr, sizeof(currentTimeStr), "%I:%M %p", &localTime);

	struct tm estimatedTime = localTime;

	// Calculate estimated arrival time (current time + estimation time)
	estimatedTime.tm_min += path.total_distance;

	mktime(&estimatedTime); // Normalize the time structure after adding minutes


	// Format and output the estimated time as "hh:mm am/pm"
	strftime(currentTimeStr, sizeof(currentTimeStr), "%I:%M %p", &estimatedTime);
	cout << GREY << "\n\t\t>   Estimation Time (mins): " << MAGENTA << path.total_distance << RESET << endl;
	cout << GREY << "\t\t>   Estimated Arrival Time: " << MAGENTA << currentTimeStr << RESET << endl;
	cout << "\t\t+=================================================================================\n";

	pathTime = path.total_distance;

	return pathTime;

}

// Recursive function to find all paths between two vertices
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

void printAllPaths(int starting_location, int destination_location, const std::vector<Path>& all_paths, const string locationName[], int bestTime) {
	int count = 1;
	Path bestpath;
	int pathTime;
	int tempPathTime = 10000000;
	int bestPathTime{};
	if (all_paths.empty()) {
		cout << RED << "No paths found." << RESET << endl;
	}
	else {
		cout << RESET << "\n\t\t+=================================================================================" << endl;
		cout << RESET << "\t\t" << AQUA << "    All possible paths from " << locationName[starting_location] << " to " << locationName[destination_location] << ":" << RESET << endl;
		cout << "\t\t+=================================================================================\n";
		for (const Path& path : all_paths) {
			cout << YELLOW << "\t\t>   Path " << count << ":" << RESET << endl;
			pathTime = printPath(path, locationName); // Pass locationName to printPath
			if (pathTime < tempPathTime) {
				tempPathTime = pathTime;
				bestpath = path;
			}

			//cout << "\n";
			count++;
			bestPathTime = tempPathTime;
		}

		if (bestPathTime == bestTime) {
			cout << "\n\n\t\t+=================================================================================\n";
			cout << RESET << "\t\t" << GREEN << "    Shortest Path from " << locationName[starting_location] << " to " << locationName[destination_location] << "(Suggests):" << RESET << endl;
			cout << "\t\t+=================================================================================\n";
			printPath(bestpath, locationName);
		}

	}
}

// print update function to shows the algo working
void printUpdateMessage(int my_rank, int v, int new_dis) {
	cout << GREEN << "\t\tProcess " << my_rank << RESET << ": Updating dist[" << v << "] = " << new_dis << endl;
	//printf("\t\t\t\t%sUpdating dist[%d] = %d%s\n", GREEN, v, new_dist, RESET);
}

// to output current time
void printCurrentTime() {
	// Get the current time
	time_t currentTime = time(nullptr);
	struct tm localTime;
	if (localtime_s(&localTime, &currentTime) != 0) {
		std::cerr << "\t\tError getting local time." << endl;
		return;
	}

	char currentTimeStr[10];
	strftime(currentTimeStr, sizeof(currentTimeStr), "%I:%M %p", &localTime);

	cout << "\t\t" << GREEN << "+---------------------------+";
	cout << "\n\t\t" << GREEN << "|  Current time: " << currentTimeStr << "   |" << RESET << endl;
	cout << "\t\t" << GREEN << "+---------------------------+\n";
}

/// <summary>
/// for performance 
/// </summary>
/// to get memory usage
double GetMemoryUsageMB() {
	PROCESS_MEMORY_COUNTERS pmc;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(PROCESS_MEMORY_COUNTERS))) {
		return static_cast<double>(pmc.WorkingSetSize) / (1024 * 1024);
	}
	else {
		std::cerr << "Error getting process memory info: " << GetLastError() << std::endl;
		return -1.0; // Return a negative value to indicate an error
	}
}

void getStartPointEndPoint(int& start_point, int& end_point) {
	cout << AQUA;

	cout << "\t\t\t   +===========================================================+\n";
	cout << "\t\t\t   ||                                                         ||\n";
	cout << "\t\t\t   || [ Bellman-Ford Algorithm for Campus Navigation System ] ||\n";
	cout << "\t\t\t   ||                                                         ||\n";
	cout << "\t\t\t   +===========================================================+\n\n";
	cout << CYAN;
	cout << "\t\t\t\t\t+..............................+\n";
	cout << "\t\t\t\t\t|  Parallel Method Used: CUDA  |\n";
	cout << "\t\t\t\t\t+..............................+\n\n";
	cout << YELLOW;
	cout << "\t\t\t\t+-----+------------------------------------------+\n";
	cout << "\t\t\t\t| No  |	Location Name                            |\n";
	cout << "\t\t\t\t+-----+------------------------------------------+\n";
	cout << RESET;

	for (int i = 0; i < sizeof(locationName) / sizeof(locationName[0]); i++) {
		if (i < 9) {
			std::cout << MAGENTA << "\t\t\t\t|  " << i + 1 << "  | " << RESET << locationName[i] << std::endl;
			cout << YELLOW;
			cout << "\t\t\t\t+-----+------------------------------------------+\n";
		}
		else {
			std::cout << MAGENTA << "\t\t\t\t| " << i + 1 << "  | " << RESET << locationName[i] << std::endl;
			cout << YELLOW;
			cout << "\t\t\t\t+-----+------------------------------------------+\n";
			cout << RESET;
		}
	}

	cout << RESET;

	// Get user input only from rank 0
	cout << "\n\t\t\t\tEnter the starting point (1-17): ";
	cin >> start_point;
	cout << "\n\t\t\t\tEnter the destination point (1-17): ";
	cin >> end_point;

	// Validate input and handle any errors
	while (start_point < 1 || start_point > 17 || end_point < 1 || end_point > 17) {
		cout << RED;
		cout << "\n\t\t\t\tInvalid start or end point. Please enter valid points.\n" << endl;
		cout << RESET;
		cout << "\n\t\t\t\tEnter the starting point (1-17): ";
		cin >> start_point;
		cout << "\n\t\t\t\tEnter the destination point (1-17): ";
		cin >> end_point;
	}
}

template <typename T>
// to output the variable value in each rank for debug purpose
void debugoutput(const std::string& variableName, const T& value, int rank) {
	//std::cout << variableName << " : " << value << " Rank: " << rank << std::endl;
	cout << "Rank: " << rank << " " << variableName << " : " << value << endl;
}

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
			outputf << "\t\tFOUND NEGATIVE CYCLE!" << endl;
		}
		outputf.close();
		return 0;
	}

	void print_serial_process(bool has_negative_cycle, int* dist) {
		if (!has_negative_cycle) {
			for (int i = 0; i < N; i++) {
				if (dist[i] > INF)
					dist[i] = INF;
				cout << GREEN << "Process " << i << RESET << ": Updating dist[" << i << "] = " << dist[i] << endl;
			}
		}
		else {
			cout << RED << "FOUND NEGATIVE CYCLE!" << RESET << endl;
		}
	}


	void print_result_to_console(bool has_negative_cycle, int* dist) {
		if (!has_negative_cycle) {
			for (int i = 0; i < N; i++) {
				if (dist[i] > INF)
					dist[i] = INF;
				std::cout << "Vertex " << i << ": " << dist[i] << '\n';
			}
			cout << std::flush;
		}
		else {
			cout << "FOUND NEGATIVE CYCLE!" << endl;
		}
	}

}

void bellman_ford_serial(int n, int* mat, int* dist, bool* has_negative_cycle, int start_point, int end_point, std::vector<int>& prev) {

	// to check if negative cycle exists in the graph
	*has_negative_cycle = false;

	for (int i = 0; i < n; i++) {
		// initialize the distance to all vertices as infinity
		dist[i] = INF;
		// Initialize previous node to -1 for all locations
		prev[i] = -1; 
	}

	// set distance to ssource point to 0
	dist[start_point] = 0;

	bool has_change;

	// worst case scenario is n-1 to get the shortest path
	for (int i = 0; i < n - 1; i++) {

		// to keep track any update on the distence in current iteration
		has_change = false;

		// to loop through all pair of vertices (u,v)
		for (int u = 0; u < n; u++) {
			for (int v = 0; v < n; v++) {
				// get the weight for the pair
				int weight = mat[utils::convert_dimension_2D_1D(u, v, n)];
				if (weight < INF && weight >0) {
					// check if distance to v through u is shorter than 
					// current distance to v, if yes then update
					if (dist[u] + weight < dist[v]) {
						has_change = true;
						dist[v] = dist[u] + weight;
						prev[v] = u; // Update previous node for location v
						std::cout << GREEN << "\t\t\t\tUpdating dist[" << v << "] = " << dist[v] << std::endl;
					}
				}
			}
		}
		// if no update, break the loop early
		// shortest path have been found
		if (!has_change) {
			return;
		}
	}

	// Check negative cycle, if negative cycle, terminate the program
	for (int u = 0; u < n; u++) {
		for (int v = 0; v < n; v++) {
			int weight = mat[utils::convert_dimension_2D_1D(u, v, n)];
			if (weight < INF) {
				if (dist[u] + weight < dist[v]) {
					*has_negative_cycle = true;
					return;
				}
			}
		}
	}
}

void bellman_ford_parallel(int rank, int size, int n, int* mat, int* dist, bool* has_negative_cycle, int start_point, int end_point) {

	int loc_n; // need a local copy for N
	int loc_start, loc_end;
	int* loc_mat; //local matrix
	int* loc_dist; //local distance

	//step 1: broadcast N for all processes to know total number of vertices
	if (rank == 0) {
		loc_n = n;
	}
	MPI_Bcast(&loc_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//step 2: find local task range, each rank is handling own dataset
	int ave = loc_n / size;
	loc_start = ave * rank;
	loc_end = ave * (rank + 1);
	if (rank == size - 1) {
		loc_end = loc_n;
	}

	//step 3: memory allocation for local data
	loc_mat = (int*)malloc(loc_n * loc_n * sizeof(int));
	loc_dist = (int*)malloc(loc_n * sizeof(int));
	int* loc_pred = (int*)malloc(loc_n * sizeof(int));

	//step 4: broadcast matrix mat, show same copy of graph structure
	if (rank == 0)
		memcpy(loc_mat, mat, sizeof(int) * loc_n * loc_n);
	MPI_Bcast(loc_mat, loc_n * loc_n, MPI_INT, 0, MPI_COMM_WORLD);

	//step 5: bellman-ford algorithm
	for (int i = 0; i < loc_n; i++) {
		// initialize the distance to all vertices as infinity
		loc_dist[i] = INF;
		// Initialize predecessor node to -1 for all locations
		loc_pred[i] = -1; 
	}

	loc_dist[start_point] = 0;  // Set the starting vertex's distance to 0
	MPI_Barrier(MPI_COMM_WORLD);

	bool loc_has_change;
	int loc_iter_num = 0;

	// worst case scenario is n-1 to get the shortest path
	for (int iter = 0; iter < loc_n - 1; iter++) {
		loc_has_change = false;
		loc_iter_num++;
		for (int u = loc_start; u < loc_end; u++) {
			for (int v = 0; v < loc_n; v++) {
				int weight = loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
				if (weight < INF && weight > 0) {
					if (loc_dist[u] + weight < loc_dist[v] ) {
						loc_has_change = true;
						loc_dist[v] = loc_dist[u] + weight;
						loc_pred[v] = u;  
						std::cout << GREEN << "\t\t\t\tUpdating dist[" << v << "] = " << loc_pred[v] << std::endl;
					}
				}
			}
		}

		// check any of the proces has change the flag value, if all no, means shortest path been found
		MPI_Allreduce(MPI_IN_PLACE, &loc_has_change, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
		if (!loc_has_change)
			break;
		// find the minimum distance among all process and update loc_dist
		MPI_Allreduce(MPI_IN_PLACE, loc_dist, loc_n, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	}

	//check negative cycle
	if (loc_iter_num == loc_n - 1 || loc_iter_num == end_point) {
		loc_has_change = false;
		for (int u = loc_start; u < loc_end; u++) {
			for (int v = 0; v < loc_n; v++) {
				int weight = loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
				if (weight < INF && weight > 0) {
					if (loc_dist[u] + weight < loc_dist[v]) {
						loc_dist[v] = loc_dist[u] + weight;
						loc_has_change = true;
						break;
					}
				}
			}
		}
		MPI_Allreduce(&loc_has_change, has_negative_cycle, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
	}

	//step 6: retrieve results back
	if (rank == 0)
		memcpy(dist, loc_dist, loc_n * sizeof(int));

	//step 7: remember to free memory
	free(loc_mat);
	free(loc_dist);
}

int main(int argc, char** argv) {
	string filename = "C:\\Users\\mingf\\Desktop\\dspcAssignment\\dspcAssignment\\dataset3.txt";
	
	// MPI initialization
	MPI_Init(&argc, &argv);

	// declare all require variable
	int* serial_dist = nullptr; // Initialize dist pointer to nullptr
	int* parallel_dist = nullptr; // Initialize dist pointer to nullptr
	bool has_negative_cycle = false;
	int start_point;
	int end_point;
	std::vector<int> prev;


	// ========== for performance comparison =======================
	// Get a handle to the current process
	HANDLE hProcess = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	memset(&pmc, 0, sizeof(PROCESS_MEMORY_COUNTERS));
	// Store memory usage for all processes
	double serial_memory_usages_before;
	double serial_memory_usages_after;
	double parallel_memory_usages_before;
	double parallel_memory_usages_after;

	// Time counter
	double serial_execution_time = 0; // Store serial execution time
	double parallel_execution_time = 0; // Store parallel execution time

	// MPI
	int size;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get process id
	MPI_Comm_size(MPI_COMM_WORLD, &size);	// get num of processor

	// Only rank 0 process perform the task
	if (rank == 0) {

		start_point = 17;
		end_point = 1;
		assert(utils::read_file(filename) == 0);
		serial_dist = (int*)malloc(sizeof(int) * utils::N);
		parallel_dist = (int*)malloc(sizeof(int) * utils::N);
		getStartPointEndPoint(start_point, end_point);

		start_point--;
		end_point--;
		if (start_point < 0 || end_point < 0) {
			utils::abort_with_error_message("starting_vertex or destination_vertex WAS NOT FOUND!");
		}


		// Measure memory usage before running the serial algorithm
		serial_memory_usages_before = GetMemoryUsageMB();

		cout << CYAN << "\n\n\n\t\t\t\tSerial Processing...\n\n";

		std::vector<int> prev(utils::N);
		auto t1 = std::chrono::steady_clock::now();
		bellman_ford_serial(utils::N, utils::mat, serial_dist, &has_negative_cycle, start_point, end_point, prev);
		auto t2 = std::chrono::steady_clock::now();
		// getting the time execution
		serial_execution_time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1e3;

		// Measure memory usage after running the serial algorithm
		serial_memory_usages_after = GetMemoryUsageMB();

		cout << YELLOW;
		cout << "\t\t\t\t+------------------------------------------------+\n";
		cout << RESET;
		std::cerr << AQUA << std::fixed << "\t\t\t\t|   Serial Execution Time   : " << RESET << serial_execution_time << " sec       |" << endl;
		cout << "\t\t\t\t+------------------------------------------------+\n";


		cout << CYAN << "\n\n\n\t\t\t\tParallel Processing...\n\n";
	}

	//// Broadcast input values to all processes
	MPI_Bcast(&start_point, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&end_point, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//// Time counter
	double t1, t2;

	// Measure memory usage before running the serial algorithm
	parallel_memory_usages_before = GetMemoryUsageMB();

	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	bellman_ford_parallel(rank, size,  utils::N, utils::mat, parallel_dist, &has_negative_cycle, start_point, end_point);

	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();

	// Measure memory usage after running the serial algorithm
	parallel_memory_usages_after = GetMemoryUsageMB();

	if (rank == 0) {

		cout << CYAN << "\n\n\t\t+=================================================================================+";
		cout << CYAN << "\n\t\t|                                   RESULTS                                       |\n";
		cout << "\t\t+=================================================================================+\n";
		printCurrentTime();

		// Find all paths from initial_node to final_node, and print the result
		std::vector<int> current_path;
		std::vector<Path> all_paths;
		findAllPaths(start_point, end_point, current_path, all_paths, 0, utils::mat, utils::N);
		int bestTime = parallel_dist[end_point];
		printAllPaths(start_point, end_point, all_paths, locationName, bestTime);

		// ==================  PERFORMANCE EVALUATION SECTION ===========================
		cout << YELLOW << "\n\n\t\t+=================================================================================+\n";
		cout << YELLOW << "\t\t\t\t            PERFORMANCE EVALUATION            \n";
		cout << "\t\t+=================================================================================+\n";
		cout << "\n\t\t" << CYAN << "+===================================SPEED UP======================================+\n";


		double parallel_execution_time = t2 - t1;
		double speedup = serial_execution_time - parallel_execution_time;

		// display the execution time and speed up
		cout << AQUA << std::setprecision(6) << "\t\tSerial Execution Time 	: " << RESET << serial_execution_time << " seconds" << endl;
		cout << AQUA << std::setprecision(6) << "\t\tParallel Execution Time	: " << RESET << parallel_execution_time << " seconds" << endl;
		cout << AQUA << "\t\tSpeedup			: " << RESET << speedup << endl;
		cout << endl;

		cout << "\n\t\t" << CYAN << "+===================================MEMORY USAGE==================================+\n";
		// Display memory usage
		cout << AQUA << "\t\tMU Before Serial Execution		: " << RESET << std::fixed << std::setprecision(4) << serial_memory_usages_before << " MB" << endl;
		cout << AQUA << "\t\tMU After Serial Execution		: " << RESET << std::fixed << std::setprecision(4) << serial_memory_usages_after << " MB" << endl;
		cout << AQUA << "\t\tMU Before Parallel Execution		: " << RESET << std::fixed << std::setprecision(4) << parallel_memory_usages_before << " MB" << endl;
		cout << AQUA << "\t\tMU After Parallel Execution		: " << RESET << std::fixed << std::setprecision(4) << parallel_memory_usages_after << " MB" << endl;

		// Calculate memory usage difference for serial and parallel execution
		double serial_memory_usage_diff = serial_memory_usages_after - serial_memory_usages_before;
		double parallel_memory_usage_diff = parallel_memory_usages_after - parallel_memory_usages_before;

		cout << AQUA << "\n\t\tMU Difference for Serial Execution	: " << RESET << std::fixed << std::setprecision(4) << serial_memory_usage_diff << " MB" << endl;
		cout << AQUA << "\t\tMU Difference for Parallel Execution	: " << RESET << std::fixed << std::setprecision(4) << parallel_memory_usage_diff << " MB" << endl;

		cout << YELLOW << "\t\t===================================================================================" << RESET << endl;


		cout << AQUA << "\n\t\t+===================================== END =======================================+\n" << RESET << endl;

		// free all the variable or pointer
		free(parallel_dist);
		free(serial_dist);
		free(utils::mat);
	}

	// finalize or end the mpi environment
	MPI_Finalize();
	return 0;
}

