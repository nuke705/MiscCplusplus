/* Zhi Ji, IE 523 hw 1, N_queens problem*/
#ifndef N_queens
#define N_queens
using namespace std;
class Board {
	int size;											//size of chessboard
	int count = 1;										//count of solution
	int **chessboard;									//setup board using p2p
	//check if the position is safe to put a queen
	bool is_this_position_safe(int row, int col) {
		int n = size;
		/*check column*/
		for (int i = 0; i < n; i++) {
			if (chessboard[i][col] == 1) 
				return false;		
		}
		/*check row*/
		for (int i = 0; i < n; i++) {
			if (chessboard[row][i] == 1) 
				return false;
		}
		/* upper diag on left side*/
		for (int i = row, j = col; i >= 0 && j >= 0; i--, j--) {
			if (chessboard[i][j] == 1) 
				return false;
		}
		/* Check lower diagonal on left side */
		for (int i = row, j = col; j >= 0 && i<n; i++, j--) {
			if (chessboard[i][j] == 1) 
				return false;
		}
		//only left side needs to be checked because the right side is
		//empty during the process
		return true;
	}

	//set up an empty chessboard using p2p
	void initialize(int n) {
		size = n;
		chessboard = new int*[n];
		for (int i = 0; i < n; i++) {
			chessboard[i] = new int[n];
			for (int j = 0; j< n; j++) {
				chessboard[i][j] = 0;			//put 0 for any empty position
			}
		}
	}

	//print the board
	void print_board() {
		int n = size;
		std::cout << n << "-Quenn Problem Solution " /*<< count */<< std::endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (chessboard[i][j] == 1) {
					cout << "Q ";
				}
				else {
					cout << "+ ";
				}
			}
			cout << endl;
		}
		cout << endl;
	}
	//solve the puzzle recursively for one solution
	bool solve(int col) {	
		if (col >= size && size > 1) {		
			//print_board();    reserved for all solution
			//count++;
			return true;
		} else {
			for (int i = 0; i< size; i++) {
				if (is_this_position_safe(i, col)) {
					chessboard[i][col] = 1;           //put queen
					if (solve(col + 1)) {			//recursion
						return true;
					}
					else {             
						chessboard[i][col] = 0;   //remove if not safe
					}
				}
			}
		}
		return false;
	}

// public for the main program to use
public:
	void nQueens(int n) {
		initialize(n);
		if (solve(0)) {     //find a solution
			print_board();
		}
		else {				//cannot find a solution given size n
			std::cout << "There is no solution to the " << n << "-Queen Problem" << std::endl;
		}
	}

};

#endif 
