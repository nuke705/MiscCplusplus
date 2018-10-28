/* Zhi Ji, IE 523 hw 1, N_queens problem all solution*/
#ifndef N_queens
#define N_queens
using namespace std;
class Board {
	int size;											//size of chessboard
	int count = 0;										//count of solution
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

	//print the board for a solution
	void print_board() {
		std::cout << size << " -Quenn Problem Solution " << count << std::endl;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
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

	//solve for all solution recursively
	bool solve(int col) {
		int n = size;
		if (col >= n && n>1) {
			count++;					//find a solution, add a count
			print_board();				//print every solution
			//return true;
		}
		else {
			for (int i = 0; i<n; i++) {
				if (is_this_position_safe(i, col)) {
					chessboard[i][col] = 1;
					if (solve(col + 1)) {
						return true;
					}
					else {
						chessboard[i][col] = 0;
					}
				}
			}
		}
		return false;
	}

// method for main program
public:
	void nQueens(int n) {
		initialize(n);
		solve(0);            //find all solution
		if (count > 0) {     //print the summary of solutoins when there exists solutions
			std::cout << "There are " << count << " solutions to " << n << "-Queen Problem" << std::endl;
		}
		else {               //print this when there is no solution
			std::cout << "There is no solution to the " << n << "-Queen Problem" << std::endl;
		}
	}

};

#endif 