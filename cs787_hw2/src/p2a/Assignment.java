package p2a;

import java.util.Random;

public class Assignment {
	
	public static int randInt(int min, int max) {

		// NOTE: Usually this should be a field rather than a method
		// variable so that it is not re-seeded every call.
		Random rand = new Random();

		// nextInt is normally exclusive of the top value,
		// so add 1 to make it inclusive
		int randomNum = rand.nextInt((max - min) + 1) + min;

		return randomNum;
	}
	
	
	public static void main(){

		int n;
		int cost_range = 10000;
		for(n = 10; n <= 100; n += 10) {
//		gen a random n x n matrix
		int [][] cost_array = new int[n][n];
		for(int i = 0; i < n; ++i) {
			for(int j = 0; j < n; ++j) {
				cost_array[i][j] = randInt(0, cost_range);
			}	
		}
		
		
			
//		assignment by GREEDY, get total cost
//		assignment by OPT / Hungarian, get total cost
//		show the ratio			
		}

	}
}
