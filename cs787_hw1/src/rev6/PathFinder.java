
package rev6;
// gw: note, corrected find-solution, trace back according to M[i][v] = M[i-1][w] + cost_vw
// gw: note, the graph framework is modified slightly by adding adj_rev[], for tracking all edges incident TO an vertex.
// gw: note, memoization is impl'd, it is crucial to running time.
//Test: m = 15, n = 300, result as below:
//Min value is:-38765.0
//Source at :75
//
//i: 1, index_v: 75, height_v: 75
//i: 2, index_v: 438, height_v: 137
//i: 3, index_v: 791, height_v: 189
//i: 4, index_v: 1134, height_v: 231
//i: 5, index_v: 1467, height_v: 263
//i: 6, index_v: 1790, height_v: 285
//i: 7, index_v: 2103, height_v: 297
//i: 8, index_v: 2407, height_v: 300
//i: 9, index_v: 2705, height_v: 297
//i: 10, index_v: 2993, height_v: 284
//i: 11, index_v: 3271, height_v: 261
//i: 12, index_v: 3539, height_v: 228
//i: 13, index_v: 3798, height_v: 186
//i: 14, index_v: 4047, height_v: 134
//i: 15, index_v: 4286, height_v: 72
//Sink at: 4515

import java.util.Arrays;
import java.util.*;

public class PathFinder {
	public DoubleLinkedEdgeWeightedDigraph G;
	public double [][] M;
	public boolean [][] isEmpty;
	
	public PathFinder(DoubleLinkedEdgeWeightedDigraph G, double [][] M, boolean [][] isEmpty) {
		this.G = G;
		this.M = M;
		this.isEmpty = isEmpty;
	}
	

	public static void main(String[] args) {

		
		// Graph set-up 2:
		int m = 15, n = 300;
		double mass = 2.0;// g = 10.0 * m * m / n;
		double g = 10.0;

		DoubleLinkedEdgeWeightedDigraph G = new DoubleLinkedEdgeWeightedDigraph((m + 1) * (n + 1));
		
//		Assumption 1: can only go forward in time
//		< m, not <= m
		// current time = i for index_j		
		for(int i = 0; i < m; ++i) {
			// current height = j for index_j
			for(int j = 0; j <= n; ++j) {
				// -------
				// next  height = k for index_k
				for (int k = 0; k <= n; ++k) {
					
					// next time = l for index_k
					for(int l = i + 1; l <= m; ++l ) {

						// gw: correction made, index_k should be at time = i + 1
						int index_j = i * (n + 1) + j,
								index_k = l * (n + 1) + k;

						double y_prev = j,
								y_curr = k,
								dy_curr = y_curr - y_prev,
								velocity_curr = dy_curr / (l - i), // dt = l - i 
								kinetic_curr = (mass / 2.0) * velocity_curr * velocity_curr,
								potential_curr = mass * g * y_curr,
								integrand_curr = kinetic_curr - potential_curr; //
								//integrand_curr = Math.abs(kinetic_curr - potential_curr); //
//						edge_weight should consider both initial and final point
//						edge_weight = delta_integrand
//									= 0.5 * mass * (dy_curr* dy_curr) / (l - i) - 0.5 * g * mass * (l-i) * (y_prev + y_curr)
						//double edge_weight = integrand_curr;
						//TODO: why this expression works?
						double edge_weight = 0.5 * mass * (dy_curr* dy_curr) / (l - i) - 0.5 * g * mass * (l-i) * (y_prev + y_curr);
						G.addEdge(new DirectedEdge(index_j, index_k, edge_weight));


						
					}
				}
			}
		}
		
		// => set s and t here
		// Shortest-Path (G = G, s = index_(0,n), t = index_(m,0))		
		int n_node = G.V(), 
			num_index = (m + 1) * (n + 1), 
			index_s = 0 * (n + 1) + n/4,
			//index_s = 0 * (n + 1) + n, 
			index_t = m * (n + 1) + 0;
		// M[edge count][vertices]		
			// edge count need to include 0 edge, so +1
		// TODO: init to infinity?
		double [][] M = new double[(n_node - 1) + 1][num_index];
		boolean [][] isEmpty = new boolean[(n_node - 1) + 1][num_index];
		for(int i  = 0; i != n_node; ++i) {
			for(int j = 0; j != num_index; ++j) {
				M[i][j] = Double.POSITIVE_INFINITY;
				isEmpty[i][j] = true;
			}
		}

		PathFinder PF = new PathFinder(G, M, isEmpty);	
		for(int i = 0; i != num_index; ++ i) {
			if (i == index_s) { 
				PF.M[0][i] = 0;
				PF.isEmpty[0][i] = false;
			}
			else {
				PF.M[0][i] = Double.POSITIVE_INFINITY;
				PF.isEmpty[0][i] = false;
			}
		}
		

		
		for(int i = 1; i != n_node; ++i) {	
			for(int j = 0; j != num_index; ++j) {
				// gw: currection, should use i, (not i - 1)				
				PF.M[i][j] = PF.memoized_opt_recurr(index_s, i, j);
				PF.isEmpty[i][j] = false;
//				TODO: implement first pointer : find when M is updated
			}
		}
		
//		now get the path out of M
				
		int edge_min = Integer.MAX_VALUE;
		double min = Double.POSITIVE_INFINITY;
		for (int k = 0; k <= m; ++k) {
			// not necessarily i-1, can be from 0 ..  i-1
			double value = M[k][index_t];	
			if (value < min) {
				min = value;
				edge_min = k;
			}	
		}
		System.out.println("Min value is:" + min);

//		(done!)old_todo: unfinished, think about what i to pass in
		
		
		
		System.out.print(PF.find_solution(index_s, m, n, edge_min, index_t) + '\n');
		System.out.println("Sink at: " + index_t);
		

	}
	
//	G, index_t, M are fixed params while recurring
	private double memoized_opt_recurr(int index_s, int i, int index_w) {
//		TODO: think more about base case
//		base case
		if (i == 0) 
			return M[0][index_w];
		else if(!isEmpty[i][index_w]){
			// TODO: tbc, likely no update is involved in the full table of M[][], either it is empty or filled with a final value.
			return M[i][index_w];
		}
		else
			return min(
					memoized_opt_recurr(index_s, i-1, index_w),
					//gw: finding index_w is inside the function, here we pass in simply index_v	
					// gw: note the (i-1) is executed inside the function, so here we simply put i.			
					min_among_connection(index_s, i, index_w )
					);
	} 
	
	private static double min(double a, double b) {
		return a < b ? a : b;
	}
	
	
//	TODO: memoize
//	TODO: recovering solution
	
	
	private double min_among_connection(int index_s, int i, int index_w) {
//		(done!)old_todo: verify that no base case is needed (probably because its called proc has base case.)
//		ans: no need base case, will reduce to opt_recurr's base case
		double min = Double.POSITIVE_INFINITY;

		for(DirectedEdge e : G.adj_rev(index_w)) {

			int temp_index_v = e.from();

			double value = memoized_opt_recurr(index_s, i-1, temp_index_v) + e.weight();
			if(value < min) {
				min = value;
			}
							
			
		}		

		return min;
	}
	
// 	TODO: double check
//	G, index_s, M, m, n are fixed
//	i is number of edges to index_t
	private String find_solution(int index_s, int m, int n, int i, int index_w){

		if(i == 0)
		//if(index_w == index_s)
			return (new String("Source at :" + index_w + '\n'));

		else {
			int index_v = index_w;
			
			// For a given index_w, find value = (c_vw + M[i-1][index_v]) and the index_v s.t. value is min.
			// It is guaranteed there exists such a index_v
			double min = Double.POSITIVE_INFINITY;
			
			// calculate the minimum of all c_vw where vw is an edge.
			// TODO: index_w is ending point, you cannot get all edge(index_v, index_w) using adj(w)
			
			double epsilon = 1E-3;
			for(DirectedEdge e : G.adj_rev(index_w)) {
				int temp_index_v = e.from();
				
				//if(temp_index_v == index_s) continue;

				double value = M[i-1][temp_index_v] + e.weight();	
				if (Math.abs(value - M[i][index_w]) < epsilon) {
					min = value;
					index_v = temp_index_v;
				}	

			}			
			
			return (new String(find_solution(index_s, m, n, i-1, index_v)+ '\n'
					+ "i: " + i + ", index_v: " + index_v + ", height_v: " + (index_v % (n + 1))));


		}
			
	}
	
	

}
