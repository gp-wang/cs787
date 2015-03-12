package p2a;
import java.util.*;

public class PathFinder_rev2 {

	public static void main(String[] args) {
	// graph setup 1:
//		
//		int m = 3, n = 3;
//		int mass = 2, g = 10 * m * m / n, dt = 1; 
//
//		EdgeWeightedDigraph G = new EdgeWeightedDigraph((m + 1) * (n + 1));
//		
////		< m, not <= m 
//		for(int i = 0; i < m; ++i) {
//			// current point at time = i
//			for(int j = 0; j <= n; ++j) {
//				// next point at time = (i + 1)
//				for (int k = 0; k <= n; ++k) {
//					// gw: correction made, index_k should be at time = i + 1
//					int index_j = i * (n + 1) + j,
//						index_k = (i + 1) * (n + 1) + k;
//						
//					int y_prev = j,
//						y_curr = k,
//						dy_curr = y_curr - y_prev,
//						velocity_curr = dy_curr / dt, 
//						kinetic_curr = (mass / 2) * velocity_curr * velocity_curr,
//						potential_curr = mass * g * y_curr,
//						integrand_curr = kinetic_curr - potential_curr; //
//					
//					//TODO: udpate edge_weight formula
//					double edge_weight = integrand_curr;
//					G.addEdge(new DirectedEdge(index_j, index_k, edge_weight));
//				}
//			}
//		}
		//Graph Setup 2:
		//int m = 15, n = 300;
		int m = 4, n = 4;
		double mass = 2.0;// g = 10.0 * m * m / n;
		double g = 10.0;

		EdgeWeightedDigraph G = new EdgeWeightedDigraph((m + 1) * (n + 1));
		
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
			index_s = 0 * (n + 1) + n / 4, 
			index_t = (m) * (n + 1) + 0;

		int s = index_s;
        

        BellmanFordSP sp = new BellmanFordSP(G, s);

        // print negative cycle
        if (sp.hasNegativeCycle()) {
            for (DirectedEdge e : sp.negativeCycle())
                StdOut.println(e);
        }

        // print shortest paths
        else {
            for (int v = 0; v < G.V(); v++) {
//            	gw: only print the corner point
            	if(v != index_t) continue;
                if (sp.hasPathTo(v)) {
                    StdOut.printf("%d to %d (%5.2f)  ", s, v, sp.distTo(v));
                    for (DirectedEdge e : sp.pathTo(v)) {
                        StdOut.print(e + "   ");
                    }
                    StdOut.println();
                }
                else {
                    StdOut.printf("%d to %d           no path\n", s, v);
                }
            }
        }
		System.out.println("index_s: " + index_s + ", index_t: " + index_t);
		

	}
	
//	G, index_t, M are fixed params while recurring
	private static double opt_recurr(EdgeWeightedDigraph G, int index_t, double [][] M, int i, int index_v) {
//		TODO: think more about base case
//		base case
		if (i == 0) 
			return M[0][index_v];
		else
			return min(
					opt_recurr(G, index_t, M, i-1, index_v),
					//gw: finding index_w is inside the function, here we pass in simply index_v	
					// gw: note the (i-1) is executed inside the function, so here we simply put i.			
					min_among_connection(G, index_t, M, i, index_v )
					);
	} 
	
	private static double min(double a, double b) {
		return a < b ? a : b;
	}
	
	
//	TODO: memoize
//	TODO: recovering solution
	
	
	private static double min_among_connection(EdgeWeightedDigraph G, int index_t, double [][] M, int i, int index_v) {
//		(done!)old_todo: verify that no base case is needed (probably because its called proc has base case.)
//		ans: no need base case, will reduce to opt_recurr's base case
		double min = Double.MAX_VALUE;
		
//		(done!)old_todo: verify that e is LEAVING v
		for(DirectedEdge e : G.adj(index_v)) {
			int index_w = e.to();
			double value = opt_recurr(G, index_t, M, i-1, index_w) + e.weight();
			if(value < min) {
				min = value;
			}
				
		}
		return min;
	}
	
// 	TODO: double check
//	G, index_s, M, m, n are fixed
//	i is number of edges to index_t
	private static String find_solution(EdgeWeightedDigraph G, int index_s, double [][] M, int m, int n, int i, int index_w){
//		TODO: verify base case logic
		//if( index_w == index_s)
//		if(i == 0)
//		TODO: verify, 
		if(i == m)
			return (new String("Note: source reached!" + '\n'));
		else {
			int index_v = -1;
			
			// For a given index_w, find value = (c_vw + M[i-1][index_v]) and the index_v s.t. value is min.
			// It is guaranteed there exists such a index_v
			double min = Double.MAX_VALUE;
			
			// calculate the minimum of all c_vw where vw is an edge.
			// TODO: index_w is ending point, you cannot get all edge(index_v, index_w) using adj(w)
			// gw: correction, build a new iterable
			ArrayList<DirectedEdge> edges_to_w = new ArrayList<DirectedEdge>();
			for(DirectedEdge e : G.edges()) {
				if(e.to() == index_w) {
					edges_to_w.add(e);
				}
			}

			for(DirectedEdge e : edges_to_w) {
				int temp_index_v = e.from();
				//double value = M[i-1][temp_index_v] + e.weight();
				//TODO gw: tentative correction: i-1 => m-(i - 1)
				//double value = M[i-1][temp_index_v] + e.weight();
				//gw: correction from i -1 to i + 1, because we are moving backward to temp_index_v from index_w
				//double value = M[i+1][temp_index_v] + e.weight();
				//gw: likely no need to add e.weight(), 
				double value = M[i+1][temp_index_v] + e.weight();	
				if (value < min) {
					min = value;
					index_v = temp_index_v;
				}
			}
			
			//gw: correction, i - 1 => i + 1
			if (min < M[i+1][index_w]) {
				// index_v is in the opt Path
				//				x = (m - i) // m +1 time instants, m intervals, i remaining edges to t, 
				return (new String("i: " + i + ", index_v: " + index_v + ", height_v: " + (index_v % (n + 1)) + '\n' 
						+ find_solution(G, index_s, M, m, n, i+1, index_v) 
						));
			}
			else {
				// index_v is not in the opt Path
				// gw: correction i - 1 => i + 1
				return find_solution(G, index_s, M, m, n, i+1, index_w) + '\n';
			}
			
			
		}
			
	}
	
	

}
