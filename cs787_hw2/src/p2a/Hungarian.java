package p2a;

/*************************************************************************
 *  Compilation:  javac Hungarian.java
 *  Execution:    java Hungarian N
 *  Dependencies: FordFulkerson.java FlowNetwork.java FlowEdge.java 
 *
 *  Solve an N-by-N assignment problem. Bare-bones implementation:
 *     - takes N^5 time in worst case.
 *     - assumes weights are >= 0  (add a large constant if not)
 *  
 *  For N^4 version: http://pages.cs.wisc.edu/~m/cs787/hungarian.txt
 *  Can be improved to N^3
 *
 *********************************************************************/

public class Hungarian {
    private int N;              // number of rows and columns
    private double[][] weight;  // the N-by-N weight matrix
    private double[] x;         // dual variables for rows
    private double[] y;         // dual variables for columns
    private int[] xy;           // xy[i] = j means i-j is a match
    private int[] yx;           // yx[j] = i means i-j is a match
 
    public Hungarian(double[][] weight) {
        this.weight = weight.clone();
        N = weight.length;
        x = new double[N];
        y = new double[N];
        xy = new int[N];
        yx = new int[N];
        for (int i = 0; i < N; i++) xy[i] = -1;
        for (int j = 0; j < N; j++) yx[j] = -1;

        while (true) {

            // build graph of 0-reduced cost edges
            FlowNetwork G = new FlowNetwork(2*N+2);
            int s = 2*N, t = 2*N+1;
            for (int i = 0; i < N; i++) {
                if (xy[i] == -1) G.addEdge(new FlowEdge(s, i, 1.0));
                else             G.addEdge(new FlowEdge(s, i, 1.0, 1.0));
            }
            for (int j = 0; j < N; j++) {
                if (yx[j] == -1) G.addEdge(new FlowEdge(N+j, t, 1.0));
                else             G.addEdge(new FlowEdge(N+j, t, 1.0, 1.0));
            }
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (reduced(i, j) == 0) {
                        if (xy[i] != j) G.addEdge(new FlowEdge(i, N+j, 1.0));
                        else            G.addEdge(new FlowEdge(i, N+j, 1.0, 1.0));
                    }
                }
            }

            // to make N^4, start from previous solution
            FordFulkerson ff = new FordFulkerson(G, s, t);

            // current matching
            for (int i = 0; i < N; i++) xy[i] = -1;
            for (int j = 0; j < N; j++) yx[j] = -1;
            for (int i = 0; i < N; i++) {
                for (FlowEdge e : G.adj(i)) {
                    if ((e.from() == i) && (e.flow() > 0)) {
                        xy[i] = e.to() - N;
                        yx[e.to() - N] = i;
                    }
                }
            }

            // perfect matching
            if (ff.value() == N) break;

            // find bottleneck weight
            double max = Double.POSITIVE_INFINITY;
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    if (ff.inCut(i) && !ff.inCut(N+j) && (reduced(i, j) < max))
                        max = reduced(i, j);

            // update dual variables
            for (int i = 0; i < N; i++)
                if (!ff.inCut(i))   x[i] -= max;
            for (int j = 0; j < N; j++)
                if (!ff.inCut(N+j)) y[j] += max;

    //        StdOut.println("value = " + ff.value());
        }
        assert check();
    }

    // reduced cost of i-j
    private double reduced(int i, int j) {
        return weight[i][j] - x[i] - y[j];
    }

    private double weight() {
        double totalWeight = 0.0;
        for (int i = 0; i < N; i++) totalWeight += weight[i][xy[i]];
        return totalWeight;
    }

    private int sol(int i) {
        return xy[i];
    }


    // check optimality conditions
    private boolean check() {
        // check that xy[] is a permutation
        boolean[] perm = new boolean[N];
        for (int i = 0; i < N; i++) {
            if (perm[xy[i]]) {
                StdOut.println("Not a perfect matching");
                return false;
            }
            perm[xy[i]] = true;
        }

        // check that all edges in xy[] have 0-reduced cost
        for (int i = 0; i < N; i++) {
            if (reduced(i, xy[i]) != 0) {
                StdOut.println("Solution does not have 0 reduced cost");
                return false;
            }
        }

        // check that all edges have >= 0 reduced cost
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (reduced(i, j) < 0) {
                    StdOut.println("Some edges have negative reduced cost");
                    return false;
                }
            }
        }
        return true;
    }

    public static void main(String[] args) {

        //int N = Integer.parseInt(args[0]);
    	int N = 10;
    	while (N <= 2000) {
    		double[][] weight = new double[N][N];
    		//        gw: edited loop
 //   		System.out.println();
    		for (int i = 0; i < N; i++) {
 //   			System.out.println();	
    			for (int j = 0; j < N; j++){
    				weight[i][j] = StdRandom.random();
 //   				System.out.print(weight[i][j] + ", ");
    			}            
    		}
 //   		System.out.println();
 //   		System.out.println();        
    		Hungarian assignment = new Hungarian(weight);
    		//StdOut.println("weight = " + assignment.weight());
    		//        for (int i = 0; i < N; i++)
    		//            StdOut.println(i + "-" + assignment.sol(i));

    		//Greedy
    		Boolean[] row_marker = new Boolean[N];
    		Boolean[] col_marker = new Boolean[N];
    		for(int i = 0; i < N; ++i) {
    			row_marker[i] = false;
    			col_marker[i] = false;
    		}

    		int counter = N;
    		double greedy_total = 0.0;
    		while(counter > 0) {
    			double min = Double.MAX_VALUE;
    			int min_row = -1, min_col = -1;
    			for(int i = 0; i < N; ++i) {
    				// if row is marked, it is crossed out, skip it
    				if(row_marker[i])
    					continue;
    				for(int j = 0; j < N; ++j) {
    					if(col_marker[j])
    						continue;
    					if(weight[i][j] < min) {
    						min = weight[i][j];
    						min_row = i;
    						min_col = j;
    					}
    				}
    			}
    			row_marker[min_row] = true;
    			col_marker[min_col] = true;
    			greedy_total += min;	        		

    			-- counter;
    		}
    		System.out.println("N = " + N + ", Ratio = " + greedy_total / assignment.weight());
    		
    		if(N < 100) N += 10;
    		else if(N < 1000) N += 100;
    		else N += 1000;
    	}
    }

}
