  
hw Problem description  
--  
2. These questions deal with the assignment problem. Use the minimizing version, as presented  
in class.  
a) Consider the following heuristic. Applied to an n × n array of cost values, GREEDY  
chooses an entry of minimum cost, and crosses off its row and column. This is done  
recursively, on the remaining (n−1)×(n−1) matrix, stopping when a complete assignment  
has been made.  
Show by examples that GREEDY has no approximation factor. To do this, you should  
exhibit a sequence of cost arrays P1, P2, P3, . . . for which  
	limit_infinity ( GREEDY(Pn)/OPT(Pn) ) = infinity  
b) For n = 2, 3, . . ., let Mn be the matrix with entries  
	mij = gcd(i, j), 2 ≤ i, j ≤ n.  
Compute An, the value of an optimal assignment for Mn. (The answers are different  
depending on the parity of n. Consider odd n first.)
