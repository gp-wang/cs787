package p2a;

public class GCD_Matrix_Assignment {
	
	static int gcd(int a, int b)
	{
		if(a == 0 || b == 0) return a+b; // base case
		return gcd(b,a%b);
	}
	public static void main(String[] args) {
		for (int i = 0; i < 20; i++) {
			System.out.println();
			for (int j = 0; j < 20; j++) {
				System.out.print(gcd(i + 1,j + 1) + " ");
			}
			
		}
	}
}
