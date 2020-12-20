
public class Test2 {

	public static void main(String[] args) {
		double b[] = { 1, 1 };
		double d[] = {1, 1, 1};
        double C[][] = { { 1, 0 }, { 0, 1 } };
        double D[][] = { { 1, 2, 3 }, { 4, 5, 6 }, {7, 8, 9} };
        double E[][] = { { 7, 3, -5 }, { -1, -2, 4 }, {-4, 1, -3} };
        double F[][] = { { 1, 0, -5 }, { 0, 0, 4 }, {0, 0, -3} };
        double e[] = {-12, 5, 1};
        
        //printIt(Gauss.solve(F, e));
        printIt(Gauss.solveSing(D));
	}
	
	public static void printIt(double[] x) {
		System.out.print("[");
		for(int i = 0; i < x.length; i++) {
			System.out.print(x[i]);
			if(i < x.length - 1) {
				System.out.print(", ");
			}
		}
		System.out.println("]");
	}
	
	public static void test1(double[][] R) {
		swap(R, 0, 1);
	}
	
	private static void swap(double[][] R, int i, int j) {
    	double[] tmp = R[i];
    	R[i] = R[j];
    	R[j] = tmp;
    }
}
