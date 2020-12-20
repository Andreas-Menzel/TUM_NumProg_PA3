
public class Gauss {

    /**
     * Diese Methode soll die Loesung x des LGS R*x=b durch
     * Rueckwaertssubstitution ermitteln.
     * PARAMETER:
     * R: Eine obere Dreiecksmatrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] backSubst(double[][] R, double[] b) {
    	int n = R.length;						// R ist n x n Matrix.
    	double[] x = new double[n];				// x ist Vektor der Länge n.
    	int m = n - 1;							// m als Indexwert von n.
    	for(int i = m; i >= 0; i--) {			// Matrix von unten nach oben durchgehen.
    		double r = 0.0;						// Zwischenergebnis aller bereits berechneten Koeffizienten.
			for(int j = i + 1; j < n; j++) {	// Alle bereits berechneten x benutzen.
				r += R[i][j] * x[j];			// Koeffizient mal korrespondierendes x wird aufaddiert.
			}
			x[i] = (b[i] - r) / R[i][i];		// Die einfache Gleichung (R[i][i] * x[i] + r = b[i]) wird nach x aufgelöst.
    	}
        return x;
    }

    /**
     * Diese Methode soll die Loesung x des LGS A*x=b durch Gauss-Elimination mit
     * Spaltenpivotisierung ermitteln. A und b sollen dabei nicht veraendert werden.
     * PARAMETER: A:
     * Eine regulaere Matrix der Groesse n x n
     * b: Ein Vektor der Laenge n
     */
    public static double[] solve(double[][] A, double[] b) {
    	int n = A.length;												// Dimension n.
    	double[][] R = new double[n][n];								// R als Kopie von A.
    	double[] r = new double[n];										// r als Kopie von b.
    	for(int i = 0; i < n; i++) {									// R und r initialisieren.
    		r[i] = b[i];
    		for(int j = 0; j < n; j++) {
    			R[i][j] = A[i][j];
    		}
    	}
    	for(int i = 0; i < n; i++) {									// Pivot-Suche.
    		int pivotIndex = i;											// Pivot-Index.
        	double pivot = R[i][i];										// Pivot-Wert (i-tes Element in Spalte als Init-Wert).
        	for(int j = i; j < n; j++) {								// Pivot finden.
        		if(Math.abs(R[j][i]) > Math.abs(pivot)) {				// Falls neues maximales Element gefunden:
        			pivot = R[j][i];									// aktualisiere pivot.
        			pivotIndex = j;										// aktualisiere pivotIndex.
        		}
        	}
        	/*
        	 * Falls der Betrag von pivot kleiner als 1E-10 ist, so kann kein Pivot-Element gefunden werden.
        	 * In diesem Fall machen weitere Berechnungen keinen Sinn, da die Hauptdiagonale nun ab hier ausschließlich 0er enthält, durch
        	 * welche wir im Verlauf dieses Algorithmus teilen würden, wobei dann als Ergebnis dieser Berechnungen Infinity übernommen wird.
        	 * Im Anschluss rechnet dann bakSubst mit Infinity, was ebenfalls nicht definiert ist.
        	 * Da die Methode solveSing() nun jedoch darauf abzielt, genau dieses erste |pivot| == 0 zu finden, verwenden wir hier diese undefinierten
        	 * Fälle als nützliche Ausgabe für solveSing() und sparen uns dort kopierten Code.
        	 * Es wird nun für einen solchen Fall ein neuer, spezieller Vektor zurückgegeben, welcher als ersten Eintrag NaN als Schlüssel enthält (dadurch können
        	 * valide Berechnungen, also falls die Matrix doch invertierbar war, von diesem Fall unterschieden werden). Außerdem enthält der Vektor noch
        	 * die Indexkoordinaten des Pivot-Elements. Dieses steht an Stelle (i, i), weshalb nur i zurückgegeben werden muss.
        	 */
        	if(Math.abs(pivot) < 0.0000000001) {
        		return new double[]{Double.NaN, i};
        	}
        	double[] tmp1 = R[i];										// pivotIndex an oberste Stelle tauschen.
        	R[i] = R[pivotIndex];
        	R[pivotIndex] = tmp1;
        	double tmp2 = r[i];											// Ergebnisvektor r muss ebenfalls getauscht werden.
        	r[i] = r[pivotIndex];
        	r[pivotIndex] = tmp2;
        	for(int j = i + 1; j < n; j++) {							// Alle folgenden Spaltenelemente unter Zeile i durchgehen.
        		double times = R[j][i];									// Wie oft man die i-te Spalte abziehen muss.
        		for(int k = i; k < n; k++) {							// Jede j-te Zeile (von i + 1 nach unten) durchgehen.
        			R[j][k] = R[j][k] - ((R[i][k] / pivot) * times);	// Vorherige Zahl - ((Korrespondierende Zahl in Zeile i / pivot) * times)
        		}
        		r[j] = r[j] - ((r[i] / pivot) * times);					// Selbige Berechnung für den Ergebnisvektor durchführen.
        	}
    	}
        return backSubst(R, r);											// Verwendung der Rücksubstitution zum Lösen des LGS.
    }

    /**
     * Diese Methode soll eine Loesung p!=0 des LGS A*p=0 ermitteln. A ist dabei
     * eine nicht invertierbare Matrix. A soll dabei nicht veraendert werden.
     *
     * Gehen Sie dazu folgendermassen vor (vgl.Aufgabenblatt):
     * -Fuehren Sie zunaechst den Gauss-Algorithmus mit Spaltenpivotisierung
     *  solange durch, bis in einem Schritt alle moeglichen Pivotelemente
     *  numerisch gleich 0 sind (d.h. <1E-10)
     * -Betrachten Sie die bis jetzt entstandene obere Dreiecksmatrix T und
     *  loesen Sie Tx = -v durch Rueckwaertssubstitution
     * -Geben Sie den Vektor (x,1,0,...,0) zurueck
     *
     * Sollte A doch intvertierbar sein, kann immer ein Pivot-Element gefunden werden(>=1E-10).
     * In diesem Fall soll der 0-Vektor zurueckgegeben werden.
     * PARAMETER:
     * A: Eine singulaere Matrix der Groesse n x n
     */
    public static double[] solveSing(double[][] A) {
    	int n = A.length;										// Dimension von A.
    	double[] b = new double[n];								// Lösungsvektor eines homogenen LGS (Nullvektor).
    	double[] f = solve(A, b);								// Zwischenergebnis. Findet heraus, ob das LGS invertierbar ist oder nicht.
    	if(Double.isNaN(f[0])) {								// Nicht invertierbar.
    		int firstZero = (int) f[1];							// Index der ersten 0 in der Hauptdiagonalen.
    		double[][] T = new double[firstZero][firstZero];	// Teilmatrix.
    		double[] v = new double[firstZero];					// Teilvektor nach Aufgabenbeschreibung.
    		for(int i = 0; i < firstZero; i++) {				// Beide initialisieren.
    			v[i] = -A[i][firstZero];
    			for(int j = 0; j < firstZero; j++) {
        			T[i][j] = A[i][j];
        		}
    		}
    		double[] sol = solve(T, v);							// Neues LGS lösen.
    		double[] newSol = new double[n];					// Dimension der Lösung auf n erweitern.
    		for(int i = 0; i < firstZero; i++) {				// Kopieren der ersten firstZero Stellen.
    			newSol[i] = sol[i];
    		}
    		newSol[firstZero] = 1;								// An der Stelle der ersten 0 kommt eine 1 (danach nur noch 0er).
    		return newSol;
    	}
    	else {
    		return b;											// Matrix invertierbar (und homogenes LGS) bedeutet Nullvektor als Ergebnis.
    	}
    }

    /**
     * Diese Methode berechnet das Matrix-Vektor-Produkt A*x mit A einer nxm
     * Matrix und x einem Vektor der Laenge m. Sie eignet sich zum Testen der
     * Gauss-Loesung
     */
    public static double[] matrixVectorMult(double[][] A, double[] x) {
        int n = A.length;
        int m = x.length;

        double[] y = new double[n];

        for (int i = 0; i < n; i++) {
            y[i] = 0;
            for (int j = 0; j < m; j++) {
                y[i] += A[i][j] * x[j];
            }
        }

        return y;
    }
}
