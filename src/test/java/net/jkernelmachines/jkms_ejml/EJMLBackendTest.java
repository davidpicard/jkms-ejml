/**
 * 
 */
package net.jkernelmachines.jkms_ejml;

import static org.junit.Assert.*;

import org.junit.Test;

import fr.lip6.jkernelmachines.util.algebra.MatrixOperations;

/**
 * @author picard
 *
 */
public class EJMLBackendTest {

	/**
	 * Test method for
	 * {@link net.jkernelmachines.jkms_ejml.EJMLBackend#inv(double[][])}.
	 */
	@Test
	public final void testInv() {
		EJMLBackend ejml = new EJMLBackend();

		double[][] A = { { 4, 0, 0 }, { 0, 3, 0 }, { 0, 0, 2 } };

		double[][] Ainv = ejml.inv(A);

		// A*(invA) = I
		double[][] I = MatrixOperations.mul(A, Ainv);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == j) {
					assertEquals(1.0, I[i][j], 1e-10);
				} else {
					assertEquals(0.0, I[i][j], 1e-10);
				}
			}
		}

		// test larger matrix
		int n = 128;
		double[][] X = new double[n][n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				X[i][j] = Math.random() * 2 - 1.0;

		double[][] G = MatrixOperations.transMul(X, X);
		for (int i = 0; i < n; i++) {
			G[i][i] += 1.0;
		}
		// A*(invA) = I
		double[][] Ginv = ejml.inv(G);

		I = MatrixOperations.mul(G, Ginv);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					assertEquals(1.0, I[i][j], 1e-10);
				} else {
					assertEquals(0.0, I[i][j], 1e-10);
				}
			}
		}

	}

	/**
	 * Test method for
	 * {@link net.jkernelmachines.jkms_ejml.EJMLBackend#eig(double[][])}.
	 */
	@Test
	public final void testEig() {

		EJMLBackend ejml = new EJMLBackend();
		
		double[][] A = { { 4, 3, 1 }, { 3, -3, -2 }, { 1, -2, 2 } };

		double[][][] eig = ejml.eig(A);

		// is U orthogonal?
		double[][] UtU = MatrixOperations.transMul(eig[0], eig[0]);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i == j) {
					assertEquals(1.0, UtU[i][j], 1e-10);
				} else {
					assertEquals(0, UtU[i][j], 1e-10);
				}
			}
		}

		// is L diagonal
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				if (i != j) {
					assertEquals(0.0, eig[1][i][j], 1e-10);
				}
			}
		}

//		// Does Q*L*Q' reconstruct A
//		double[][] C = MatrixOperations.mul(eig[0],
//				MatrixOperations.mul(eig[1], MatrixOperations.trans(eig[0])));
//		for (int i = 0; i < 3; i++) {
//			for (int j = 0; j < 3; j++) {
//				assertEquals(A[i][j], C[i][j], 1e-10);
//			}
//		}

		// test larger matrix
		int n = 256;
		double[][] X = new double[n][n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				X[i][j] = Math.random() * 2 - 0.5;

		double[][] G = MatrixOperations.transMul(X, X);
		double[][][] ei = ejml.eig(G);
		// is U orthogonal?
		UtU = MatrixOperations.transMul(ei[0], ei[0]);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					assertEquals(1.0, UtU[i][j], 1e-10);
				} else {
					assertEquals(0, UtU[i][j], 1e-10);
				}
			}
		}
		// is L diagonal
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					assertEquals(0.0, ei[1][i][j], 1e-10);
				} else {
					assertTrue((ei[1][i][i] + 1e-10) >= 0);
				}
			}
		}
		double[][] rec = MatrixOperations.mul(ei[0],
				MatrixOperations.mul(ei[1], MatrixOperations.trans(ei[0])));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				assertEquals(G[i][j], rec[i][j], 1e-10);
			}
		}
	}

}
