/**
 * 
 */
package net.jkernelmachines.jkms_ejml;

import java.util.Arrays;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;

import fr.lip6.jkernelmachines.util.algebra.AlgebraBackend;
import fr.lip6.jkernelmachines.util.algebra.ThreadedMatrixOperations;

/**
 * @author picard
 *
 */
public class EJMLBackend extends AlgebraBackend {

	/* (non-Javadoc)
	 * @see fr.lip6.jkernelmachines.util.algebra.AlgebraBackend#inv(double[][])
	 */
	@Override
	public double[][] inv(double[][] A) {
		int n = A.length;
		DenseMatrix64F A_ejml = new DenseMatrix64F(A);
		
		if(!CommonOps.invert(A_ejml))
			return null;
		double[][] A_inv = new double[n][n];
		for(int i = 0 ; i < n ; i++) {
			for(int j = 0 ; j < n ; j++) {
				A_inv[i][j] = A_ejml.get(i, j);
			}
		}
		return A_inv;
	}

	/* (non-Javadoc)
	 * @see fr.lip6.jkernelmachines.util.algebra.AlgebraBackend#eig(double[][])
	 */
	@Override
	public double[][][] eig(double[][] A) {
		int n = A.length;
		EigenDecomposition<DenseMatrix64F> dec = DecompositionFactory.eig(n, true, true);
		DenseMatrix64F G_ejml = new DenseMatrix64F(A);
		dec.decompose(G_ejml);
		double[][] lambda = new double[n][n];
		double[][] U = new double[n][];
		for(int i = 0 ; i < n ; i++) {
			lambda[i][i] = dec.getEigenvalue(i).getMagnitude();
			DenseMatrix64F u = dec.getEigenVector(i);
			U[i] = Arrays.copyOf(u.data, n);
			
		}
		return new double[][][]{ThreadedMatrixOperations.transi(U), lambda};
	}

}
