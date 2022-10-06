package test.beast.evolution.substmodel;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.util.Randomizer;

public class ApproximateSubstModelTest {

	double a = 1;
	double b = 1;
	
	double [] times2;
	double [][] ratesThroughTime2;

	private void run() {
		Randomizer.setSeed(127196);
        Double [] f = new Double[]{Randomizer.nextDouble(),Randomizer.nextDouble(),Randomizer.nextDouble(),Randomizer.nextDouble()};
        double sum = 0;
        for (double d : f) {
        	sum += d;
        }
        for (int i = 0; i < 4; i++) {
        	f[i] = f[i] / sum;
        }
		String pi = f[0] + " " + f[1] + " " + f[2] + " " + f[3];
		RealParameter f2 = new RealParameter(pi);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", f2, "estimate", false);

        GTR gtr = new GTR();
        Double [] rates = new Double[]{Randomizer.nextDouble(),Randomizer.nextDouble(),Randomizer.nextDouble(),Randomizer.nextDouble(),Randomizer.nextDouble(),Randomizer.nextDouble()};
        gtr.initByName("rateAC", new RealParameter(rates[0]+""),
        		"rateAG", new RealParameter(rates[1]+""),
        		"rateAT", new RealParameter(rates[2]+""),
        		"rateCG", new RealParameter(rates[3]+""),
        		"rateCT", new RealParameter(rates[4]+""),
        		"rateGT", new RealParameter(rates[5]+""),
        		"frequencies", freqs);

        double t = 0;
        int N = 101;
        double [][] ratesThroughTime = new double[N][16];
        double [] times = new double[N];

        ratesThroughTime2 = new double[N][16];
        times2 = new double[N];
        
        int k = 0;
        while (t < 10) {
        	double [] transprobs = new double[16];
        	gtr.getTransitionProbabilities(null, t, 0, 1.0, transprobs);
    		System.out.print(t+ "\t");
    		System.arraycopy(transprobs, 0, ratesThroughTime[k], 0, 16);
        	for (int i = 0; i < transprobs.length; i++) {
        		System.out.print(transprobs[i] + "\t");
        	}
    		System.out.println();
    		times[k] = t;
    		
    		
        	gtr.getTransitionProbabilities(null, t/3.0, 0, 1.0, transprobs);
    		System.arraycopy(transprobs, 0, ratesThroughTime2[k], 0, 16);
    		times2[k] = t/3.0;
    		
    		
    		
    		if (t == 0) {
    			t += 1e-4;
    		} else {
    			t += 0.1;
    		}
    		k++;
        }

        for (int i = 0; i < 16; i++) {
        	approximate(ratesThroughTime, times, i);
        }
        
	}

	
	private void approximate(double[][] ratesThroughTime, double[] times, int index) {
		SplineInterpolator interpolator = new SplineInterpolator();
		double [] x = new double[16];
		double [] y = new double[16];
		int i = 0;
		set(i++, x, y, 0, times, ratesThroughTime, index);
		set(i++, x, y, 1, times, ratesThroughTime, index);
		set(i++, x, y, 2, times, ratesThroughTime, index);
		set(i++, x, y, 4, times, ratesThroughTime, index);
		set(i++, x, y, 6, times, ratesThroughTime, index);
		
		set(i++, x, y, 8, times, ratesThroughTime, index);
		set(i++, x, y, 10, times, ratesThroughTime, index);
		set(i++, x, y, 13, times, ratesThroughTime, index);
		set(i++, x, y, 16, times, ratesThroughTime, index);
		set(i++, x, y, 20, times, ratesThroughTime, index);

		set(i++, x, y, 25, times, ratesThroughTime, index);
		set(i++, x, y, 32, times, ratesThroughTime, index);
		set(i++, x, y, 42, times, ratesThroughTime, index);		
		set(i++, x, y, 57, times, ratesThroughTime, index);
		set(i++, x, y, 77, times, ratesThroughTime, index);		
		set(i++, x, y, 100, times, ratesThroughTime, index);
		
		
		PolynomialSplineFunction f = interpolator.interpolate(x, y);
		double error = calcError(f, times, ratesThroughTime, index);
		System.out.println(index  + "\t" + error + "\t");
		//System.out.println(Arrays.toString(diffs));

		
	}

	private double calcError(PolynomialSplineFunction poly, double[] times, double[][] ratesThroughTime, int index) {
		double error1 = 0;
		double error2 = 0;
		int N = times.length;
		diffs = new double[N];
		double maxDiff = 0;
		for (int i = 0; i < N; i++) {
			double f = poly.value(times[i]);
			double diff = (f - ratesThroughTime[i][index]);
			diffs[i] = diff;
			error2 += diff * diff;
			error1 += Math.abs(diff);
			maxDiff = Math.max(maxDiff, Math.abs(diff));
		}
		// return maxDiff;
		// return error1 / N;
		return error2 / N;
	}


	private void set(int i, double[] x, double[] y, int j, double[] times, double[][] ratesThroughTime, int index) {
		x[i] = times[j];
		y[i] = ratesThroughTime[j][index];
	}


	private void approximate2(double[][] ratesThroughTime, double[] times, int index) {
		final int N = times.length;
		double f0 = ratesThroughTime[N-1][index];
		a = ratesThroughTime[0][index]-f0;
		
//		double minE = 1e100;
//		int minK = -1;
//		for (int k = 0; k < 20; k++) {
//			b = -Math.log(Math.abs((ratesThroughTime[k][index] -f0) / a))/times[k];
//			double error = calcError(f0, a, b, index, ratesThroughTime, times);
//			if (error < minE) {
//				minE = error;
//				minK = k;
//			}
//		}
//		int k = minK;
		int k = 10;
		b = -Math.log(Math.abs((ratesThroughTime[k][index] -f0) / a))/times[k];
		double error = calcError(f0, a, b, index, ratesThroughTime, times);
		System.out.println(index + " " + k + "\t" + a + "\t" + b + "\t" + error + "\t");
		System.out.println(Arrays.toString(diffs));

        UnivariateFunction fa = new UnivariateFunction() {    		
			@Override
			public double value(double a) {
				double error = calcError(f0, a, b, index, ratesThroughTime, times);
				return error;
			}        	
        };
        
        UnivariateFunction fb = new UnivariateFunction() {    		
			@Override
			public double value(double b) {
				double error = calcError(f0, a, b, index, ratesThroughTime, times);
				return error;
			}
        };

	
        
        double EPS_ABS = 1e-6;
        double EPS_REL = 1e-6;
        
        for (int i = 0; i < 10; i++) {
	        BrentOptimizer optimizerA = new BrentOptimizer(EPS_ABS, EPS_REL);
	        UnivariatePointValuePair pa =
			optimizerA.optimize(new MaxEval(200),
	                new UnivariateObjectiveFunction(fa),
	                GoalType.MINIMIZE,
	                new InitialGuess(new double [] {a}),
	                new SearchInterval(-5,5));
	        a = pa.getPoint();
	        
	        BrentOptimizer optimizerB = new BrentOptimizer(EPS_ABS, EPS_REL);
	        UnivariatePointValuePair pb =
			optimizerB.optimize(new MaxEval(200),
	                new UnivariateObjectiveFunction(fb),
	                GoalType.MINIMIZE,
	                new InitialGuess(new double [] {b}),
	                new SearchInterval(-5,5));
	        b = pb.getPoint();
        }
        
		error = calcError(f0, a, b, index, ratesThroughTime, times);
		System.out.println(index + "\t" + a + "\t" + b + "\t" + error + "\t");
		//System.out.println(index + " " + error);

	}

	double [] diffs;

	private double calcError(double f0, double a, double b, int index, double[][] ratesThroughTime, double[] times) {
		double error1 = 0;
		double error2 = 0;
		int N = times.length;
		diffs = new double[N];
		double maxDiff = 0;
		for (int i = 0; i < N; i++) {
			double f = f0 + a * Math.exp(-b * times[i]);
			double diff = (f - ratesThroughTime[i][index]);
			diffs[i] = diff;
			error2 += diff * diff;
			error1 += Math.abs(diff);
			maxDiff = Math.max(maxDiff, Math.abs(diff));
		}
		// return maxDiff;
		// return error1 / N;
		return error2 / N;
	}        	


	public static void main(String[] args) {
		ApproximateSubstModelTest s = new ApproximateSubstModelTest();
		s.run();
	}


}
