package sample;
/**
 * Created by Mark on 18.10.2014.
 */
import Jama.*;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.Math;
public class Conjugate_gradient_method {
    int dim1;
    int dim2;

    Matrix A;
    //   Matrix ATA=A.transpose().times(A);
    Matrix b;
    Matrix x;

    double[] X;

    public Conjugate_gradient_method(double[][] F, double[] l, int dim1, int dim2) {
        A = new Matrix(F, dim1, dim2);
        b = new Matrix(l, dim1);
        x = new Matrix(dim2,1);

        double rsnew = 0;
        Matrix ATA = (A.transpose().times(A)).copy();
        Matrix ATb=(A.transpose().times(b)).copy();

        Matrix r = ATb.minus(ATA.times(x)).copy();
        Matrix p = r.copy();
        double rsold = (r.transpose().times(r)).get(0, 0);

        for (int i = 0; i < 10e4; i++) {
            Matrix Ap = ATA.times(p);
            double alpha = rsold / (p.transpose().times(Ap).get(0, 0));
            x = p.times(alpha).plus(x);
            r = r.minus(Ap.times(alpha));
            rsnew = (r.transpose().times(r)).get(0, 0);
            if (Math.sqrt(rsnew) < 1e-10) {

                break;
            }
            p = p.times(rsnew / rsold).plus(r);
            rsold = rsnew;
        }



        X = new double[dim2];
        for (int i = 0; i < dim2; i++)
            X[i] = x.get(i, 0);

    }

    public Conjugate_gradient_method(double[][] F, double[] l, double[][] Regularization, int dim1, int dim2) {
        A = new Matrix(F, dim1, dim2);
        b = new Matrix(l, dim1);
        x = new Matrix(dim2,1);
        Matrix RegMat = new Matrix(Regularization,dim2, dim2);

        double rsnew = 0;
        Matrix ATA = (A.transpose().times(A)).copy();
        ATA = (ATA.plus(RegMat));
        Matrix ATb=(A.transpose().times(b)).copy();

        Matrix r = ATb.minus(ATA.times(x)).copy();
        Matrix p = r.copy();
        double rsold = (r.transpose().times(r)).get(0, 0);

        for (int i = 0; i < 10e4; i++) {
            Matrix Ap = ATA.times(p);
            double alpha = rsold / (p.transpose().times(Ap).get(0, 0));
            x = p.times(alpha).plus(x);
            r = r.minus(Ap.times(alpha));
            rsnew = (r.transpose().times(r)).get(0, 0);
            if (Math.sqrt(rsnew) < 1e-10) {

                break;
            }
            p = p.times(rsnew / rsold).plus(r);
            rsold = rsnew;
        }



        X = new double[dim2];
        for (int i = 0; i < dim2; i++)
            X[i] = x.get(i, 0);

    }

    public static void WriteToFile(String name, double[] X) throws FileNotFoundException, IOException {
        PrintWriter fout = new PrintWriter(new FileWriter(name));
        //fout.write(this.getSize());
        int ColumnLength = X.length;


        for (int j = 0; j < ColumnLength; ++j) {
            fout.print(String.format("%.4f", X[j]));
            fout.write("  ");
        }
        fout.write("\n");

        fout.close();

    }
}