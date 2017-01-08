package sample;
/**
 * Created by Mark on 18.10.2014.
 */
import Jama.Matrix;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
public class Conjugate_gradient_method1 {
    int dim1;
    int dim2;

    Matrix A;
    Matrix b;
    Matrix x;

    double[] X;

    public Conjugate_gradient_method1(double[][] F, double[] l, int dim1, int dim2) {
        A = new Matrix(F, dim1, dim2);
        b = new Matrix(l, dim1);
        x = new Matrix(dim2,1);
        X = new double[101];
        Matrix r = b.minus(A.times(x)).copy();
        Matrix r1=r.copy();

        Matrix v = r.copy();
        Matrix p = r.copy(); Matrix s = r.copy(); Matrix t = r.copy();
        double alpha=1,p0=1,p0_new=1,w=1,b1=1;

        for (int i = 1; i < 10e3; i++) {
            p0_new=r1.transpose().times(r).get(0,0);
            b1=(p0_new/p0)*(alpha/w);
            p=r.plus(p.minus(v.times(w))).times(b1);
            v=A.times(p);
            alpha=p0_new/(r1.transpose().times(v).get(0,0));
            s=r.minus(v.times(alpha));

            t=A.times(s);
            w=t.times(s).get(0,0)/(t.times(t).get(0,0));
            x=x.plus(p.times(alpha)).plus(s.times(w));
            X[i]=A.times(x).minus(b).normInf();
            if (A.times(x).minus(b).normInf()<1) {
                break;
            }

            r=s.minus(t.times(w));
            p0=p0_new;
        }
        try{
            WriteToFile("NormDataX.txt", X);
        }
        catch (IOException e) {
            e.printStackTrace();
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