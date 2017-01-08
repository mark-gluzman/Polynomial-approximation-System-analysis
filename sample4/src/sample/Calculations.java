package sample;

import Jama.Matrix;
import com.sun.rowset.internal.Row;

import javax.swing.*;
import java.io.IOException;

/**
 * Created by Gregory on 17.10.2014.
 */
enum Polynom {Chebyshev,Lejandr, Lagger,Hermit}

public class Calculations {
//--------- Variables---------------------
    // dim = 3 array
    int[] PolyOrderJ;

    // for Y recovering
    double MaxY;
    double MinY;
    double coef;
    Conjugate_gradient_method Method;


    double[] Ax;
    // row number in NormalizedData
    int NumberOfSamples;

    Polynom type;
    // Here is cofficients for a_i
    // dimension n1, n2, n3
    double[]A1;
    double[]A2;
    double[]A3;
    double[]AincludeAll;
    double[] residualBigF_System;
    double residualBigF_System_max;
    double []  ResidVecLambda;
    double ResidVecLambda_max;
double[] C;
    // the third level of the hierarchy:
        double[] BigF1;
        double c1;
        double[] BigF2;
        double c2;
        double[] BigF3;
        double c3;
        // BigF for system;
        double[][] BigF_System;




    //for the second method of lambda calculation
    double [][] A_Lambda1;
    double [][] A_Lambda2;
    double [][] A_Lambda3;

    // Psi vector
    // dimention n1,n2,n3
    double[] Psi1, Psi2, Psi3;

    //PsiMatrices
    double[][] PsMat1, PsMat2, PsMat3;
    double[][] PsMat;

    // arrays of lambda in each of the three double sum
	// n1*polyOrderJi
    double[][] Lambda1, Lambda2, Lambda3;
    double[] LambdaVec;

    // Matrix for lambda Linear System A*lambda=b
    double[][] A;

    // arrays of polynom value
    // n1*polyOrderJi
    double[][] Tp1, Tp2, Tp3;

    // recovering
    // Y_model - model recovering
    // Y_true - true recovering
    double[] Y_model;
    double[] Y_true;


    // The time is an index in massive now

    int indexPredict;

    // ou indexes for training sets are :[time - 40, time].
    // after it can be the real time.




// -----------------------------------end of variables------------------------------------



   // elements of this Rezult will be the arguments of the double sum
    public double[][] ElemByElemMatrixMultiplication(double[][] Mat1, double[][] Mat2) {
        int RowSize = Mat1.length;
        int ColumnSize = Mat1[0].length;
        if (RowSize == Mat2.length && ColumnSize == Mat2[0].length) {
            double[][] Result = new double[RowSize][ColumnSize];
            for (int i = 0; i < RowSize; ++i) {
                for (int j = 0; j < ColumnSize; ++j) {
                    Result[i][j] = Mat1[i][j] * Mat2[i][j];
                }

            }
            return Result;
        } else return null;
    }

    // calculation of double sum
 
    public double ElementSum(double[][] Mat) {
        int RowSize = Mat.length;
        int ColumnSize = Mat[0].length;
        double sum = 0;
        for (int i = 0; i < RowSize; ++i) {
            for (int j = 0; j < ColumnSize; ++j) {
                sum = sum + Mat[i][j];
            }
        }
        return sum;
    }

    //-----!!! polynom generation
    public double Chebyshev(double x, int rang) {
        if (rang < 0) return 0;
        if (rang == 0) return 1;
        if (rang == 1) return -1 + 2 * x;

        return 2 * (-1 + 2 * x) * Chebyshev(x, rang - 1) - Chebyshev(x, rang - 2);
    }
    public double Lejandr(double x, int rang) {
        if (rang < 0) return 0;
        if (rang == 0) return 1;
        if (rang == 1) return x;

        return (1.0/(rang))*((2*rang - 1)*x*Lejandr(x,rang-1) -
                (rang-1)*Lejandr(x,rang - 2));
    }

    public double Lagger(double x, int rang) {
        if (rang < 0) return 0;
        if (rang == 0) return 1;
        if (rang == 1) return -x + 1;

        return (2 * rang - 1 - x) * Lagger(x, rang - 1) -
                (rang - 1) * (rang - 1) * Lagger(x, rang - 2);
    }

    public double Hermit(double x, int rang) {
        if (rang < 1) return 0;
        if (rang == 1) return 2*x;
        if (rang == 2) return 4*x*x - 2;

        return 2*x*Hermit(x,rang - 1)  - 2*(rang - 1)*Hermit(x, rang - 2);
    }

    public double [][] Transpose(double [][] Mat){
        int RowSize = Mat.length;
        int ColumnSize = Mat[0].length;
        double[][] MatTrans = new double[ColumnSize][RowSize];
        double temp;
        for(int i = 0; i < RowSize; ++i){
            for(int j = 0; j < ColumnSize; ++j){
                MatTrans[j][i] = Mat[i][j];
            }
        }
        return MatTrans;
    }

    // !!!
    // Make PolyFeatures from one vector;
    double[] PolynomVector(double [] x){

        int RowSize = x.length;
        double[][] PolyVec = new double[RowSize][PolyOrderJ[0]+1];
        switch (type) {
            case Chebyshev:
                for (int i = 0; i < RowSize; ++i) {
                    for (int j = 0; j <= PolyOrderJ[0]; ++j) {
                        PolyVec[i][j] = Chebyshev(x[i], j);
                    }
                }
                break;
            case Lagger:
                for (int i = 0; i < RowSize; ++i) {
                    for (int j = 0; j <= PolyOrderJ[0]; ++j) {
                        PolyVec[i][j] = Lagger(x[i], j);
                    }
                }
                break;
            case Lejandr:
                for (int i = 0; i <RowSize; ++i) {
                    for (int j = 0; j <= PolyOrderJ[0]; ++j) {
                        PolyVec[i][j] = Lejandr(x[i], j);
                    }
                }
                break;
            case Hermit:
                for (int i = 0; i < RowSize; ++i) {
                    for (int j = 0; j <= PolyOrderJ[0]; ++j) {
                        PolyVec[i][j] = Hermit(x[i], j);
                    }
                }
                break;
        }
        return MatrixShapeToVector(PolyVec);

    }

    // ------ matrix computation using polynomials

    double[][] PolynomMatrix(double [][] Xj,int indexQ, int PolyOrder, Polynom type) {
        // Here Xj must be transposed, for example Xj has 2*45 dimensionality
        // indexQ - q1 or q2 or q3 or q0 in big formula
        int N1 = Xj.length;
        double[][] PolyMat = new double[N1][PolyOrder + 1];
        switch (type) {
            case Chebyshev:
                for (int i = 0; i < N1; ++i) {
                    for (int j = 0; j <= PolyOrder; ++j) {
                        PolyMat[i][j] = Chebyshev(Xj[i][indexQ], j);
                    }
                }
                break;
            case Lagger:
                for (int i = 0; i < N1; ++i) {
                    for (int j = 0; j <= PolyOrder; ++j) {
                        PolyMat[i][j] = Lagger(Xj[i][indexQ], j);
                    }
                }
                break;
            case Lejandr:
                for (int i = 0; i < N1; ++i) {
                    for (int j = 0; j <= PolyOrder; ++j) {
                        PolyMat[i][j] = Lejandr(Xj[i][indexQ], j);
                    }
                }
                break;
            case Hermit:
                for (int i = 0; i < N1; ++i) {
                    for (int j = 0; j <= PolyOrder; ++j) {
                        PolyMat[i][j] = Hermit(Xj[i][indexQ], j);
                    }
                }
                break;
        }
        return PolyMat;

    }

    double[] MatrixShapeToVector(double[][]  Mat){
        int RowSize = Mat.length;
        int ColumnSize = Mat[0].length;
        int dim = Mat.length * Mat[0].length;
        double[] Vector = new double[dim];
        int ind = 0;

        for(int i = 0; i < RowSize; ++i){
            for(int j = 0; j < ColumnSize; ++j){
                Vector[ind] = Mat[i][j];
                ++ind;
            }
        }
        return Vector;
    }

    // !!!
    // return 1*3 array
    // ExDim  = {Rows, Columns, Rows * Columns}
    int[] ExtractDim(double[][] Mat){
        int[] ret = new int[3];
        ret[0] = Mat.length;
        ret[1] = Mat[0].length;
        ret[2] = ret[0]*ret[1];

        return ret;
    }

    // This is for new class overloading
    void MatrixForLinearEquation(double [][]X /*int PolyOrder, Polynom type, double lambda*/){
        Tp1 = PolynomMatrix(X, 0, PolyOrderJ[0], type);
        int[] DimTp1 = ExtractDim(Tp1);

        int RowSizeA = NumberOfSamples;
        int sum;
        int ColumnSizeA = DimTp1[2];

        A = new double[RowSizeA][ColumnSizeA];

        for(int i = 0; i < RowSizeA; ++i ) {
            Tp1 = PolynomMatrix(X, i, PolyOrderJ[0], type);
            sum = 0;
            for (int j = 0; j < DimTp1[2]; ++j) {
                A[i][sum] = MatrixShapeToVector(Tp1)[j];
                ++sum;

                //A_Lambda1[i][j] = MatrixShapeToVector(Tp1)[j];
            }
        }
    }

    void MatrixForLinearEquation(double[][] X1,double[][] X2,double[][] X3){

        // DimTpI is required this matrices by my initialization :)
        Tp1 = PolynomMatrix(X1, 0, PolyOrderJ[0], type);
        Tp2 = PolynomMatrix(X2,0, PolyOrderJ[1],type);
        Tp3 = PolynomMatrix(X3,0, PolyOrderJ[2],type);

        int[] DimTp1 = ExtractDim(Tp1);
        int[] DimTp2 = ExtractDim(Tp2);
        int[] DimTp3 = ExtractDim(Tp3);

        int RowSizeA = NumberOfSamples;
        int ColumnSizeA = DimTp1[2] + DimTp2[2] + DimTp3[2];


       A_Lambda1= new double [RowSizeA][DimTp1[2]];
        A_Lambda2= new double [RowSizeA][DimTp2[2]];
        A_Lambda3= new double [RowSizeA][DimTp3[2]];


        // !!!
        // Initialization of Matrix A
        A = new double[RowSizeA][ColumnSizeA];

        // first row initialization
        int sum = 0;

        for(int i = 0; i < RowSizeA; ++i ){
            Tp1 = PolynomMatrix(X1,i, PolyOrderJ[0],type);
            Tp2 = PolynomMatrix(X2,i, PolyOrderJ[1],type);
            Tp3 = PolynomMatrix(X3,i, PolyOrderJ[2],type);
            sum = 0;
            for(int j = 0; j < DimTp1[2]; ++j){
                A[i][sum] = MatrixShapeToVector(Tp1)[j];
                ++sum;

                A_Lambda1[i][j]=MatrixShapeToVector(Tp1)[j];
                }
            for(int j = 0; j < DimTp2[2]; ++j){
                A[i][sum] = MatrixShapeToVector(Tp2)[j];
                ++sum;

                A_Lambda2[i][j]=MatrixShapeToVector(Tp2)[j];
            }
            for(int j = 0; j < DimTp3[2]; ++j){
                A[i][sum] = MatrixShapeToVector(Tp3)[j];
                ++sum;

                A_Lambda3[i][j]=MatrixShapeToVector(Tp3)[j];
            }
        }
    }



    double ScalarF(double [][] X1, double[][] X2, double[][] X3, int indexQ){
        // X1, X2, X3 must be transposed, for example X1 has 2*45 dimensionality

        Tp1 = PolynomMatrix(X1,indexQ, PolyOrderJ[0],type);
        Tp2 = PolynomMatrix(X2,indexQ, PolyOrderJ[1],type);
        Tp3 = PolynomMatrix(X3,indexQ, PolyOrderJ[2],type);
        return ElementSum(ElemByElemMatrixMultiplication(Lambda1, Tp1)) +
               ElementSum(ElemByElemMatrixMultiplication(Lambda2, Tp2)) +
               ElementSum(ElemByElemMatrixMultiplication(Lambda3, Tp3));
    }
    private void LambdaMatrixCreation(double[] LambdaVec){
        int RowSize = Lambda1.length;
        int ColumnSize = Lambda1[0].length;
        int sum =0;

        for(int i = 0; i < RowSize; ++i){
            for(int j = 0; j < ColumnSize; ++j){
                Lambda1[i][j] = LambdaVec[sum];
                ++sum;
            }
        }

        RowSize = Lambda2.length;
        ColumnSize = Lambda2[0].length;

        for(int i = 0; i < RowSize; ++i){
            for(int j = 0; j < ColumnSize; ++j){
                Lambda2[i][j] = LambdaVec[sum];
                ++sum;
            }
        }

        RowSize = Lambda3.length;
        ColumnSize = Lambda3[0].length;

        for(int i = 0; i < RowSize; ++i){
            for(int j = 0; j < ColumnSize; ++j){
                Lambda3[i][j] = LambdaVec[sum];
                ++sum;
            }
        }

    }

    private double[] CalculateF_iUsing_ScalarF_Function(double [][] X1, double[][] X2, double[][] X3){
        double[] yCalc = new double[NumberOfSamples];
        for(int i = 0; i < NumberOfSamples; ++i){
            yCalc[i] = ScalarF(X1, X2, X3,i);
        }
        return yCalc;
    }

    private double[] PsiInitialization(double[][] X1, double[][] X2, double[][] X3,
                                       double[][] Lambda1, double[][] Lambda2, double [][] Lambda3,
                                       int indexQ,int NumberOfPsi){
        Psi1 = new double[X1.length];
        Psi2 = new double[X2.length];
        Psi3 = new double[X3.length];
        switch(NumberOfPsi) {
            case 0:
                Tp1 = PolynomMatrix(X1, indexQ, PolyOrderJ[0], type);
                for(int i = 0; i < Psi1.length;++i){
                    Psi1[i] = PsiElementCalc(Lambda1[i],Tp1[i]);
                }
                return Psi1;

            case 1:
                Tp2 = PolynomMatrix(X2, indexQ, PolyOrderJ[1], type);
                for(int i = 0; i < Psi2.length;++i){
                    Psi2[i] = PsiElementCalc(Lambda2[i],Tp2[i]);
                }
                return Psi2;

            case 2:
                Tp3 = PolynomMatrix(X3, indexQ, PolyOrderJ[2], type);
                for(int i = 0; i < Psi3.length;++i){
                    Psi3[i] = PsiElementCalc(Lambda3[i],Tp3[i]);
                }
                return Psi3;


        }
       return null;
    }

    public double[][] PsiMatrixCreation(double[][] X1,double[][] X2,double[][] X3,
                                        double[][] Lambda1, double[][] Lambda2, double [][] Lambda3,
                                        int NumberOfPsiMatrix){
        // we must solve three different systems on the second level
        // it return one of three matrix - defined by NumberOfPsiMatrix
        double[][] PsiMatrix;
        switch(NumberOfPsiMatrix){
            case 0:
                PsiMatrix = new double[NumberOfSamples][];
                for (int i = 0; i < NumberOfSamples; ++i){
                    PsiMatrix[i] = PsiInitialization(X1, X2, X3, Lambda1, Lambda2, Lambda3, i, 0);
                }
                return PsiMatrix;
            case 1:
                PsiMatrix = new double[NumberOfSamples][];
                for (int i=0; i < NumberOfSamples; ++i)
                    PsiMatrix[i] = PsiInitialization(X1,X2,X3,Lambda1,Lambda2,Lambda3,i,1);
                return PsiMatrix;
            case 2:
                PsiMatrix = new double[NumberOfSamples][];
                for(int i = 0; i < NumberOfSamples; ++i){
                    PsiMatrix[i] = PsiInitialization(X1, X2, X3, Lambda1, Lambda2, Lambda3, i, 2);
                }
                return PsiMatrix;
        }
        return null;
    }

    private double PsiElementCalc(double[] Lambda, double[] TpLambda){
        // TpLambda - Cut Tp
       
        double sum = 0;
        for(int i = 0; i < Lambda.length; ++i){
            sum = sum + Lambda[i]*TpLambda[i];
        }
        return sum;

    }
    private double[][] BigF_ForSystem(){
        double[][] Mat = new double[3][];

        Mat[0] = BigF1;
        Mat[1] = BigF2;
        Mat[2] = BigF3;
        return Transpose(Mat);
    }

    void A1A2A3_Extraction_from_AIncludeAll(){
        A1 = new double[PsMat1[0].length];
        A2 = new double[PsMat2[0].length];
        A3 = new double[PsMat3[0].length];

        int sum = 0;
        for(int i=0; i < PsMat1[0].length; ++i){
            A1[i] = AincludeAll[sum];
            ++sum;
        }
        for(int i=0; i < PsMat2[0].length; ++i){
            A2[i] = AincludeAll[sum];
            ++sum;
        }
        for(int i=0; i < PsMat3[0].length;++i){
            A3[i] = AincludeAll[sum];
            ++sum;
        }

    }
    double[][]  MakeOnePsiMatrix(){
        int ColumnSizePsMat1 = PsMat1[0].length;
        int ColumnSizePsMat2 = PsMat2[0].length;
        int ColumnSizePsMat3 = PsMat3[0].length;

        double[][] PsiBigMatrix = new double[NumberOfSamples][ColumnSizePsMat1 + ColumnSizePsMat2 + ColumnSizePsMat3];

        for(int i = 0; i < NumberOfSamples; ++i) {
            int sum = 0;
            for (int j = 0; j < ColumnSizePsMat1; ++j){
                PsiBigMatrix[i][sum] = PsMat1[i][j];
                ++sum;
            }
            for(int j =0; j < ColumnSizePsMat2; ++j){
                PsiBigMatrix[i][sum] = PsMat2[i][j];
                ++sum;
            }
            for(int j = 0; j < ColumnSizePsMat3; ++j){
                PsiBigMatrix[i][sum] = PsMat3[i][j];
                ++sum;
            }

        }
        return PsiBigMatrix;
    }

    double[] Recovered_Y(double[] Y_par){
        double[] Recov = new double[Y_par.length];
        coef=(MaxY - MinY);
        for(int i=0; i < Y_par.length; ++i){
            Recov[i] = Y_par[i]*coef+MinY;
        }
        return Recov;
    }

    private void BigF_initialization(){
        BigF1 = new double[NumberOfSamples];
        BigF2 = new double[NumberOfSamples];
        BigF3 = new double[NumberOfSamples];
        int ColumnSizePsiMat1 = PsMat1[0].length;
        int ColumnSizePsiMat2 = PsMat2[0].length;
        int ColumnSizePsiMat3 = PsMat3[0].length;
        // the sum from j1 = 0 to j1 = n1 is writing here from
        // from lecture 11 page 15
        double sum = 0;

        for(int i = 0; i < NumberOfSamples; ++i ){
            for(int j = 0; j < ColumnSizePsiMat1; ++j){
                sum = sum + A1[j]*PsMat1[i][j];
            }
            BigF1[i] = sum;
            sum = 0;

            for(int j = 0; j < ColumnSizePsiMat2; ++j){
                sum = sum + A2[j]*PsMat2[i][j];
            }
            BigF2[i] = sum;
            sum =0;

            for(int j = 0; j < ColumnSizePsiMat3; ++j){
                sum = sum + A3[j]*PsMat3[i][j];
            }
            BigF3[i] = sum;
            sum =0;
        }
    }

    public double[] SolveEquation(double[][] A, double []Y, double[] residual, double[] Ax){

        double[] Solution = new double[A[0].length];
        Method = new Conjugate_gradient_method(A, Y, A.length, A[0].length);

        Matrix AA = new Matrix(A, A.length, A[0].length);
        Matrix x = new Matrix(Method.X,Method.X.length);
        Matrix b = new Matrix(Y,Y.length);
        Matrix resid = b.minus(AA.times(x));
        Matrix AxComputed = AA.times(x);

        for(int i = 0; i < A[0].length; ++i){
            Solution[i] = x.get(i,0);
        }

        for (int i = 0; i < Y.length; ++i) {
            residual[i] = resid.get(i, 0);
            Ax[i] = AxComputed.get(i,0);
        }
        return Solution;
    }

    public double[] SolveEquation(double[][] A, double []Y, double[] residual, double[] Ax, double [] Y_real){

        double[] Solution = new double[A[0].length];
        Method = new Conjugate_gradient_method(A, Y, A.length, A[0].length);

        Matrix AA = new Matrix(A, A.length, A[0].length);
        Matrix x = new Matrix(Method.X,Method.X.length);
        Matrix b = new Matrix(Y,Y.length);
        Matrix resid = b.minus(AA.times(x));
        Matrix AxComputed = AA.times(x);

        for(int i = 0; i < A[0].length; ++i){
            Solution[i] = x.get(i,0);
        }

        for (int i = 0; i < Y.length; ++i) {
            residual[i] = resid.get(i, 0);
            Ax[i] = AxComputed.get(i,0);
        }
        return Solution;
    }

    // !!!
    // for regularization - we will solve equation (A'A+lambda*I)*x = A'y
    public double[] SolveEquation(double[][] A, double []Y, double lambda, double[] residual, double[] Ax){

        double[] Solution = new double[A[0].length];
        double[][] RegMat = new double[A[0].length][A[0].length];

        // Make Lambda in the each diagonal element of RegMat

        for(int i = 0; i < A[0].length; ++i ){
            RegMat[i][i] = lambda;
        }

        Method = new Conjugate_gradient_method(A, Y, RegMat, A.length, A[0].length);

        Matrix AA = new Matrix(A, A.length, A[0].length);
        Matrix x = new Matrix(Method.X,Method.X.length);
        Matrix b = new Matrix(Y,Y.length);
        Matrix resid = b.minus(AA.times(x));
        Matrix AxComputed = AA.times(x);

        for(int i = 0; i < A[0].length; ++i){
            Solution[i] = x.get(i,0);
        }

        for (int i = 0; i < Y.length; ++i) {
            residual[i] = resid.get(i, 0);
            Ax[i] = AxComputed.get(i,0);
        }
        return Solution;
    }


    // Predict function
    // Have an x  as a parameter - one row,
    // make the polynomial features, and after that will make scalar multiplication on vector
    // Lambda
    public double PredictNewX(double[] x){
        double[] PolyX = PolynomVector(x);
        double sum = 0;
        for(int i = 0; i < LambdaVec.length; ++i){
            sum = sum + LambdaVec[i]*PolyX[i];
        }
        return sum;
    }

    // overloading for next class
    // lambda is for regularization
    public Calculations(double[][] X, double[] Y, double lambda, Polynom type, int PolyOrder,
                        int NumSamples, double MaxY, double MinY)
    {
        PolyOrderJ = new int[1];
        PolyOrderJ[0] = PolyOrder;
        this.type = type;
        NumberOfSamples = NumSamples;
        this.MaxY = MaxY;
        this.MinY = MinY;
        Y_true = Recovered_Y(Y);
        Y_model = new double[NumSamples];
        // здесь формирует матрицу A
        MatrixForLinearEquation(X);
        ResidVecLambda = new double[Y.length];
        Ax = new double[Y.length];

        LambdaVec = SolveEquation(A,Y, lambda, ResidVecLambda,Ax);
        ResidVecLambda_max = 0;
        for (int i=0;i < ResidVecLambda.length;i++)
        {
            if (Math.abs(ResidVecLambda[i])>ResidVecLambda_max)
                ResidVecLambda_max=Math.abs(ResidVecLambda[i]);

        }
       // double[][] TransposedX = Transpose(X);
       // for(int i = 0; i < NumberOfSamples; ++i)
          //  Y_model[i] = PredictNewX(TransposedX[i]);

        residualBigF_System = ResidVecLambda;
        //residualBigF_System_max =
        Y_model = Recovered_Y(Ax);

    }

    public Calculations(double[][] X1, double [][] X2, double[][] X3, double[] Y, Polynom type, int[] PolyOrder,
                        int NumSamples, double MaxY, double MinY, boolean LambdaTogether, boolean SecondLevelTogether)
    {
        // X1, X2, X3, and maybe Y is transposed, for example X1 has 2*45 dimensionality
        PolyOrderJ = PolyOrder;
        this.type = type;
        NumberOfSamples = NumSamples;
        this.MaxY = MaxY;
        this.MinY = MinY;
        Y_true = Recovered_Y(Y);

        // matrix A formation
        MatrixForLinearEquation(X1, X2, X3);

        Ax = new double[Y.length];
        ResidVecLambda = new double[Y.length];

        if(!LambdaTogether){
            Lambda_one_matrix(X1,X2,X3,Y,SecondLevelTogether);
        }
        else{
            LambdaThreeMatrix(X1,X2,X3,Y,SecondLevelTogether);
        }

    }

    //public Calculation(double )
    public Calculations(double[][] X1, double [][] X2, double[][] X3, double[] Y, Polynom type, int[] PolyOrder,
                        int NumSamples, double MaxY, double MinY, boolean LambdaTogether, boolean SecondLevelTogether, double [] Y_real)
    {
        // X1, X2, X3, and maybe Y is transposed, for example X1 has 2*45 dimensionality
        PolyOrderJ = PolyOrder;
        this.type = type;
        NumberOfSamples = NumSamples;
        this.MaxY = MaxY;
        this.MinY = MinY;
        Y_true = Recovered_Y(Y_real);

        // matrix A formation
        MatrixForLinearEquation(X1, X2, X3);

        Ax = new double[Y.length];
        ResidVecLambda = new double[Y.length];

        if(!LambdaTogether){
            Lambda_one_matrix(X1,X2,X3,Y,SecondLevelTogether, Y_real);
        }
        else{
            LambdaThreeMatrix(X1,X2,X3,Y,SecondLevelTogether, Y_real);
        }

    }
    public void Lambda_one_matrix (double[][] X1, double [][] X2, double[][] X3, double[] Y, boolean SecondLevelTogether)
       {
           double[] LambdaVec = SolveEquation(A,Y,ResidVecLambda,Ax);
           ResidVecLambda_max=0;
           for (int i=0;i < ResidVecLambda.length;i++)
           {
               if (Math.abs(ResidVecLambda[i])>ResidVecLambda_max)
                   ResidVecLambda_max=Math.abs(ResidVecLambda[i]);
           }

           double[][] ResidVecA = new double[3][Y.length];
           double[] ResidForBigPsi = new double[Y.length];

           Lambda1 = new double[X1.length][PolyOrderJ[0]+1];
           Lambda2 = new double[X2.length][PolyOrderJ[1]+1];
           Lambda3 = new double[X3.length][PolyOrderJ[2]+1];

           LambdaMatrixCreation(LambdaVec);

           PsMat1 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,0);
           PsMat2 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,1);
           PsMat3 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,2);



           if(SecondLevelTogether) {
               PsMat = MakeOnePsiMatrix();
               AincludeAll = SolveEquation(PsMat, Y, ResidForBigPsi,Ax);
               A1A2A3_Extraction_from_AIncludeAll();
           }
           else{
               A1 = SolveEquation(PsMat1,Y, ResidVecA[0],Ax);
               A2 = SolveEquation(PsMat2,Y, ResidVecA[1],Ax);
               A3 = SolveEquation(PsMat3,Y, ResidVecA[2],Ax);
           }

           BigF_initialization();
           BigF_System = BigF_ForSystem();

           residualBigF_System = new double[Y.length];
           // double [3] C  = {c1, c2, c3}
           C = SolveEquation(BigF_System, Y, residualBigF_System,Ax);

           residualBigF_System_max=0;
           for (int i=0;i<residualBigF_System.length;i++)
           {
               if (Math.abs(residualBigF_System[i])>(residualBigF_System_max))
                   residualBigF_System_max=Math.abs(residualBigF_System[i]);
           }

           Y_model = Recovered_Y(Ax);

           // Here is situated the reshaping of lambda

           try{
               DataMapping.WriteToFile("Lambda1.txt", Lambda1);
               DataMapping.WriteToFile("Lambda2.txt", Lambda2);
               DataMapping.WriteToFile("Lambda3.txt", Lambda3);
               Conjugate_gradient_method.WriteToFile("ResidForBigPsi.txt", ResidForBigPsi);
              // DataMapping.WriteToFile("A.txt", A);
               DataMapping.WriteToFile("ResidVecASecondLevel.txt", ResidVecA);
               Conjugate_gradient_method.WriteToFile("Lambda_vec.txt", LambdaVec);
               Conjugate_gradient_method.WriteToFile("C_vec_2_level.txt", C);
               Conjugate_gradient_method.WriteToFile("Resid_Vec_Lambda.txt", ResidVecLambda);
               Conjugate_gradient_method.WriteToFile("ResidBig_Fsystem.txt", residualBigF_System);
           }
           catch (IOException e) {
               e.printStackTrace();
           }


       }

    public void LambdaThreeMatrix(double[][] X1, double [][] X2, double[][] X3, double[] Y, boolean SecondLevelTogether){

        for(int i = 1; i <= 3;++i){
            Lambda_three_matrix_initialize(X1,X2,X3,Y,i);
        }

        PsMat1 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,0);
        PsMat2 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,1);
        PsMat3 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,2);

        double[] ResidForBigPsi = new double[Y.length];
        double[][] ResidVecA = new double[3][Y.length];

        if(SecondLevelTogether) {
            try {
                PsMat = MakeOnePsiMatrix();
                AincludeAll = SolveEquation(PsMat, Y, ResidForBigPsi,Ax);
                A1A2A3_Extraction_from_AIncludeAll();
            }
            catch(java.lang.NullPointerException ex){
                JOptionPane.showMessageDialog(null,"Nullpointer");
            }

        }
        else{
            A1 = SolveEquation(PsMat1,Y, ResidVecA[0],Ax);
            A2 = SolveEquation(PsMat2,Y, ResidVecA[1],Ax);
            A3 = SolveEquation(PsMat3,Y, ResidVecA[2],Ax);
        }

        BigF_initialization();
        BigF_System = BigF_ForSystem();

        residualBigF_System = new double[Y.length];
        // double [3] C  = {c1, c2, c3}
        C = SolveEquation(BigF_System, Y, residualBigF_System,Ax);

        residualBigF_System_max=0;
        for (int i=0;i<residualBigF_System.length;i++)
        {
            if (Math.abs(residualBigF_System[i])>(residualBigF_System_max))
                residualBigF_System_max=Math.abs(residualBigF_System[i]);
        }

        Y_model = Recovered_Y(Ax);

        // Here is situated the reshaping of lambda

        try{
            Conjugate_gradient_method.WriteToFile("ResidForBigPsi.txt", ResidForBigPsi);
            // DataMapping.WriteToFile("A.txt", A);
            DataMapping.WriteToFile("ResidVecASecondLevel.txt", ResidVecA);
          //  Conjugate_gradient_method.WriteToFile("Lambda_vec.txt", LambdaVec);
            Conjugate_gradient_method.WriteToFile("C_vec_2_level.txt", C);
           // Conjugate_gradient_method.WriteToFile("Resid_Vec_Lambda.txt", ResidVecLambda);
            Conjugate_gradient_method.WriteToFile("ResidBig_Fsystem.txt", residualBigF_System);
        }
        catch (IOException e) {
            e.printStackTrace();
        }

    }



    public void Lambda_one_matrix (double[][] X1, double [][] X2, double[][] X3, double[] Y, boolean SecondLevelTogether, double [] Y_real)
    {
        double[] LambdaVec = SolveEquation(A,Y,ResidVecLambda,Ax, Y_real);
        ResidVecLambda_max=0;
        for (int i=0;i < ResidVecLambda.length;i++)
        {
            if (Math.abs(ResidVecLambda[i])>ResidVecLambda_max)
                ResidVecLambda_max=Math.abs(ResidVecLambda[i]);
        }

        double[][] ResidVecA = new double[3][Y.length];
        double[] ResidForBigPsi = new double[Y.length];

        Lambda1 = new double[X1.length][PolyOrderJ[0]+1];
        Lambda2 = new double[X2.length][PolyOrderJ[1]+1];
        Lambda3 = new double[X3.length][PolyOrderJ[2]+1];

        LambdaMatrixCreation(LambdaVec);

        PsMat1 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,0);
        PsMat2 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,1);
        PsMat3 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,2);



        if(SecondLevelTogether) {
            PsMat = MakeOnePsiMatrix();
            AincludeAll = SolveEquation(PsMat, Y_real, ResidForBigPsi,Ax);
            A1A2A3_Extraction_from_AIncludeAll();
        }
        else{
            A1 = SolveEquation(PsMat1,Y_real, ResidVecA[0],Ax);
            A2 = SolveEquation(PsMat2,Y_real, ResidVecA[1],Ax);
            A3 = SolveEquation(PsMat3,Y_real, ResidVecA[2],Ax);
        }

        BigF_initialization();
        BigF_System = BigF_ForSystem();

        residualBigF_System = new double[Y.length];
        // double [3] C  = {c1, c2, c3}
        C = SolveEquation(BigF_System, Y_real, residualBigF_System,Ax);

        residualBigF_System_max=0;
        for (int i=0;i<residualBigF_System.length;i++)
        {
            if (Math.abs(residualBigF_System[i])>(residualBigF_System_max))
                residualBigF_System_max=Math.abs(residualBigF_System[i]);
        }

        Y_model = Recovered_Y(Ax);

        // Here is situated the reshaping of lambda

        try{
            DataMapping.WriteToFile("Lambda1.txt", Lambda1);
            DataMapping.WriteToFile("Lambda2.txt", Lambda2);
            DataMapping.WriteToFile("Lambda3.txt", Lambda3);
            Conjugate_gradient_method.WriteToFile("ResidForBigPsi.txt", ResidForBigPsi);
            // DataMapping.WriteToFile("A.txt", A);
            DataMapping.WriteToFile("ResidVecASecondLevel.txt", ResidVecA);
            Conjugate_gradient_method.WriteToFile("Lambda_vec.txt", LambdaVec);
            Conjugate_gradient_method.WriteToFile("C_vec_2_level.txt", C);
            Conjugate_gradient_method.WriteToFile("Resid_Vec_Lambda.txt", ResidVecLambda);
            Conjugate_gradient_method.WriteToFile("ResidBig_Fsystem.txt", residualBigF_System);
        }
        catch (IOException e) {
            e.printStackTrace();
        }


    }

    public void LambdaThreeMatrix(double[][] X1, double [][] X2, double[][] X3, double[] Y, boolean SecondLevelTogether, double[] Y_real){

        for(int i = 1; i <= 3;++i){
            Lambda_three_matrix_initialize(X1,X2,X3,Y,i);
        }

        PsMat1 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,0);
        PsMat2 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,1);
        PsMat3 = PsiMatrixCreation(X1,X2,X3,Lambda1, Lambda2,Lambda3,2);

        double[] ResidForBigPsi = new double[Y.length];
        double[][] ResidVecA = new double[3][Y.length];

        if(SecondLevelTogether) {
            try {
                PsMat = MakeOnePsiMatrix();
                AincludeAll = SolveEquation(PsMat, Y_real, ResidForBigPsi,Ax);
                A1A2A3_Extraction_from_AIncludeAll();
            }
            catch(java.lang.NullPointerException ex){
                JOptionPane.showMessageDialog(null,"Nullpointer");
            }

        }
        else{
            A1 = SolveEquation(PsMat1,Y_real, ResidVecA[0],Ax);
            A2 = SolveEquation(PsMat2,Y_real, ResidVecA[1],Ax);
            A3 = SolveEquation(PsMat3,Y_real, ResidVecA[2],Ax);
        }

        BigF_initialization();
        BigF_System = BigF_ForSystem();

        residualBigF_System = new double[Y.length];
        // double [3] C  = {c1, c2, c3}
        C = SolveEquation(BigF_System, Y_real, residualBigF_System,Ax);

        residualBigF_System_max=0;
        for (int i=0;i<residualBigF_System.length;i++)
        {
            if (Math.abs(residualBigF_System[i])>(residualBigF_System_max))
                residualBigF_System_max=Math.abs(residualBigF_System[i]);
        }

        Y_model = Recovered_Y(Ax);

        // Here is situated the reshaping of lambda

        try{
            Conjugate_gradient_method.WriteToFile("ResidForBigPsi.txt", ResidForBigPsi);
            // DataMapping.WriteToFile("A.txt", A);
            DataMapping.WriteToFile("ResidVecASecondLevel.txt", ResidVecA);
            //  Conjugate_gradient_method.WriteToFile("Lambda_vec.txt", LambdaVec);
            Conjugate_gradient_method.WriteToFile("C_vec_2_level.txt", C);
            // Conjugate_gradient_method.WriteToFile("Resid_Vec_Lambda.txt", ResidVecLambda);
            Conjugate_gradient_method.WriteToFile("ResidBig_Fsystem.txt", residualBigF_System);
        }
        catch (IOException e) {
            e.printStackTrace();
        }

    }






    public void Lambda_three_matrix_initialize (double[][] X1, double [][] X2, double[][] X3, double[] Y, int i)
    {
        double [][]A_Lambda_i={{0}};
        switch(i){
            case 1:
                A_Lambda_i = A_Lambda1.clone();
                break;

            case 2:
                A_Lambda_i = A_Lambda2.clone();
                break;

            case 3:
                A_Lambda_i = A_Lambda3.clone();
                break;

        }




        Method = new Conjugate_gradient_method(A_Lambda_i, Y, A_Lambda_i.length, A_Lambda_i[0].length);

        Matrix AA = new Matrix(A_Lambda_i, A_Lambda_i.length, A_Lambda_i[0].length);
        Matrix x = new Matrix(Method.X, Method.X.length);
        Matrix b = new Matrix(Y, Y.length);
        Matrix resid = b.minus(AA.times(x));

        double[] ResidVecLambda = new double[Y.length];
        //  double[] ResidMatLambda = new double[Y.length];
        double[] LambdaVec = Method.X;
        for (int j = 0; j < Y.length; ++j) {
            ResidVecLambda[j] = resid.get(j, 0);
        }

        // Here is situated the reshaping of lambda

        int RowSize = 0;
        int ColumnSize = 0;

        switch(i){
            case 1:
                RowSize = X1.length;
                ColumnSize = PolyOrderJ[0] + 1;
                break;

            case 2:
                RowSize = X2.length;
                ColumnSize = PolyOrderJ[1] + 1;
                break;
            case 3:
                RowSize = X3.length;
                ColumnSize = PolyOrderJ[2] + 1;
                break;
        }
        double [][] Lambda = new double[RowSize][ColumnSize];

        int sum =0;

        for(int k = 0; k < RowSize; ++k){
            for(int j = 0; j < ColumnSize; ++j){
                Lambda[k][j] = LambdaVec[sum];
                ++sum;
            }
        }


        switch(i){
            case 1:
                Lambda1 = Lambda.clone();
                break;

            case 2:
                Lambda2 = Lambda.clone();
                break;

            case 3:
                Lambda3 = Lambda.clone();
                break;

        }


        try {
            DataMapping.WriteToFile("Lambda"+Integer.toString(i)+"_three_matrix.txt", Lambda);
            Conjugate_gradient_method.WriteToFile("Lambda"+Integer.toString(i)+"_three_vec.txt", LambdaVec);
            Conjugate_gradient_method.WriteToFile("Resid_Vec_Lambda"+Integer.toString(i)+"_three.txt", ResidVecLambda);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


        // ------------- This Code is for Checking the function that the bottom
        // They works well
        // please, dont remove this comments ;)

        /*
        double[] Fi = CalculateF_iUsing_ScalarF_Function(X1,X2,X3);
        for(int i = 0; i < Y.length; ++i){
            ResidMatLambda[i] = Y[i] - Fi[i];
        }
        // first row by LambdaVec
        // Second By LambdaMat
        double[][] Residual = new double[2][Y.length];
        Residual[0] = ResidVecLambda;
        Residual[1] = ResidMatLambda;
        try{
            DataMapping.WriteToFile("Residual_within_lambda.txt", Residual);
            DataMapping.WriteToFile("A.txt",A);
        }
        catch (IOException e) {
            e.printStackTrace();
        }*/





    }









