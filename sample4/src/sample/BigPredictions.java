package sample;
import java.io.*;
import java.util.Scanner;
/**
 * Created by Gregory on 23.11.2014.
 */
public class BigPredictions {
    private
        Calculations OnlyForPrediction;
        int NumberOfSamples;
        int WindowSize;
        int PredictLength;
        double lambda = 0;
        int RealIndex;

        // Sizes X1, X2, X3, Y
        int[] sizeVars;
        int PolyOrder;
        double [][][] VarWindows;
       /* double[][] X1Window;
        double[][] X2Window;
        double[][] X3Window;
        double[][] YWindow;*/

    public BigPredictions(double[][] X, int[][] IndexesOfAllVar, double Lambda, int NumberSam, int PolyOrder,
                          int WindowSize, int LastRealIndex, int PredictLength)
    {
        this.PredictLength = PredictLength;
        NumberOfSamples = NumberSam;
        sizeVars = new int[IndexesOfAllVar.length];
        for(int i = 0; i < sizeVars.length; ++i){
            sizeVars[i] = IndexesOfAllVar[i].length;
        }
        this.WindowSize = WindowSize;
        this.PolyOrder = PolyOrder;

        try {
            DataMapping.WriteToFile("X.txt",DataMapping.Transpose(X));
        }
        catch (IOException e) {
            e.printStackTrace();
        }

        /*for(int i = 0; i < SizesOfAllVar.length; ++i){
            sizeVars[i] = SizesOfAllVar[i];
        }
        /*X1Window = new double[sizeVars[0]][WindowSize];
        X2Window = new double[sizeVars[1]][WindowSize];
        X3Window = new double[sizeVars[2]][WindowSize];
        YWindow = new double[sizeVars[3]][WindowSize];*/
        VarWindows = new double[sizeVars.length][][];
       // for(int i=0; i < sizeVars.length; ++i){
       //     VarWindows[i] = new double[sizeVars[i]][WindowSize + PredictLength];
       // }
        if(Lambda >=0){
            this.lambda = Lambda;
            SetFirstWindows(X, LastRealIndex, IndexesOfAllVar);
        }


    }
    public double[] MaxVersusRows(double[][]Y){
        int RowSize = Y.length;
        double[] Max = new double[RowSize];
        for(int i = 0; i < RowSize; ++i ){
            Max[i] = DataMapping.Max(Y[i]);
        }
        return Max;
    }
    public double[] MinVersusRows(double[][]Y){
        int RowSize = Y.length;
        double[] Min = new double[RowSize];
        for(int i = 0; i < RowSize; ++i ){
            Min[i] = DataMapping.Min(Y[i]);
        }
        return Min;
    }
    public double[][] ReturnArrayWithSpecializedIndexes(double[][] X, int[] Ind){
        double[][] Dst = new double[Ind.length][WindowSize + PredictLength];
        int sum = 0;
        for(int i = 0; i < X.length; ++i){
            for(int j = 0; j < Ind.length; ++j){
                if(i == Ind[j]){
                    int IndexConstraint = Dst[0].length;
                    for(int d = 0; d < IndexConstraint; ++d){
                        Dst[sum][d]= X[i][RealIndex + PredictLength - IndexConstraint + 1 + d];
                    }
                    sum++;
                }
            }
        }
        return Dst;
    }
    public  void SetFirstWindows(double [][]X, int LastRealIndex, int[][] IndexesOfAllVar/*, double [][]*/){
        if (LastRealIndex > WindowSize && LastRealIndex <= NumberOfSamples - PredictLength){
            int FirstInd = LastRealIndex - WindowSize + 1;
            RealIndex = LastRealIndex;
            int HelpLastIndex = LastRealIndex + PredictLength;
            int sum = 0;
            for(int l =0; l < VarWindows.length; ++l){
                VarWindows[l] = ReturnArrayWithSpecializedIndexes(X,IndexesOfAllVar[l]);
                String s = "VarWindow_" + String.valueOf(l)+".txt";
                double[][] Temp = DataMapping.Transpose(VarWindows[l]);
                try {
                    DataMapping.WriteToFile(s, Temp);
                }
                catch (IOException e) {
                    e.printStackTrace();
                }
            }

            /*for(int l =0; l < VarWindows.length; ++l) {
                for (int j = 0; j < X.length; ++j) {
                    if(j == IndexesOfAllVar[l][j]) {
                        for (int i = FirstInd; i <= HelpLastIndex; ++i) {
                            VarWindows[l][sum][i - FirstInd] = X[j][i];
                        }
                        ++sum;
                    }
                }
            }*/

           /* for(int j = 0; j < sizeVars[1]; ++j) {
                for (int i = FirstInd; i <= LastRealIndex; ++i) {
                    VarWindows[1][j][i] = X[sizeVars[0]+j][i];
                }
            }

            for(int j = 0; j < sizeVars[2]; ++j) {
                for (int i = FirstInd; i <= LastRealIndex; ++i) {
                    VarWindows[2][j][i] = X[sizeVars[0]+sizeVars[1]+j][i];
                }
            }

            for(int j = 0; j < sizeVars[3]; ++j) {
                for (int i = FirstInd; i <= LastRealIndex; ++i) {
                    VarWindows[3][j][i] = X[sizeVars[0]+sizeVars[1]+sizeVars[2]+j][i];
                }
            }*/

        }
    }

    public double[][] CopyLeftUpperPartOfMatrix(double[][] Mat, int RowSize, int ColumnSize){
        if(Mat.length >= RowSize && Mat[0].length >=ColumnSize) {
            double[][] Src = new double[RowSize][ColumnSize];
            for (int i = 0; i < RowSize; ++i) {
                for (int j = 0; j < ColumnSize; ++j) {
                    Src[i][j] = Mat[i][j];
                }
            }
            return Src;
        }
        else return null;
    }

    public double[] CopyUpperPartOfVector(double[] Vec, int Size){
        if(Vec.length >= Size) {
            double[] Src = new double[Size];
            for (int i = 0; i < Size; ++i) {
                Src[i] = Vec[i];
            }
            return Src;
        }
        return null;
    }

    public double[][] getPredictedY_i()
    {
        if (RealIndex != 0) {
            double[][] PredictY = new double[sizeVars[3]][PredictLength];
            int SizeVarsLength = sizeVars.length;
            int XVariablesNumber = SizeVarsLength-1;
            int YVariableIndex = SizeVarsLength-1;

            //OnlyForPrediction = new Calculations[sizeVars[3]];
            for (int i = 0; i < XVariablesNumber; ++i){
                double[][] CurrentMat = CopyLeftUpperPartOfMatrix(VarWindows[i],VarWindows[i].length, WindowSize);
                double[] CurrentVec = CopyUpperPartOfVector(VarWindows[YVariableIndex][i], WindowSize);
                double MaxOfCurrentVec = DataMapping.Max(CurrentVec);
                double MinOfCurrentVec = DataMapping.Min(CurrentVec);

                OnlyForPrediction = new Calculations(CurrentMat, CurrentVec,
                        lambda,Polynom.Chebyshev, PolyOrder, WindowSize, MaxOfCurrentVec,
                        MinOfCurrentVec);

                double[][] Rows_With_X = DataMapping.Transpose(VarWindows[i]);
                for(int j = 0; j < PredictLength; ++j ){
                    PredictY[i][j] = OnlyForPrediction.PredictNewX(Rows_With_X[WindowSize+j]);
                }
            }
                return PredictY;
        }
        return null;
    }
}
