package sample; /**
 * Created by Gregory on 12.10.2014.
 */
import java.io.*;
import java.util.Scanner;

public class DataMapping {
    double [][] NormalizedData;
    double [][] RealData;
    int DimX1;
    int DimX2;
    int DimX3;
    int DimY;
    int NumSamples;
    double[] MaxArray;
    double[] MinArray;

    public static double Max(double[] array){
       int len = array.length;
       double max = array[0];
       for(int i = 1; i < len; ++i ){
           if (max < array[i]) max = array[i];
       }
       return  max;
    }
    public static double Min(double[] array){
        int len = array.length;
        double min = array[0];
        for(int i = 1; i < len; ++i ){
            if (min > array[i]) min = array[i];
        }
        return  min;
    }
    public static double [][] Transpose(double [][] Mat){
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
    public double[][] Normalize(double[][] Mat){
        int RowSize = Mat.length;
        int ColumnSize = Mat[0].length;
        MaxArray = new double[ColumnSize];
        MinArray = new double[ColumnSize];

        //We will return this matrix, but it will not be transposed
        // we want maximum within fixed column
        double[][] HelpMat = Transpose(Mat);

        // Process of normalization
        for(int i = 0; i < ColumnSize; ++i){
            // The max min versus columns is needed here.
            MaxArray[i] = Max(HelpMat[i]);
            MinArray[i] = Min(HelpMat[i]);
            for(int j = 0; j < RowSize; ++j ){
                HelpMat[i][j] = (HelpMat[i][j] - MinArray[i])/ (MaxArray[i] - MinArray[i]);
            }
        }
        return Transpose(HelpMat);
    }
    public static void WriteToFile(String name, double[][] Mat)throws FileNotFoundException, IOException {
        PrintWriter fout = new PrintWriter(new FileWriter(name));
        //fout.write(this.getSize());
        int ColumnLength =  Mat[0].length;
        int RowLength = Mat.length;

        for ( int i = 0; i < RowLength; ++i ) {
            for ( int j = 0; j < ColumnLength; ++j ) {
                fout.print(String.format("%.4f", Mat[i][j]));
                fout.write("  ");
            }
            fout.write("\n");
        }
        fout.close();
    }
    public DataMapping(String FileName, int DimX1, int DimX2, int DimX3, int DimY, int FileLength ){
        this.NumSamples = FileLength;
        this.DimX1 = DimX1;
        this.DimX2 = DimX2;
        this.DimX3 = DimX3;
        this.DimY = DimY;
        int ColumnSize = DimX1 + DimX2 + DimX3 + DimY;
        int RowSize = FileLength;
        String str;
        NormalizedData = new double[RowSize][ColumnSize];
        RealData = new double[RowSize][ColumnSize];
        try{
            File InFile = new File(FileName);
            Scanner s = new Scanner(InFile);
            for(int i=0; i < RowSize; ++i) {
                for (int j = 0; j < ColumnSize; ++j) {
                    //str = s.next();
                   // NormalizedData[i][j] = Double.parseDouble(str);
                    NormalizedData[i][j] = s.nextDouble();
                    RealData[i][j]=NormalizedData[i][j];
                    //s.next().toDouble
                }
                //a1[i] = s.nextInt();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        NormalizedData = Normalize(NormalizedData);
        try{
            WriteToFile("NormData.txt", NormalizedData);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }




    public DataMapping(String FileName1, String FileName2,int DimX1, int DimX2, int DimX3, int DimY, int FileLength ){
        this.NumSamples = FileLength;
        this.DimX1 = DimX1;
        this.DimX2 = DimX2;
        this.DimX3 = DimX3;
        this.DimY = DimY;
        int ColumnSize = DimX1 + DimX2 + DimX3 + DimY;
        int RowSize = FileLength;
        String str;
        NormalizedData = new double[RowSize][ColumnSize];
        RealData = new double[RowSize][ColumnSize];
        try{
            File InFile1 = new File(FileName1);
            Scanner s1 = new Scanner(InFile1);
            File InFile2 = new File(FileName2);
            Scanner s2 = new Scanner(InFile2);
            for(int i=0; i < RowSize; ++i) {
                for (int j = 0; j < ColumnSize; ++j) {

                     if (j<ColumnSize-1) {
                         NormalizedData[i][j] = s1.nextDouble();
                         RealData[i][j] = NormalizedData[i][j];
                     }
                    else
                     {
                         NormalizedData[i][j] = s2.nextDouble();
                         RealData[i][j] = NormalizedData[i][j];
                     }
                }
                //a1[i] = s.nextInt();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        NormalizedData = Normalize(NormalizedData);
        try{
            WriteToFile("NormData.txt", NormalizedData);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }



}
