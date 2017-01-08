package sample;

/**
 * Created by Mark on 23.11.2014.
 */
public class Robust {
    int[][] condition;
    int length_of_window;
    int step_now;
    int new_length_of_sample;
    double [][] X_new;
    int length_of_xy;
    int[][]fail;
String S;
    Robust(int length_of_xy1, int length_of_window1, int step_now1) {

        length_of_xy=length_of_xy1;
        if (length_of_window1 < step_now1) {
            condition = new int[length_of_window1][length_of_xy];
            for (int i = 0; i < length_of_window1; i++)
                for (int j = 0; j < length_of_xy; j++)
                    condition[i][j] = 0;
        }
        else {
            condition = new int[step_now1+1][length_of_xy];
            for (int i = 0; i < step_now1; i++)
                for (int j = 0; j < length_of_xy; j++)
                    condition[i][j] = 0;
        }
        length_of_window = length_of_window1;
        step_now = step_now1;
        new_length_of_sample=length_of_window1;

        S="";
    }

    double [][] calculation_of_condition(double[][] X, double[][] X_norm) {



        if (length_of_window > step_now) {
            for (int i = 1; i < step_now - 1; i++) {
                //                   if ((Math.abs(X[i - 1][0] + X[i + 1][0]) / 2 - X[i][0]) > 1) condition[i][0] = 11;
                if (X[i][2] > 450) condition[i][2] = 21;
                if ((X[i][2] > 1) && (X[i][2] < 50)) condition[i][1] = 22;
                if ((X[i][3] > 0) && (X[i][3] < 1100)) condition[i][2] = 31;


                if (Math.abs((X[i - 1][4] + X[i + 1][4]) / 2 - X[i][4]) > 300) condition[i][3] = 41;

                if (X[i][5] > 12.8) condition[i][5] = 51;
                if (X[i][5] < 11.4) condition[i][5] = 52;

                if (Math.abs((X[i - 1][6] + X[i + 1][6]) / 2 - X[i][6]) > 1) condition[i][5] = 61;
                if (X[i][7] > 12.5) condition[i][7] = 71;
                if (X[i][7] < 11.4) condition[i][7] = 72;
            }

        } else {
            for (int i = step_now - length_of_window + 1; i < step_now - 1; i++) {
                //                   if ((Math.abs(X[i - 1][0] + X[i + 1][0]) / 2 - X[i][0]) > 1) condition[i][0] = 11;
                if (X[i][2] > 450)
                    condition[i-step_now + length_of_window-1][1] = 21;
                if ((X[i][2] > 1) && (X[i][2] < 50)) condition[i-step_now + length_of_window-1][1] = 22;
                if ((X[i][3] > 0) && (X[i][3] < 1100)) condition[i-step_now + length_of_window-1][2] = 31;


                if (Math.abs((X[i - 1][4] + X[i + 1][4]) / 2 - X[i][4]) > 300) condition[i-step_now + length_of_window][3] = 41;

                if (X[i][5] > 13) condition[i-step_now + length_of_window-1][4] = 51;
                if (X[i][5] < 10) condition[i-step_now + length_of_window-1][4] = 52;

                if (Math.abs((X[i - 1][6] + X[i + 1][6]) / 2 - X[i][6]) > 1) condition[i-step_now + length_of_window-1][5] = 61;
                if (X[i][7] > 13) condition[i-step_now + length_of_window-1][6] = 71;
                if (X[i][7] < 10) condition[i-step_now + length_of_window-1][6] = 72;
            }
        }

        int sum = 0;
/*
for (int i=0;i<condition.length;i++) {
    sum = 0;
    for (int j = 0; j < condition[0].length; j++) {
        sum = sum + condition[i][j];}
        if (sum != 0)
            new_length_of_sample--;


}
        int z=0;
        X_new= new double[new_length_of_sample][length_of_xy];
        for (int i=0;i<condition.length;i++) {
            sum = 0;
            for (int j = 0; j < condition[0].length; j++) {
                sum = sum + condition[i][j];}
                if (sum == 0)
                {
                  for (int k=0;k<X[0].length;k++) {
                    X_new[z][k] = X_norm[i+step_now][k];}
                    z++;


                }

            }


*/

        int leg=(int)(Math.ceil(this.length_of_window / 10));
        double[] test = new double[this.length_of_window + leg*X[0].length];
        fail = new int[X[0].length][leg];
//double [] test= new double[]{2, 3, 4, 4, 5, 5, 6, 6, 6, 80};
        for (int num=1;num<X[0].length;num++ )
        {
        for (int i = 1; i < test.length; i++)
            test[i] = X[step_now - test.length + i][num];


        double max_x = 0, min_x = 10000000, SPWE = 0, sum1 = 0, sum2 = 0, dist = 0;
        for (int i = 0; i < test.length; i++) {
            if (test[i] < min_x)
                min_x = test[i];
            if (test[i] > max_x)
                max_x = test[i];
        }

        for (int i = 0; i < test.length; i++)
            for (int j = test.length-1; j > i; j--) {
                sum1 += (test[i] + test[j]) * Math.cos((Math.PI / 2) * (test[i] - test[j]) / (max_x - min_x));
                sum2 += Math.cos((Math.PI / 2) * (test[i] - test[j]) / (max_x - min_x));
            }

        SPWE = sum1 / sum2/2;
        boolean p=false;


        for (int i = 0; i != fail[0].length; i++) {
            dist = 0;


            for (int j = 0; j < test.length; j++) {
                p=false;
                if (dist < Math.abs(test[j] - SPWE)) {
      //              double o=Math.abs(test[j] - SPWE);
                    for (int k=0;k<i;k++) {
                        if (j==fail[num][k])
                            p=true;}
                        if (p==false) {
                    dist = Math.abs(test[j] - SPWE);
                    fail[num][i] = j;

                    }
                }
        }
        }
    }

            return X_new;
        }

    String Massage ()
    {
        S="";
        for (int i=0;i<condition.length;i++)
            for (int j=0;j<condition[0].length;j++)
            {
                switch (condition[i][j]) {
                    case 11:
                        S=S+"Problem with Acc. sensor#1\n";
                        break;
                    case 21:
                        S=S+"Problem with crankshaft sensor#1\n";
                        break;
                    case 22:
                        S=S+"Problem with crankshaft sensor#2\n";
                        break;
                    case 31:
                        S=S+"Problem with extra generator voltage sensor#1\n";
                        break;
                    case 41:
                        S=S+"Problem total voltage sensor#1\n";
                        break;

                    case 51:
                        S=S+"Problem grid voltage sensor#1\n";
                        break;
                    case 52:
                        S=S+"Problem grid voltage sensor#2\n";
                        break;
                    case 61:
                        S=S+"Problem fuel sensor#1\n";
                        break;
                    case 71:
                        S=S+"Problem acc. voltage sensor#1\n";
                        break;
                    case 72:
                        S=S+"Problem acc. voltage sensor#2\n";
                        break;


                }
            }
        if (S=="") S="Data is fine\n";
        return S;
    }

}