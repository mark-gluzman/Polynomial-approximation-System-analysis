package sample;

/**
 * Created by Mark on 25.11.2014.
 */

public class Estimation_of_situation {

    double[][] persent;
    double[]  function_percent;
    double level_crush;
    double[] level_est_contingency;
    double [] level_est_critical;

Estimation_of_situation(int step) {
    level_est_contingency = new double[]{11.7, step * 0.003, 11.5};
    level_est_critical = new double[]{10.5, 1, 10};
}
    double [] crucial (double Y[][], int prediction_length) {
        persent = new double[3][prediction_length];
        function_percent=new double[prediction_length];

        double max_f=0;
        for (int i = 0; i < prediction_length; i++){
            for (int j = 0; j < 3; j++)
            {
                if ((Y[j][i] < level_est_contingency[j])&&(Y[j][i]>level_est_critical[j]))
                    persent[j][i] = Math.abs(Y[j][i] - level_est_contingency[j]) / Math.abs(level_est_contingency[j] - level_est_critical[j]);
                else persent[j][i] = 0;
            }
                function_percent[i] = 1 - (1 - persent[0][i]) * (1 - persent[1][i]) * (1 - persent[2][i]);
            if (max_f<function_percent[i]) max_f=function_percent[i];
            }
        level_crush=max_f;

        return  function_percent;
    }

}

