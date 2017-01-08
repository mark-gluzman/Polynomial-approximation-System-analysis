package sample;


import javafx.beans.property.SimpleStringProperty;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.scene.chart.AreaChart;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.*;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.util.Callback;

import javax.swing.*;
import java.io.*;
import java.util.Arrays;

public class ControllerSample {


    public Polynom type_now;

    public TextField text_file_input;
    public TableView table_data_input;
    public TextField box_y;
    public TextField box_x2;
    public TextField box_x3;
    public TextField box_x1;
    public TextField box_length;
    public TextField Text_res;

    public double[][] X1;
    public double[][] X2;
    public double[][] X3;
    public double[][] Y;
    public double[] MaxY;
    public double[] MinY;
    public RadioButton Radio1;
    public RadioButton Radio2;
    public RadioButton Radio3;
    public RadioButton Radio4;

    public TextField box_p1;
    public TextField box_p2;
    public TextField box_p3;
    public TableView data_real;
    public LineChart y0_first_level;
    public LineChart Y0_last_level;
    public TextField Y_order;
    public ToggleButton SecondLevelTogether;
    public LineChart Graph_Y;
    public TextArea Text_results;
    public TextField Text_Y;
    public Initial_data data_in;
    public Button box_down_main;
    public Button box_down_y;
    public DataMapping DaMa;
    public Calculations LambdaDefinition;
    public TextField Text_err1;
    public TextField Text_err_res;
    public TextField Text_err_res_recover;
    public LineChart Graph_Y_norm;
    public ToggleButton lambda3Systems;
    public ToggleButton ButArith_mean;
    public TextArea Text_robust;

    public Button But_fore;
    public TextField Text_time;
    public LineChart Graph_grid;
    public LineChart Graph_acc;
    public LineChart Graph_fuel;
    public Slider SliderTime;
    public TextField Text_window;
    public TextField Text_prediction;
    public TextField Text_fuel;
    public TextField Text_grid;
    public TextField Text_acc;
    public ProgressBar Progress3;
    public ProgressBar Progress1;
    public ProgressBar Progress2;
    public AreaChart Graph_level;
    public TextArea Text_foR_for;

    public void open_file_input(javafx.event.ActionEvent actionEvent) {


        final JFileChooser fileopen = new JFileChooser();
        int ret = fileopen.showOpenDialog(null);
        if (ret == JFileChooser.APPROVE_OPTION) {
            File file_input = fileopen.getSelectedFile();
            text_file_input.setText(file_input.getPath());
        }

    }


    public void test(ActionEvent actionEvent) {

        if (Radio1.selectedProperty().getValue()) {
            type_now = Polynom.Chebyshev;
        }

        if (Radio2.selectedProperty().getValue()) {
            type_now = Polynom.Lejandr;
        }
        if (Radio3.selectedProperty().getValue()) {
            type_now = Polynom.Lagger;
        }
        if (Radio4.selectedProperty().getValue()) {
            type_now = Polynom.Hermit;
        }
        data_in.add(type_now, Integer.parseInt(box_p1.getText()), Integer.parseInt(box_p2.getText()), Integer.parseInt(box_p3.getText()));

        int Y_ord = Integer.parseInt(Y_order.getText()) - 1;
        int[] Order = {data_in.p1, data_in.p2, data_in.p3};


        double[] AverageY = new double[data_in.length_of_sample];
        double help_max = 0, help_min = 0;
        for (int i = 0; i < data_in.length_of_sample; i++) {
            help_max = help_min = Y[0][i];

            for (int j = 0; j < data_in.length_of_y; j++) {
                if (Y[j][i] > help_max)
                    help_max = Y[j][i];
                if (Y[j][i] < help_min)
                    help_min = Y[j][i];

            }
            AverageY[i] = (help_max + help_min) / 2;
        }


        if (ButArith_mean.selectedProperty().getValue())
            LambdaDefinition = new Calculations(X1, X2, X3, AverageY, data_in.type, Order, Y[Y_ord].length, MaxY[Y_ord], MinY[Y_ord], lambda3Systems.selectedProperty().getValue(), SecondLevelTogether.selectedProperty().getValue(), Y[Y_ord]);
        else
            LambdaDefinition = new Calculations(X1, X2, X3, Y[Y_ord], data_in.type, Order, Y[Y_ord].length, MaxY[Y_ord], MinY[Y_ord], lambda3Systems.selectedProperty().getValue(), SecondLevelTogether.selectedProperty().getValue());


        double Y_error_rec = 0;
        for (int i = 0; i < LambdaDefinition.Y_model.length; i++) {
            if (Math.abs(LambdaDefinition.Y_model[i] - LambdaDefinition.Y_true[i]) > Y_error_rec)
                Y_error_rec = Math.abs(LambdaDefinition.Y_model[i] - LambdaDefinition.Y_true[i]);
        }

        Text_err_res_recover.setText(Double.toString(Y_error_rec));
        Text_err1.setText(Double.toString(LambdaDefinition.ResidVecLambda_max));

        Text_err_res.setText(Double.toString(LambdaDefinition.residualBigF_System_max));


        Write_Results(Text_res.getText(), LambdaDefinition, data_in, Y_error_rec);
        Write_to_text(Text_res.getText());

    }


    public void DistributeData(double[][] Mat, int DimX1, int DimX2, int DimX3, int DimY, int NumSamples) {

        X1 = new double[DimX1][NumSamples];
        for (int i = 0; i < DimX1; ++i) {
            for (int j = 0; j < NumSamples; ++j) {
                X1[i][j] = Mat[j][i];
            }
        }
        X2 = new double[DimX2][NumSamples];
        for (int i = 0; i < DimX2; ++i) {
            for (int j = 0; j < NumSamples; ++j) {
                X2[i][j] = Mat[j][DimX1 + i];
            }
        }
        X3 = new double[DimX3][NumSamples];
        for (int i = 0; i < DimX3; ++i) {
            for (int j = 0; j < NumSamples; ++j) {
                X3[i][j] = Mat[j][DimX1 + DimX2 + i];
            }
        }
        Y = new double[DimY][NumSamples];
        for (int i = 0; i < DimY; ++i) {
            for (int j = 0; j < NumSamples; ++j) {
                Y[i][j] = Mat[j][DimX1 + DimX2 + DimX3 + i];
            }
        }

    }

    public void OnRadio1(ActionEvent actionEvent) {
        if (!Radio1.selectedProperty().getValue()) {
            Radio1.setSelected(true);
            type_now = Polynom.Chebyshev;
        } else {
            Radio2.setSelected(false);
            Radio3.setSelected(false);
            Radio4.setSelected(false);
        }

    }

    public void OnRadio2(ActionEvent actionEvent) {
        if (!Radio2.selectedProperty().getValue()) {
            Radio2.setSelected(true);
            type_now = Polynom.Lejandr;
        } else {
            Radio1.setSelected(false);
            Radio3.setSelected(false);
            Radio4.setSelected(false);
        }
    }

    public void OnRadio3(ActionEvent actionEvent) {
        if (!Radio3.selectedProperty().getValue()) {
            Radio3.setSelected(true);
            type_now = Polynom.Lagger;
        } else {
            Radio2.setSelected(false);
            Radio1.setSelected(false);
            Radio4.setSelected(false);
        }
    }

    public void OnRadio4(ActionEvent actionEvent) {
        if (!Radio4.selectedProperty().getValue()) {
            Radio4.setSelected(true);
            type_now = Polynom.Hermit;
        } else {
            Radio2.setSelected(false);
            Radio3.setSelected(false);
            Radio1.setSelected(false);
        }
    }

    public void OnButtonRes(ActionEvent actionEvent) {
        JFileChooser fileopen = new JFileChooser();
        int ret = fileopen.showDialog(null, "Cоxранить");
        if (ret == JFileChooser.APPROVE_OPTION) {
            File file_input = fileopen.getSelectedFile();
            Text_res.setText(file_input.getPath());
        }
    }

    public void OnClean(ActionEvent actionEvent) {
        y0_first_level.getData().removeAll();
        y0_first_level.getData().clear();
        Graph_Y.getData().removeAll();
        Graph_Y.getData().clear();
        Y0_last_level.getData().removeAll();
        Y0_last_level.getData().clear();
        Graph_Y_norm.getData().removeAll();
        Graph_Y_norm.getData().clear();

    }

    public void Write_to_text(String file_name) {
        try {
            File file = new File(file_name);

            BufferedReader br = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(file), "UTF-8"
                    )
            );
            String line = null;
            while ((line = br.readLine()) != null) {
                Text_results.setText(Text_results.getText() + line + "\n");
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void Write_Results(String file_rez_name, Calculations LambdaDefinition, Initial_data data_in, double err) {
        try {

            PrintWriter fout = new PrintWriter(new FileWriter(file_rez_name));
            fout.write("The coefficients which were obtained\n");
            fout.write("The first level (lambda1)\n");


            File file = new File("Lambda1.txt");

            BufferedReader br = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(file), "UTF-8"
                    )
            );
            String line = null;
            while ((line = br.readLine()) != null) {
                fout.println(line);
            }
            br.close();


            fout.write("The first level (Lambda2)\n");


            file = new File("Lambda2.txt");

            br = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(file), "UTF-8"
                    )
            );
            line = null;
            while ((line = br.readLine()) != null) {
                fout.println(line);
            }
            br.close();


            fout.write("The first level (Lambda3)\n");


            file = new File("Lambda3.txt");

            br = new BufferedReader(
                    new InputStreamReader(
                            new FileInputStream(file), "UTF-8"
                    )
            );
            line = null;
            while ((line = br.readLine()) != null) {
                fout.println(line);
            }
            br.close();


            int ColumnLength = LambdaDefinition.A1.length;

            fout.write("The second level(a1)\n");
            for (int j = 0; j < ColumnLength; ++j) {
                fout.print(String.format("%.4f", LambdaDefinition.A1[j]));
                fout.write("  ");
            }
            fout.write("\n");

            ColumnLength = LambdaDefinition.A2.length;

            fout.write("The second level (a2)\n");
            for (int j = 0; j < ColumnLength; ++j) {
                fout.print(String.format("%.4f", LambdaDefinition.A2[j]));
                fout.write("  ");
            }
            fout.write("\n");

            ColumnLength = LambdaDefinition.A3.length;

            fout.write("The second level (a3)\n");
            for (int j = 0; j < ColumnLength; ++j) {
                fout.print(String.format("%.4f", LambdaDefinition.A3[j]));
                fout.write("  ");
            }
            fout.write("\n");


            ColumnLength = LambdaDefinition.C.length;

            fout.write("The third level (C)\n");
            for (int j = 0; j < ColumnLength; ++j) {
                fout.print(String.format("%.4f", LambdaDefinition.C[j]));
                fout.write("  ");
            }
            fout.write("\n");


            String s1 = "", s2 = "", s3 = "";

            for (int i = 0; i < data_in.length_of_x1; i++)
                for (int j = 0; j < data_in.p1; j++)
                    s1 = s1 + "(" + Double.toString(LambdaDefinition.Lambda1[i][j]) + ")" + "*T" + Integer.toString(j) + "(x1" + Integer.toString(i) + ")+";
            for (int i = 0; i < data_in.length_of_x2; i++)
                for (int j = 0; j < data_in.p2; j++)
                    s2 = s2 + "(" + Double.toString(LambdaDefinition.Lambda2[i][j]) + ")" + "*T" + Integer.toString(j) + "(x2" + Integer.toString(i) + ")+";
            for (int i = 0; i < data_in.length_of_x3; i++)
                for (int j = 0; j < data_in.p3; j++)
                    s3 = s3 + "(" + Double.toString(LambdaDefinition.Lambda3[i][j]) + ")" + "*T" + Integer.toString(j) + "(x3" + Integer.toString(i) + ")+";
            s3 = s3.substring(0, s3.length() - 1);
            String s1z = s1 + s2 + s3;
            fout.write("psi functions\n");
            fout.write(s1z);


            String s12 = "", s22 = "", s32 = "";

            for (int i = 0; i < data_in.length_of_x1; i++)

                s12 = s12 + "(" + Double.toString(LambdaDefinition.A1[i]) + ")" + "*psi1(x1" + Integer.toString(i) + ")+";
            for (int i = 0; i < data_in.length_of_x2; i++)
                s22 = s22 + "(" + Double.toString(LambdaDefinition.A2[i]) + ")" + "*psi2(x2" + Integer.toString(i) + ")+";
            for (int i = 0; i < data_in.length_of_x3; i++)
                s32 = s32 + "(" + Double.toString(LambdaDefinition.A3[i]) + ")" + "*psi3(x3" + Integer.toString(i) + ")+";

            s32 = s32.substring(0, s32.length() - 1);
            //   String s12z=s12+s22+s32;
            //   fout.write("\n");
            fout.write("The second level\n");
            fout.write("F1=" + s12 + "\n");
            fout.write("F2=" + s22 + "\n");
            fout.write("F3=" + s32 + "\n");


            fout.write("\n");
            String s14 = "";
            fout.write("The third level\n");
            for (int i = 0; i < 3; i++)

                s14 = s14 + Double.toString(LambdaDefinition.C[i]) + "*F(x1)+";
            s14 = s14.substring(0, s14.length() - 1);
            fout.write(s14);

            fout.write("\n");

            fout.write("Function which was obtained\n");


            String s10 = "", s20 = "", s30 = "";


            for (int i = 0; i < data_in.length_of_x1; i++)

                s10 = s10 + Double.toString(LambdaDefinition.A1[i]) + "*(" + s1 + ")+";
            for (int i = 0; i < data_in.length_of_x2; i++)
                s20 = s20 + Double.toString(LambdaDefinition.A2[i]) + "*(" + s2 + ")+";

            for (int i = 0; i < data_in.length_of_x3; i++)
                s30 = s30 + Double.toString(LambdaDefinition.A3[i]) + "*(" + s3 + ")+";

            s30 = s30.substring(0, s30.length() - 1);


            s10 = Double.toString(LambdaDefinition.C[0]) + "*(" + s10 + ")" + Double.toString(LambdaDefinition.C[0]) + "*(" + s20 + ")" + Double.toString(LambdaDefinition.C[0]) + "*(" + s30 + ")";
            fout.write(s10 + "\n");
            fout.write("Function which was obtained\n");
            s10 = Double.toString(LambdaDefinition.MinY) + "+" + Double.toString(LambdaDefinition.coef) + "*(" + s10 + ")";
            fout.write(s10 + "\n");

            fout.write("Error of function recovering \n");
            s10 = Double.toString(err);
            fout.write(s10 + "\n");

            fout.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void open_file_input_y(ActionEvent actionEvent) {
        JFileChooser fileopen = new JFileChooser();
        int ret = fileopen.showDialog(null, "Открыть Fайл");
        if (ret == JFileChooser.APPROVE_OPTION) {
            File file_input = fileopen.getSelectedFile();
            Text_Y.setText(file_input.getPath());
        }

    }


    public void ONbox_down_main(ActionEvent actionEvent) {


        data_real.getColumns().removeAll();
        table_data_input.getColumns().removeAll();
        data_real.getColumns().clear();
        table_data_input.getColumns().clear();
        ObservableList<String[]> data1 = FXCollections.observableArrayList();
        ObservableList<String[]> data = FXCollections.observableArrayList();
        data1.removeAll();
        table_data_input.setItems(data1);

        data.removeAll();
        data_real.setItems(data);
        data_real.getColumns().removeAll();
        table_data_input.getColumns().removeAll();


        data_in = new Initial_data(Integer.parseInt(box_length.getText()), Integer.parseInt(box_x1.getText()), Integer.parseInt(box_x2.getText()), Integer.parseInt(box_x3.getText()), Integer.parseInt(box_y.getText()));
        int RowSize = data_in.length_of_sample;
        int ColumnSize = data_in.column_sum;

        DaMa = new DataMapping(text_file_input.getText(), data_in.length_of_x1, data_in.length_of_x2, data_in.length_of_x3, data_in.length_of_y, data_in.length_of_sample);
        MaxY = new double[data_in.length_of_y];
        MinY = new double[data_in.length_of_y];
        int step1 = data_in.length_of_x1 + data_in.length_of_x2 + data_in.length_of_x3;
        int step2 = step1 + data_in.length_of_y;
        for (int i = step1; i < step2; ++i) {
            MaxY[i - step1] = DaMa.MaxArray[i];
            MinY[i - step1] = DaMa.MinArray[i];
        }


        String[][] staffArray = new String[RowSize + 1][ColumnSize];

        for (int j = 0; j < data_in.length_of_x1; ++j)
            staffArray[0][j] = "X1" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x2; ++j)
            staffArray[0][data_in.length_of_x1 + j] = "X2" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x3; ++j)
            staffArray[0][data_in.length_of_x1 + data_in.length_of_x2 + j] = "X3" + Integer.toString(1 + j);

        for (int j = 0; j < data_in.length_of_y; ++j)
            staffArray[0][data_in.length_of_x1 + data_in.length_of_x3 + data_in.length_of_x2 + j] = "Y" + Integer.toString(j + 1);


        for (int i = 0; i < RowSize; ++i) {
            for (int j = 0; j < ColumnSize; ++j)
                staffArray[i + 1][j] = String.valueOf(DaMa.RealData[i][j]);
        }
        //    StackPane root = new StackPane();


        data = FXCollections.observableArrayList();
        data.addAll(Arrays.asList(staffArray));
        data.remove(0);//remove titles from data
        for (int i = 0; i < staffArray[0].length; i++) {
            TableColumn tc = new TableColumn(staffArray[0][i]);
            final int colNo = i;
            tc.setCellValueFactory(new Callback<CellDataFeatures<String[], String>, ObservableValue<String>>() {
                @Override
                public ObservableValue<String> call(CellDataFeatures<String[], String> p) {
                    return new SimpleStringProperty((p.getValue()[colNo]));
                }
            });

            data_real.getColumns().add(tc);
        }
        data_real.setItems(data);


        String[][] staffArray1 = new String[RowSize + 1][ColumnSize];

        for (int j = 0; j < data_in.length_of_x1; ++j)
            staffArray1[0][j] = "X1" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x2; ++j)
            staffArray1[0][data_in.length_of_x1 + j] = "X2" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x3; ++j)
            staffArray1[0][data_in.length_of_x1 + data_in.length_of_x2 + j] = "X3" + Integer.toString(1 + j);

        for (int j = 0; j < data_in.length_of_y; ++j)
            staffArray1[0][data_in.length_of_x1 + data_in.length_of_x3 + data_in.length_of_x2 + j] = "Y" + Integer.toString(j + 1);


        for (int i = 0; i < RowSize; ++i) {
            for (int j = 0; j < ColumnSize; ++j)
                staffArray1[i + 1][j] = String.valueOf(DaMa.NormalizedData[i][j]);
        }


        data1 = FXCollections.observableArrayList();
        data1.addAll(Arrays.asList(staffArray1));
        data1.remove(0);//remove titles from data
        for (int i = 0; i < staffArray1[0].length; i++) {
            TableColumn tc = new TableColumn(staffArray1[0][i]);
            final int colNo = i;
            tc.setCellValueFactory(new Callback<CellDataFeatures<String[], String>, ObservableValue<String>>() {
                @Override
                public ObservableValue<String> call(CellDataFeatures<String[], String> p) {
                    return new SimpleStringProperty((p.getValue()[colNo]));
                }
            });

            table_data_input.getColumns().add(tc);
        }
        table_data_input.setItems(data1);


        // X_I and Y initialization
        DistributeData(DaMa.NormalizedData, DaMa.DimX1, DaMa.DimX2, DaMa.DimX3, DaMa.DimY, DaMa.NumSamples);
        try {
            DataMapping.WriteToFile("X1.txt", X1);
            DataMapping.WriteToFile("X2.txt", X2);
            DataMapping.WriteToFile("X3.txt", X3);
            DataMapping.WriteToFile("Y.txt", Y);
        } catch (IOException e) {
            e.printStackTrace();
        }

        //   box_down_main.setDisable(true);
        box_down_y.setDisable(false);
    }

    public void Onbox_down_y(ActionEvent actionEvent) {
        data_real.getColumns().removeAll();
        table_data_input.getColumns().removeAll();
        data_real.getColumns().clear();
        table_data_input.getColumns().clear();
        ObservableList<String[]> data1 = FXCollections.observableArrayList();
        ObservableList<String[]> data = FXCollections.observableArrayList();
        data1.removeAll();
        table_data_input.setItems(data1);

        data.removeAll();
        data_real.setItems(data);
        data_real.getColumns().removeAll();
        table_data_input.getColumns().removeAll();


        box_y.setText(Integer.toString(Integer.parseInt(box_y.getText()) + 1));
        data_in = new Initial_data(Integer.parseInt(box_length.getText()), Integer.parseInt(box_x1.getText()), Integer.parseInt(box_x2.getText()), Integer.parseInt(box_x3.getText()), Integer.parseInt(box_y.getText()));
        int RowSize = data_in.length_of_sample;
        int ColumnSize = data_in.column_sum;

        DataMapping DaMa = new DataMapping(text_file_input.getText(), Text_Y.getText(), data_in.length_of_x1, data_in.length_of_x2, data_in.length_of_x3, data_in.length_of_y, data_in.length_of_sample);
        MaxY = new double[data_in.length_of_y];
        MinY = new double[data_in.length_of_y];
        int step1 = data_in.length_of_x1 + data_in.length_of_x2 + data_in.length_of_x3;
        int step2 = step1 + data_in.length_of_y;
        for (int i = step1; i < step2; ++i) {
            MaxY[i - step1] = DaMa.MaxArray[i];
            MinY[i - step1] = DaMa.MinArray[i];
        }


        String[][] staffArray = new String[RowSize + 1][ColumnSize];

        for (int j = 0; j < data_in.length_of_x1; ++j)
            staffArray[0][j] = "X1" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x2; ++j)
            staffArray[0][data_in.length_of_x1 + j] = "X2" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x3; ++j)
            staffArray[0][data_in.length_of_x1 + data_in.length_of_x2 + j] = "X3" + Integer.toString(1 + j);

        for (int j = 0; j < data_in.length_of_y; ++j)
            staffArray[0][data_in.length_of_x1 + data_in.length_of_x3 + data_in.length_of_x2 + j] = "Y" + Integer.toString(j + 1);


        for (int i = 0; i < RowSize; ++i) {
            for (int j = 0; j < ColumnSize; ++j)
                staffArray[i + 1][j] = String.valueOf(DaMa.RealData[i][j]);
        }
        //    StackPane root = new StackPane();


        //  ObservableList<String[]>
        data = FXCollections.observableArrayList();
        data.removeAll();
        data.addAll(Arrays.asList(staffArray));
        data.remove(0);//remove titles from data
        for (int i = 0; i < staffArray[0].length; i++) {
            TableColumn tc = new TableColumn(staffArray[0][i]);
            final int colNo = i;
            tc.setCellValueFactory(new Callback<CellDataFeatures<String[], String>, ObservableValue<String>>() {
                @Override
                public ObservableValue<String> call(CellDataFeatures<String[], String> p) {
                    return new SimpleStringProperty((p.getValue()[colNo]));
                }
            });

            data_real.getColumns().add(tc);

        }
        data_real.setItems(data);


        String[][] staffArray1 = new String[RowSize + 1][ColumnSize];

        for (int j = 0; j < data_in.length_of_x1; ++j)
            staffArray1[0][j] = "X1" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x2; ++j)
            staffArray1[0][data_in.length_of_x1 + j] = "X2" + Integer.toString(j + 1);

        for (int j = 0; j < data_in.length_of_x3; ++j)
            staffArray1[0][data_in.length_of_x1 + data_in.length_of_x2 + j] = "X3" + Integer.toString(1 + j);

        for (int j = 0; j < data_in.length_of_y; ++j)
            staffArray1[0][data_in.length_of_x1 + data_in.length_of_x3 + data_in.length_of_x2 + j] = "Y" + Integer.toString(j + 1);


        for (int i = 0; i < RowSize; ++i) {
            for (int j = 0; j < ColumnSize; ++j)
                staffArray1[i + 1][j] = String.valueOf(DaMa.NormalizedData[i][j]);
        }


        //   ObservableList<String[]>
        data1 = FXCollections.observableArrayList();
        data1.removeAll();
        data1.addAll(Arrays.asList(staffArray1));

        data1.remove(0);//remove titles from data
        for (int i = 0; i < staffArray1[0].length; i++) {
            TableColumn tc = new TableColumn(staffArray1[0][i]);
            final int colNo = i;
            tc.setCellValueFactory(new Callback<CellDataFeatures<String[], String>, ObservableValue<String>>() {
                @Override
                public ObservableValue<String> call(CellDataFeatures<String[], String> p) {
                    return new SimpleStringProperty((p.getValue()[colNo]));
                }
            });

            table_data_input.getColumns().add(tc);
        }
        table_data_input.setItems(data1);


        // X_I and Y initialization
        DistributeData(DaMa.NormalizedData, DaMa.DimX1, DaMa.DimX2, DaMa.DimX3, DaMa.DimY, DaMa.NumSamples);
        try {
            DataMapping.WriteToFile("X1.txt", X1);
            DataMapping.WriteToFile("X2.txt", X2);
            DataMapping.WriteToFile("X3.txt", X3);
            DataMapping.WriteToFile("Y.txt", Y);
        } catch (IOException e) {
            e.printStackTrace();
        }


        box_down_y.setDisable(true);


    }

    public void Onbox_clean_data(ActionEvent actionEvent) {
        // table_data_input = new TableView();
        data_real.getColumns().removeAll();
        table_data_input.getColumns().removeAll();
        data_real.getColumns().clear();
        table_data_input.getColumns().clear();
        ObservableList<String[]> data1 = FXCollections.observableArrayList();
        ObservableList<String[]> data = FXCollections.observableArrayList();
        data1.removeAll();
        table_data_input.setItems(data1);

        box_down_main.setDisable(false);
        box_down_y.setDisable(true);
        data.removeAll();
        data_real.setItems(data);
        data_real.getColumns().removeAll();
        table_data_input.getColumns().removeAll();
    }

    public void OnButton_Draw_real(ActionEvent actionEvent) {
        XYChart.Series series = new XYChart.Series();
        series = new XYChart.Series();
        series.setName("Real_data");
        for (int i = 0; i < Y[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i + 1), LambdaDefinition.Y_true[i]));

        Graph_Y.getData().add(series);

        series = new XYChart.Series();
        series = new XYChart.Series();
        series.setName("Real_data");
        for (int i = 0; i < Y[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i + 1), Y[Integer.parseInt(Y_order.getText()) - 1][i]));

        Graph_Y_norm.getData().add(series);


    }

    public void OnButton_Draw_model(ActionEvent actionEvent) {

        XYChart.Series series = new XYChart.Series();
        series.setName(data_in.name_of_graph);
        for (int i = 0; i < Y[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i + 1), LambdaDefinition.ResidVecLambda[i]));

        y0_first_level.getData().add(series);

        series = new XYChart.Series();
        series.setName(data_in.name_of_graph);
        for (int i = 0; i < Y[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i + 1), LambdaDefinition.residualBigF_System[i]));

        Y0_last_level.getData().add(series);


        series = new XYChart.Series();
        series.setName(data_in.name_of_graph);
        for (int i = 0; i < Y[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i + 1), LambdaDefinition.Y_model[i]));

        Graph_Y.getData().add(series);


        series = new XYChart.Series();
        series.setName(data_in.name_of_graph);
        for (int i = 0; i < Y[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i + 1), LambdaDefinition.Ax[i]));

        Graph_Y_norm.getData().add(series);
    }


    public void On_button_forecast(ActionEvent actionEvent) {
        int step_now = Integer.parseInt(Text_time.getText());
        int window = Integer.parseInt(Text_window.getText());
        int predict_length=Integer.parseInt(Text_prediction.getText());
      Robust robust_1 = new Robust(DaMa.RealData[0].length, window, step_now);

        double[][] Forecast_array=robust_1.calculation_of_condition(DaMa.RealData, DaMa.NormalizedData);











        Text_robust.setText(robust_1.Massage());

        int[][]  ArrayWithIndexes = {{1,2,3, 4},{0,2,3},{2,3,4},{5,6, 7}};

        BigPredictions BigPredictions_1 = new BigPredictions(DaMa.Transpose(DaMa.NormalizedData), ArrayWithIndexes, 0, DaMa.NormalizedData.length, Integer.parseInt(box_p1.getText()), window, step_now, predict_length);
     //   BigPredictions BigPredictions_1 = new BigPredictions(Forecast_array, ArrayWithIndexes, 0, Forecast_array.length, data_in.p1, window, step_now, predict_length);
        double[][] Y_prediction = BigPredictions_1.getPredictedY_i();
        double[][] Y_prediction_recover=new double[Y_prediction.length][Y_prediction[0].length];

        for (int i=0;i<Y_prediction.length;i++)
            for (int j=0;j<Y_prediction[0].length;j++)
                Y_prediction_recover[i][j]=Y_prediction[i][j]*(MaxY[i]-MinY[i])+MinY[i];

        Estimation_of_situation estimationOfSituation =new Estimation_of_situation(step_now);
        double[] funk = estimationOfSituation.crucial(Y_prediction_recover, predict_length);
        if(estimationOfSituation.level_crush<0.15) {
            Progress1.setVisible(true);
            Progress2.setVisible(false);
            Progress3.setVisible(false);
            Progress1.setProgress(estimationOfSituation.level_crush);
        }
        if((estimationOfSituation.level_crush>0.15)&&(estimationOfSituation.level_crush<0.7)) {
            Progress2.setVisible(true);
            Progress1.setVisible(false);
            Progress3.setVisible(false);
            Progress2.setProgress(estimationOfSituation.level_crush);
        }
        if(estimationOfSituation.level_crush>0.7) {
            Progress3.setVisible(true);
            Progress2.setVisible(false);
            Progress1.setVisible(false);
            Progress3.setProgress(estimationOfSituation.level_crush);
        }



        String string="";
        string="Time="+Integer.toString(step_now*10);
        string=string+"  Total level of fatality="+Integer.toString((int)(estimationOfSituation.level_crush*8))+"\n";
        string=string+"Probability of steps:\n";
        for (int j=0;j<Y_prediction.length;j++)
            for (int i=0;i<Y_prediction[0].length;i++)
                string=string+" For y="+Integer.toString(j+1)+" For step="+Integer.toString(i+1)+":"+Double.toString(estimationOfSituation.persent[j][i]);


        Text_foR_for.setText(string);








        Graph_fuel.getData().removeAll();
        Graph_fuel.getData().clear();
        Graph_acc.getData().removeAll();
        Graph_acc.getData().clear();
        Graph_grid.getData().removeAll();
        Graph_grid.getData().clear();
        Graph_level.getData().removeAll();
        Graph_level.getData().clear();


double [][] Y_recover = new double[Y.length][Y[0].length];
        for (int i=0;i<Y.length;i++)
            for(int j=0;j<Y[0].length;j++)
                Y_recover[i][j]=Y[i][j]*(MaxY[i]-MinY[i])+MinY[i];

        XYChart.Series series = new XYChart.Series();
        series = new XYChart.Series();
        series.setName("Real value");
        if (step_now > window) {
            for (int i = 0; i < window+Y_prediction[0].length+1; i++)

                series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), Y_recover[0][step_now - window + i]));
        } else {
            for (int i = 0; i < step_now; i++)
                series.getData().add(new XYChart.Data(Integer.toString(10 * i), Y_recover[0][i]));
        }
        Graph_grid.getData().add(series);






        series = new XYChart.Series();
        series.setName("Forecast");
        for (int i = 0; i < Y_prediction[0].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now + 1 + i)), Y_prediction_recover[0][i]));

        Graph_grid.getData().add(series);


   //     Graph_grid.setCreateSymbols(false);


        series = new XYChart.Series();
        series.setName("level of contingency");
        for (int i = 0; i < window+Y_prediction[0].length+1; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), estimationOfSituation.level_est_contingency[0]));

        Graph_grid.getData().add(series);


        series = new XYChart.Series();
        series.setName("level of critical");
        for (int i = 0; i < window+Y_prediction[0].length+1; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), estimationOfSituation.level_est_critical[0]));

        Graph_grid.getData().add(series);



        series = new XYChart.Series();
        series.setName("Real value");
        if (step_now > window) {
            for (int i = 0; i < window+Y_prediction[1].length+1; i++)

                series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), Y_recover[1][step_now - window + i]));
        } else {
            for (int i = 0; i < step_now; i++)
                series.getData().add(new XYChart.Data(Integer.toString(10 * i), Y_recover[1][i]));
        }
        Graph_fuel.getData().add(series);



        series = new XYChart.Series();
        series.setName("Forecast");
        for (int i = 0; i < Y_prediction[1].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now + 1 + i)), Y_prediction_recover[1][i]));

        Graph_fuel.getData().add(series);



        series = new XYChart.Series();
        series.setName("level of contingency");
        for (int i = 0; i < window+Y_prediction[0].length+1; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), estimationOfSituation.level_est_contingency[1]));

        Graph_fuel.getData().add(series);


        series = new XYChart.Series();
        series.setName("level of critical");
        for (int i = 0; i < window+Y_prediction[0].length+1; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), estimationOfSituation.level_est_critical[1]));

        Graph_fuel.getData().add(series);









        series = new XYChart.Series();
        series.setName("Real value");
        if (step_now > window) {
            for (int i = 0; i < window+Y_prediction[2].length+1; i++)

                series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), Y_recover[2][step_now - window + i]));
        } else {
            for (int i = 0; i < step_now; i++)
                series.getData().add(new XYChart.Data(Integer.toString(10 * i), Y_recover[2][i]));
        }
        Graph_acc.getData().add(series);



        series = new XYChart.Series();
        series.setName("Forecast");
        for (int i = 0; i < Y_prediction[2].length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now + 1 + i)), Y_prediction_recover[2][i]));

        Graph_acc.getData().add(series);





        series = new XYChart.Series();
        series.setName("level of contingency");
        for (int i = 0; i < window+Y_prediction[0].length+1; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), estimationOfSituation.level_est_contingency[2]));

        Graph_acc.getData().add(series);


        series = new XYChart.Series();
        series.setName("level of critical");
        for (int i = 0; i < window+Y_prediction[0].length+1; i++)
            series.getData().add(new XYChart.Data(Integer.toString(10 * (step_now - window + i)), estimationOfSituation.level_est_critical[2]));

        Graph_acc.getData().add(series);






        series = new XYChart.Series();
        series.setName("Critical level");
        for (int i = 0; i < funk.length; i++)
            series.getData().add(new XYChart.Data(Integer.toString(i+1), funk[i]));

        Graph_level.getData().add(series);


     //   int step_now=Integer.parseInt(Text_time.getText());
        SliderTime.setValue(step_now);
        Text_grid.setText(Double.toString(Y[0][step_now]*(MaxY[0]-MinY[0])+MinY[0]));
        Text_fuel.setText(Double.toString(Y[1][step_now]*(MaxY[1]-MinY[1])+MinY[1]));
        Text_acc.setText(Double.toString(Y[2][step_now]*(MaxY[2]-MinY[2])+MinY[2]));
    }



}