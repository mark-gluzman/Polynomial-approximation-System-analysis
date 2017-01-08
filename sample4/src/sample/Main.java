package sample;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;

import java.util.Timer;
import java.util.TimerTask;


public class Main extends Application {

    @Override
    public void start(final Stage primaryStage) throws Exception{
        Parent root1 = FXMLLoader.load(getClass().getResource("sample.fxml"));
        Parent root = FXMLLoader.load(getClass().getResource("sample.fxml"));
        primaryStage.setTitle("SA lab4");
        Scene scene1=new Scene(root1);
        final Scene scene=new Scene(root);
        scene.getStylesheets().add(getClass().getResource("test.css").toExternalForm());
        primaryStage.setScene(scene1);
        ServiceSample newr =new ServiceSample();






       primaryStage.setScene(scene);
       newr.start(primaryStage);


        Timer timer = new Timer();
        timer.schedule(new TimerTask() {
            @Override
            public void run() {
                Platform.runLater(new Runnable() {
                    public void run() {

                        primaryStage.setScene(scene);
                        primaryStage.setMaximized(true);
       //                 primaryStage.show();
                    }
                });
            }
        }, 3000);



     //   primaryStage.setScene(scene);
     //
    }
    public static void main(String[] args) {       launch(args);


    }


}
