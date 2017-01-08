package sample;

/**
 * Created by Mark on 14.10.2014.
 */

public class Initial_data {

    int length_of_sample=45;
    int length_of_x1=2;
    int length_of_x2=2;
    int length_of_x3=3;
    int length_of_y=4;
    int column_sum=11;
    int [] length_of_all;
    public String name_of_graph;
    Polynom type;

    int p1=3;
    int p2=3;
    int p3=3;

    public Initial_data(int ls, int l1, int l2, int l3, int ly,Polynom type1, int p11, int p12, int p13 )
    {
        length_of_sample=ls;
        length_of_x1=l1;
        length_of_x2=l2;
        length_of_x3=l3;
        length_of_y=ly;
        column_sum=length_of_x1+length_of_x2+length_of_x3+length_of_y;
        type=type1;
        p1=p11;
        p2=p12;
        p3=p13;


        if ( type1==Polynom.Chebyshev)
        {
            name_of_graph="Chebyshev "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }

        if ( type1==Polynom.Lejandr)
        {
            name_of_graph="Lejandr "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }
        if (type1==Polynom.Lagger)
        {
            name_of_graph="Lagger "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }
        if (type1==Polynom.Hermit)
        {
            name_of_graph="Hermit "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }

      }
    public void add(Polynom type1, int p11, int p12, int p13 )
    {

        column_sum=length_of_x1+length_of_x2+length_of_x3+length_of_y;
        type=type1;
        p1=p11;
        p2=p12;
        p3=p13;


        if ( type1==Polynom.Chebyshev)
        {
            name_of_graph="Chebyshev "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }

        if ( type1==Polynom.Lejandr)
        {
            name_of_graph="Lejandr "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }
        if (type1==Polynom.Lagger)
        {
            name_of_graph="Lagger "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }
        if (type1==Polynom.Hermit)
        {
            name_of_graph="Hermit "+Integer.toString(p11)+" "+Integer.toString(p12)+" "+Integer.toString(p13);
        }

    }

    public Initial_data(int ls, int l1, int l2, int l3, int ly)
    {
        length_of_sample=ls;
        length_of_x1=l1;
        length_of_x2=l2;
        length_of_x3=l3;
        length_of_y=ly;
        column_sum=length_of_x1+length_of_x2+length_of_x3+length_of_y;
        length_of_all= new int[4];
        length_of_all[0]=length_of_x1;
        length_of_all[1]=length_of_x2;
        length_of_all[2]=length_of_x3;
        length_of_all[3]=length_of_y;
    }

    public void set_sample(int i)

    {
        length_of_sample=i;
    }

    public void set_x3(int i)

    {
        length_of_x3=i;
    }
    public void set_x1(int i)

    {
        length_of_x1=i;
    }
    public void set_x2(int i)

    {
        length_of_x2=i;
    }
    public void set_y(int i)

    {
        length_of_y=i;
    }
}
