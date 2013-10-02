package edu.gmu.cds.sim;

public class Parameters
{
    public int potential_type;
    public int ndim = 3;

    public double mass1;
    public double epsilon1;
    public double rin1;
    public double rout1;
    public double rscale1[] = new double[3];
    public double theta1;
    public double phi1;
    public double heat1;
    public int    opt1;

    public double mass2;
    public double epsilon2;
    public double rin2;
    public double rout2;
    public double rscale2[] = new double[3];
    public double theta2;
    public double phi2;
    public double heat2;
    public int    opt2;

    public double eps1;
    public double eps2;

    public double tcurrent;


//  real (kind=8), dimension(:,:), allocatable :: x0, xout

    public int n, n1, n2;

    public double time, tstart, tend;
    public boolean tIsSet;
    public double inclination_degree;
    public double omega_degree;
    public double rmin;
    public double velocity_factor;
    public double mratio;
    public double secondary_size;

    public double sec_vec[] = new double[6];
    public boolean use_sec_vec;

    public double h;
    public static final double hbase = 0.001;

    public int nstep;
    public int nout;
    public int iout;
    public int unit;
    public int istep;

    public String fname;
    public int iostat;
    public boolean header_on;

    // info about time of closest approach
    public double tmin;
    public double vrmin;
    public double beta;

    public Parameters()
    {
        defaultParameters();
    }

    public void defaultParameters()
    {
        standardGalaxy1();
        standardGalaxy2();
        testCollision();

        customCollision();
    }

    public void standardGalaxy1()
    {
        // disk profile #1
        mass1 = 1.0;
        epsilon1  = 0.3;
        rin1 = 0.05;
        rout1 = 1.0;
        rscale1 = new double[]{3.0,3.0,3.0};
        theta1 = 0.0;
        phi1 = 0.0;
        opt1 = 1;
        heat1 = 0.0;
        eps1 = epsilon1*epsilon1;
    }

    public void standardGalaxy2()
    {
        // disk profile #2
        mass2 = 1.0;
        epsilon2  = 0.3;
        rin2 = 0.05;
        rout2 = 1.0;
        rscale2 = new double[]{3.0,3.0,3.0};
        theta2 = 0.0;
        phi2 = 0.0;
        opt2 = 1;
        heat2 = 0.0;
        eps2 = epsilon2*epsilon2;
    }

    public void testCollision()
    {
        // collision parameters
        inclination_degree = 90.0d;
        omega_degree = 0.0;
        rmin = 1.0;
        velocity_factor = 1.0;
        time = -3.0;

        // time step
        h = hbase;
        nout = 5;
        nstep = 500;

        // particle numbers
        n1 = 1000;
        n2 = 1000;
        n = n1 + n2;
    }

    public void customCollision()
    {
        phi1 = 5.0;
        theta1 = 5.0;
        rscale1 = new double[]{1.0,1.0,1.0};
        rout1 = 1.0;
        mass1 = 1.0;
        epsilon1 = 0.3;
        eps1 = epsilon1*epsilon1;
        n1 = 1000;
        heat1 = 0.0;
        opt1 = 1;

        phi2 = 0.0;
        theta2 = 0.0;
        rscale2 = new double[]{0.3,0.3,0.3};
        rout2 = 0.5;
        mass2 = 0.5;
        epsilon2 = 0.3;
        eps2 = epsilon2*epsilon2;
        n2 = 500;
        heat2 = 0.0;
        opt2 = 1;

        inclination_degree = 20.0;
        omega_degree = 0.0;
        rmin = 0.9;
        velocity_factor = 0.9;

        h = hbase;
        time = -5;

        n = n1 + n2;
    }
}
