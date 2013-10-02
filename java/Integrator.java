package edu.gmu.cds.sim;

/**
 *
 */
public class Integrator
{
    public static final double od6 = 1.0d/6.0d;
    public static final double od3 = 1.0d/3.0d;

    public ForceModel force;

    /** temporary array used in integration */
    public double x[][] = null;

    /** temporary array used in integration */
    public double xe[][] = null;

    /** internal array to hold force,acceleration calculations */
    public double f[][] = null;
 
    public Integrator(ForceModel force)
    {
        this.force = force;
    }

    /**
     * Initialize temporary storage arrays
     */
    public void initRKVar(int n)
    {
        x  = new double[n+1][6];
        xe = new double[n+1][6];
        f  = new double[n+1][6];

        if(force != null)
        {
            force.initVars(n+1);
        }
    }

    /**
     * Cleanup temporary storage.
     */
    public void deallocateRKVar()
    {
        x = null;
        xe = null;
        f = null;

        if(force != null)
        {
            force.deallocateVars();
        }
    }

    /**
     * Using the provided ForceModel, compute the updated positions and velocities
     * by taking a step of size h.  These are computed for all particles.
     */
    public void rk4(double x0[][], double xout[][], double h)
    {
        int n = x0.length;
        int i=0;

        // pre-compute products
        double hod6 = h*od6;
        double h0p5 = h*0.5d;
        double hod3 = h*od3;
     
        // we unroll the second dimension of 6 elements to speed computation       
//Prof.startCall("rk4_1");
        for(i=0;i<n;i++)
        {
            x[i][0] = x0[i][0];
            x[i][1] = x0[i][1];
            x[i][2] = x0[i][2];
            x[i][3] = x0[i][3];
            x[i][4] = x0[i][4];
            x[i][5] = x0[i][5];
        }
//Prof.endCall("rk4_1");
          
        force.diffeq(x,f);

//Prof.startCall("rk4_2");
        for(i=0;i<n;i++)
        {
            xe[i][0] = x0[i][0] + f[i][0] * hod6;
            x[i][0]  = x0[i][0] + f[i][0] * h0p5;
            xe[i][1] = x0[i][1] + f[i][1] * hod6;
            x[i][1]  = x0[i][1] + f[i][1] * h0p5;
            xe[i][2] = x0[i][2] + f[i][2] * hod6;
            x[i][2]  = x0[i][2] + f[i][2] * h0p5;
            xe[i][3] = x0[i][3] + f[i][3] * hod6;
            x[i][3]  = x0[i][3] + f[i][3] * h0p5;
            xe[i][4] = x0[i][4] + f[i][4] * hod6;
            x[i][4]  = x0[i][4] + f[i][4] * h0p5;
            xe[i][5] = x0[i][5] + f[i][5] * hod6;
            x[i][5]  = x0[i][5] + f[i][5] * h0p5;
        }
//Prof.endCall("rk4_2");

        force.diffeq(x,f);

//Prof.startCall("rk4_3");
        for(i=0;i<n;i++)
        {
            xe[i][0] = xe[i][0] + f[i][0] * hod3;
            x[i][0]  = x0[i][0] + f[i][0] * h0p5;
            xe[i][1] = xe[i][1] + f[i][1] * hod3;
            x[i][1]  = x0[i][1] + f[i][1] * h0p5;
            xe[i][2] = xe[i][2] + f[i][2] * hod3;
            x[i][2]  = x0[i][2] + f[i][2] * h0p5;
            xe[i][3] = xe[i][3] + f[i][3] * hod3;
            x[i][3]  = x0[i][3] + f[i][3] * h0p5;
            xe[i][4] = xe[i][4] + f[i][4] * hod3;
            x[i][4]  = x0[i][4] + f[i][4] * h0p5;
            xe[i][5] = xe[i][5] + f[i][5] * hod3;
            x[i][5]  = x0[i][5] + f[i][5] * h0p5;
        }
//Prof.endCall("rk4_3");

        force.diffeq(x,f);

//Prof.startCall("rk4_4");
        for(i=0;i<n;i++)
        {
            xe[i][0] = xe[i][0] + f[i][0] * hod3;
            x[i][0]  = x0[i][0] + h * f[i][0];
            xe[i][1] = xe[i][1] + f[i][1] * hod3;
            x[i][1]  = x0[i][1] + h * f[i][1];
            xe[i][2] = xe[i][2] + f[i][2] * hod3;
            x[i][2]  = x0[i][2] + h * f[i][2];
            xe[i][3] = xe[i][3] + f[i][3] * hod3;
            x[i][3]  = x0[i][3] + h * f[i][3];
            xe[i][4] = xe[i][4] + f[i][4] * hod3;
            x[i][4]  = x0[i][4] + h * f[i][4];
            xe[i][5] = xe[i][5] + f[i][5] * hod3;
            x[i][5]  = x0[i][5] + h * f[i][5];
        }
//Prof.endCall("rk4_4");

        force.diffeq(x,f);

//Prof.startCall("rk4_5");
        // combine last sub-step with copying to output array
        for(i=0;i<n;i++)
        {
            xe[i][0] = xe[i][0] + f[i][0] * hod6;
            xe[i][1] = xe[i][1] + f[i][1] * hod6;
            xe[i][2] = xe[i][2] + f[i][2] * hod6;
            xe[i][3] = xe[i][3] + f[i][3] * hod6;
            xe[i][4] = xe[i][4] + f[i][4] * hod6;
            xe[i][5] = xe[i][5] + f[i][5] * hod6;
              
            xout[i][0] = xe[i][0];
            xout[i][1] = xe[i][1];
            xout[i][2] = xe[i][2];
            xout[i][3] = xe[i][3];
            xout[i][4] = xe[i][4];
            xout[i][5] = xe[i][5];
        }
//Prof.endCall("rk4_5");

        return;
    } // end rk4

    /**
     * Using the provided ForceModel, compute the updated position and velocity of
     * the secondary galaxy by taking a step of size h.
     */
    public void rk41(double xx0[], double xxe[], double h)
    {
        double x[] = new double[7];
        double f[] = new double[7];

        int i=0;
        int n=6;

        // pre-compute products
        double hod6 = h*od6;
        double h0p5 = h*0.5d;
        double hod3 = h*od3;

        //for(i=0;i<n;i++)
        //{
                //x[i] = xx0[i];
        //}

        x[0] = xx0[0];
        x[1] = xx0[1];
        x[2] = xx0[2];
        x[3] = xx0[3];
        x[4] = xx0[4];
        x[5] = xx0[5];

        force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
            //xxe[i] = xx0[i] + h * f[i] * od6;
            //x[i]   = xx0[i] + h * f[i] * 0.5d;
            xxe[0] = xx0[0] + f[0] * hod6;
            x[0]   = xx0[0] + f[0] * h0p5;
            xxe[1] = xx0[1] + f[1] * hod6;
            x[1]   = xx0[1] + f[1] * h0p5;
            xxe[2] = xx0[2] + f[2] * hod6;
            x[2]   = xx0[2] + f[2] * h0p5;
            xxe[3] = xx0[3] + f[3] * hod6;
            x[3]   = xx0[3] + f[3] * h0p5;
            xxe[4] = xx0[4] + f[4] * hod6;
            x[4]   = xx0[4] + f[4] * h0p5;
            xxe[5] = xx0[5] + f[5] * hod6;
            x[5]   = xx0[5] + f[5] * h0p5;
        //}

        force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
                //xxe[i] = xxe[i] + h * f[i] * od3;
                //x[i]   = xx0[i] + h * f[i] *0.5d;
                xxe[0] = xxe[0] + f[0] * hod3;
                x[0]   = xx0[0] + f[0] * h0p5;
                xxe[1] = xxe[1] + f[1] * hod3;
                x[1]   = xx0[1] + f[1] * h0p5;
                xxe[2] = xxe[2] + f[2] * hod3;
                x[2]   = xx0[2] + f[2] * h0p5;
                xxe[3] = xxe[3] + f[3] * hod3;
                x[3]   = xx0[3] + f[3] * h0p5;
                xxe[4] = xxe[4] + f[4] * hod3;
                x[4]   = xx0[4] + f[4] * h0p5;
                xxe[5] = xxe[5] + f[5] * hod3;
                x[5]   = xx0[5] + f[5] * h0p5;
        //}

        force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
                //xxe[i] = xxe[i] + h * f[i] * od3;
                //x[i]   = xx0[i] + h*f[i];
                xxe[0] = xxe[0] + f[0] * hod3;
                x[0]   = xx0[0] + h*f[0];
                xxe[1] = xxe[1] + f[1] * hod3;
                x[1]   = xx0[1] + h*f[1];
                xxe[2] = xxe[2] + f[2] * hod3;
                x[2]   = xx0[2] + h*f[2];
                xxe[3] = xxe[3] + f[3] * hod3;
                x[3]   = xx0[3] + h*f[3];
                xxe[4] = xxe[4] + f[4] * hod3;
                x[4]   = xx0[4] + h*f[4];
                xxe[5] = xxe[5] + f[5] * hod3;
                x[5]   = xx0[5] + h*f[5];
        //}

        force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
                //xxe[i] = xxe[i] + h * f[i] * od6;
                xxe[0] = xxe[0] + f[0] * hod6;
                xxe[1] = xxe[1] + f[1] * hod6;
                xxe[2] = xxe[2] + f[2] * hod6;
                xxe[3] = xxe[3] + f[3] * hod6;
                xxe[4] = xxe[4] + f[4] * hod6;
                xxe[5] = xxe[5] + f[5] * hod6;
        //}

        return;
    } // end rk41
}
