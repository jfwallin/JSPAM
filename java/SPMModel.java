package edu.gmu.cds.sim;

/**
 * This ForceModel represents a softened point-mass.
 *
 */
public class SPMModel implements ForceModel
{
    // temporary storage arrays dependent on n
    public double[] r22, r21;
    public double[] r1, r2, rn;
    public double[] a1, a2, a3;

    // temporary storage independent of n
    public double xn[] = new double[6];
    public double xp = 0.0d;
    public double yp = 0.0d;
    public double zp = 0.0d;

    // force parameters
    public Parameters params = null;
    public double m1, m2, m3;
    public double eps1, eps2;


    public SPMModel()
    {
    }

    /**
     * Sets the simulation parameters to this model.
     */
    public void setParameters(Parameters p)
    {
        params = p;

        if(params != null)
        {
            m1 = params.mass1;
            m2 = params.mass2;
            m3 = params.mass2;
            eps1 = params.epsilon1 * params.epsilon1;
            eps2 = params.epsilon2 * params.epsilon2;
        }
    }


    /** 
     * Initialize temporary storage.
     */
    public void initVars(int n)
    {
        r22 = new double[n];
        r21 = new double[n];
        r1 = new double[n];
        r2 = new double[n];
        rn = new double[n];
        a1 = new double[n];
        a2 = new double[n];
        a3 = new double[n];
    }

    /** 
     * Cleanup temporary storage.
     */
    public void deallocateVars()
    {
        r22 = null;
        r21 = null;
        r1 = null;
        r2 = null;
        rn = null;
        a1 = null;
        a2 = null;
        a3 = null;
    }

    /**
     * Compute the circular velocity for a particle at a distance r from the specified mass.
     * The rout scale of the disk and softening length, eps, are provided.
     */
    public double circularVelocity(double mass, double r, double rout, double eps)
    {
        double ftotal = mass / ( r*r + eps );
        double v = Math.sqrt(ftotal*r);

        return v;
    }

    /**
     * For the given particle positions and velocities, calculate the accelerations.
     */
    public void diffeq(double x[][], double f[][])
    {
        double tmp1, tmp2, tmp3, tmp4, tmp5;
        int n = x.length;
        int i;

        for(i=0;i<6;i++)
        {
            xn[i] = x[n-1][i];
        }

        // make temps to handle perturber galaxy
        xp = xn[0];
        yp = xn[1];
        zp = xn[2];

        tmp4 = xp*xp+yp*yp+zp*zp;   // r2n
        tmp5 = Math.sqrt(tmp4);     // rn
        tmp4 = -m3/(tmp4+eps2);

        for(i=0;i<n;i++)
        {
            tmp1 = (x[i][0]-xp);
            tmp2 = (x[i][1]-yp);
            tmp3 = (x[i][2]-zp);

            r22[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;

            tmp1 = x[i][0];
            tmp2 = x[i][1];
            tmp3 = x[i][2];

            r21[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;

            r2[i]  = Math.sqrt(r22[i]);
            r1[i]  = Math.sqrt(r21[i]);
            rn[i]  = tmp5;

            a1[i] = -m1 / (r21[i] + eps1);
            a2[i] = -m2 / (r22[i] + eps2);
            a3[i] = tmp4;
        }

        // this is a correction to prevent NaN errors in the vectorized
        // function evalution at the location of the second mass
        r2[n-1] = 1.0d;

        // calculate the RHS of the diffeq

        // f(:,1) = x(:,4)
        // f(:,2) = x(:,5)
        // f(:,3) = x(:,6)

        // f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn
        // f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn
        // f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn

        for(i=0; i<n; i++)
        {
            tmp1 = a1[i]/r1[i];
            tmp2 = a2[i]/r2[i];
            tmp3 = a3[i]/rn[i];

            f[i][0] = x[i][3];
            f[i][1] = x[i][4];
            f[i][2] = x[i][5];

            f[i][3] = tmp1*x[i][0] + tmp2*(x[i][0]-xn[0]) + tmp3*xn[0];
            f[i][4] = tmp1*x[i][1] + tmp2*(x[i][1]-xn[1]) + tmp3*xn[1];
            f[i][5] = tmp1*x[i][2] + tmp2*(x[i][2]-xn[2]) + tmp3*xn[2];
        }
    }

    /**
     * Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
     */
    public void diffq1(double x[], double f[])
    {
        double r21, r1, a1, a2;

        r21 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        r1  = Math.sqrt(r21);

        a1 = -m1 / (r21 + eps1);
        a2 = -m2 / (r21 + eps2);

        a1 = a1 + a2;

        f[0] = x[3];
        f[1] = x[4];
        f[2] = x[5];
        f[3] = a1 * x[0] / r1;
        f[4] = a1 * x[1] / r1;
        f[5] = a1 * x[2] / r1;
    }
}
