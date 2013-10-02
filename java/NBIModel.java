package edu.gmu.cds.sim;

/**
 * This ForceModel represents an n-body halo/disk/bulge potential that is sampled
 * and then interpolated to calculate the force.
 *
 */
public class NBIModel implements ForceModel
{
    // force parameters
    public Parameters params = null;
    public double m1, m2, m3;
    public double mm1rs2, mm2rs2, mm3rs2;
    public double eps1, eps2;

    protected static boolean isInit = false;
    protected static int nnn = 10000;
    protected static int nnm1 = nnn-1;
    protected static double rad[];
    protected static double[] rho_halo, mass_halo;
    protected static double[] rho_disk, mass_disk;
    protected static double[] rho_bulge, mass_bulge;
    protected static double[] rho_total, mass_total;
    protected static double[] masses, radius, density;
    protected static double[] vr2, vr, new_vr2, new_vr;
    protected static double[] acceleration, acceleration_particle;
    protected static double[] new_mass, new_rho, phi;

    protected static double rrout1 = 0;
    protected static double rrout2 = 0;

    protected int pn = 0;
    protected int pn1 = 0;
    protected int pn2 = 0;


    protected static double rs_internal = 10.0d;
    protected static double rs_internal2 = rs_internal*rs_internal;
    protected static double rs_internal3 = rs_internal*rs_internal*rs_internal;
    protected static double rmax_scale = 100.0d;
    protected static final double sqrtpi = Math.sqrt(Math.PI);

    protected static double pscale;
    protected static double tmpdfind;

    protected static double lnl;

    // temporary storage arrays dependent on n
    public double[] r22, r21;
    public double[] r1, r2, rn;
    public double[] a1, a2, a3;
    public int ival11[] = null;
    public int ival22[] = null;
    public int ivaln[] = null;
    public double df_force11[] = null;
    public double df_force22[] = null;
    public double df_forcen[] = null;
    public double c3n[] = null;

    // temporary storage independent of n
    public double xn[] = new double[6];
    public double xp = 0.0d;
    public double yp = 0.0d;
    public double zp = 0.0d;
    public double vxp = 0.0d;
    public double vyp = 0.0d;
    public double vzp = 0.0d;


    public NBIModel()
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
            mm1rs2 = -m1*rs_internal2;
            mm2rs2 = -m2*rs_internal2;
            mm3rs2 = -m3*rs_internal2;
            eps1 = params.epsilon1 * params.epsilon1;
            eps2 = params.epsilon2 * params.epsilon2;
            pn = params.n;
            pn1 = params.n1;
            pn2 = params.n2;
            rrout1 = params.rout1;
            rrout2 = params.rout2;
            if(rad != null)
            {
                //initDistribution();
            }
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
 
        ival11 = new int[n];
        ival22 = new int[n];
        ivaln = new int[n];

        df_force11 = new double[n];
        df_force22 = new double[n];
        df_forcen = new double[n];

        c3n = new double[n];

        initDistribution();
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
 
        ival11 = null;
        ival22 = null;
        ivaln = null;

        df_force11 = null;
        df_force22 = null;
        df_forcen = null;

        c3n = null;
    }

    /**
     * Build an n-body halo/disk/bulge
     */
    protected synchronized void initDistribution()
    {
        if(isInit) return;
        isInit = true;

        double rmax;
        double mold, dmold, mtotal, mtot;
        double rscale;
        double dx, x;
        double alphahalo, qhalo, gammahalo, mhalo, rchalo, rhalo, epsilon_halo;
        double zdisk, zdiskmass, hdisk, zdiskmax;
        double hbulge, mbulge;
        double rho_tmp;
        double G, factor;
        double r, m,  pi;
        double p1, rd, rho_local;
        double p, rr, dr, rh, dp, mnew, dm;
        double acc_merge, rad_merge, acc;
        double xmax;

        int j, nmax, k, nmerge, ntotal, jj;

        rad = new double[nnn];
        rho_halo = new double[nnn];
        mass_halo = new double[nnn];
        rho_disk = new double[nnn];
        mass_disk = new double[nnn];
        rho_bulge = new double[nnn];
        mass_bulge = new double[nnn];
        rho_total = new double[nnn];
        mass_total = new double[nnn];
        masses = new double[nnn];
        radius = new double[nnn];
        density = new double[nnn];
        vr2 = new double[nnn];
        vr = new double[nnn];
        new_vr2 = new double[nnn];
        new_vr = new double[nnn];
        acceleration = new double[nnn];
        acceleration_particle = new double[nnn];
        new_mass = new double[nnn];
        new_rho = new double[nnn];
        phi = new double[nnn];

        // set the constant for dynamical friction
        lnl = 0.00d;

        // set up the parameters for the halo
        mhalo = 5.8d;
        rhalo = 10.0d;
        rchalo = 10.0d;
        gammahalo = 1.0d;
        epsilon_halo = 0.4d;
        pi = Math.PI;

        // derive additional constants
        qhalo = gammahalo / rchalo;
        alphahalo = 1.0d / ( 1.0d - sqrtpi * qhalo * Math.exp(qhalo*qhalo) * (1.0d - erf(qhalo)) );

        // set the integration limits and zero integration constants
        rmax = 20;
        nmax = 2000;
        dr = rmax / nmax;
        mold = 0;

        rscale = 5;
        ntotal = nnn;

        // set the limits for integration, and zero integration constants
        k = nmax / 2;
        dx = 1.0  / k;
        x = 0.0d;
        dmold = 0.0d;
        mtot = 0.0d;
        m = 0.0d;
        G = 1;

        /////
        // set the fundamental disk parameters
        zdisk = 0.2;
        zdiskmax = 3.5;
        hdisk = 1.0;


        /////
        // set the fundamental bulge parameters
        hbulge = 0.2;
        mbulge = 0.3333;


        ///////////////////////////////////////////////3
        ///// set up the radius array
        for(j = 0; j<nmax; j++)
        {
          x = x + dx;
          rad[j]= x * rchalo;
        }

        ///////////////////////////////////////////////3
        /////

        dr = rad[1] - rad[0];
        dx = dr / rchalo;

        for(j =0; j<nmax; j++)
        {
            // set the position
            r = rad[j];
            x = r / rchalo;

            // calculate the local rho based
            rho_tmp = alphahalo / (2*sqrtpi*sqrtpi*sqrtpi ) * (Math.exp(-x*x) / (x*x + qhalo*qhalo));

            // renormalize for the new halo size
            rho_tmp = rho_tmp / ( rchalo * rchalo * rchalo);

            // calculate mass in local shell, and update total mass
            dm = rho_tmp * 4 * pi * r * r *dr;
            mtot = mtot + dm;

            // store values in an array
            rho_halo[j] = rho_tmp * mhalo;
            mass_halo[j] = mtot * mhalo;
        }

        /////
        // now calculate the potential
        for(j = 0; j<nmax; j++)
        {
            r = rad[j];
            m = mass_halo[j];
            p1 = -G * m / r;
            phi[j] = p1;
        }


        ///////////////////////////////////////////////3
        // disk model

        /////
        // loop over the distribution

        for(j=0; j<nmax; j++)
        {
            // set the radius
            rd = rad[j];

            // find the local density in the disk
            rho_local  = Math.exp(-rd/hdisk)/ (8*pi*hdisk*hdisk);
            rho_disk[j] = rho_local;

            // find the mass in the spherical shell
            mnew = 4 * pi * rho_local * rd *rd * dr;

            mass_disk[j] = mnew + mold;
            mold = mass_disk[j];
        }

        ///////////////////////////////////////////////3
        // bulge model

        /////
        // loop over the distribution
        mold = 0.0;
        for(j=0; j<nmax; j++)
        {
            // set the radius
            rd = rad[j];

            // find the local density in the disk
            rho_local  = Math.exp(-rd*rd/(hbulge*hbulge));
            rho_bulge[j] = rho_local;

            // find the mass in the spherical shell
            mnew = 4 * pi * rho_local * rd *rd * dr;

            mass_bulge[j] = mnew + mold;
            mold = mass_bulge[j];
        }

        // renormalize distribution
        factor = mbulge / mass_bulge[nmax-1];
        for(j=0; j<nmax; j++)
        {
            mass_bulge[j] = mass_bulge[j] * factor;
            rho_bulge[j]  = rho_bulge[j]  * factor;
        }


        dr = rad[1] - rad[0];

        //////////////////////////////////////////////////
        j = 0;
        mass_total[j]=  (mass_halo[j] + mass_disk[j] + mass_bulge[j]);
        r = rad[j];
        rho_total[j] = mass_total[j] /  (4.0d/3.0d * pi * r * r * r);
        dr = rad[1] - rad[0];

        for(j = 1; j<nmax; j++)
        {
            r = rad[j];
            mass_total[j]=  (mass_halo[j] + mass_disk[j] + mass_bulge[j]);

            dm = mass_total[j] - mass_total[j-1];
            rho_total[j] = dm / (4 * pi * r * r * dr);
        }

        //////////////////////////////////////////////
        // find the velocity dispersion v_r**2

        masses = mass_total;
        radius = rad;
        density = rho_total;


        for(j=0; j<nmax; j++)
        {
            p = 0.0d;
            rr = radius[j];
            dr = radius[nmax-1] / nmax;

            for(jj = j; jj<nmax; jj++)
            {
                m  = masses[jj];
                rh = density[jj];
                rr = rr + dr;

                dp = rh * G * m / (rr*rr) * dr;
                p = p + dp;
            }

            vr2[j] = 1/density[j] * p;
            vr[j] = Math.sqrt(vr2[j]);
        }

        //////////////////////////////////////////////
        // find the velocity dispersion v_r**2
        masses = mass_total;
        radius = rad;
        density = rho_total;

        for(j=0; j<nmax; j++)
        {
            p = 0.0d;
            rr = radius[j];
            dr = radius[nmax-1] / nmax;
            for(jj = j; jj<nmax; jj++)
            {
                m  = masses[jj];
                rh = density[jj];
                rr = rr + dr;
                dp = rh * G * m / (rr*rr) * dr;
                p = p + dp;
            }

            vr2[j] = 1/density[j] * p;
            vr[j] = Math.sqrt(vr2[j]);
        }

        //////////////////////////////////////////////
        // find the accelerations felt by the particles and center of mass
        masses = mass_total;
        radius = rad;
        density = rho_total;

        for(j=0; j<nmax; j++)
        {
            rr = radius[j];
            m  = masses[j];
            acceleration[j] = G * m / (rr*rr);
            acceleration_particle[j] = acceleration[j];
        }

        //acceleration_particle = acceleration;
        nmerge = 50;
        acc_merge = acceleration[nmerge];
        rad_merge = rad[nmerge];

        for(j=0; j<nmerge; j++)
        {
            rr = radius[j];
            m  = masses[j];

            // smoothed acceleration
            acc = G * m / (rr*rr + .1* (rad_merge -rr));
            acceleration_particle[j] = acc;
        }

        //////////////////////////////////////////////
        // rederive the masses from the new particle acceleration
        radius = rad;
        dr = rad[1] - rad[0];

        // find the accelerations felt by the particles and center of mass
        radius = rad;

        for(j = 0; j<nmax; j++)
        {
            rr = radius[j];
            new_mass[j] = rr*rr * acceleration_particle[j] / G;
            new_rho[j]  = new_mass[j] / (4 * pi * rr * rr * dr);
        }


        //////////////////////////////////////////////
        // find the velocity dispersion v_r**2 using the new density and masses

        masses = new_mass;
        radius = rad;
        density = new_rho;


        for(j=0; j<nmax; j++)
        {
            p = 0.0d;
            rr = radius[j];
            dr = radius[nmax-1] / nmax;

            for(jj=j; jj<nmax; jj++)
            {
                m  = masses[jj];
                rh = density[jj];
                rr = rr + dr;
                dp = rh * G * m / (rr*rr) * dr;
                p = p + dp;
            }

            new_vr2[j] = 1/density[j] * p;
            new_vr[j] = Math.sqrt(new_vr2[j]);
        }


        //////////////////////////////////////////////
        // extend the values to large rmax
        for(j=nmax; j<ntotal; j++)
        {
            mass_total[j] = mass_total[nmax-1];
            mass_halo[j] = mass_halo[nmax-1];
            mass_disk[j] = mass_disk[nmax-1];
            mass_bulge[j] = mass_bulge[nmax-1];
            new_mass[j] = new_mass[nmax-1];

            rho_total[j] = 0.0d;
            new_rho[j]   = 0.0d;

            vr[j]      = 1d-6;
            vr2[j]     = 1d-6;
            new_vr[j]  = 1d-6;
            new_vr2[j] = 1d-6;

            m = mass_total[nmax-1];
            rr = rad[nmax-1] + dr*(j - nmax);
            rad[j] = rr;
            acc = G * m / (rr*rr);
            acceleration_particle[j] = acc;
            acceleration[j] = acc;
        }

        //////////////////////////////////////////////
        // normalize to the unit mass

        for(j=0; j<ntotal; j++)
        {
            mass_total[j]  = mass_total[j] / 7.13333d;
            mass_halo[j]   = mass_halo[j]  / 7.13333d;
            mass_disk[j]   = mass_disk[j]  / 7.13333d;
            mass_bulge[j]  = mass_bulge[j] / 7.13333d;
            new_mass[j]    = new_mass[j]   / 7.13333d;

            rho_total[j]   = rho_total[j]  / 7.13333d;
            new_rho[j]     = new_rho[j]    / 7.13333d;

            vr[j]          = vr[j]      / 7.13333d;
            vr2[j]         = vr2[j]     / 7.13333d;
            new_vr[j]      = new_vr[j]  / 7.13333d;
            new_vr2[j]     = new_vr2[j] / 7.13333d;

            rad[j]         = rad[j];

            acceleration_particle[j] = acceleration_particle[j] / 7.13333d;
            acceleration[j]          = acceleration[j]  / 7.13333d;
        }

        pscale = 1.0d;

        tmpdfind = pscale*rs_internal/rmax_scale*nnn;
    } // end initDistribution

    /**
     * Determine the index to use for interpolating the force.
     */
    protected int dfIndex(double rin, double rs)
    {
        int ival = 0;
        //double tmp = rin*pscale*rs_internal/rmax_scale;
        double tmp = rin*tmpdfind;
        //tmp = tmp*nnn + 1;
        ival = (int)tmp;
        ival = Math.min(ival, nnn-1);
        //ival = min(int(  (rin * pscale * rs_internal/ rmax_scale) * nnn + 1), nnn)
/*
        ival--; // java is 0-based
        if(ival < 0)
        {
            ival = 0;
        }
*/

        return ival;
    }

    /**
     * Compute the circular velocity for a particle at a distance r from the specified mass.
     * The rout scale of the disk and softening length, eps, are provided.
     */
    public double circularVelocity(double mass, double r, double rout, double eps)
    {
        double ftotal = 0;

        int ival = dfIndex(r, rout);

        ftotal = mass * acceleration_particle[ival] * rs_internal * rs_internal;

        double v = Math.sqrt(ftotal*r);

        return v;
    }

    /**
     * For the given particle positions and velocities, calculate the accelerations.
     */
    public void diffeq(double x[][], double f[][])
    {
//Prof.startCall("diffeq");
        double tmp1, tmp2, tmp3, tmp4, tmp5;
        tmp1=tmp2=tmp3=0;
        int n = x.length;
        int i = 0;

        double df_sigma = 0;
        double df_rho = 0;
        double c1 = 0;
        double c2 = 0;
        double xvalue = 0;
        double v1 = 0;
        double v21 = 0;

        for(i=0; i<6; i++)
        {
            xn[i] = x[n-1][i];
        }

        xp  = xn[0]; yp  = xn[1]; zp  = xn[2];
        vxp = xn[3]; vyp = xn[4]; vzp = xn[5];

        tmp4 = xp*xp+yp*yp+zp*zp;

        // distance between the two galaxies - the tidal force
        tmp5 = Math.sqrt(tmp4);

int iv1 = 0;
int iv2 = 0;
int ivn = 0;
double df1 = 0;
double df2 = 0;
double dfn = 0;
double tx[] = null;
double r22d = 0;
double r21d = 0;
double aa1 = 0;
double aa2 = 0;
double aa3 = 0;




        // get the forces, sigma and rho, and rescale them
        //ivaln[0] = dfIndex(rn[0],rrout2);
        ivn = dfIndex(tmp5,rrout2);

        df_sigma  =  new_vr2[ivn] * rs_internal2;
        df_rho    =  new_rho[ivn] * rs_internal3;

        // df
        v21 = vxp*vxp + vyp*vyp + vzp*vzp;
        v1  = Math.sqrt(v21);

        xvalue = v1 / df_sigma;
        c1 = erf(xvalue) - 2.0d * xvalue / sqrtpi * Math.exp(-xvalue*xvalue);

        // df formula with G=1
        c2  = 4.0d * Math.PI * m2 * lnl / v21;
        for(i=0; i<n; i++)
        {
            c3n[i] = 0;
        }


        for(i=pn1; i<n; i++)
        {
            c3n[i]  = c1 * c2 * df_rho;
        }

double tf[] = null;
double tv1 = 1.0d/v1;
double c3tv = 0;
double dv = 0;
        ivn  = dfIndex(tmp5, rrout2);

        for(i=0; i<n; i++)
        {
            // distance between the companion and the particle
            tx = x[i];

            //tmp1 = (x[i][0]-xp);
            //tmp2 = (x[i][1]-yp);
            //tmp3 = (x[i][2]-zp);
            tmp1 = (tx[0]-xp);
            tmp2 = (tx[1]-yp);
            tmp3 = (tx[2]-zp);
            //r22[i] = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
            r22d = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

            // distance between the main galaxy and the particle
            //tmp1 = x[i][0];
            //tmp2 = x[i][1];
            //tmp3 = x[i][2];
            tmp1 = tx[0];
            tmp2 = tx[1];
            tmp3 = tx[2];
            //r21[i] = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
            r21d = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

            //r2[i] = Math.sqrt(r22[i]);
            //r1[i] =  Math.sqrt(r21[i]);
            //r2[i] = Math.sqrt(r22d);
            //r1[i] =  Math.sqrt(r21d);
            //rn[i] = tmp5;

            r22d = Math.sqrt(r22d);
            r21d = Math.sqrt(r21d);
            //r2[i] = r22d;
            //r1[i] = r21d;
            //rn[i] = tmp5;


            //ival11[i]  = dfIndex(r1[i], rrout1);
            //ival22[i] = dfIndex(r2[i], rrout2);
            //ivaln[i]  = dfIndex(rn[i], rrout2);
            //iv1  = dfIndex(r1[i], rrout1);
            //iv2 = dfIndex(r2[i], rrout2);
            //ivn  = dfIndex(rn[i], rrout2);
            //iv1  = dfIndex(r21d, rrout1);
            //iv2 = dfIndex(r22d, rrout2);
            //ivn  = dfIndex(tmp5, rrout2);

            dv = r21d*tmpdfind;
            iv1 = (int)dv; if(iv1>nnm1)iv1=nnm1;
            dv = r22d*tmpdfind;
            iv2 = (int)dv; if(iv2>nnm1)iv2=nnm1;

            //df_force11[i] = acceleration_particle[ival11[i]] * rs_internal2;
            //df_force22[i] = acceleration_particle[ival22[i]] * rs_internal2;
            //df_forcen[i]  = acceleration_particle[ivaln[i]]  * rs_internal2;
            //df1 = acceleration_particle[iv1] * rs_internal2;
            //df2 = acceleration_particle[iv2] * rs_internal2;
            //dfn  = acceleration_particle[ivn]  * rs_internal2;

            //a1[i] = -m1 * df_force11[i];
            //a2[i] = -m2 * df_force22[i];
            //a3[i] = -m3 * df_forcen[i];
            //a1[i] = -m1 * df1;
            //a2[i] = -m2 * df2;
            //a3[i] = -m3 * dfn;
            //a1[i] = acceleration_particle[iv1]*mm1rs2;
            //a2[i] = acceleration_particle[iv2]*mm2rs2;
            //a3[i] = acceleration_particle[ivn]*mm3rs2;
            aa1 = acceleration_particle[iv1]*mm1rs2;
            aa2 = acceleration_particle[iv2]*mm2rs2;
            aa3 = acceleration_particle[ivn]*mm3rs2;

            //tmp1 = a1[i]/r1[i];
            //tmp2 = a2[i]/r2[i];
            //tmp3 = a3[i]/rn[i];
            tmp1 = aa1/r21d;
            tmp2 = aa2/r22d;
            tmp3 = aa3/tmp5;
            c3tv = c3n[i]*tv1;

            //f[i][0] = x[i][3];
            //f[i][1] = x[i][4];
            //f[i][2] = x[i][5];

            //f[i][3] = tmp1*x[i][0] + tmp2*(x[i][0]-xn[0]) + tmp3*xn[0] - c3n[i] * xn[3] / v1;
            //f[i][4] = tmp1*x[i][1] + tmp2*(x[i][1]-xn[1]) + tmp3*xn[1] - c3n[i] * xn[4] / v1;
            //f[i][5] = tmp1*x[i][2] + tmp2*(x[i][2]-xn[2]) + tmp3*xn[2] - c3n[i] * xn[5] / v1;

            tf = f[i];
            //tx = x[i];
            tf[0] = tx[3];
            tf[1] = tx[4];
            tf[2] = tx[5];

            tf[3] = tmp1*tx[0] + tmp2*(tx[0]-xn[0]) + tmp3*xn[0] - xn[3] * c3tv;
            tf[4] = tmp1*tx[1] + tmp2*(tx[1]-xn[1]) + tmp3*xn[1] - xn[4] * c3tv;
            tf[5] = tmp1*tx[2] + tmp2*(tx[2]-xn[2]) + tmp3*xn[2] - xn[5] * c3tv;
        }

// redo acceleration for last particle
r22d=1.0d;
tmp2 = aa2/r22d;
            tf[3] = tmp1*tx[0] + tmp2*(tx[0]-xn[0]) + tmp3*xn[0] - xn[3] * c3tv;
            tf[4] = tmp1*tx[1] + tmp2*(tx[1]-xn[1]) + tmp3*xn[1] - xn[4] * c3tv;
            tf[5] = tmp1*tx[2] + tmp2*(tx[2]-xn[2]) + tmp3*xn[2] - xn[5] * c3tv;


//Prof.endCall("diffeq");
    }

    /**
     * Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
     */
    public void diffq1(double x[], double f[])
    {
        double r21, r1, a1, a2, at;
        double c1, c2, c3, v21, v1, xvalue;

        int ival, ival2;
        double df_force1, df_force2;
        double df_sigma, df_rho;

        int i;
        double rr, rlocal, ivalue, dr, mmm, dm;

        r21 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        r1  = Math.sqrt(r21);

        // get the index for the calculations
        ival  =  dfIndex(r1, rrout1   );
        ival2 =  dfIndex(r1, rrout2   );

        // get the forces, sigma and rho, and rescale them
        df_force1 = acceleration_particle[ival] * rs_internal * rs_internal;
        df_force2 = acceleration_particle[ival2]* rs_internal * rs_internal;

        df_sigma  = new_vr2[ival] * rs_internal * rs_internal;
        df_rho    = new_rho[ival] * ( rs_internal * rs_internal * rs_internal );

        // interpolated forces
        a1 = -m1 * df_force1;
        a2 = -m2 * df_force2;
        at = a1 + a2;

        // df
        v21 = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
        v1  = Math.sqrt(v21);

        xvalue = v1 / df_sigma;
        c1 = erf(xvalue) - 2.0d * xvalue / sqrtpi * Math.exp(-xvalue*xvalue);

        // df formula with G=1
        c2 = -4.0d * Math.PI * m2 * lnl / v21;
        c3 = c1 * c2 * df_rho;

        f[0] = x[3];
        f[1] = x[4];
        f[2] = x[5];

        f[3] = at * x[0]/r1 - c3 * x[3]/v1;
        f[4] = at * x[1]/r1 - c3 * x[4]/v1;
        f[5] = at * x[2]/r1 - c3 * x[5]/v1;
    }
   
    /**
     * Computes the error function.  Based on code from Numerical Recipies and the following URL: 
     * http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
     * fractional error in math formula less than 1.2 * 10 ^ -7.
     * although subject to catastrophic cancellation when z in very close to 0
     * from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
     */
    public double erf(double z) 
    {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        double ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 + 
                                            t * ( 0.09678418 + 
                                            t * (-0.18628806 + 
                                            t * ( 0.27886807 + 
                                            t * (-1.13520398 + 
                                            t * ( 1.48851587 + 
                                            t * (-0.82215223 + 
                                            t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }
}
