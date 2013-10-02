package edu.gmu.cds.sim;

import edu.gmu.cds.util.MathUtil;

import java.util.Random;

public class SetupUtil
{
    public ForceModel forceModel;
    public Integrator integrator;
    public Parameters params;
    public Random rand = null;

    public SetupUtil()
    {
        rand = new Random();
    }

    public ForceModel createForceModel(int potential_type, boolean apply)
    {
        ForceModel force = null;
        switch(potential_type)
        {
            case 0:
            default:
                force = new SPMModel();
                break;
            case 1:
                force = new NBIModel();
                break;
            case 2:
                force = new MONDModel();
                break;
        }

        if(apply)
        {
            forceModel = force;
            if(integrator != null) integrator.force=force;
        }

        return force;
    }

    public void setHelpers(ForceModel forceModel, Integrator integrator, Parameters params)
    {
        this.forceModel = forceModel;
        this.params = params;
        this.integrator = integrator;
    }

    /**
     * Determine if the caller provides a parameter file or a parameter string
     */
    public void customCollision(String args[])
    {
        params.tIsSet = false;

        if(args.length > 0)
        {
            if(args.length > 1 && args[0].toLowerCase().equals("-f"))
            {
                // parse the file
                IOUtil.readParameterFile(params,args[1]);
            }
            else
            {
                // parse the state string
                IOUtil.parseStateInfoString(params,args[0]);
                params.potential_type=1;
                params.h = Parameters.hbase;
                params.tstart = -5;
                params.tend = 0;
                params.time = -5;

                if(args.length > 1)
                { 
                    double t = IOUtil.parseDouble(args[1]);
                    if(t != 0)
                    {
                        params.tstart = t;
                        params.time = t; 
                        params.tIsSet = true;
                    }
                }
            }
        }
        else
        {
            params.phi1   = 5.0d;
            params.theta1 = 5.0d;
            params.rscale1 = new double[]{1.0d,1.0d,1.0d};
            params.rout1   = 1.0d;
            params.mass1   = 1.0d;
            params.epsilon1 = 0.3d;
            params.eps1 = params.epsilon1*params.epsilon1;
            params.n1 = 1000;
            params.heat1 = 0.0;
            params.opt1 = 1;

            params.phi2   = 0.0d;
            params.theta2 = 0.0d;
            params.rscale2 = new double[]{0.30d,0.30d,0.30d};
            params.rout2   = .50d;
            params.mass2   = .50d;
            params.epsilon2 = 0.3d;
            params.eps2 = params.epsilon2*params.epsilon2;
            params.n2 = 500;
            params.heat2 = 0.0;
            params.opt2 = 1;

            params.inclination_degree =20.0d;
            params.omega_degree = 0.0d;
            params.rmin = 0.90d;
            params.velocity_factor = 0.90d;

            params.h = Parameters.hbase;
            params.time = -5;
            params.tstart = -5;
            params.tIsSet = true;
        }


        params.n = params.n1 + params.n2;
        params.eps1 = params.epsilon1*params.epsilon1;
        params.eps2 = params.epsilon2*params.epsilon2;

    }

    /**
     * Initialize the disks and set the initial positions.
     */
    public void createCollision()
    {
        profile(params.rin1, params.rout1, params.rscale1, 0, params.n1, params.mass1,
                params.eps1, params.theta1, params.phi1, params.opt1, params.heat1,
                integrator.x);        

        profile(params.rin2, params.rout2, params.rscale2, params.n1, params.n, params.mass2,
                params.eps2, params.theta2, params.phi2, params.opt2, params.heat2,
                integrator.x);        

        // determine if we need to calculate tStart
        if( !params.tIsSet ) 
        {
            double sec_vec[] = params.sec_vec;
            double rv4min[] = {sec_vec[0],sec_vec[1],sec_vec[2],-sec_vec[3],-sec_vec[4],-sec_vec[5],0.0d};
            double tminVals[] = getTStart(rv4min, params.mass1, params.mass2, params.eps1, params.eps2 , params.h,-30.0d ,10.0d*params.rout1, params.rout1, params.rout2);

            double tmpT = tminVals[0];
            if ( tmpT < -12.0 ) 
            {
                tmpT = -5;
            }

            if ( Math.abs(tmpT) < params.h)
            {
                tmpT = -5;
            }
            params.tstart = tmpT;
            params.time = tmpT;
            params.tIsSet  = true;
        }

        double mins[] = null;
        if(params.use_sec_vec)
        {
            mins = perturberPositionVec(params.sec_vec, params.mass1, params.mass2, 
                                        params.eps1, params.eps2,
                                        params.h, params.n, params.n1, params.time, integrator.x);
        }
        else
        {
            mins = perturberPosition(params.inclination_degree, params.omega_degree, 
                                     params.rmin, params.velocity_factor, params.rout1,
                                     params.mass1, params.mass2, params.eps1, params.eps2,
                                     params.h, params.n, params.n1, params.time, 
                                     integrator.x, params.sec_vec);
        }

        params.tmin = mins[0];
        params.rmin = mins[1];
        params.vrmin = mins[2];
        params.beta = mins[3];
    }

    public double randm()
    {
        return rand.nextDouble();
    }

    public double distrb(double r1, int opt, double rscale[])
    {
    	  double distrb = 0;
    	  
    	  if (opt == 1)
    	  {
    	    distrb = 1.0/r1;
    	  }
    	  else if (opt == 2)
    	  {
    	    distrb = Math.exp(-r1/rscale[0]);
    	  }
    	  else if (opt == 3)
    	  {
    	    distrb = Math.exp(-r1*r1*rscale[0] - rscale[1]*r1 - rscale[2] );
    	  }
    	  
    	  return distrb;
    }//	end function distrb


    public void profile(double rin, double rout, double rscale[], int nstart, int ntot, 
                        double mass, double eps, double theta, double phi, int opt, 
                        double heat, double x0[][])
    {
        double stheta,ctheta,sphi,cphi,pi;
        double x3,y3,z3,xv3,yv3,zv3,x2,y2,z2,xv2,yv2,zv2;
        double x,y,z,xv,yv,zv;

        int i, j, n, nprof;

        double rnorm;
        double[] rp, r, angle, v, p_ring, cp_ring;
        double st, ct, dr, ran, r1, r2, ptot;
        int[] n_ring;
        int nring, dnring, is, ie, iring, tct;

        double xp, yp, zp;
        double fx, fy, fz, ftotal;
        double tmp;

        int ival;

        n = x0.length;
if(ntot > n)
{
ntot = n;
System.err.println("ntot exceeds n");
}
        r = new double[n];
        angle = new double[n];
        v = new double[n];

        pi = Math.PI;

        stheta = Math.sin(theta*pi/180.0);    
        ctheta = Math.cos(theta*pi/180.0);    
        sphi   = Math.sin(phi*pi/180.0);    
        cphi   = Math.cos(phi*pi/180.0);    

        // set up the probability distribution for the disk
        nprof = 1000;
        nring = nprof / 10;

        dnring = nprof/nring;
        rp = new double[nprof];
        n_ring = new int[nprof];
        p_ring = new double[nprof];
        cp_ring = new double[nprof];

        // set the differential sum of the probability funtion into a vector
        rnorm = 0.0;
        dr = (rout - rin)/((double)nprof);
        for(i=0; i<nprof; i++)
        {
            r1 = i*dr + rin;
            rp[i] = distrb(r1, opt, rscale) * r1 * dr * 2.0 * pi;
            rnorm = rnorm + rp[i];
        }

        // normalize the vector
        for(i=0; i<nprof; i++)
        {
            rp[i] /= rnorm;
        }

        // take the fine bins and put them into the selection bins
        tct = 0;
        for(iring = 0; iring < nring; iring++)
        {
        //do iring =  1, nring
          is = (iring) * dnring + 1;
          ie = (iring+1) * dnring ;

          ptot = 0.0d;
          for(i=is;i<ie;i++)
          {
          //do i = is, ie
            ptot = ptot + rp[i];
          } //enddo
          p_ring[iring] = ptot;
        } //enddo

    // formulative cumulative distribution function
    //
        cp_ring[0] = p_ring[0];
        for(iring=1;iring<nring;iring++)
        {
        //do iring = 2, nring
          cp_ring[iring] = cp_ring[iring -1] + p_ring[iring];

        } //enddo


    // find the number of particles in each bin
    //
        n_ring = new int[nprof];
        for(i=nstart-1;i<ntot;i++)
        {

            // find the radial position bin
            ran = randm();
            j = 0;
           // while (j < nprof && ran > rp[j]/rnorm)
                //while (j < nprof && ran > rp[j])
            while(j<nring && ran > cp_ring[j])
            {
              j = j + 1;
            } //enddo
            j--;
            if(j<0)
           {
                j=0;
            }
            n_ring[j] = n_ring[j] + 1;
        } //enddo

        tct = 0;
        i = nstart;
        for(iring=0;iring<nring&&i<n;iring++)
        {
        //do iring =  1, nring
          is = (iring) * dnring + 1;
          ie = (iring+1) * dnring;

          r1 = (is)*dr + rin;
          r2 = (ie)*dr + rin;

          for(j=0; j<n_ring[iring]&&i<n; j++)
          {
          //do j = 1, n_ring(iring)
            ran = randm();
            r[i] = r1 + ran * (r2 - r1);
            i = i + 1;
          } //enddo
        } //enddo

    //   set the angular positions and orbital velocities
    //
        for(i=nstart;i<ntot;i++)
        {
        //do i=nstart,ntot
          angle[i] = 2.0d * pi * randm();
            v[i] = forceModel.circularVelocity(mass,r[i],rout,eps);
        }

        // set position and velocity based on the distribution parameters
        for(i= nstart; i<ntot; i++)
        {
            st  =  Math.sin(angle[i]);
            ct  =  Math.cos(angle[i]);

            x   =  ct*r[i];
            y   =  st*r[i];
            z   =  0.0;

            xv  = -v[i]*st;
            yv  =  v[i]*ct;
            zv  =  0.0;

            x2  =   x * ctheta +  z * stheta;
            y2  =   y;
            z2  =  -x * stheta +  z * ctheta;
            xv2 =  xv * ctheta + zv * stheta;
            yv2 =  yv;
            zv2 = -xv * stheta + zv * ctheta;

            x3  =  x2  * cphi -  y2 * sphi;
            y3  =  x2  * sphi +  y2 * cphi;
            z3  =  z2;
            xv3 =  xv2 * cphi - yv2 * sphi;
            yv3 =  xv2 * sphi + yv2 * cphi;
            zv3 =  zv2;
      

            x0[i][0] = x3;
            x0[i][1] = y3;
            x0[i][2] = z3;
            x0[i][3] = xv3  + randm()*heat;
            x0[i][4] = yv3  + randm()*heat;
            x0[i][5] = zv3  + randm()*heat;

        }
    }

    /**
     *
     */
    public double[] perturberPosition(double inclinationDegree, double omegaDegree, 
                                      double rMin, double velocityFactor, double rout1,
                                      double mass1, double mass2, double eps1, double eps2,
                                      double h, int n, int n1, double t0, double x0[][], double sec_vec[])
    {
        double xx0[] = new double[6];

        double omega = Math.toRadians(omegaDegree);
        double inc = Math.toRadians(inclinationDegree);

	double v = Math.sqrt(2.0d)*forceModel.circularVelocity(mass1+mass2,rMin,rout1,eps1);

        v = -v * velocityFactor;

        //      setup the transformaton based on the matrix in
        //      fundamentals of astrodynamics p. 82 by
        //      bate, mueller, and white (1971)
        //
        xx0[0] = Math.cos(omega) * rMin;
        xx0[1] = Math.sin(omega) * Math.cos(inc) * rMin;
        xx0[2] = Math.sin(omega) * Math.sin(inc) * rMin;

        xx0[3] = -Math.sin(omega) * v;
        xx0[4] =  Math.cos(omega) * Math.cos(inc) * v;
        xx0[5] =  Math.cos(omega) * Math.sin(inc) * v;

        sec_vec[0] = xx0[0];
        sec_vec[1] = xx0[1];
        sec_vec[2] = xx0[2];
        sec_vec[3] = -xx0[3];
        sec_vec[4] = -xx0[4];
        sec_vec[5] = -xx0[5];

        return perturberPositionVec(sec_vec, mass1, mass2, eps1, eps2, h, n, n1, t0, x0);
    }

    public double[] perturberPositionVec(double xx0[], double mass1, double mass2, 
                                         double eps1, double eps2, 
                                         double h, int n, int n1, double t0, double x0[][])
    {
        double xxe[] = new double[6];
        double tcurrent = 0;
        int i = 0;

        // reverse the velocity for backward integration
        xx0[3] = -xx0[3];
        xx0[4] = -xx0[4];
        xx0[5] = -xx0[5];

        double tmin = tcurrent;
        double tmpr = 0;
        double tmpv = 0;

        // avoid multiple calls to sqrt, save it until the end
        double rmin = xx0[0]*xx0[0]+xx0[1]*xx0[1]+xx0[2]*xx0[2];
        double vrmin = xx0[3]*xx0[3]+xx0[4]*xx0[4]+xx0[5]*xx0[5];

        // now move position back to t0 from t=0.0
        while( t0 < tcurrent )
        {
            integrator.rk41(xx0, xxe, h);
            xx0[0] = xxe[0];  xx0[1] = xxe[1]; xx0[2] = xxe[2];
            xx0[3] = xxe[3];  xx0[4] = xxe[4]; xx0[5] = xxe[5];

            tcurrent = tcurrent - h;

            tmpr = xx0[0]*xx0[0]+xx0[1]*xx0[1]+xx0[2]*xx0[2];
            if(tmpr < rmin)
            {
                rmin = tmpr;
                vrmin = xx0[3]*xx0[3]+xx0[4]*xx0[4]+xx0[5]*xx0[5];
                tmin = tcurrent;
            }
        }

        // reverse the velocity for forward integration
        xx0[3] = -xx0[3];
        xx0[4] = -xx0[4];
        xx0[5] = -xx0[5];

        // now adjust the test particles from the second disk
        // to the proper velocity and positions
        if(params.n > params.n1)
        {
            for(i=params.n1; i<params.n; i++)
            {
                x0[i][0] += xx0[0];
                x0[i][1] += xx0[1];
                x0[i][2] += xx0[2];
                x0[i][3] += xx0[3];
                x0[i][4] += xx0[4];
                x0[i][5] += xx0[5];
            }
        }

        // include the perturbing galaxy
        i = params.n;
        x0[i][0] += xx0[0];
        x0[i][1] += xx0[1];
        x0[i][2] += xx0[2];
        x0[i][3] += xx0[3];
        x0[i][4] += xx0[4];
        x0[i][5] += xx0[5];

        vrmin = Math.sqrt(vrmin);
        // beta = (m1+m2)/(rmin*rmin + vrmin)
        double beta = (mass1+mass2)/(rmin*vrmin);
        rmin = Math.sqrt(rmin);
        return new double[]{tmin,rmin,vrmin,beta};
    }

    /**
     * Convert the position and velocity to classical orbital
     * elements.
     * 
     *         coe[0] = p;
        coe[1] = ecc;
        coe[2] = inc;
        coe[3] = LAN;
        coe[4] = w;
        coe[5] = v;
	coe[6] = u;        

     * @param r
     * @param v
     * @param mu
     * @return
     */
    public static double[] rvToCoe(double r[], double v[], double mu)
    {
        double muInv = 1.0/mu;

        double rmag = MathUtil.mag(r);
        double vmag = MathUtil.mag(v);
    
        double[] h = MathUtil.cross(r,v);
        double hmag = MathUtil.mag(h);
   
        double K[] = {0.0,0.0,1.0}; 
        double[] n = MathUtil.cross(K,h);
        double nmag = MathUtil.mag(n);
    
        double tmp1 = vmag*vmag - mu/rmag;
        double tmp2 = MathUtil.dot(r,v);
    
        double[] v1 = MathUtil.scale(tmp1,r);
        double[] v2 = MathUtil.scale(tmp2,v);
    
        double[] ev = MathUtil.sub(v1,v2);
        ev = MathUtil.scale(muInv,ev);
    
        double p = hmag*hmag*muInv;
        double ecc = MathUtil.mag(ev);
        double cosi = h[2]/hmag;
        double cosO = n[0]/nmag;
        double cosw = MathUtil.dot(n,ev)/(nmag*ecc);
        double cosv = MathUtil.dot(ev,r)/(ecc*rmag);
        double cosu = MathUtil.dot(n,r)/(nmag*rmag);

        double outv[] = new double[7];
        outv[0] = p;
        outv[1] = ecc;
        outv[2] = Math.acos(cosi);
    
        tmp1 = Math.acos(cosO);
            
        if (n[0] < 0)
        {
            tmp1 = 2.0*Math.PI-tmp1;
        }
        outv[3] = tmp1; 
        
        tmp1 = Math.acos(cosw);
    
        if (ev[1] < 0 ) 
        {
            tmp1 = 2.0*Math.PI-tmp1;
        }
        outv[4] = tmp1; 
    
        tmp1 = Math.acos(cosv);
    
        if (MathUtil.dot(r,v) < 0)
        {
            tmp1 = 2.0*Math.PI-tmp1;
        }
        outv[5] = tmp1; 
        
        if(cosu > 1.0 || cosu < -1.0) 
        {
            outv[6] = Double.NaN;
        }
        else
        {
            tmp1 = Math.acos(cosu);
            if(r[1]>0) 
            {
                tmp1 = 2.0*Math.PI-tmp1;
            }
            outv[6] = tmp1;
        }

        return outv;
    }    

   /**
    
	! * Find the time of rmin, assuming earlier than now, given
	! * the r and v values.  Returns r and v at time of rmin
	! * by replacing r and v.  r and v are given as
	! * {rx,ry,rz,vx,vy,vz}.
*/	
    public double[] getTStart(double rv[], double mass1, double mass2, double eps1, double eps2, double h, double tmin, double mind, double rout1, double rout2)
    {
        int i = 0;

        double mu = mass1+mass2;
        double t = 0.0;
		
        double []r = {rv[0], rv[1], rv[2]};
        double []v = {-rv[3],-rv[4], -rv[5]};
        double []coe = rvToCoe(r,v,mu);
        double ecc = coe[1];
        double a = coe[0]/(1.0-ecc*ecc);
        double period = 0.0d;
        double apocenter = a*(1.0+ecc);
        double a2 = apocenter*apocenter;
        double tApp = 0.0;
    	
        boolean isEllipse = false;
    	    	
        if (ecc < 1.0)
        {
            isEllipse = true;
            period = 2.0*Math.PI/Math.sqrt(mu)*Math.pow(a,1.5d);
            period = period * 1.0;
        }
    
    	
    	double xxe[] = new double[7];
    	double rvmin[] = new double[7];
    	
    	for(i=0;i<7;i++)
	{
	    rvmin[i] = rv[i];
	}	
      
        double distNew = rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2];
        double distOld = 2.0*distNew;
    	
        double distNearApp = -1e30;
   
        // keep looping as long as distance is decreasing
        while(tmin < t) 
        {
            r[0] = rv[0]; r[1] = rv[1]; r[2] = rv[2];
            v[0] = rv[3]; v[1] = rv[4]; v[2] = rv[5];
            coe = rvToCoe(r,v,mu);
            xxe[6] =t +h;
            //call wrap_rk41(rv, h, mass1, mass2, eps1, eps2, xxe)
            integrator.rk41(rv, xxe, h);
    	   
            distNew = xxe[0]*xxe[0] + xxe[1]*xxe[1] + xxe[2]*xxe[2];
    	    
            // if it's ellipse and it's near apocenter, take this time
            if ( isEllipse && (Math.abs(distNew-a2)/a2 < 0.05d) )
            {
              if(distNew > distNearApp)
              {
                distNearApp = distNew;
                tApp = t;
              }
            }
    	    
            if (distNew < distOld)
            {
              distOld = distNew;
    	      for(i=0;i<7;i++)
	      {
	          rvmin[i] = xxe[i];
	      }	
              rvmin[6] = rvmin[6] - h;
            }

    	    for(i=0;i<7;i++)
	    {
	      rv[i] = xxe[i];
	    }	
            rv[6] = xxe[6] - h * 2.0d;
            t = t - h;
        } // end loop

    	    for(i=0;i<7;i++)
	    {
	      rv[i] = rvmin[i];
	    }	
    	
      double minDist = Math.sqrt(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2]);
      double minVel = Math.sqrt(rv[3]*rv[3] + rv[4]*rv[4] + rv[5]*rv[5]);
    	
      t = rv[6];
    	
      if(isEllipse && tApp < 0.0d)
      {
        t = tApp;
      }
      else
      {
        t = t - mind/minVel;
      }
    	
        double outStuff[] = {t,minDist,minVel,rv[6]};
    	
    	
        return outStuff;
    }
    
}
