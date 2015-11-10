/*
Fortran Version
  Copyright 2008 by John Wallin
  Copyright 1987, 1989 by Starchild Software

Java Version 
  Copyright 2008 by Anthony Holincheck
  
Licensed under the Academic Free License version 3.0
See the file "LICENSE" for more information
*/
package edu.gmu.cds.sim.trim;

import java.util.Random;

/**
* Java version of SPAM written in FORTRAN by John Wallin.
* 
* @author Anthony Holincheck
*
*/
public class JSpam
{
	public static final double pi = Math.PI;
	public static final double pid2 = pi/2.0;
	public static final double sqrtpi = Math.sqrt(Math.PI);
	public static final double twosqrtpicubed = 2*sqrtpi*sqrtpi*sqrtpi;
  public static final double od6 = 1.0d/6.0d;
  public static final double od3 = 1.0d/3.0d;
      
  public static final int nnn = 10000;
	public static final double rs_internal = 10.0d;
	public static final double rs_internal2 = rs_internal*rs_internal;
	public static final double rs_internal3 = rs_internal*rs_internal*rs_internal;
      
	public Random rand = null;
	
	public JSpam()
	{
		rand = new Random();
	}
	
	/**
	 * Based on findRmin - computes start time.
	 * 
	 * Find the time of rmin, assuming earlier than now, given
	 * the r and v values.  Returns r and v at time of rmin
	 * by replacing r and v.  r and v are given as
	 * {0,rx,-vx,ry,-vy,rz,-vz}.  It's important that
	 * velocity be in reverse time direction.
	 * 
	 * @param rv
	 * @param mass1
	 * @param mass2
	 * @param eps1
	 * @param eps2
	 * @param h
	 * @param tmin
	 * @return
	 */
  public double[] getTStart(double rv[], double mass1, double mass2, 
  		                  double eps1, double eps2, double h, 
  		                  double tmin, double mind, double rout1, double rout2)
  {
  	int i;
  	double t = 0;
  	double distOld = 0.0;
  	double distNew = 0.0;
  	/*
  	
  	System.out.println("rv[0] = " + rv[0]);
  	System.out.println("rv[1] = " + rv[1]);
  	System.out.println("rv[2] = " + rv[2]);
  	System.out.println("rv[3] = " + rv[3]);
  	System.out.println("rv[4] = " + rv[4]);
  	System.out.println("rv[5] = " + rv[5]);
  	System.out.println("rv[6] = " + rv[6]);
  	System.out.println("mass1 = " + mass1);
  	System.out.println("mass2 = " + mass2);
  	System.out.println("eps1 = " + eps1);
  	System.out.println("eps2 = " + eps2);
  	System.out.println("h = " + h);
  	System.out.println("tmin = " + tmin);
  	System.out.println("mind = " + mind);
  	System.out.println("rout1 = " + rout1);
  	System.out.println("rout2 = " + rout2);
  	*/
  	double mu = mass1+mass2;
  	
  	double r[]={rv[1],rv[3],rv[5]};
  	double v[]={-rv[2],-rv[4],-rv[6]};
  	double coe[] = OrbitUtil.rvToCoe(r,v,mu);
  	
  	
  	double ecc = coe[1];
  	double a = coe[0]/(1.0d-ecc*ecc);
  	double period = 0.0d;
  	double apocenter = a*(1+ecc);
  	double a2 = apocenter*apocenter;
  	double tApp = 0.0d;
  	
  	boolean isEllipse = false;
  	
  	if(ecc < 1.0d)
  	{
  		isEllipse = true;
  		period = 2*Math.PI/Math.sqrt(mu)*Math.pow(a,1.5);
  		period *= 1.0d;
  	}
  	
  	double xxe[] = new double[7];
  	double rvmin[] = new double[7];
  	
  	for(i=0;i<7;i++)
	    {
	        rvmin[i] = rv[i];
	    }
  	
  	//distNew = Math.sqrt(rv[1]*rv[1]+ rv[3]*rv[3] + rv[5]*rv[5] );
  	distNew = rv[1]*rv[1]+ rv[3]*rv[3] + rv[5]*rv[5];
  	distOld = 2.0*distNew;
  	
  	double distNearApp = -Double.MAX_VALUE;
  	
  	// keep looping as long as distance is decreasing
  	while(tmin < t) 
  	{
  		coe = OrbitUtil.rvToCoe(rv[1],rv[3],rv[5],rv[2],rv[4],rv[6],mu);
  		rk41(rv, h, mass1, mass2, eps1, eps2, xxe, rout1, rout2);
  		distNew = xxe[1]*xxe[1]+ xxe[3]*xxe[3] + xxe[5]*xxe[5];
  	    
  	    // if it's ellipse and it's near apocenter, take this time
  	    if(isEllipse && Math.abs(distNew-a2)/a2 < 0.05)
  	    {
  	    	if(distNew > distNearApp)
  	    	{
  	    		distNearApp = distNew;
  	    		tApp = t;
  	    	}
  	    }
  	    
  	    if(distNew < distOld)
  	    {
  	    	distOld = distNew;
  	    	
      	    for(i=0;i<7;i++)
      	    {
      	        rvmin[i] = xxe[i];
      	    }    	 
      	    rvmin[0] -= h;
  	    }

  	    for(i=0;i<7;i++)
  	    {
  	        rv[i] = xxe[i];
  	    }
  	    
  	    rv[0] = xxe[0] - h * 2.0d;
  	    t -= h;
      }

  	for(i=0;i<7;i++)
	    {
	        rv[i] = rvmin[i];
	    }
  	
  	//distNew = Math.sqrt(rv[1]*rv[1]+ rv[3]*rv[3] + rv[5]*rv[5] );
  	double minDist = Math.sqrt(rv[1]*rv[1] + rv[3]*rv[3] + rv[5]*rv[5]);
  	// velocity at min distance
  	double minVel = Math.sqrt(rv[2]*rv[2] + rv[4]*rv[4] + rv[6]*rv[6]);
  	
  	t = rv[0];
  	
  	// if it's an ellipse take the time of apogee
  	if(isEllipse && tApp < 0.0d)
  	{
  		t = tApp;
  	}
  	else
  	{
  	    // take the velocity at time of rmin and subtract
  		// time it would take to go another mind away from rmin
  		t -= mind/minVel;
  	}
  	
  	//System.out.println("\n\ntmin = " + t + "\trmin = " + distNew + "\n\n");
  	
  	double out[] = {t,minDist,minVel,rv[0]};
  	
  	return out;
  }
	

  	
	/**
	 * Find the time of rmin, assuming earlier than now, given
	 * the r and v values.  Returns r and v at time of rmin
	 * by replacing r and v.  r and v are given as
	 * {0,rx,-vx,ry,-vy,rz,-vz}.  It's important that
	 * velocity be in reverse time direction.
	 * 
	 * @param rv
	 * @param mass1
	 * @param mass2
	 * @param eps1
	 * @param eps2
	 * @param h
	 * @param tmin
	 * @return
	 */
  public double[] findRmin(double rv[], double mass1, double mass2, 
  		                   double eps1, double eps2, double h, double tmin,
  		                   double rout1, double rout2)
  {
  	int i;
  	double t = 0;
  	double distOld = 0.0;
  	double distNew = 0.0;
  	double xxe[] = new double[7];
  	double rvmin[] = new double[7];
  	
  	for(i=0;i<7;i++)
	    {
	        rvmin[i] = rv[i];
	    }
  	
  	//distNew = Math.sqrt(rv[1]*rv[1]+ rv[3]*rv[3] + rv[5]*rv[5] );
  	distNew = rv[1]*rv[1]+ rv[3]*rv[3] + rv[5]*rv[5];
  	distOld = 2.0*distNew;
  	
  	// keep looping as long as distance is decreasing
  	while(tmin < t) 
  	{
  		rk41(rv, h, mass1, mass2, eps1, eps2, xxe, rout1, rout2);
  	    //distNew = Math.sqrt(xxe[1]*xxe[1]+ xxe[3]*xxe[3] + xxe[5]*xxe[5] );
  	    distNew = xxe[1]*xxe[1]+ xxe[3]*xxe[3] + xxe[5]*xxe[5];
//System.out.println(t + "\t" + Math.sqrt(distNew) + "\t" + Math.sqrt(xxe[2]*xxe[2]+ xxe[4]*xxe[4] + xxe[6]*xxe[6]));
  	    if(distNew < distOld)
  	    {
  	    	distOld = distNew;
  	    	
      	    for(i=0;i<7;i++)
      	    {
      	        rvmin[i] = xxe[i];
      	    }    	 
      	    rvmin[0] -= h;
  	    }

  	    for(i=0;i<7;i++)
  	    {
  	        rv[i] = xxe[i];
  	    }
  	    
  	    rv[0] = xxe[0] - h * 2.0d;
  	    t -= h;
      }

  	for(i=0;i<7;i++)
	    {
	        rv[i] = rvmin[i];
	    }
  	
  	//distNew = Math.sqrt(rv[1]*rv[1]+ rv[3]*rv[3] + rv[5]*rv[5] );
  	
  	t = rv[0];
  	
  	//System.out.println("\n\ntmin = " + t + "\trmin = " + distNew + "\n\n");
  	
  	double minDist = Math.sqrt(rv[1]*rv[1] + rv[3]*rv[3] + rv[5]*rv[5]);
  	// velocity at min distance
  	double minVel = Math.sqrt(rv[2]*rv[2] + rv[4]*rv[4] + rv[6]*rv[6]);
  	
  	double out[] = {t,minDist,minVel};
  	
  	return out;
  }
  
  public void perturberPosition(double inclination_degree, double omega_degree, 
  		                      double rmin, double velocity_factor, double mass1,
  		                      double mass2, double eps1, double eps2, double h,
  		                      int n, int n1, double t0, double x0[][],
  		                      double rout1, double rout2)
  {
  	double en, v1;
  	double xx0[] = new double[7];
     	double omega, incl;
  	
  	// change inclination and omega into radians
  	incl = inclination_degree * pi/180.0d;
  	omega = omega_degree * pi/180.0d;

  	// energy from mass1
  	if ( eps2 > 0.0d)
  	{
  		en = mass1/eps1 * (pid2 - Math.atan(rmin/eps1 ) ); 
  	}
  	else
  	{
  		en = mass1 / rmin;
  	}

  	// energy from mass2
  	if ( eps1 > 0)
  	{
  		en = en + mass2/eps1 * (pid2 - Math.atan(rmin/eps2 ) );
  	}
  	else
  	{
  	 	en = en +  eps2 / rmin;
  	}

  	// calculate escape velocity and velocity at rmin
  	v1 = Math.sqrt(2.0d * en);
  	v1= -v1 * velocity_factor;


  	//	setup the transformaton based on the matrix in
  	//	fundamentals of astrodynamics p. 82 by
  	//	bate, mueller, and white (1971)
  	//
  	xx0[0] = 0.0d;
  	xx0[1] = Math.cos(omega) * rmin;
  	xx0[3] = Math.sin(omega) * Math.cos(incl) * rmin;
  	xx0[5] = Math.sin(omega) * Math.sin(incl) * rmin;

  	xx0[2] = -Math.sin(omega) * v1;
  	xx0[4] =  Math.cos(omega) * Math.cos(incl) * v1;
  	xx0[6] =  Math.cos(omega) * Math.sin(incl) * v1;
  	
  	perturberPosition(xx0,mass1,mass2,eps1,eps2,h,n,n1,t0,x0,rout1,rout2);
  }
  
  /**
   * Update the position of the stars in the secondary galaxy.
   * 
   * xx0 should be 0,rx,-vx,ry,-vy,rz,-vz.  That is rmin and
   * v at rmin (just in opposite direction).
   * 
   * @param xx0
   * @param mass1
   * @param mass2
   * @param eps1
   * @param eps2
   * @param h
   * @param n
   * @param n1
   * @param t0
   * @param x0
   */
  public void perturberPosition(double xx0[], double mass1,
          double mass2, double eps1, double eps2, double h,
          int n, int n1, double t0, double x0[][],
          double rout1, double rout2)
  {
  	int i = 0;
  	int j = 0;
  	double xxe[] = new double[7];
  	//double dist1;
  	
  	//
  	//       now move position back to t0 from t=0.0
  	//
  	/**/
  	while(t0 < xx0[0]) 
  	{
  		rk41(xx0, h, mass1, mass2, eps1, eps2, xxe, rout1, rout2);
  	    //dist1 = Math.sqrt(xx0[1]*xx0[1]+ xx0[3]*xx0[3] + xx0[5]*xx0[5] );

  	    for(i=0;i<7;i++)
  	    {
  	        xx0[i] = xxe[i];
  	    }
  	    
  	    xx0[0] = xxe[0] - h * 2.0d;
      }
  	/**/
  	//
  	//       set the time to t0
  	xx0[0] = t0;
  	xx0[2] = -xx0[2];
  	xx0[4] = -xx0[4];
  	xx0[6] = -xx0[6];
  	
  	//	now move the test particles from the 
  	//	second disk to the proper velocity and positions
  	//
  	if (n > n1)
  	{
  		for(i=n1;i<n;i++)
  		{
  	    	for(j=0;j<7;j++)
  	    	{
  	    	    x0[i][j] = x0[i][j] + xx0[j];
  	    	}
  	    } 
  	}
  	
  	//
  	// 	include the perturbing galaxy
  	//
  	n = n + 1;
  	
  	for(i=0; i<7; i++)
  	{
  	    x0[n-1][i] = xx0[i];
  	}

  	        
  	// set the particle times
  	// todo figure out length of x0
  	int np = x0.length;
  	for(i=0; i<np; i++)
  	{
  	    x0[i][0]=t0;
  	}

      return;

  } //end subroutine perturber_position


  //---------------------------------------------------------
  // subroutine rk41(xx0, h, mass1, mass2, eps1, xxe)
  public void rk41(double xx0[], double h, double mass1, double mass2,
  		         double eps1, double eps2, double xxe[], double rout1, double rout2)	  
  {
  	double x[] = new double[7];
  	double f[] = new double[7];
      
      int i=0;
      int n=7;
      
      for(i=0;i<n;i++)
      {
      	x[i] = xx0[i];
      }
      
  	diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	
  	for(i=0;i<n;i++)
  	{
          xxe[i] = xx0[i] + h * f[i] * od6;
          x[i]   = xx0[i] + h * f[i] * 0.5d;
  	}
  	diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	  
  	for(i=0;i<n;i++)
  	{
  		xxe[i] = xxe[i] + h * f[i] * od3;
  		x[i]   = xx0[i] + h * f[i] *0.5d;
  	}
  	diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	  
      for(i=0;i<n;i++)
  	{
      	xxe[i] = xxe[i] + h * f[i] * od3;
      	x[i]   = xx0[i] + h*f[i];
  	}
      diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	  
  	for(i=0;i<n;i++)
  	{
  		xxe[i] = xxe[i] + h * f[i] * od6;
  	}  
  	
      /*
      x[0] = xx0[0];
      x[1] = xx0[1];
      x[2] = xx0[2];
      x[3] = xx0[3];
      x[4] = xx0[4];
      x[5] = xx0[5];
      x[6] = xx0[6];
      
      
  	diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);

      xxe[0] = xx0[0] + h * f[0] * od6;
      x[0]   = xx0[0] + h * f[0] * 0.5d;
      xxe[1] = xx0[1] + h * f[1] * od6;
      x[1]   = xx0[1] + h * f[1] * 0.5d;
      xxe[2] = xx0[2] + h * f[2] * od6;
      x[2]   = xx0[2] + h * f[2] * 0.5d;
      xxe[3] = xx0[3] + h * f[3] * od6;
      x[3]   = xx0[3] + h * f[3] * 0.5d;
      xxe[4] = xx0[4] + h * f[4] * od6;
      x[4]   = xx0[4] + h * f[4] * 0.5d;
      xxe[5] = xx0[5] + h * f[5] * od6;
      x[5]   = xx0[5] + h * f[5] * 0.5d;
      xxe[6] = xx0[6] + h * f[6] * od6;
      x[6]   = xx0[6] + h * f[6] * 0.5d;

  	diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	  
  	for(i=0;i<n;i++)
  	{
  		xxe[i] = xxe[i] + h * f[i] * od3;
  		x[i]   = xx0[i] + h * f[i] *0.5d;
  	}
  	diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	  
      for(i=0;i<n;i++)
  	{
      	xxe[i] = xxe[i] + h * f[i] * od3;
      	x[i]   = xx0[i] + h*f[i];
  	}
      diffq1(x, mass1, mass2, eps1, eps2, f, rout1, rout2);
  	  
  	for(i=0;i<n;i++)
  	{
  		xxe[i] = xxe[i] + h * f[i] * od6;
  	}  
  	*/
  	return;
  }//	end subroutine rk41

 	//-----------------------------------------------------------
  //	subroutine diffq1(x, mass1, mass2, eps1, f)
  public void diffq1(double x[], double mass1, double mass2, double eps1, double eps2, double f[], double rout1, double rout2)
  {
  	double r21, r1, a1, m1, invr1, a2, v21, invv1;
  	r21=x[1]*x[1]+x[3]*x[3]+x[5]*x[5];
  	r1=Math.sqrt(r21);
  	invr1 = 1.0/r1;
  	
      a1 = -mass1 / (r21+ eps1*eps1) - mass2 / (r21 + eps2*eps2);


      
  	f[0] = 1.0d;
  	f[1] = x[2];
  	//f[2] = a1 * x[1] / r1;
  	f[2] = a1 * x[1]*invr1;
  	f[3] = x[4];
  	//f[4] = a1 * x[3]/r1;
  	f[4] = a1 * x[3]*invr1;
  	f[5] = x[6];
  	//f[6] = a1 * x[5]/r1;
  	f[6] = a1 * x[5]*invr1;


  	
  	return;
  } //  	end subroutine diffq1
  	       

  	  
  	//------------------------------------------------------
  	//
  	//
  	//
 // 	subroutine profile(rin, rout, rscale, nstart, ntot, mass, eps, &
 // 	     theta, phi, opt, heat, seed, t0, x0)

  	//
  	//
  	//       variables-
  	//       opt - option for the distribution
  	//       rin - inner radius
  	//       rout - outer radius
  	//       rscale - scale of brightness drop
  	//	     nstart - start number for placement of particles
  	//       ntot - number of particles to be placed
  	//       heat - heat parameter
  	//       m - mass of galaxy
  	//       sl - softening length
  	//       nring - number of rings
  	//       npart - number of particle per ring (opt)
  	//       x0 - position of center of mass
  	//
  // 	subroutine profile(rin, rout, rscale, nstart, ntot, mass, eps, &
  // 	     theta, phi, opt, heat, seed, t0, x0)
  public double[] profile(double rin, double rout, double rscale[], int nstart, int ntot,
  		                double mass, double eps, double theta, double phi, int opt,
  		                double heat, double seed, double t0, double x0[][])
  {
  	int len = ntot-nstart+1;
  	double rVals[] = new double[len];
  	double vVals[] = new double[len];
  	
  	double stheta,ctheta,sphi,cphi;
  	double x3,y3,z3,xv3,yv3,zv3,x2,y2,z2,xv2,yv2,zv2;
  	double x,y,z,xv,yv,zv;
  	
  	int i,j,n,nprof;
  	
  	double rnorm;
  	double[] rp, r, angle, v;
  	double st, ct, dr, ran, r1;


  	n = x0.length;
  	r = new double[n];
  	angle = new double[n];
  	v = new double[n];

  	  
      stheta = Math.sin(theta*pi/180.0d);
  	ctheta = Math.cos(theta*pi/180.0d);
  	sphi   = Math.sin(phi*pi/180.0d);
  	cphi   = Math.cos(phi*pi/180.0d);

  	//
  	//       set up the probablity distribution for the disk

  	nprof = 250;
  	rp = new double[nprof];

  	rnorm = 0.0d;
  	
  	dr = (rout - rin)/(double)(nprof);
  	
  	// TODO check limits
  	for(i=0;i<nprof;i++)
  	{
  	    r1 = i*dr + rin;
  	    rp[i] =  distrb(r1, opt, rscale ) * r1 * dr * 2.0d* pi;
  	    rnorm = rnorm + rp[i];
  	    rp[i] = rnorm;
  	}


  	//	set up radial and angular positions
  	//
//  	 TODO check limits
  	for(i=nstart-1;i<ntot;i++)
  	{
  	
  	    // find the radial position bin
  	    ran = randm(seed);
  	    j = 0;
  	    while (ran > rp[j]/rnorm && j < nprof)
  	    {
  	      j = j + 1;
  	    }
  	    
  	    // set the position and angle
  	    r[i] = j*dr + rin;
  	    angle[i] = 2.0d* pi * randm(seed);

  	    // set the orbital velocity for this particle
  	    v[i] = Math.sqrt( mass * r[i]/ ( r[i]*r[i] + eps*eps ) );
  	}

  	// sort r and v in order of ascending r
  	
  	System.arraycopy(r,nstart-1,rVals,0,len);
  	System.arraycopy(v,nstart-1,vVals,0,len);
  	shell2(rVals,vVals);
  	System.arraycopy(rVals,0,r,nstart-1,len);
  	System.arraycopy(vVals,0,v,nstart-1,len);
  	
  	
  	// set position and velocity based on the distribution parameters
  	for(i=nstart-1;i<ntot;i++)
  	{
  	    
  	    st=Math.sin(angle[i]);
  	    ct=Math.cos(angle[i]);
  	    
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
  	    
  	    
  	    x0[i][1] = x3; 
  	    x0[i][3] = y3; 
  	    x0[i][5] = z3;
  	    
  	    //x0[i][2] = xv3; 
  	    //x0[i][4] = yv3;
  	    //x0[i][6] = zv3;
  	    
  	    x0[i][2] = xv3  + randm(seed)*heat*v[i];
  	    x0[i][4] = yv3  + randm(seed)*heat*v[i];
  	    x0[i][6] = zv3  + randm(seed)*heat*v[i];
  	    x0[i][0] = t0;    
  	}

  	  //deallocate(r)
  	  //deallocate(angle)
  	  //deallocate(v)
  	  //deallocate(rp)
  	
  	return rVals;
  }// 	end subroutine profile

  public double circularVelocity(double rout, double mass, double eps)
  {
  	return Math.sqrt( mass * rout/ ( rout*rout + eps*eps ) );
  }
  
  	//
  	//
  	//-----------------------------------------------------------------
  	//
  	//
  	//function distrb(r1,opt,rscale)
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

	//
	//
	//-----------------------------------------------------------------
	//
	//
	public double randm(double seed)
	{
	  /*
	  real (kind=8) :: randm
	  real (kind=8) seed

	  call random_number(randm)
	  return
	
	//end function randm
	  */
	  return rand.nextDouble();
  }
	
	/**
	 * From numerical recipes.
	 * Sorts an array a[1..n] into ascending numerical order by Shell's method
	 * (diminishing increment sort). n is input; a is replaced on output by its
	 * sorted rearrangement.
	 */
	public void shell2(double a[], double b[])
	{
		int n = a.length;
	    int i,j,inc;
	    double v, v2;
	    inc=1;          //Determine the starting increment.
	    do 
	    {
	    	inc *= 3;
			inc++;
		} while (inc <= n);
	    
	    do 
	    { 				// Loop over the partial sorts.
			inc /= 3;
			for (i=inc+1;i<=n;i++) 
			{           // Outer loop of straight insertion.
				v=a[i-1];
				v2=b[i-1];
				j=i;
				while (a[j-inc-1] > v) 
				{ 		// Inner loop of straight insertion.
					a[j-1]=a[j-inc-1];
					b[j-1]=b[j-inc-1];
					j -= inc;
					if (j <= inc) break;
				}
				a[j-1]=v;
				b[j-1]=v2;
			}
		} while (inc > 1);
	}
	
  /**
   * Calculate the intensity map for the specified output index
   * frame.
   *  
   * @param n
   * @param ind
   * @param limits
   * @return
   */
  public double[][] getIntensityMap(int n, double x0[][], double alpha, double beta, double min, double max, double mass1, double mass2, int n1, int n2, double r1Vals[], double r2Vals[], double r1, double r2)
  {
  	int offset = n1;
  	int nn1 = 0;
  	int nn2 = 0;
  	
  	// find cutoff based upon target r1
  	nn1 = n1;
  	for(int i=0; i<n1; i++)
  	{
  		if(r1Vals[i]>r1)
  		{
  			nn1 = i;
  			break;
  		}
  	}
  	
  	// find cutoff based upon target r2
  	nn2 = n2;
  	for(int i=0; i<n2; i++)
  	{
  		if(r2Vals[i]>r2)
  		{
  			nn2 = i;
  			break;
  		}
  	}
  	
  	return getIntensityMap(n, x0, alpha, beta, min, max, mass1, mass2, nn1, nn2, offset);
  }
  
  /**
   * Calculate the intensity map for the specified output index
   * frame.
   *  
   * @param n
   * @param ind
   * @param limits
   * @return
   */
  public double[][] getIntensityMap(int n, double x0[][], double alpha, double beta, double min, double max, double mass1, double mass2, int n1, int n2, int offset)
  {
      double map[][] = null;
      double tmpMap[][] = new double[n][n];
      
  	double ca, sa, cb, sb;
		ca = Math.cos(Math.toRadians(alpha));
		sa = Math.sin(Math.toRadians(alpha));
		cb = Math.cos(Math.toRadians(beta));
		sb = Math.sin(Math.toRadians(beta));
		
      double x,y;
      double xt,yt,zt,xp,yp;

      double rng = (max-min);
      double div = rng/n;
      double ninmap = 0;

      int size = n1;
      
      for(int i=0; i<size; i++)
      {
  		xp=x0[i][1];
  		yp=x0[i][3];
  		zt=x0[i][5];
  		xt = xp*cb - yp*sb;
  		yt = xp*sb + yp*cb;
  		x = xt*ca - zt*sa;
  		y = yt;
  		ninmap += addToMap(tmpMap,x,y,n,rng,div,min,max);
      }
  	
      ninmap *= mass1/((double)n1);
      for(int i=0; i<n; i++)
      {
      	for(int j=0; j<n; j++)
      	{
      		tmpMap[i][j] /= ninmap;
      	}
      }

      map = tmpMap;
      
      // now do second galaxy
      tmpMap = new double[n][n];
      ninmap = 0.0;
      
      size = offset+n2; 
      for(int i=offset; i<size; i++)
      {
  		xp=x0[i][1];
  		yp=x0[i][3];
  		zt=x0[i][5];
  		xt = xp*cb - yp*sb;
  		yt = xp*sb + yp*cb;
  		x = xt*ca - zt*sa;
  		y = yt;
  		ninmap += addToMap(tmpMap,x,y,n,rng,div,min,max);
      }
  	
      ninmap *= mass2/((double)n2);
      for(int i=0; i<n; i++)
      {
      	for(int j=0; j<n; j++)
      	{
      		tmpMap[i][j] /= ninmap;
      	}
      }
      
      // add maps and get total
      ninmap = 0.0d;
      
      for(int i=0; i<n; i++)
      {
      	for(int j=0; j<n; j++)
      	{
      		map[i][j] += tmpMap[i][j];
      		ninmap += map[i][j];
      	}
      }
      
      for(int i=0; i<n; i++)
      {
      	for(int j=0; j<n; j++)
      	{
      		map[i][j] /= ninmap;
      	}
      }
      /*
      int size = x0.length;
      
      double ninmap = 0;
      for(int i=0; i<size; i++)
      {
  		xp=x0[i][1];
  		yp=x0[i][3];
  		zt=x0[i][5];
  		xt = xp*cb - yp*sb;
  		yt = xp*sb + yp*cb;
  		x = xt*ca - zt*sa;
  		y = yt;
  		ninmap += addToMap(map,x,y,n,rng,div,min,max);
      }
  	
      for(int i=0; i<n; i++)
      {
      	for(int j=0; j<n; j++)
      	{
      		map[i][j] /= ninmap;
      	}
      }
      */
      return map;
  }
  
  public int addToMap(double map[][], double x, double y, int n, double rng, double div, double min, double max)
  {
      int xind = (int)((x-min)/div);	
      int yind = (int)((y-min)/div);
      
      if(xind < 0)
      {
      	//xind = 0;
      	return 0;
      }
      
      if(xind >= n)
      {
      	//xind = n-1;
      	return 0;
      }
      
      if(yind < 0)
      {
      	//yind = 0;
      	return 0;
      }
      
      if(yind >= n)
      {
      	//yind = n-1;
      	return 0;
      }
      
      map[xind][yind]++;
      
      return 1;
  }
  
  public integrator getIntegrator()
  {
  	return new integrator();
  }
  
	// RK4 integrator
	public static class integrator
	{
	 
      double xe[][] = null;
		double x[][] = null;
		double f[][] = null;
		
		double[] r22, r21, r2n;
		double[] r1, r2, rn;
		double[] a1, a2, a3;

		double m1, m2, m3;
		double eps1, eps2;
	double xn[] = new double[7];
              double xp = 0.0d;
              double yp = 0.0d;
              double zp = 0.0d;
              
      protected integrator()
      {
      	
      }
      
		//-----------------------------------------------------------------
		//
		//public void init_rkvar(double x0[][], double mass1, double mass2, double epsilon1, double epsilon2)
		//subroutine init_rkvar(x0, mass1, mass2, epsilon1, epsilon2, t1, p1, t2, p2, rs1, rs2, ro1, ro2, nn, nn1, nn2)
		public void init_rkvar(double x0[][], double mass1, double mass2, double epsilon1, double epsilon2, double t1, double p1, double t2, double p2, double rs1[], double rs2[], double ro1, double ro2, int nn, int n1, int n2)
		{
			int n = x0.length;
			x = new double[n][7];
			f = new double[n][7];
			xe = new double[n][7];

			r22 = new double[n];
			r21 = new double[n];
			r2n = new double[n];
			r1 = new double[n];
			r2 = new double[n];
			rn = new double[n];
			a1 = new double[n];
			a2 = new double[n];
			a3 = new double[n];


		  m1 = mass1;
		  m2 = mass2;
		  m3 = mass2;
		  eps1 = epsilon1*epsilon1;
		  eps2 = epsilon2*epsilon2;
		}// end subroutine init_rkvar

		//-----------------------------------------------------------------
		//
		//subroutine rk4(x0, h, xout)
		public void rk4(double x0[][], double h, double xout[][])
		{
		
			int n = x0.length;
			int i=0;
			double hod6 = h*od6;
			double h0p5 = h*0.5d;
			double hod3 = h*od3;
			
			for(i=0;i<n;i++)
			{
                          /*
				for(j=0;j<d;j++)
				{
					x[i][j] = x0[i][j];		
				}
                           */
                          x[i][0] = x0[i][0];
                          x[i][1] = x0[i][1];
                          x[i][2] = x0[i][2];
                          x[i][3] = x0[i][3];
                          x[i][4] = x0[i][4];
                          x[i][5] = x0[i][5];
                          x[i][6] = x0[i][6];
                          
			}
		  
		  diffeq(x);
		  for(i=0;i<n;i++)
			{
                    /*
				for(j=0;j<d;j++)
				{
			xe[i][j] = x0[i][j] + h * f[i][j] * od6;
		  x[i][j]  = x0[i][j] + h * f[i][j] * 0.5d;
				}
                     **/
			  /*
	  			xe[i][0] = x0[i][0] + h * f[i][0] * od6;
	            x[i][0]  = x0[i][0] + h * f[i][0] * 0.5d;
	            xe[i][1] = x0[i][1] + h * f[i][1] * od6;
	            x[i][1]  = x0[i][1] + h * f[i][1] * 0.5d;
	            xe[i][2] = x0[i][2] + h * f[i][2] * od6;
	            x[i][2]  = x0[i][2] + h * f[i][2] * 0.5d;
	            xe[i][3] = x0[i][3] + h * f[i][3] * od6;
	            x[i][3]  = x0[i][3] + h * f[i][3] * 0.5d;
	            xe[i][4] = x0[i][4] + h * f[i][4] * od6;
	            x[i][4]  = x0[i][4] + h * f[i][4] * 0.5d;
	            xe[i][5] = x0[i][5] + h * f[i][5] * od6;
	            x[i][5]  = x0[i][5] + h * f[i][5] * 0.5d;
	            xe[i][6] = x0[i][6] + h * f[i][6] * od6;
	            x[i][6]  = x0[i][6] + h * f[i][6] * 0.5d;
             */       
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
	            xe[i][6] = x0[i][6] + f[i][6] * hod6;
	            x[i][6]  = x0[i][6] + f[i][6] * h0p5;
			}

		  diffeq(x);
		  for(i=0;i<n;i++)
			{
                    /*
				for(j=0;j<d;j++)
				{
			xe[i][j] = xe[i][j] + h * f[i][j] * od3;
		  x[i][j]  = x0[i][j] + h * f[i][j] * 0.5d;
				}
                     */
			  /*
			  	xe[i][0] = xe[i][0] + h * f[i][0] * od3;
		        x[i][0]  = x0[i][0] + h * f[i][0] * 0.5d;
		        xe[i][1] = xe[i][1] + h * f[i][1] * od3;
		        x[i][1]  = x0[i][1] + h * f[i][1] * 0.5d;
		        xe[i][2] = xe[i][2] + h * f[i][2] * od3;
		        x[i][2]  = x0[i][2] + h * f[i][2] * 0.5d;
		        xe[i][3] = xe[i][3] + h * f[i][3] * od3;
		        x[i][3]  = x0[i][3] + h * f[i][3] * 0.5d;
		        xe[i][4] = xe[i][4] + h * f[i][4] * od3;
		        x[i][4]  = x0[i][4] + h * f[i][4] * 0.5d;
		        xe[i][5] = xe[i][5] + h * f[i][5] * od3;
		        x[i][5]  = x0[i][5] + h * f[i][5] * 0.5d;
		        xe[i][6] = xe[i][6] + h * f[i][6] * od3;
		        x[i][6]  = x0[i][6] + h * f[i][6] * 0.5d;
                  */
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
		        xe[i][6] = xe[i][6] + f[i][6] * hod3;
		        x[i][6]  = x0[i][6] + f[i][6] * h0p5;
			  
			}

		  diffeq(x);
		  for(i=0;i<n;i++)
			{
                    /*
				for(j=0;j<d;j++)
				{
			xe[i][j] = xe[i][j] + h * f[i][j] * od3;
		  x[i][j]  = x0[i][j] + h * f[i][j];
				}
                     */
			  /*
			  	xe[i][0] = xe[i][0] + h * f[i][0] * od3;
		        x[i][0]  = x0[i][0] + h * f[i][0];
		        xe[i][1] = xe[i][1] + h * f[i][1] * od3;
		        x[i][1]  = x0[i][1] + h * f[i][1];
		        xe[i][2] = xe[i][2] + h * f[i][2] * od3;
		        x[i][2]  = x0[i][2] + h * f[i][2];
		        xe[i][3] = xe[i][3] + h * f[i][3] * od3;
		        x[i][3]  = x0[i][3] + h * f[i][3];
		        xe[i][4] = xe[i][4] + h * f[i][4] * od3;
		        x[i][4]  = x0[i][4] + h * f[i][4];
		        xe[i][5] = xe[i][5] + h * f[i][5] * od3;
		        x[i][5]  = x0[i][5] + h * f[i][5];
		        xe[i][6] = xe[i][6] + h * f[i][6] * od3;
		        x[i][6]  = x0[i][6] + h * f[i][6];
		        */
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
		        xe[i][6] = xe[i][6] + f[i][6] * hod3;
		        x[i][6]  = x0[i][6] + h * f[i][6];
			}

		  diffeq(x);
		  for(i=0;i<n;i++)
			{
                    /*
				for(j=0;j<d;j++)
				{
			xe[i][j] = xe[i][j] + h * f[i][j] * od6;
				}
                     **/
			  /*
              xe[i][0] = xe[i][0] + h * f[i][0] * od6;
              xe[i][1] = xe[i][1] + h * f[i][1] * od6;
              xe[i][2] = xe[i][2] + h * f[i][2] * od6;
              xe[i][3] = xe[i][3] + h * f[i][3] * od6;
              xe[i][4] = xe[i][4] + h * f[i][4] * od6;
              xe[i][5] = xe[i][5] + h * f[i][5] * od6;
              xe[i][6] = xe[i][6] + h * f[i][6] * od6;
			*/
            xe[i][0] = xe[i][0] + f[i][0] * hod6;
            xe[i][1] = xe[i][1] + f[i][1] * hod6;
            xe[i][2] = xe[i][2] + f[i][2] * hod6;
            xe[i][3] = xe[i][3] + f[i][3] * hod6;
            xe[i][4] = xe[i][4] + f[i][4] * hod6;
            xe[i][5] = xe[i][5] + f[i][5] * hod6;
            xe[i][6] = xe[i][6] + f[i][6] * hod6;
            
            xout[i][0] = xe[i][0];
            xout[i][1] = xe[i][1];
            xout[i][2] = xe[i][2];
            xout[i][3] = xe[i][3];
            xout[i][4] = xe[i][4];
            xout[i][5] = xe[i][5];
            xout[i][6] = xe[i][6];
			}

		  /*
			for(i=0;i<n;i++)
			{
                   
				//for(j=0;j<d;j++)
				//{
			 //xout[i][j] = xe[i][j];
				//}
                     
                          xout[i][0] = xe[i][0];
                          xout[i][1] = xe[i][1];
                          xout[i][2] = xe[i][2];
                          xout[i][3] = xe[i][3];
                          xout[i][4] = xe[i][4];
                          xout[i][5] = xe[i][5];
                          xout[i][6] = xe[i][6];
                          
			}
		*/
		return;
		}
		//end subroutine rk4
		//
		//
		//----------------------------------------------------------------
		//subroutine diffeq(x)
		public void diffeq(double x[][])
		{
			
			int n = x.length;
			int i;
			
			for(i=0;i<7;i++)
			{
				xn[i] = x[n-1][i];
			}
	
                      // make temps to handle perturber galaxy
                      xp = xn[1];
                      yp = xn[3];
                      zp = xn[5];
                      double tmp4 = xp*xp+yp*yp+zp*zp;
                      double tmp5 = Math.sqrt(tmp4);
                      tmp4 = -m3/(tmp4+eps2);
                      
			double tmp1, tmp2, tmp3;
			for(i=0;i<n;i++)
			{
                          /*
				tmp1 = (x[i][1]-xn[1]);
				tmp2 = (x[i][3]-xn[3]);
				tmp3 = (x[i][5]-xn[5]);
                           **/
                              tmp1 = (x[i][1]-xp);
                              tmp2 = (x[i][3]-yp);
                              tmp3 = (x[i][5]-zp);
                              
				r22[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;
				
				tmp1 = x[i][1];
				tmp2 = x[i][3];
				tmp3 = x[i][5];
				r21[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;
				/*
				tmp1 = xn[1];
				tmp2 = xn[3];
				tmp3 = xn[5];
				r2n[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;
				*/
				  r2[i] = Math.sqrt(r22[i]);
				  r1[i] =  Math.sqrt(r21[i]);
				  //rn[i] =  Math.sqrt(r2n[i]);
                                rn[i] = tmp5;
                                
				  a1[i] = -m1 / (r21[i] + eps1);
				  a2[i] = -m2 / (r22[i] + eps2);
				  //a3[i] = -m3 / (r2n[i] + eps2);
                                a3[i] = tmp4;
			}


		  // this is a correction to prevent NaN errors in the vectorized
		  // function evalution at the location of the second mass
		  r2[n-1] = 1.0d;

		  // calculate the RHS of the diffeq
		  for(i=0; i<n; i++)
		  {
                    tmp1 = a1[i]/r1[i];
                    tmp2 = a2[i]/r2[i];
                    tmp3 = a3[i]/rn[i];
                    
		  f[i][0] = 1.0d;
		  f[i][1] = x[i][2];
                f[i][2] = tmp1*x[i][1] + tmp2*(x[i][1]-xn[1]) + tmp3*xn[1];
		  //f[i][2] = a1[i]*x[i][1]/r1[i] + a2[i]*(x[i][1]-xn[1])/r2[i] + a3[i]*xn[1]/rn[i];
		  f[i][3] = x[i][4];
                f[i][4] = tmp1*x[i][3] + tmp2*(x[i][3]-xn[3]) + tmp3*xn[3];
		  //f[i][4] = a1[i]*x[i][3]/r1[i] + a2[i]*(x[i][3]-xn[3])/r2[i] + a3[i]*xn[3]/rn[i];
		  f[i][5] = x[i][6];
                f[i][6] = tmp1*x[i][5] + tmp2*(x[i][5]-xn[5]) + tmp3*xn[5];
		  //f[i][6] = a1[i]*x[i][5]/r1[i] + a2[i]*(x[i][5]-xn[5])/r2[i] + a3[i]*xn[5]/rn[i];
		  }
		  return;
		}//end subroutine diffeq

		//----------------------------------------------------------------



		//end module integrator

	}
	
	public static void main(String args[])
	{
		JSpam jspam = new JSpam();

		double mass1, epsilon1, rin1, rout1, theta1, phi1;
		double mass2, epsilon2, rin2, rout2, theta2, phi2;
      double heat, seed, inclination_degree, omega_degree, rmin, velocity_factor, time, h;
		int opt1, opt2;
		int nout, nstep, n1, n2, n, istep;
		  int i,j,d;
		  int iout = 0;
		
		//double rscale1[] = new double[3];
		//double rscale2[] = new double[3];
		double rscale1[] = {3.0,3.0,3.0};
		double rscale2[] = {3.0,3.0,3.0};
		
		  // disk profile #1
		  mass1 = 1.0d;
		  epsilon1  = 0.3d;
		  rin1 = 0.05d;
		  rout1 = 1.0d;
		  //rscale1 = 3.0d;
		  theta1 = 0.0d;
		  phi1 = 0.0d;
		  opt1 = 1;

		  // disk profile #2
		  mass2 = 1.0d;
		  epsilon2  = 0.3d;
		  rin2 = 0.05d;
		  rout2 = 1.0d;
		  //rscale2 = 3.0d;
		  theta2 = 0.0d;
		  phi2 = 0.0d;
		  opt2 = 1;
		  
		  // generic disk parameters
		  heat = 0.0d;
		  seed = 0.0d;

		  // collision parameters
		  inclination_degree = 10.0d;
		  omega_degree = 0.0d;
		  rmin = 1.0d;
		  velocity_factor = 1.0d;
		  time = -3.0;
		  
		  // time step parameters
		  h = 0.02d;
		  nout = 5;
		  nstep = 500;

		  // particle numbers for the total and each disk
		  n1 = 1000;
		  n2 = 1000;
		  n = n1 + n2;

		  // allocate the space for the particles
		  double x0[][] = new double[n+1][7];
		  double xout[][] = new double[n+1][7];
		  
		  
		  // create the disks
		  jspam.profile(rin1, rout1, rscale1, 1, n1, mass1, epsilon1,
		       theta1, phi1, opt1, heat, seed, time, x0);
		 
		  jspam.profile(rin2, rout2, rscale2, n1+1, n, mass2, epsilon2,
		       theta2, phi2, opt2, heat, seed, time, x0);

		  // set the perturber galaxy position
		  jspam.perturberPosition(inclination_degree, omega_degree, rmin,
		       velocity_factor, mass1, mass2, epsilon1, epsilon2, h, n, n1, time, x0, rout1, rout2);
		  
		  
		  integrator intg = new integrator();

		  intg.init_rkvar(x0,mass1,mass2,epsilon1,epsilon2,theta1,phi1,theta2,phi2,rscale1,rscale2,rout1,rout2,n,n1,n2);
		  
		  d = 7;
		  for(istep = 1;istep<=nstep;istep++)
		  {

		    // take a step
		    intg.rk4(x0, h, xout);
		    for(i=0;i<n;i++)
		    {
		    	for(j=0;j<d;j++)
		    	{
		    		x0[i][j] = xout[i][j];		
		    	}
		    }
		    

		    // io routine
		    if( istep % nout == 0)
		    {
		      iout++; 
		    }
		  }
		  for(i=0; i<n; i++)
		  {
			  System.out.println(x0[i][1]+"\t"+x0[i][3]+"\t"+x0[i][5]+"\t"+x0[i][2]+"\t"+x0[i][4]+"\t"+x0[i][6]);
		  }
	}
}
