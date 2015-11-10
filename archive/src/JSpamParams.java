/*
Copyright 2008 by Anthony Holincheck
Licensed under the Academic Free License version 3.0
See the file "LICENSE" for more information
*/
package edu.gmu.cds.sim.trim;


public class JSpamParams
{
	public double mass1, epsilon1, rin1, rout1, theta1, phi1;
	public double mass2, epsilon2, rin2, rout2, theta2, phi2;
  public double heat, seed, inclination_degree, omega_degree, rmin, velocity_factor, time, h;
  public int opt1, opt2;
	public int nout, nstep, n1, n2, n, istep;

	public double alpha, beta, gamma;
	public double ca, sa, cb, sb, cg, sg;
	
	public double scalex;
	public double scaley;
	
	public double r1min = 0;
	public double r2min = 0;
	public double r1max = 0;
	public double r2max = 0;
	public double a1 = 0.0d;
	public double b1 = 0.0d;
	public double a2 = 0.0d;
	public double b2 = 0.0d;
	public double x1 = 0.0d;
	public double y1 = 0.0d;
	public double x2 = 0.0d;
	public double y2 = 0.0d;
	public double ang1 = 0.0d;
	public double ang2 = 0.0d;
	public double rx = 0.0d;
	public double ry = 0.0d;
	public double rz = 0.0d;
	public double vx = 0.0d;
	public double vy = 0.0d;
	public double vz = 0.0d;
	public double initVz = 0.0d;
	public double sf = 0.0d;
	public double minx = 0.0d;
	public double maxx = 0.0d;
	public double miny = 0.0d;
	public double maxy = 0.0d;
	
	public double simValues[][] = null;
	
	public JSpamParams()
	{
		setDefaults();
	}
	
	public Object clone()
	{
		JSpamParams jsp = new JSpamParams();
		
	    jsp.mass1 = mass1;
	    jsp.epsilon1 = epsilon1;
	    jsp.rin1 = rin1;
	    jsp.rout1 = rout1;
	    jsp.theta1 = theta1;
	    jsp.phi1 = phi1;
	    
	    jsp.mass2 = mass2;
	    jsp.epsilon2 = epsilon2;
	    jsp.rin2 = rin2;
	    jsp.rout2 = rout2;
	    jsp.theta2 = theta2;
	    jsp.phi2 = phi2;
	    
		jsp.heat = heat;
		jsp.seed = seed;
		jsp.inclination_degree = inclination_degree;
		jsp.omega_degree = omega_degree;
		jsp.rmin = rmin;
		jsp.velocity_factor = velocity_factor;
		jsp.time = time;
		jsp.h = h;
		
		jsp.opt1 = opt1;
		jsp.opt2 = opt2;
		
	    jsp.nout = nout;
	    jsp.nstep = nstep;
	    jsp.n1 = n1;
	    jsp.n2 = n2;
	    jsp.n = n;
	    jsp.istep = istep;

		jsp.setAlpha(alpha);
		jsp.setBeta(beta);
		jsp.setGamma(gamma);
		
		jsp.scalex = scalex;
		jsp.scaley = scaley;
		
		jsp.r1min = r1min;
		jsp.r2min = r2min;
		jsp.r1max = r1max;
		jsp.r2max = r2max;

		jsp.a1 = a1;
		jsp.b1 = b1;
		jsp.a2 = a2;
		jsp.b2 = b2;
		jsp.x1 = x1;
		jsp.y1 = y1;
		jsp.x2 = x2;
		jsp.y2 = y2;
		jsp.ang1 = ang1;
		jsp.ang2 = ang2;
		jsp.rx = rx;
		jsp.ry = ry;
		jsp.rz = rz;
		jsp.vx = vx;
		jsp.vy = vy;
		jsp.vz = vz;
		jsp.initVz = initVz;
		
		jsp.sf = sf;
		jsp.minx = minx;
		jsp.maxx = maxx;
		jsp.miny = miny;
		jsp.maxy = maxy;
		
		return jsp;
	}
	
	public void setAlpha(double alpha)
	{
		this.alpha = alpha;
		ca = Math.cos(Math.toRadians(alpha));
		sa = Math.sin(Math.toRadians(alpha));
	}
	
	public void setBeta(double beta)
	{
		this.beta = beta;
		cb = Math.cos(Math.toRadians(beta));
		sb = Math.sin(Math.toRadians(beta));		
	}
	
	public void setGamma(double gamma)
	{
		this.gamma = gamma;
		cg = Math.cos(Math.toRadians(gamma));
		sg = Math.sin(Math.toRadians(gamma));		
	}
	
	public void setDefaults()
	{
		minx=-2;
		maxx=2;
		miny=-2;
		maxy=2;
		scalex = 4.0;
		scaley = 4.0;
		
		// disk profile #1
		mass1 = 1.0d;
		epsilon1  = 0.1d;
		rin1 = 0.05d;
		rout1 = 1.0d;
		//rscale1 = 3.0d;
		theta1 = 0.0d;
		phi1 = 0.0d;
		opt1 = 1;
		
		// disk profile #2
		mass2 = 1.0d;
		epsilon2  = 0.1d;
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
		n2 = 500;
		n = n1 + n2;
		
		setAlpha(0);
		setBeta(0);
		setGamma(0);
	}
	
	public void adjustTotalN(int n)
	{
		this.n1 = (int)(mass1/(mass1+mass2)*((double)n));
		this.n2 = (int)(mass2/(mass1+mass2)*((double)n));
		this.n = this.n1 + this.n2;
	}
}
