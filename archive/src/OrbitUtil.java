/*
Copyright 2008 by Anthony Holincheck
Licensed under the Academic Free License version 3.0
See the file "LICENSE" for more information
*/
package edu.gmu.cds.sim.trim;


import java.util.ArrayList;
import java.util.List;


/**
* Utilities used for converting between r, v and
* orbital elements.
* 
* @author Anthony Holincheck
*
*/
public class OrbitUtil
{
  public static final double[] I = {1.0,0.0,0.0};
  public static final double[] J = {0.0,1.0,0.0};
  public static final double[] K = {0.0,0.0,1.0};
  
	public static final double TWOPI = 2.0d*Math.PI;

  public static final double G = 6.673e-11;
  public static final double solarMassToKg = 1.98892e30;
  public static final double metersPerMpc = 3.08568025e22;
  public static final double kmPerMpc = 3.08568025e19; // or m per kpc
  public static final double degToRad = Math.PI/180.0d;
  
  public static final double MU_kg = 1.0e11*solarMassToKg;
  public static final double DU_m = 15.0*kmPerMpc; // meters in a kpc
  public static final double TU_s = Math.sqrt(DU_m*DU_m*DU_m/(G*MU_kg));

  public static final double MU_sm = 1.0e11;
  public static final double DU_mpc = 15.0d/1000.0d; 
  
  public static final double VEL_KMS = DU_m/1000.0d/TU_s;
  public static final double A_MSS = DU_m/TU_s/TU_s;

  public static final double A0_MKS = 1.2e-10;
  public static final double A0 = A0_MKS/A_MSS;
  
  public static final double DEFAULT_EPS = Math.sqrt(0.1d);
  
  /**
   * Computes a family of orbits that fit the specified rx, ry, and vz.
   * It uses the mass value to calculate an escape velocity at the rx+ry distance.
   * It uses the specified diameter to determine rz, which is + or - 5 diameter.
   * 
   * @param rx distance in Mpc
   * @param ry distance in Mpc
   * @param vz velocity in km/s
   * @param m1 mass in solar masses
   * @param m2 mass in solar masses
   * @param d1 diameter in Mpc
   * @return rvs
   */
public static double[][] computeOrbits(double rx, double ry, double vz, double m1, double m2, double d1)
  {
  	int i=0;
  	int j=0;
  	int k=0;
  	
  	int ind = 0;
  	
  	int nr = 51;
  	int nvm = 51;
  	int nva = 36;
  	
  	//double rvs[][] = new double[nr*nvm*nva][6];
  	
  	List validRVs = new ArrayList((int)(nr*nvm*nva/2));
  	double rvs[][] = null;
  	double tmpRv[] = null;
  	
  	double angInc = (int)(360.0d/(double)(nva-1));
  	angInc *= degToRad;
  	
  	double vmag = 0.0;
  	double vmag1 = 0.0;
  	double vesc = 0.0;
  	double vinc = 0.0;
  	double vang = 0.0;
  	double vx,vy;
  	double coe[] = null;
  	
  	// calculate gravitational parameter
  	double mu = G*(m1+m2)*solarMassToKg;

  	// convert to units of meters
  	double r[] = {rx*metersPerMpc,ry*metersPerMpc,0};
  	
  	// convert to units of meters/s
  	double v[] = {0,0,vz*1000.0d};
  	
  	// calculate escape velocity
  	double rmag = Math.sqrt(r[0]*r[0] + r[1]*r[1]);
  	vesc = Math.sqrt(2.0d*mu/rmag);
  	
  	vmag = Math.sqrt(4*vesc*vesc-v[2]*v[2]);
  	vmag1 = vmag;
  	vinc = vmag/((double)nvm-1);
  	
  	// start rz at nd times diameter
  	double nd = 10.0d;
  	double rinc = 2.0d*nd/((double)nr-1)*d1;
  	
      double rz = -nd*d1;

      double minE = 0.0;
      double maxE = 4.0;
      
      try
      {
	        // loop over possible values of rz
	        for(i=0; i<nr; i++)
	        {
	            r[2] = rz*metersPerMpc;
	            vmag = vmag1;
	            // loop over velocity magnitude increments
	            for(j=0;j<nvm;j++)
	            {
	            	// loop over velocity angles
	            	for(k=0;k<=nva;k++)
	            	{
	            		//vang = k*5.0d*degToRad;
	            		vang = k*angInc;
	            		vx = vmag*Math.cos(vang);
	            		vy = vmag*Math.sin(vang);
	            		v[0] = vx;
	            		v[1] = vy;
	            		coe = rvToCoe(r,v,mu);
	            		// check eccentricity
	            		if(coe[1] <= maxE && coe[1]>=minE)
	            		{
	            			tmpRv = new double[6];
	            			tmpRv[0] = r[0];
	            			tmpRv[1] = r[1];
	            			tmpRv[2] = r[2];
	            			tmpRv[3] = v[0];
	            			tmpRv[4] = v[1];
	            			tmpRv[5] = v[2];
	            			validRVs.add(tmpRv);
	            		}
	            		ind++;
	            	}
	            	vmag -= vinc;
	            }
	            
	        	// increment rz
	        	//rz += 0.04*d1;
	            rz += rinc;
	        }
	        int size = validRVs.size();
	        rvs = new double[size][6];
  		
	        for(i = 0; i<size; i++)
	        {
	        	tmpRv = (double[])validRVs.get(i);
		        rvs[i][0] = tmpRv[0];
	    		rvs[i][1] = tmpRv[1];
	    		rvs[i][2] = tmpRv[2];
	    		rvs[i][3] = tmpRv[3];
	    		rvs[i][4] = tmpRv[4];
	    		rvs[i][5] = tmpRv[5];
	        }
      }
      catch(Exception ex)
      {
      	ex.printStackTrace();
      }
      finally
      {
      	//if(fw != null)try{fw.close();}catch(Exception ex){};
      }
      
      System.out.println(ind + "\t" + rvs.length);
      return rvs;
  }
  /*
  public static double[][] computeOrbits(double rx, double ry, double vz, double m1, double m2, double d1)
  {
  	int i=0;
  	int j=0;
  	int k=0;
  	
  	int ind = 0;
  	
  	int nr = 51;
  	int nvm = 20;
  	int nva = 71;
  	
  	double rvs[][] = new double[nr*nvm*nva][6];
  	
  	double angInc = (int)(360.0d/(double)nva);
  	angInc *= degToRad;
  	
  	double vmag = 0.0;
  	double vesc = 0.0;
  	double vinc = 0.0;
  	double vang = 0.0;
  	double vx,vy;
  	double coe[] = null;
  	
  	// calculate gravitational parameter
  	double mu = G*(m1+m2)*solarMassToKg;

  	// convert to units of meters
  	double r[] = {rx*metersPerMpc,ry*metersPerMpc,0};
  	
  	// convert to units of meters/s
  	double v[] = {0,0,vz*1000.0d};
  	
  	// calculate escape velocity
  	double rmag = Math.sqrt(r[0]*r[0] + r[1]*r[1]);
  	vesc = Math.sqrt(2.0d*mu/rmag);
  	
  	vmag = Math.sqrt(1*vesc*vesc-v[2]*v[2]);
  	vinc = vmag/((double)nvm);
  	
  	// start rz at -3 times diameter
  	double rinc = 6.0d/((double)nr)*d1;
  	
      double rz = -3.0d*d1;
      //FileWriter fw = null;
      
      try
      {
      	/*
      	fw = new FileWriter("orbitfamily.txt");
      	fw.write("variables: 5\r\n");
      	fw.write("labels:\r\n");
      	fw.write("P\r\n");
      	fw.write("e\r\n");
      	fw.write("inc\r\n");
      	fw.write("omega\r\n");
      	fw.write("w\r\n");
      
      // loop over possible values of rz
      for(i=0; i<nr; i++)
      {
          r[2] = rz*metersPerMpc;
          
          // loop over velocity magnitude increments
          for(j=0;j<nvm;j++)
          {
          	// loop over velocity angles
          	for(k=0;k<nva;k++)
          	{
          		//vang = k*5.0d*degToRad;
          		vang = k*angInc;
          		vx = vmag*Math.cos(vang);
          		vy = vmag*Math.sin(vang);
          		v[0] = vx;
          		v[1] = vy;
          		//coe = rvToCoe(r,v,mu);
          		/*
          		for(int q=0;q<5;q++)
          		{
          			if(q >0)
          			{
          				fw.write("\t");
          			}
          			fw.write(String.valueOf(coe[q]));
          		}
          		fw.write("\r\n");
          		dumpCoe(coe);
          		
          		rvs[ind][0] = r[0];
          		rvs[ind][1] = r[1];
          		rvs[ind][2] = r[2];
          		rvs[ind][3] = v[0];
          		rvs[ind][4] = v[1];
          		rvs[ind][5] = v[2];
          		
          		ind++;
          	}
          	vmag -= vinc;
          }
          
      	// increment rz
      	//rz += 0.04*d1;
          rz += rinc;
      }
      }
      catch(Exception ex)
      {
      	ex.printStackTrace();
      }
      finally
      {
      	//if(fw != null)try{fw.close();}catch(Exception ex){};
      }
      
      return rvs;
  }
  */
  /**
   * Converts r and v in m and m/s to canonical units.  The
   * diameter of the primary galaxy, in Mpc is given as well.
   * The simulation assumes radius of galaxy as 1 and mass of 
   * galaxy as 1.  So we know all we need to make canonical units.
   * The output is a vector ready to use in JSpam to find the 
   * rmin.  The velocity has opposite direction of what is given here.
   * 
   * @param r
   * @param v
   * @param d1
   * @return
   */
  public static double[] toCanonical(double r[], double v[], double d1, double mu)
  {
  	double rv[] = new double[9];
  	
  	// calculating r is easy
  	rv[1] = r[0]/DU_m;
  	rv[3] = r[1]/DU_m;
  	rv[5] = r[2]/DU_m;
  	
  	
  	rv[2] = -v[0]*TU_s/DU_m;
  	rv[4] = -v[1]*TU_s/DU_m;
  	rv[6] = -v[2]*TU_s/DU_m;
  	
  	rv[7] = DU_m;
  	rv[8] = TU_s;

  	return rv;
  }
  /*
  public static double[] toCanonical(double r[], double v[], double d1, double mu)
  {
  	double rv[] = new double[9];
  	
  	double DU = d1*0.5*metersPerMpc;
  	
  	// calculating r is easy
  	rv[1] = r[0]/DU;
  	rv[3] = r[1]/DU;
  	rv[5] = r[2]/DU;
  	
  	// to get v we need to get a TU
  	double TU = Math.sqrt(DU*DU*DU/mu);
  	
  	rv[2] = -v[0]*TU/DU;
  	rv[4] = -v[1]*TU/DU;
  	rv[6] = -v[2]*TU/DU;
  	
  	rv[7] = DU;
  	rv[8] = TU;
  	return rv;
  }
  */
  
  public static double toDegreesInRange(double rad, double min, double max)
  {
  	double deg = 0;
  
  	double range = max-min;
  	
  	deg = Math.toDegrees(rad);
  	
  	while(deg > max)
  	{
  		deg-=range;
  	}
  	
  	while(deg < min)
  	{
  		deg += range;
  	}
  	
  	return deg;
  }
  
  /**
   * Convert coe to r and v vector.
   * @param coe
   * @return
   */
  public static double[] coeToRv(double coe[], double mu)
  {
  	double rv[] = new double[6];
  	
  	double rcoe[] = new double[3];
  	double r[] = new double[3];
  	double vcoe[] = new double[3];
  	double v[] = new double[3];
  	
  	double p = coe[0];
  	double e = coe[1];
  	double nu = coe[5];
  	
  	double rmag = p/(1.0+e*Math.cos(nu));
  	rcoe[0]=rmag*Math.cos(nu);
  	rcoe[1]=rmag*Math.sin(nu);
  	
  	double vmag = Math.sqrt(mu/p);
  	vcoe[0]=-vmag*Math.sin(nu);
  	vcoe[1]=vmag*(e+Math.cos(nu));

  	double cosi = Math.cos(coe[2]);
  	double sini = Math.sin(coe[2]);
  	double cosO = Math.cos(coe[3]);
  	double sinO = Math.sin(coe[3]);
  	double cosw = Math.cos(coe[4]);
  	double sinw = Math.sin(coe[4]);
  	
  	double R[][] = new double[3][3];
  	
  	R[0][0] = cosO*cosw-sinO*sinw*cosi;
  	R[0][1] = -cosO*sinw-sinO*cosw*cosi;
  	R[0][2] = sinO*sini;
  	R[1][0] = sinO*cosw+cosO*sinw*cosi;
  	R[1][1] = -sinO*sinw+cosO*cosw*cosi;
  	R[1][2] = -cosO*sini;
  	R[2][0] = sinw*sini;
  	R[2][1] = cosw*sini;
  	R[2][2] = cosi;
  	
  	for(int i=0; i<3; i++)
  	{
  		r[i] = 0;
  		v[i] = 0;
  		for(int j=0; j<3; j++)
  		{
  			r[i] += R[i][j]*rcoe[j];
  			v[i] += R[i][j]*vcoe[j];
  		}
  	}
  	
  	rv[0]=r[0];
  	rv[1]=r[1];
  	rv[2]=r[2];
  	rv[3]=v[0];
  	rv[4]=v[1];
  	rv[5]=v[2];
  	//System.out.println(rv[0]+","+rv[1]+","+rv[2]+","+rv[3]+","+rv[4]+","+rv[5]);
  	return rv;
  }
  
  /**
   * Convert the position and velocity to classical orbital
   * elements.
   * 
   * @param r
   * @param v
   * @param mu
   * @return
   */
  public static double[] rvToCoe(double rx, double ry, double rz, double vx, double vy, double vz, double mu)
  {
  	double r[] = new double[]{rx,ry,rz};
  	double v[] = new double[]{vx,vy,vz};
  	//System.out.println(rx+","+ry+","+rz+","+vx+","+vy+","+vz);
  	return rvToCoe(r,v,mu);
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
  	double muInv = 1.0d/mu;
  	//System.out.println("mu\t"+mu + "\t" + muInv);
  	//System.out.println("r\t"+r[0]+"\t"+r[1]+"\t"+r[2]);
  	//System.out.println("v\t"+v[0]+"\t"+v[1]+"\t"+v[2]);
      double coe[] = new double[7];
     
      double rmag = mag(r);
      double vmag = mag(v);
      //System.out.println(rmag + "\t" + vmag);
      // angular momentum
      double h[] = cross(r,v);
      double hmag = mag(h);
      
      // node vector
      double n[] = cross(K,h);
      double nmag = mag(n);
      
      // eccentricity vector
      double tmp1 = 0.0;
      double tmp2 = 0.0;
      
      tmp1 = vmag*vmag - mu/rmag;
      tmp2 = dot(r,v);
      
      double v1[] = null;
      double v2[] = null;
      
      v1 = scale(tmp1,r);
      v2 = scale(tmp2,v);
      double ev[] = sub(v1,v2);
      ev = scale(muInv,ev);
      
      // p
      double p = hmag*hmag*muInv;
      
      // ecc
      double ecc = mag(ev);
      
      // cosi
      double cosi = h[2]/hmag;
      
      // cosO
      double cosO = n[0]/nmag;
      
      // cosw
      double cosw = dot(n,ev)/(nmag*ecc);
      
      // cosv
      double cosv = dot(ev,r)/(ecc*rmag);
      
      // cosu
      double cosu = dot(n,r)/(nmag*rmag);
      
      //System.out.println(p + "\t" + ecc + "\t" + cosi + "\t" + cosO + "\t" + cosw + "\t" + cosv + "\t" + cosu);
      coe[0] = p;
      coe[1] = ecc;
      coe[2] = Math.acos(cosi);
      
      tmp1 = Math.acos(cosO);
      if(n[1] < 0)
      {
      	// should be less than pi
      	tmp1 = TWOPI-tmp1;
      }
      coe[3] = tmp1;
      
      tmp1 = Math.acos(cosw);
      if(ev[2] < 0)
      {
      	// should be less than pi
      	tmp1 = TWOPI-tmp1;
      }
      coe[4] = tmp1;
      
      tmp1 = Math.acos(cosv);
      if(dot(r,v) < 0)
      {
      	// should be less than pi
      	tmp1 = TWOPI-tmp1;
      }
      coe[5] = tmp1;
      
      if(cosu > 1.0 || cosu < -1.0)
      {
      	coe[6] = Double.NaN;
      }
      else
      {
      	tmp1 = Math.acos(cosu);
      	if(r[2]>0)
      	{
      		// should be less than pi
      	}
      	coe[6] = tmp1;
      }
      
      return coe;
  }
  
  /**
   * Calculates cross product of two vectors.
   * 
   * @param v1
   * @param v2
   * @return
   */
  public static double[] cross(double v1[], double v2[])
  {
  	double v[] = new double[3];
  	
  	v[0] = v1[1]*v2[2]-v1[2]*v2[1];
  	v[1] = v1[2]*v2[0]-v1[0]*v2[2];
  	v[2] = v1[0]*v2[1]-v1[1]*v2[0];
  	
  	return v;
  }
  
  /**
   * Calculates the dot product of two vectors.
   * 
   * @param v1
   * @param v2
   * @return
   */
  public static double dot(double v1[], double v2[])
  {
  	int size = v1.length;
  	double cp = 0.0;
  	
  	for(int i=0; i<size; i++)
  	{
  		cp += v1[i]*v2[i];
  	}
  	
  	return cp;
  }
  
  /**
   * Calculates the magnitude of the vector.
   * 
   * @param v
   * @return
   */
  public static double mag(double v[])
  {
  	return Math.sqrt(dot(v,v));
  }
  
  /**
   * Add two vectors.
   * 
   * @param v1
   * @param v2
   * @return
   */
  public static double[] add(double v1[], double v2[])
  {
  	int size = v1.length;
  	double v[] = new double[size];
  	
  	for(int i=0; i<size; i++)
  	{
  		v[i] = v1[i]+v2[i];
  	}
  	
  	return v;
  }
  
  /**
   * Subtract v2 from v1;
   * 
   * @param v1
   * @param v2
   * @return
   */
  public static double[] sub(double v1[], double v2[])
  {
  	int size = v1.length;
  	double v[] = new double[size];
  	
  	for(int i=0; i<size; i++)
  	{
  		v[i] = v1[i] - v2[i];
  	}
  	return v;
  }
  
  /**
   * Scale the vector by the specified multiplier.
   * 
   * @param sc
   * @param v1
   * @return
   */
  public static double[] scale(double sc, double v1[])
  {
  	int size = v1.length;
      double v[] = new double[size];	
      
      for(int i=0; i<size; i++)
      {
      	v[i] = sc*v1[i];
      }
      
      return v;
  }
  
  public static void dumpCoe(double coe[])
  {
      System.out.println("p = " + coe[0]);
      System.out.println("e = " + coe[1]);
      System.out.println("i = " + coe[2] + "\t" + toDegreesInRange(coe[2],0,180));
      System.out.println("O = " + coe[3] + "\t" + toDegreesInRange(coe[3],0,360));
      System.out.println("w = " + coe[4] + "\t" + toDegreesInRange(coe[4],0,360));
      System.out.println("v = " + coe[5] + "\t" + toDegreesInRange(coe[5],0,360));
      System.out.println("u = " + coe[6] + "\t" + toDegreesInRange(coe[6],0,360));
  }
  
  /**
   * Uses f and g series to propogate orbit.
   * 
   * @param ro
   * @param vo
   * @param dt
   * @param r
   * @param v
   */
  public static void fandgSeries(double ro[], double vo[], double dt, double r[], double v[])
  {
  	double mu = G*solarMassToKg;
  	fandgSeries(ro,vo,dt,mu,r,v);
  }
  
  /**
   * Uses f and g series to propogate orbit.
   * 
   * @param ro
   * @param vo
   * @param dt
   * @param mu
   * @param r
   * @param v
   */
  public static void fandgSeries(double ro[], double vo[], double dt, double mu, double r1[], double v1[])
  {
  	double romag = mag(ro);
  	double vomag = mag(vo);
  	
  	double rdotv = dot(ro,vo);
  	
  	double vo2 = vomag*vomag;
  	
  	//double en = vo2*0.5 - mu/romag;
  	double alpha = -vo2/mu + 2/romag;
  	
  	double psi = 0.0;
  	double xo = 0.0;
  	
  	double tmp1 = 0.0;
  	double tmp2 = 0.0;
  	
  	// elliptic
  	if(alpha > 1e-6)
  	{
  	    xo = Math.sqrt(mu)*dt*alpha;	
  	}
  	// parabolic
  	else if(alpha > -1e-6)
  	{
  		double h[] = cross(ro,vo);
  		double hmag = mag(h);
  		double p = hmag*hmag/mu;
  		double cot2s = 3.0*Math.sqrt(mu/(p*p*p))*dt;
  		
  		double s = acot(cot2s);
  		s = 0.5*s;
  		double tan3w = Math.tan(s);
  		double w  = Math.pow(tan3w,1.0d/3.0d);
  		w = Math.atan(w);
  		xo = Math.sqrt(p)*2*cot(2*w);
  	}
  	// hyperbolic
  	else
  	{
  		double a = 1.0d/alpha;
  		double signum = 1.0d;
  		if(dt < 0.0d)
  		{
  			signum = -1.0d;
  		}
  		
  		tmp1 = -2.0d*mu*alpha*dt;
  		tmp2 = signum*Math.sqrt(-mu*a)*(1.0d-romag*alpha);
  		tmp1 = tmp1/(rdotv + tmp2);
  		tmp1 = Math.log(tmp1);
  		xo = signum*Math.sqrt(-a)*tmp1;
  	}
  	
  	// Now we have xo, let's iterate to improve
  	double xn = 0.0;
  	double xn2 = 0.0;
  	double xnp1 = xo;
  	double sqmu = Math.sqrt(mu);
  	
  	double c[] = new double[2];
  	double c2 = 0;
  	double c3 = 0;
  	double r = 0;
  	
  	// iterate to get xn
  	do
  	{
  		xn = xnp1;
  		xn2 = xn*xn;
  		
  		psi = xn2*alpha;
  		getC2C3(psi,c);
  		c2 = c[0];
  		c3 = c[1];
  		
  		r = xn2*c2+rdotv/sqmu*xn*(1.0d-psi*c3)+romag*(1.0d-psi*c2);
  		tmp1 = sqmu*dt-xn2*xn*c3-rdotv/sqmu*xn2*c2-romag*xn*(1-psi*c3);
  		tmp1 = tmp1/r;
  		xnp1 = xn + tmp1;
  	}while(Math.abs(xn-xnp1) > 1e-6);
  	
  	xn = xnp1;
  	xn2 = xn*xn;
  	
  	double f = 1.0d - xn2/romag*c2;
  	double g = dt -xn2*xn/sqmu*c3;
  	double gdot = 1.0d - xn2/r*c2;
  	double fdot = sqmu/(r*romag)*xn*(psi*c3-1.0d);
  	
  	double rt1[] = add(scale(f,ro),scale(g,vo));
  	double vt1[] = add(scale(fdot,ro),scale(gdot,vo));
  	
  	r1[0] = rt1[0];
  	r1[1] = rt1[1];
  	r1[2] = rt1[2];
  	v1[0] = vt1[0];
  	v1[1] = vt1[1];
  	v1[2] = vt1[2];
  	
  	tmp1 = f*gdot-fdot*g;
  	System.out.println(tmp1);
  }
  
  /**
   * Calculates c2 and c3.  Vallado 2nd edition algorithm 1.
   * 
   * @param psi
   * @param c
   */
  public static void getC2C3(double psi, double c[])
  {
  	double tmp = 0.0;
  	
  	if(psi > 1e-6)
  	{
  		tmp = Math.sqrt(psi);
  		c[0] = (1-Math.cos(tmp))/psi;
  		c[1] = (tmp-Math.sin(tmp))/(tmp*tmp*tmp);
  	}
  	else if(psi < -1e-6)
  	{
  		tmp = Math.sqrt(-psi);
  		// Java 1.4 doesn't implement cosh and sinh in Math
  		//c[0] = (1-Math.cosh(tmp))/psi;
  	    //c[1] = (Math.sinh(tmp)-tmp)/(tmp*tmp*tmp);
  	}
  	else
  	{
  		c[0] = 0.5;
  		c[1] = 1.0d/6.0d;
  	}
  }
  
  public static double cot(double x)
  {
  	return 1.0d/Math.tan(x);
  }
  
  public static double acot(double x)
  {
  	return Math.atan(1.0d/x);
  }
}
