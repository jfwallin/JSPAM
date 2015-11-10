package edu.gmu.cds.sim.trim;

public class DFJSpam extends JSpam
{
	private static String syncObject = new String("");

	// we only need to initialize the distribution
	// once, other instances can just get a copy
	private static dfintegrator INIT_INT = null;
	
	protected dfintegrator dfint = null;
	
	public double maxDF = -1e20;
	public double timeMax = -1e20;
	public DFJSpam()
	{
		super();
		dfint = new dfintegrator();
	}
	
	// override force calculation
    public void diffq1(double x[], double mass1, double mass2, double eps1, double eps2, double f[], double rout1, double rout2)
    {
    	//super.diffq1(x,mass1,mass2,eps1,f);

    	double r21 = 0;
    	double r1 = 0;
    	double a1 = 0;
    	double a2 = 0;
    	double at = 0;
    	
    	double c1 = 0;
    	double c2 = 0;
    	double c3 = 0;
    	double v21 = 0;
    	double v1 = 0;
    	double xvalue = 0;
    	
    	int ival = 0;
    	int ival2 = 0;
    	
    	double df_force1 = 0;
    	double df_force2 = 0;
    	double df_sigma = 0;
    	double df_rho = 0;
    	
    	int i = 0;
    	double rr = 0;
    	double rlocal = 0;
    	double ivalue = 0;
    	double dr = 0;
    	double mmm = 0;
    	double dm = 0;
    	
    	//    r21 = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    	//    r1  = sqrt(r21)

    	r21=x[1]*x[1]+x[3]*x[3]+x[5]*x[5];
    	v21=x[2]*x[2]+x[4]*x[4]+x[6]*x[6];

    	r1=Math.sqrt(r21);
    	        
    	// get the index for the calculations
    	    ival  =  dfint.df_index(r1, rout1   );;
    	    ival2 =  dfint.df_index(r1, rout2   );;

    	// get the forces, sigma and rho, and rescale them

   	    //df_force1 = dfint.acceleration_particle[ival] * rs_internal * rs_internal;
   	    //df_force2 = dfint.acceleration_particle[ival2] * rs_internal * rs_internal;
   	    df_force1 = dfint.acceleration_particle[ival] * rs_internal2;
   	    df_force2 = dfint.acceleration_particle[ival2] * rs_internal2;

   	    //df_sigma  =  dfint.new_vr2[ival] * rs_internal * rs_internal;
   	    //df_rho    = dfint.new_rho[ival] * ( rs_internal * rs_internal * rs_internal );
   	    df_sigma  =  dfint.new_vr2[ival] * rs_internal2;
   	    df_rho    = dfint.new_rho[ival] * ( rs_internal3 );



    	// interpolated forces 
    	    a1 = -mass1 * df_force1;
    	    a2 = -mass2 * df_force2;
    	    at = a1 + a2;


    	// df
    	    //v21 = x(4)*x(4)+x(5)*x(5)+x(6)*x(6)
    	    v1  = Math.sqrt(v21);


    	    xvalue = v1 / df_sigma;
    	    c1 = erf(xvalue) - 2.0d * xvalue / sqrtpi * Math.exp(-xvalue*xvalue);
    	    //write(21,*) r1, c1, xvalue, v1, df_sigma, df_rho

    	// df formula with G=1
    	    c2 = -4.0d * pi * mass2 * dfint.lnl / v21;
    	    c3 = c1 * c2 * df_rho;    

    	    if(Math.abs(c3)>maxDF)
    	    {
    	    	maxDF = Math.abs(c3);
    	    	timeMax = x[0];
    	    }
    	    f[0] = 1.0d;
    	    f[1] = x[2];
    	    f[3] = x[4];
    	    f[5] = x[6];
  
    	    f[2] = at*x[1]/r1 - c3*x[2]/v1;
    	    f[4] = at*x[3]/r1 - c3*x[4]/v1;
    	    f[6] = at*x[5]/r1 - c3*x[6]/v1;
    	    //f[4] = at * x(1)/r1  - c3 * x(4)/ v1 
    	    //f[5] = at * x(2)/r1  - c3 * x(5)/ v1
    	    //f[6] = at * x(3)/r1  - c3 * x(6)/ v1
/*
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
*/
    	    return;
//    	  end subroutine diffq3

    }
    
    /*
    !------------------------------------------------------
    !
    !
    !
      subroutine profile(rin, rout, rscale, nstart, ntot, mass, eps, &
           theta, phi, opt, heat, t0, x0)

    !
    !
    !       variables-
    !       opt - option for the distribution
    !       rin - inner radius
    !       rout - outer radius
    !       rscale - scale of brightness drop
    !	nstart - start number for placement of particles
    !       ntot - number of particles to be placed
    !       heat - heat parameter
    !       m - mass of galaxy
    !       sl - softening length
    !       nring - number of rings
    !       npart - number of particle per ring (opt)
    !       x0 - position of center of mass
    !

        implicit none
        real (kind=8), intent(in) :: rin
        real (kind=8), intent(in) :: rout
        real (kind=8), dimension(3), intent(in) :: rscale
        integer (kind=4), intent(in) ::  nstart
        integer (kind=4), intent(in) ::  ntot 
        real (kind=8), intent(in) :: mass 
        real (kind=8), intent(in) :: eps
        real (kind=8), intent(in) :: theta
        real (kind=8), intent(in) :: phi
        integer(kind=4), intent(in) :: opt
        real (kind=8), intent(in) :: heat
    !    real (kind=8), intent(in) :: seed
        real (kind=8), intent(in) :: t0
        real (kind=8), intent(out), dimension(:,:) :: x0

        real (kind=8) :: stheta,ctheta,sphi,cphi,pi
        real (kind=8) :: x3,y3,z3,xv3,yv3,zv3,x2,y2,z2,xv2,yv2,zv2
        real (kind=8) :: x,y,z,xv,yv,zv

        integer (kind=4) :: i, j, n, nprof

        real (kind=8) :: rnorm
        real (kind=8), dimension(:), allocatable :: rp, r, angle, v, p_ring, cp_ring
        real (kind=8) :: st, ct, dr, ran, r1, r2, ptot
        integer (kind=4), dimension(:), allocatable :: n_ring
        integer (kind=4) :: nring, dnring, is, ie, iring, tct

        real (kind=8) :: xp, yp, zp
        real (kind=8) :: fx, fy, fz, ftotal

        integer (kind=4) :: ival
*/
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

    	int i,j,n,nprof,iring,ival;
    	int nring,dnring, is, ie, tct;
    	
    	double rnorm, ptot, r2;
    	double[] rp, r, angle, v, n_ring, p_ring, cp_ring;
    	double st, ct, dr, ran, r1;

        double xp, yp, zp;
        double fx, fy, fz, ftotal;
        
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
        //nprof = 1000;
        nring = nprof / 10;

        dnring = nprof/nring;
        rp = new double[nprof];
        n_ring = new double[nprof];
        p_ring = new double[nprof];
        cp_ring = new double[nprof];

    // set the differential sum of the probability function into a vector
        rnorm = 0.0d;
        dr = ((double)(rout - rin))/((double)nprof);
        for(i=0; i<nprof; i++)
        {
          r1 = (i+1)*dr + rin;
          rp[i] =  distrb(r1, opt, rscale ) * r1 * dr * 2.0d * pi;
          rnorm = rnorm + rp[i];
        } //enddo

    // normalize the vector
    //
        for(i=0; i<nprof;i++)
        {
        	rp[i] = rp[i] / rnorm;
        }
    //  take the fine bins and put them into the selection bins
    //
        tct = 0;
        for(iring=0;iring<nring;iring++)
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
        n_ring = new double[nprof];
    	for(i=nstart-1;i<ntot;i++)
    	{
    	
    	    // find the radial position bin
    	    ran = randm(seed);
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
        for(iring=0;iring<nring;iring++)
        {
        //do iring =  1, nring
          is = (iring) * dnring + 1;
          ie = (iring+1) * dnring; 

          r1 = (is)*dr + rin;
          r2 = (ie)*dr + rin;

          for(j=0; j<n_ring[iring]; j++)
          {
          //do j = 1, n_ring(iring)
            ran = randm(seed);
            r[i] = r1 + ran * (r2 - r1);
            i = i + 1;
          } //enddo
        } //enddo

    //   set the angular positions and orbital velocities
    //
    	for(i=nstart-1;i<ntot;i++)
    	{
        //do i=nstart,ntot
          angle[i] = 2.0d * pi * randm(seed);

          xp = r[i];
          yp = 0.0d;
          zp = 1e-3;  // included to make the log interpolation work better



    //      call INTERPOLATE_FORCE(xp, yp, zp, mass, rscale, rout, eps, fx, fy, fz, ftotal)
    //      ftotal =  mass / ( r(i)*r(i) + eps*eps )

          ival = dfint.df_index(r[i], rout);
          ftotal = mass * dfint.acceleration_particle[ival] * rs_internal * rs_internal;

          v[i] = Math.sqrt(ftotal * r[i]);
    //////      v(i) = sqrt( mass * r(i)/ ( r(i)*r(i) + eps*eps ) ) 

    //      write(*,'(i6,4f15.6)') ival, r(i), ftotal, v(i), rs_internal

        } //enddo
    //    stop

    	System.arraycopy(r,nstart-1,rVals,0,len);
    	System.arraycopy(v,nstart-1,vVals,0,len);
    	shell2(rVals,vVals);
    	System.arraycopy(rVals,0,r,nstart-1,len);
    	System.arraycopy(vVals,0,v,nstart-1,len);

    // set position and velocity based on the distribution parameters
    	for(i=nstart-1;i<ntot;i++)
    	{
        //do i= nstart,ntot

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


    ////    write(16,*) r(i)
  	    x0[i][1] = x3; 
	    x0[i][3] = y3; 
	    x0[i][5] = z3;
	    x0[i][2] = xv3  + randm(seed)*heat;
	    x0[i][4] = yv3  + randm(seed)*heat;
	    x0[i][6] = zv3  + randm(seed)*heat;
	    x0[i][0] = t0;  
	    /*
          x0(i,1) = x3; 
          x0(i,2) = y3 ;
          x0(i,3) = z3 ;
          x0(i,4) = xv3  + randm(seed)*heat;
          x0(i,5) = yv3  + randm(seed)*heat;
          x0(i,6) = zv3  + randm(seed)*heat;
    //    x0(6,i) = t0
     */

		} //enddo

    	return rVals;
      //end subroutine profile
    }
    public double circularVelocity(double rout, double mass, double eps)
    {
    	
    	int ival;
    	
    	
        double ftotal;
        double v = 0;

          ival = dfint.df_index(rout, rout);
          ftotal = mass * dfint.acceleration_particle[ival] * rs_internal * rs_internal;

          v = Math.sqrt(ftotal * rout);

    	return v;
      //end subroutine profile
    }

    
    
    
    // override to return proper integrator
    public integrator getIntegrator()
    {
    	return new dfintegrator();
    }
    
    // extend integrator to use new force
	public static class dfintegrator extends JSpam.integrator
	{
		double vxp = 0;
		double vyp = 0;
		double vzp = 0;
		/*
		real (kind=8) :: rrout1, rrout2

		  real (kind=8), dimension(:), allocatable :: df_force11, df_force22, df_forcen, c3n
		  integer (kind=4), dimension(:), allocatable :: ival11, ival22, ivaln

		  integer (kind=4) :: pn, pn1, pn2
		  */
		double rrout1 = 0;
		double rrout2 = 0;
		double df_force11[] = null;
		double df_force22[] = null;
		double df_forcen[] = null;
		double c3n[] = null;
		int ival11[] = null;
		int ival22[] = null;
		int ivaln[] = null;
		int pn = 0;
		int pn1 = 0;
		int pn2 = 0;

		double rad[] = null;
		double rho_halo[] = null;
		public double mass_halo[] = null;
		double rho_disk[] = null;
		public double mass_disk[] = null;
		double rho_bulge[] = null;
		public double mass_bulge[] = null;
		public double rho_total[] = null;
		public double mass_total[] = null;
		double masses[] = null;
		double radius[] = null;
		double density[] = null;
		double vr2[] = null;
		double vr[] = null;
		double new_vr2[] = null;
		double new_vr[] = null;
		double acceleration[] = null;
		public double acceleration_particle[] = null;
		double new_mass[] = null;
		double new_rho[] = null;
		double phi[] = null;
		
		double pscale = 0;
		double lnl = 0;
		
		double tmp1 = 0;
		//double tmp1 = pscale * rs_internal/ rmax_scale) * nnn;

		double mm1rs2 = 0;
		double mm2rs2 = 0;
		double mm3rs2 = 0;
		
		public dfintegrator()
		{
			this(true);
		}
		
		public dfintegrator(boolean initDist)
		{
			super();
			if(initDist)
			{
				init_distribution();
			}
			tmp1 = pscale * rs_internal/ 100.0d * nnn;			
		}
		
		public void init_rkvar(double x0[][], double mass1, double mass2, double epsilon1, double epsilon2, double t1, double p1, double t2, double p2, double rs1[], double rs2[], double ro1, double ro2, int nn, int n1, int n2)
		{
			super.init_rkvar(x0,mass1,mass2,epsilon1,epsilon2,t1,p1,t2,p2,rs1,rs2,ro1,ro2,nn,n1,n2);
			  pn = nn;
			  pn1 = n1;
			  pn2 = n2;
				int n = x0.length;

			  ival11 = new int[n];
			  ival22 = new int[n];
			  ivaln = new int[n];
			  
			  df_force11 = new double[n];
			  df_force22 = new double[n];
			  df_forcen = new double[n];
			  c3n = new double[n];
			  rrout1 = ro1;
			  rrout2 = ro2;
			  mm1rs2 = -m1*rs_internal2;
			  mm2rs2 = -m2*rs_internal2;
			  mm3rs2 = -m3*rs_internal2;
		}
		
		// override force calculation
		public void diffeq(double x[][])
		{
			//MonitorUtil.startCall("dfspam.diffeq");
			int n = x.length;
			int i;
			double df_sigma = 0;
			double df_rho = 0;
			double c1 = 0;
			double c2 = 0;
			double xvalue = 0;
			double v1 = 0;
			double v21 = 0;
			
			for(i=0;i<7;i++)
			{
				xn[i] = x[n-1][i];
			}
			
            // make temps to handle perturber galaxy
            xp = xn[1];
            yp = xn[3];
            zp = xn[5];
            
            vxp = xn[2];
            vyp = xn[4];
            vzp = xn[6];
            
            double tmp4 = xp*xp+yp*yp+zp*zp;
            double tmp5 = Math.sqrt(tmp4);
            tmp4 = -m3/(tmp4+eps2);
			  int iv = 0;
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
           // }

            /*
			// distance between the main galaxy and the particle
			  r21 = x(:,1)**2  + x(:,2)**2  + x(:,3)**2;
			  r1 = sqrt(r21);

			// distance between the companion and the particle
			  r22 = (x(:,1)-xn(1))**2  +(x(:,2)-xn(2))**2 + (x(:,3)-xn(3))**2;
			  r2 = sqrt(r22);


			// distance between the two galaxies - the tidal force
			  r2n = xn(1)**2 + xn(2)**2 + xn(3)**2;
			  rn = sqrt(r2n);
			*/    

        //    for(i=0; i<n; i++)
        //    {
            	//ival11[i] = df_index(r1[i],rrout1);
               	//ival22[i] = df_index(r2[i],rrout2);
               	//ivaln[i] = df_index(rn[i],rrout2);
              /* 	
  			  df_force11[i] = acceleration_particle[ival11[i]] * rs_internal2;
			  df_force22[i] = acceleration_particle[ival22[i]] * rs_internal2;
			  df_forcen[i]  = acceleration_particle[ivaln[i]]  * rs_internal2;
			  a1[i] = -m1*df_force11[i];
			  a2[i] = -m2*df_force22[i];
			  a3[i] = -m3*df_forcen[i];
*/
            	// we skip all of the extra array indexing, because no one uses ival11, ival22 or the df_forces
            	iv = df_index(r1[i],rrout1);
			  	//a1[i] = -m1*acceleration_particle[iv]*rs_internal2;
			  	a1[i] = acceleration_particle[iv]*mm1rs2;
			  	iv = df_index(r2[i],rrout2);
			  	//a2[i] = -m2*acceleration_particle[iv]*rs_internal2;
			  	a2[i] = acceleration_particle[iv]*mm2rs2;
			  	iv = df_index(rn[i],rrout2);
			  	//a3[i] = -m3*acceleration_particle[iv]*rs_internal2;
			  	a3[i] = acceleration_particle[iv]*mm3rs2;
			  	ivaln[i]=iv;
            }
            
            /*
			  do i = 1, n
			    ival11(i) =  df_index(r1(i), rrout1)
			  enddo

			  do i = 1, n
			    ival22(i)  =  df_index(r2(i), rrout2)
			  enddo

			  do i = 1, n
			    ivaln(i)  =  df_index(rn(i), rrout2)
			  enddo
*/
//combine loops?
            /*
            for(i=0; i<n; i++)
            {
			  //df_force11 = acceleration_particle[ival11] * rs_internal * rs_internal;
			  //df_force22 = acceleration_particle[ival22] * rs_internal * rs_internal;
			  //df_forcen  = acceleration_particle[ivaln]  * rs_internal * rs_internal;
			  df_force11[i] = acceleration_particle[ival11[i]] * rs_internal2;
			  df_force22[i] = acceleration_particle[ival22[i]] * rs_internal2;
			  df_forcen[i]  = acceleration_particle[ivaln[i]]  * rs_internal2;
            }
*/
			// get the forces, sigma and rho, and rescale them
			  //df_sigma  =  new_vr2[ivaln[0]] * rs_internal * rs_internal;
			  //df_rho    =  new_rho[ivaln[0]] * ( rs_internal * rs_internal * rs_internal );
			  df_sigma  =  new_vr2[ivaln[0]] * rs_internal2;
			  df_rho    =  new_rho[ivaln[0]] * ( rs_internal3 );


			  /*
			// interpolated forces 
			  for(i=0; i<n; i++)
			  {
				  a1[i] = -m1*df_force11[i];
				  a2[i] = -m2*df_force22[i];
				  a3[i] = -m3*df_forcen[i];
				  
			  //a1 = -m1 * df_force11
			  //a2 = -m2 * df_force22
			  //a3 = -m3 * df_forcen
			  }
			  */
			  
			// df
			    //v21 = xn(4)*xn(4)+xn(5)*xn(5)+xn(6)*xn(6)
			  v21 = vxp*vxp + vyp*vyp + vzp*vzp;
			    v1  = Math.sqrt(v21);
double v1Inv = 1.0d/v1;

			    xvalue = v1 / df_sigma;
			    c1 = erf(xvalue) - 2.0d * xvalue / sqrtpi * Math.exp(-xvalue*xvalue);

			// df formula with G=1
			    c2  = 4.0d * pi * m2 * lnl / v21;
			    for(i=0; i<n; i++)
			    {
			    	c3n[i] = 0;
			    }
			    
			    for(i=pn1; i<n; i++)
			    {
			    	c3n[i] = c1*c2 *df_rho;
			    }
			    //c3n(1:n-1) = 0.0d;
			    //c3n(pn1+1:n)  = c1 * c2 * df_rho;

			//    c3 = 0.0d0

			  // this is a correction to prevent NaN errors in the vectorized
			  // function evalution at the location of the second mass
			  r2[n-1] = 1.0d;


			  // calculate the RHS of the diffeq

			  for(i=0; i<n; i++)
			  {
				  tmp1 = a1[i]/r1[i];
	              tmp2 = a2[i]/r2[i];
	              tmp3 = a3[i]/rn[i];
	              //tmp4 = c3n[i]/v1;
	              tmp4 = c3n[i]*v1Inv;
	              
	              //tmp4  = 0;
	              f[i][0] = 1.0d;
	              f[i][1] = x[i][2];
	              f[i][2] = tmp1*x[i][1] + tmp2*(x[i][1]-xn[1]) + tmp3*xn[1] - tmp4*x[i][2];
	              //f[i][2] = a1[i]*x[i][1]/r1[i] + a2[i]*(x[i][1]-xn[1])/r2[i] + a3[i]*xn[1]/rn[i];
	              f[i][3] = x[i][4];
	              f[i][4] = tmp1*x[i][3] + tmp2*(x[i][3]-xn[3]) + tmp3*xn[3] - tmp4*x[i][4];
	              //f[i][4] = a1[i]*x[i][3]/r1[i] + a2[i]*(x[i][3]-xn[3])/r2[i] + a3[i]*xn[3]/rn[i];
	              f[i][5] = x[i][6];
	              f[i][6] = tmp1*x[i][5] + tmp2*(x[i][5]-xn[5]) + tmp3*xn[5] - tmp4*x[i][6];
	              //f[i][6] = a1[i]*x[i][5]/r1[i] + a2[i]*(x[i][5]-xn[5])/r2[i] + a3[i]*xn[5]/rn[i];
	              
	              // hmmmm.... it appears that the original author has changed the indices from x,vx,y,vy,z,vz 
	              // to x,y,z, vx,vy,vz
	              
	              //f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn - c3n * xn(4)/ v1 
	              //f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn - c3n * xn(5)/ v1
	              //f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn - c3n * xn(6)/ v1

			  }
		      //MonitorUtil.endCall("dfspam.diffeq");
			  
			  return;
		}
	
		public void init_distribution()
		{
			synchronized(syncObject)
			{
				if(INIT_INT == null)
				{
					INIT_INT = new dfintegrator(false);
					INIT_INT.init_distributionImpl();
					INIT_INT.tmp1 = INIT_INT.pscale * rs_internal/ 100.0d * nnn;
				}
			}
			
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
			
			// do memory copies
		    System.arraycopy(INIT_INT.rad,0,rad,0,nnn);
		    System.arraycopy(INIT_INT.rho_halo,0,rho_halo,0,nnn);
		    System.arraycopy(INIT_INT.mass_halo,0,mass_halo,0,nnn);
		    System.arraycopy(INIT_INT.rho_disk,0,rho_disk,0,nnn);
		    System.arraycopy(INIT_INT.mass_disk,0,mass_disk,0,nnn);
		    System.arraycopy(INIT_INT.rho_bulge,0,rho_bulge,0,nnn);
		    System.arraycopy(INIT_INT.mass_bulge,0,mass_bulge,0,nnn);
		    System.arraycopy(INIT_INT.rho_total,0,rho_total,0,nnn);
		    System.arraycopy(INIT_INT.mass_total,0,mass_total,0,nnn);
		    System.arraycopy(INIT_INT.masses,0,masses,0,nnn);
		    System.arraycopy(INIT_INT.radius,0,radius,0,nnn);
		    System.arraycopy(INIT_INT.density,0,density,0,nnn);
		    System.arraycopy(INIT_INT.vr2,0,vr2,0,nnn);
		    System.arraycopy(INIT_INT.vr,0,vr,0,nnn);
		    System.arraycopy(INIT_INT.new_vr2,0,new_vr2,0,nnn);
		    System.arraycopy(INIT_INT.new_vr,0,new_vr,0,nnn);
		    System.arraycopy(INIT_INT.acceleration,0,acceleration,0,nnn);
		    System.arraycopy(INIT_INT.acceleration_particle,0,acceleration_particle,0,nnn);
		    System.arraycopy(INIT_INT.new_mass,0,new_mass,0,nnn);
		    System.arraycopy(INIT_INT.new_rho,0,new_rho,0,nnn);
		    System.arraycopy(INIT_INT.phi,0,phi,0,nnn);
			
		    lnl = 0.001;
			pscale = 1.0d;
		}
		
	// df_module.f90
		//
		//     ----------------------------------------------------------------
		//     ----------------------------------------------------------------
		//
		//module df_module
		//     
		//
		//     -----Description: dynamical friction module
		//
		//
		//
		//     ----------------------------------------------------------------
		//
		//  use parameters_module
		//  implicit none
		//
		//     -----Variable declarations
		//          ---------------------
		//
		
		/*
		    integer (kind=4), parameter :: nnn = 10000
		    real (kind=8), dimension(nnn) :: rad    
		    real (kind=8), dimension(nnn) :: rho_halo, mass_halo
		    real (kind=8), dimension(nnn) :: rho_disk, mass_disk
		    real (kind=8), dimension(nnn) :: rho_bulge, mass_bulge
		    real (kind=8), dimension(nnn) :: rho_total, mass_total
		    real (kind=8), dimension(nnn) :: masses, radius, density
		    real (kind=8), dimension(nnn) :: vr2, vr, new_vr2, new_vr
		    real (kind=8), dimension(nnn) :: acceleration, acceleration_particle
		    real (kind=8), dimension(nnn) :: new_mass, new_rho, phi

		    real (kind=8), parameter :: rs_internal = 10.0d0

		    real (kind=8) :: pscale

		    real (kind=8) :: lnl
		    */
		//
		//     ----------------------------------------------------------------
		//
		//CONTAINS
		//
		//     ----------------------------------------------------------------
		//

		//
		//     ----------------------------------------------------------------
		//     ----------------------------------------------------------------
		//
		public void init_distributionImpl()
		{
			
		//  subroutine init_distribution
		//     
		//
		//     -----Description: initializes the distribution 
		//
		//
		//          on input:  
		//
		//          on output:
		//          
		//
		//     ----------------------------------------------------------------
		//
		//    implicit none
		//
		//     -----Variable declarations
		//          ---------------------
		//


		//
		//     ----------------------------------------------------------------
		//
		    double rmax;
		    double mold, dmold, mtotal, mtot;
		    double rscale;
		    double dx, x;
		    double alphahalo, qhalo, gammahalo, mhalo, rchalo, rhalo, epsilon_halo;
		    double zdisk, zdiskmass, hdisk, zdiskmax;
		    double hbulge, mbulge;
		    double rho_tmp;
		    double G, factor;
		    double r, m;
		    double p1, rd, rho_local;
		    double p, rr, dr, rh, dp, mnew, dm;
		    double acc_merge, rad_merge, acc;
		    double xmax;
		    
		    int j, nmax, k, nmerge, ntotal, jj;


		////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3

		//!!!!
		// set the constant for dynamical friction
		    lnl = 0.001;

		//////////
		// set up the parameters for the halo
		    mhalo = 5.8d;  
		    rhalo = 10.0d;
		    rchalo = 10.0d;
		    gammahalo = 1.0d;
		    epsilon_halo = 0.4d; 
		    

		//////////
		// derive additional constants
		    qhalo = gammahalo / rchalo;
		    //alphahalo = 1.0d / ( 1.0d - sqrtpi * qhalo * exp(qhalo**2) * (1.0d - erf(qhalo)) );
		    alphahalo = 1.0d / ( 1.0d - sqrtpi * qhalo * Math.exp(qhalo*qhalo) * (1.0d - erf(qhalo)) );

		//////////
		// set the integration limits and zero integration constants
		    rmax = 20;
		    nmax = 1000;
		    nmax = 2000;
		    dr = rmax / (nmax);
		    mold = 0;

		    rscale = 5;
		//    ntotal = nmax * rscale
		    ntotal = nnn;
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
			

		//////////
		// set the limits for integration, and zero integration constants
		    k = nmax / 2;
		    dx = 1.0  / k;
		    x = 0.0d;
		    dmold = 0.0d;
		    mtot = 0.0d;
		    //rad = new double[nnn];
		    m = 0.0d;
		    G = 1;
		    

		//////////
		// set the fundamental disk parameters
		    zdisk = 0.2;
		    zdiskmax = 3.5;
		    hdisk = 1.0;
		    

		//////////
		// set the fundamental bulge parameters
		    hbulge = 0.2;
		    mbulge = 0.3333;
		    
		    
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
		//!!!! set up the radius array
		    for(j=0; j<nmax; j++)
		    {
		    //do j = 1, nmax
		      x = x + dx;
		      rad[j]= x * rchalo;
		    //end do
		    }

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
		//!!!!
		 
		   dr = rad[1] - rad[0];
		   dx = dr / rchalo; 
		   
		   //do j = 1, nmax
		   for(j=0; j<nmax; j++)
		   {
		    
		// set the position
		      r = rad[j];
		      x = r / rchalo;

		// calculate the local rho based 
		      //rho_tmp = alphahalo / (2*sqrtpi**3 ) * (exp(-x**2) / (x**2 + qhalo**2))
		      rho_tmp = alphahalo / (twosqrtpicubed) * (Math.exp(-x*x) / (x*x + qhalo*qhalo));

		// renormalize for the new halo size
		      rho_tmp = rho_tmp / ( rchalo * rchalo * rchalo); 

		// calculate mass in local shell, and update total mass
		//      dm = rho_tmp * 4 * pi * x * x *dx
		      dm = rho_tmp * 4 * pi * r * r *dr;
		      mtot = mtot + dm;

		// store values in an array
		      rho_halo[j] = rho_tmp * mhalo; 
		      mass_halo[j] = mtot * mhalo;
		   }
		    //end do

		//////////
		// now calculate the potential
		   //do j = 1, nmax
		   for(j=0; j<nmax; j++)
		   {
		      r = rad[j];
		      m = mass_halo[j];
		      p1 = -G * m / r;
		      phi[j] = p1;
		   } //end do


		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
		// disk model

		////!!!
		// loop over the distribution
		   //do j = 1, nmax
		   for(j=0; j<nmax; j++)
		   {

		// set the radius
		      rd = rad[j];
		  
		// find the local density in the disk
		      //rho_local  = exp(-rd/hdisk)/ (8*pi*hdisk**2.0d0) 
		      rho_local  = Math.exp(-rd/hdisk)/ (8*pi*hdisk*hdisk);
		      rho_disk[j] = rho_local;
		      
		// find the mass in the spherical shell
		      mnew = 4 * pi * rho_local * rd *rd * dr;
		      
		      mass_disk[j] = mnew + mold;
		      mold = mass_disk[j];
		    } //end do




		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
		// bulge model

		//!!!!
		// loop over the distribution
		    mold = 0.0;

		   //do j = 1, nmax
		   for(j=0; j<nmax; j++)
		   {
		// set the radius
		      rd = rad[j];
		  
		// find the local density in the disk
		      //rho_local  = Math.exp(-rd**2/hbulge**2);
		      rho_local  = Math.exp(-rd*rd/(hbulge*hbulge));
		      rho_bulge[j] = rho_local;

		// find the mass in the spherical shell
		      mnew = 4 * pi * rho_local * rd *rd * dr;
		      
		      mass_bulge[j] = mnew + mold;
		      mold = mass_bulge[j];
		    } //end do

		// renormalize distribution
		    factor = mbulge / mass_bulge[nmax-1];

		   for(j=0; j<nmax; j++)
		   {
		      mass_bulge[j] = mass_bulge[j] * factor;
		      rho_bulge[j]  = rho_bulge[j]  * factor;
		    } //end do

		  
		    //dr = rad(2) - rad(1)      ;
		    dr = rad[1] - rad[0]      ;
		   
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
		    j = 0;
		    mass_total[j]=  (mass_halo[j] + mass_disk[j] + mass_bulge[j]); 
		    r = rad[j];
		    rho_total[j] = mass_total[j] /  (4.0d/3.0d * pi * r * r * r);
		    dr = rad[1] - rad[0];

		   for(j=1; j<nmax; j++)
		   {
		    //do j = 2,nmax
		      r = rad[j];
		      mass_total[j]=  (mass_halo[j] + mass_disk[j] + mass_bulge[j]); 

		      dm = mass_total[j] - mass_total[j-1];
		      rho_total[j] = dm / (4 * pi * r * r * dr);

		    } //end do

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// find the velocity dispersion v_r**2

		    masses = mass_total;
		    radius = rad;
		    density = rho_total;

		   for(j=0; j<nmax; j++)
		   {
//		    do j = 1,nmax

		      p = 0.0d;
		      rr = radius[j];
		      dr = radius[nmax-1] / nmax;
		      for(jj=j; jj<nmax; jj++)
		      {
		      //do jj = j,nmax
		      
		        m  = masses[jj];
		        rh = density[jj];
		        rr = rr + dr;
		        
		        dp = rh * G * m / (rr*rr) * dr;
		        p = p + dp;
		      } //end do
		      
		      vr2[j] = 1/density[j] * p;
		      vr[j] = Math.sqrt(vr2[j]);
		    } //end do


		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// find the velocity dispersion v_r**2
		    masses = mass_total;
		    radius = rad;
		    density = rho_total;
		    
		   for(j=0; j<nmax; j++)
		   {

//		    do j = 1,nmax

		      p = 0.0d;
		      rr = radius[j];
		      dr = radius[nmax-1] / nmax;
    		   for(jj=j; jj<nmax; jj++)
	    	   {

		//      do jj = j,nmax
		      
		      m  = masses[jj];
		      rh = density[jj];
		      rr = rr + dr;
		      
		      dp = rh * G * m / (rr*rr) * dr;
		      p = p + dp;
		    } //end do
		  
		    vr2[j] = 1/density[j] * p;
		    vr[j] = Math.sqrt(vr2[j]);
		   }//  enddo



		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// find the accelerations felt by the particles and center of mass
		    masses = mass_total;
		    radius = rad;
		    density = rho_total;

		   for(j=0; j<nmax; j++)
		   {

		    //do j = 1,nmax
		      rr = radius[j];
		      m  = masses[j];
		      acceleration[j] = G * m / (rr*rr);
		    } //end do


		    //acceleration_particle = acceleration;
		    System.arraycopy(acceleration,0,acceleration_particle,0,nnn);
		    nmerge = 50;
		    acc_merge = acceleration[nmerge];
		    rad_merge = rad[nmerge];
		   for(j=0; j<nmax; j++)
		   {

		    //do j = 1, nmerge
		      rr = radius[j];
		      m  = masses[j];

		// smoothed acceleration
		      acc = G * m / (rr*rr + 0.1* (rad_merge -rr)); 
		      acceleration_particle[j] = acc;
		      
		    } //end do

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// rederive the masses from the new particle acceleration
		    radius = rad;
		    dr = rad[1] - rad[0];

		// find the accelerations felt by the particles and center of mass
		    radius = rad;
		   for(j=0; j<nmax; j++)
		   {

//		    do j = 1, nmax
		      rr = radius[j];
		      new_mass[j] = rr*rr * acceleration_particle[j] / G;
		      new_rho[j]  = new_mass[j] / (4 * pi * rr * rr * dr);
		    } //end do


		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// find the velocity dispersion v_r**2 using the new density and masses
		    masses = new_mass;
		    radius = rad;
		    density = new_rho;
		    
			   for(j=0; j<nmax; j++)
			   {

//		    do j = 1, nmax

		      p = 0.0d;
		      rr = radius[j];
		      dr = radius[nmax-1] / nmax;
   		      for(jj=j; jj<nmax; jj++)
		      {
		      //do jj = j,nmax
		      
		      m  = masses[jj];
		      rh = density[jj];
		      rr = rr + dr;
		      
		      dp = rh * G * m / (rr*rr) * dr;
		      p = p + dp;
		    } //end do
		    
		    new_vr2[j] = 1/density[j] * p;
		    new_vr[j] = Math.sqrt(new_vr2[j]);
		    
		  } //end do


		////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// extend the values to large rmax
	   for(j=nmax; j<ntotal; j++)
	   {

		//do  j= nmax+1, ntotal
		  //mass_total[j] = mass_total(nmax)
		  //mass_halo[j] = mass_halo(nmax)
		  //mass_disk[j] = mass_disk(nmax)
		  //mass_bulge[j] = mass_bulge(nmax)
		 // new_mass[j] = new_mass(nmax)
		   
		  mass_total[j] = mass_total[nmax-1];
		  mass_halo[j] = mass_halo[nmax-1];
		  mass_disk[j] = mass_disk[nmax-1];
		  mass_bulge[j] = mass_bulge[nmax-1];
		  new_mass[j] = new_mass[nmax-1];

		//  rho_total[j] = 1e-3
		//  new_rho[j] = new_rho(nmax)
		  rho_total[j] = 0.0d;
		  new_rho[j]   = 0.0d;

		  vr[j]      = 1e-6;
		  vr2[j]     = 1e-6;
		  new_vr[j]  = 1e-6;
		  new_vr2[j] = 1e-6;

		  m = mass_total[nmax-1];
		  rr = rad[nmax-1] + dr*(j - nmax);
		  rad[j] = rr;
		  acc = G * m / (rr*rr);  
		  acceleration_particle[j] = acc;
		  acceleration[j] = acc;

		} //end do



		////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// normalize to the unit mass

	   for(j=0; j<ntotal; j++)
	   {
		//do  j= 1, ntotal
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

		  rad[j]         = rad[j] ;

		  acceleration_particle[j] = acceleration_particle[j] / 7.13333d;
		  acceleration[j]          = acceleration[j]  / 7.13333d;

		//!  write(11,*) rad[j], new_rho[j], new_mass[j],  new_vr[j]
		} //end do



		pscale = 1.0d;
	    }
		//end subroutine init_distribution


        public int df_index(double rin, double rs)
		{
		//
		//     ----------------------------------------------------------------
		//     ----------------------------------------------------------------
		//
		//  integer (kind=4) function df_index(rin, rs)
		//     
		//
		//     -----Description: interpolates the force for a given particle
		//
		//
		//          on input:  
		//
		//          on output:
		//          
		//
		//     ----------------------------------------------------------------
		//
		//    implicit none
		//
		//     -----Variable declarations
		//          ---------------------
		//
		//    real (kind=8), intent(in) :: rin
		//    real (kind=8), intent(in) :: rs
		//    integer (kind=4) :: ival
		//    real (kind=8), parameter :: rmax_scale = 100.0d0

		//
		//     ----------------------------------------------------------------
		//
			//MonitorUtil.startCall("dfspam.df_index");

			int ival = 0;
			//double rmax_scale = 100.0d;
		    
		//    ival = min(int(  (rin * rs_internal/ rmax_scale) * nnn + 1), nnn) 
		    //ival = Math.min((int)((rin * pscale * rs_internal/ rmax_scale) * nnn + 1), nnn); 
		    ival = Math.min((int)((rin * tmp1) + 1), nnn); 
		    ival-=1;
		    if(ival < 0)
		    {
		    	ival = 0;
		    }
			//MonitorUtil.endCall("dfspam.df_index");
		    return ival;
		}
		  //end function df_index
		//end module df_module
	}
	
	//http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
	//	 fractional error in math formula less than 1.2 * 10 ^ -7.
    // although subject to catastrophic cancellation when z in very close to 0
    // from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
    public static double erf(double z) {
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
