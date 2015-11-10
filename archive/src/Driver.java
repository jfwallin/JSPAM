package edu.gmu.cds.sim.trim;

import edu.gmu.cds.sim.trim.JSpam.integrator;

import java.io.*;

public class Driver 
{
    public static JSpamParams getParams(String str)
    {
        JSpamParams jsp = new JSpamParams();
        String sa[] = str.split(",");
        
        jsp.rx = parseDouble(sa[0]);
        jsp.ry = parseDouble(sa[1]);
        jsp.rz = parseDouble(sa[2]);
        jsp.vx = parseDouble(sa[3]);
        jsp.vy = parseDouble(sa[4]);
        jsp.vz = parseDouble(sa[5]);
        jsp.mass1 = parseDouble(sa[6]);
        jsp.mass2 = parseDouble(sa[7]);
        jsp.rout1 = parseDouble(sa[8]);
        jsp.rout2 = parseDouble(sa[9]);
        jsp.phi1 = parseDouble(sa[10]);
        jsp.phi2 = parseDouble(sa[11]);
        jsp.theta1 = parseDouble(sa[12]);
        jsp.theta2 = parseDouble(sa[13]);

        jsp.setAlpha(0); 
        jsp.setBeta(0); 
        jsp.setGamma(0); 

        return jsp;
    }
    
    public static double parseDouble(String str)
    {
        double num = 0;
        try{num = Double.parseDouble(str);}catch(Exception ex){};
        return num;
    }
    
    public static void runSim(JSpamParams jsp, double tmin, double scale, String outputDir)
    {
        //JSpam jspam = new JSpam();
        JSpam jspam = new DFJSpam();
        double rx = jsp.rx;
        double ry = jsp.ry;
        double vz = jsp.vz;
        
        // run simulation
        double rscale1[] = {3.0,3.0,3.0};
        double rscale2[] = {3.0,3.0,3.0};
        
        // allocate the space for the particles
        double x0[][] = new double[jsp.n+1][7];
        double xout[][] = new double[jsp.n+1][7];
          
        // set the initial values
        // velocities are set to negative because
        // we are going to iterate back in time first        
        double rv[] = {0.0,rx,-jsp.vx,ry,-jsp.vy,jsp.rz,-vz};
        double rv4min[] = {0.0,rx,-jsp.vx,ry,-jsp.vy,jsp.rz,-vz};
        
        // let's determine what time to use
        double tmpT = tmin;
        if(tmpT == 0)
        {
            double tminVals[] = jspam.getTStart(rv4min,jsp.mass1,jsp.mass2,jsp.epsilon1,jsp.epsilon2,jsp.h,-30.0d,10.0d*jsp.rout1,jsp.rout1,jsp.rout2);
            tmpT = tminVals[0];
            if(tmpT < -12.0d)
            {
                tmpT = -5;
            }
            
            if(Math.abs(tmpT) < jsp.h)
            {
                tmpT = -5;
            }
        }
        
        double t0 = tmpT;
        int i,j;
        int d = 7;
        
        jsp.nstep = (int)(Math.abs(t0)/jsp.h);
        
        // adjust inner radius to look pretty
        double tmprin1 = jsp.rin1;
        double tmprin2 = jsp.rin2;
        
        if(tmprin1 > 0.1*jsp.rout1)
        {
            //System.out.println("rin1 override");
            tmprin1 = 0.1*jsp.rout1;
            jsp.rin1 =tmprin1;
        }
        if(tmprin2 > 0.1*jsp.rout2)
        {
            //System.out.println("rin2 override");
            tmprin2 = 0.1*jsp.rout2;
            jsp.rin2 = tmprin2;
        }
        
        if(tmprin1 < 0.075*jsp.rout1)
        {
            //System.out.println("rin1 override");
            tmprin1 = 0.075*jsp.rout1;
            jsp.rin1 =tmprin1;
        }
        
        if(tmprin2 < 0.075*jsp.rout2)
        {
            //System.out.println("rin2 override");
            tmprin2 = 0.075*jsp.rout2;
            jsp.rin2 = tmprin2;
        }
        
        
        // create the disks
        jspam.profile(jsp.rin1, jsp.rout1, rscale1, 1, jsp.n1, jsp.mass1, jsp.epsilon1,
                jsp.theta1, jsp.phi1, jsp.opt1, jsp.heat, jsp.seed, jsp.time, x0);
         
        jspam.profile(jsp.rin2, jsp.rout2, rscale2, jsp.n1+1, jsp.n, jsp.mass2, jsp.epsilon2,
                jsp.theta2, jsp.phi2, jsp.opt2, jsp.heat, jsp.seed, jsp.time, x0);


        // set the perturber galaxy position
        jspam.perturberPosition(rv,jsp.mass1,jsp.mass2,jsp.epsilon1,jsp.epsilon2,jsp.h,jsp.n,jsp.n1,t0,x0,jsp.rout1,jsp.rout2);
        
        
        integrator intg = jspam.getIntegrator();
        intg.init_rkvar(x0,jsp.mass1,jsp.mass2,jsp.epsilon1,jsp.epsilon2,jsp.theta1,jsp.phi1,jsp.theta2,jsp.phi2,null,null,jsp.rout1,jsp.rout2,jsp.n,jsp.n1,jsp.n2);
          
        double g1[][] = new double[jsp.n1][6];
        double g2[][] = new double[jsp.n2+1][6];
        double g3[][] = new double[jsp.nstep][6];
        
        for(i=0; i<jsp.n1; i++)
        {
            g1[i][0]=x0[i][1];
            g1[i][1]=x0[i][3];
            g1[i][2]=x0[i][5];
            g1[i][3]=x0[i][2];
            g1[i][4]=x0[i][4];
            g1[i][5]=x0[i][6];
        }
        
        for(i=0; i<jsp.n2+1; i++)
        {
            g2[i][0]=x0[i+jsp.n1][1];
            g2[i][1]=x0[i+jsp.n1][3];
            g2[i][2]=x0[i+jsp.n1][5];
            g2[i][3]=x0[i+jsp.n1][2];
            g2[i][4]=x0[i+jsp.n1][4];
            g2[i][5]=x0[i+jsp.n1][6];
        }        

        outputParticles(outputDir+"/step0000.dat",x0);        

        int nv = x0.length-1;
        int ind = 0;
        String str = null;

        for(jsp.istep = 1;jsp.istep<=jsp.nstep;jsp.istep++)
        {
            intg.rk4(x0, jsp.h, xout);
            for(i=0;i<jsp.n+1;i++)
            {
                for(j=0;j<d;j++)
                {
                    x0[i][j] = xout[i][j];        
                }
            }
            for(i=0; i<jsp.n1; i++)
            {
                g1[i][0]=x0[i][1];
                g1[i][1]=x0[i][3];
                g1[i][2]=x0[i][5];
                g1[i][3]=x0[i][2];
                g1[i][4]=x0[i][4];
                g1[i][5]=x0[i][6];
            }
            
            for(i=0; i<jsp.n2+1; i++)
            {
                ind = i+jsp.n1;
                g2[i][0]=x0[ind][1];
                g2[i][1]=x0[ind][3];
                g2[i][2]=x0[ind][5];
                g2[i][3]=x0[ind][2];
                g2[i][4]=x0[ind][4];
                g2[i][5]=x0[ind][6];
            }
            
            g3[jsp.istep-1][0] = x0[nv][1];
            g3[jsp.istep-1][1] = x0[nv][3];
            g3[jsp.istep-1][2] = x0[nv][5];
            g3[jsp.istep-1][3] = x0[nv][2];
            g3[jsp.istep-1][4] = x0[nv][4];
            g3[jsp.istep-1][5] = x0[nv][6];

            if(jsp.istep%10==0)
            {
                str = String.valueOf(jsp.istep);
                while(str.length() < 4)str = "0"+str;
                outputParticles(outputDir+"/step"+str+".dat",x0);        
            }
        }
        str = String.valueOf(jsp.istep);
        while(str.length() < 4)str = "0"+str;
        outputParticles(outputDir+"/step"+str+".dat",x0);        
        
        jsp.simValues = x0;
    }

    public static void outputParticles(String file, double x0[][])
    {
        BufferedWriter osw = null; 
        OutputStream os = null;
        try
        {
            os = new FileOutputStream(file);
            osw = new BufferedWriter(new OutputStreamWriter(os));

            int size = x0.length;
            double dtmp[] = null;

            for(int i=1; i<size; i++)
            {
                dtmp = x0[i];

                for(int j=0; j<6; j++)
                {
                    osw.write(formatDouble(dtmp[j]));
                }
                osw.write("\n"); // unix style line ending
            }
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
        finally
        {
            if(osw != null)try{osw.flush();osw.close();}catch(Exception ex){};
            if(os != null)try{os.close();}catch(Exception ex){};
        }
    }

    /**
     * This method formats a double as a string to match the Fortran format.
     * For whatever reason, the DecimalFormat object does not seem to do
     * what is needed here.
     */
    public static String formatDouble(double num)
    {
        String str = "";
        String sign = "";
        if(num < 0)
        {
            sign = "-";
            num = -num;
        }

        double exp = Math.log10(num);

        exp = (int)exp; 
        int expForm = (int)exp; 
        String expSign = "+";
        if(exp < 0)
        {
            expSign = "-";
            expForm = -expForm;
        }

        String expStr = String.valueOf(expForm);
        while(expStr.length() < 2)
        {
            expStr = "0" + expStr;
        }
        expStr = "E" + expSign + expStr;

        double div = Math.pow(10.0,exp);
        num /= div;

        str = String.valueOf(num);
        if(str.length() > 10)
        {
            str = str.substring(0,10);
        }

        while(str.length() < 10)
        {
            str = str + "0";
        }
 
        str = sign+str+expStr; 

        while(str.length() < 16)
        {
            str = " " + str;
        }

        return str;
    }

    public static void main(String args[])
    {
        JSpamParams jsp = getParams(args[0]);
        jsp.adjustTotalN(4000);
        double tmin = parseDouble(args[1]);
        double scale = parseDouble(args[2]);
System.out.println(args[3] + "\t" + tmin + "\t" + scale);
        runSim(jsp,tmin,scale,args[3]);

        System.exit(0);
    }
}
