package edu.gmu.cds.sim;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import java.text.DecimalFormat;

import java.util.ArrayList;
import java.util.List;

public class IOUtil
{
    
    public static void outputParticles(String file, double x0[][])
    {
        outputParticles(new File(file), x0);
    }

    public static void outputParticles(File file, double x0[][])
    {
        FileOutputStream fos = null;

        try
        {
            fos = new FileOutputStream(file);

            outputParticles(fos, x0);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }        
        finally
        {
            if(fos != null)
            {
                try{fos.flush();}catch(Exception ex){};
                try{fos.close();}catch(Exception ex){};
            }
        }
    }

    public static void outputParticles(OutputStream os, double x0[][])
    {
        BufferedWriter osw = null;

        try
        {
            osw = new BufferedWriter(new OutputStreamWriter(os));

            int size = x0.length;
            double dtmp[] = null;

            for(int i=0; i<size; i++)
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
            if(osw != null)try{osw.flush();}catch(Exception ex){};
            if(osw != null)try{osw.close();}catch(Exception ex){};
        }
    }

    public static double[][] inputParticles(String file)
    {
        double d[][] = null;

        d = inputParticles(new File(file));

        return d;
    }

    public static double[][] inputParticles(File file)
    {
        double d[][] = null;

        long size = file.length();
        int partEst = (int)(size/97);

        FileInputStream fis = null;

        try
        {
            fis = new FileInputStream(file);

            d = inputParticles(fis,partEst);
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
        finally
        {
            if(fis != null)try{fis.close();}catch(Exception ex){};
        }

        return d;
    }

    public static double[][] inputParticles(InputStream is)
    {
        double d[][] = null;

        d = inputParticles(is, 1000);

        return d;
    }

    public static double[][] inputParticles(InputStream is, int partEst)
    {
        double d[][] = null;
        LineNumberReader lnr = null;

        try
        {
            if(partEst < 0)
            {
                partEst = 1000;
            }

            // for temporary storage of particles
            double dtmp[] = null;
            List parts = new ArrayList(partEst); 

            lnr = new LineNumberReader(new InputStreamReader(is)); 

            String line = null;
            int offset = 0;

            line = lnr.readLine();
            String sa[] = null;
            while(line != null)
            {
                line = line.trim();
                if(line.length() > 80 && !line.startsWith("#"))
                {
                    while(line.contains("  "))
                    {
                        line = line.replaceAll("  ", " ");
                    }
                    sa = line.split(" ");
                    dtmp = new double[6];
                    for(int i=0; i<6; i++)
                    {
                        dtmp[i] = parseDouble(sa[i]);
                    }
                    parts.add(dtmp);
                }
/*
                // 16 * 6
                if(line.length() == 96)
                {
                    dtmp = new double[6];

                    offset = 0;
                    for(int i=0; i<6; i++)
                    {
                        dtmp[i] = parseDouble(line.substring(offset,offset+16));
                        offset+=16;
                    }

                    parts.add(dtmp);
                }
*/
                line = lnr.readLine();
            }

            // copy to 2D array
            int size = parts.size();

            d = new double[size][6];
            for(int i=0; i<size; i++)
            {
                d[i] = (double[])parts.get(i);
            } 
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
        finally
        {
            if(lnr != null)try{lnr.close();}catch(Exception ex){};
        }
 
        return d;
    }

    /**
     * Parse a double from the string.  Catch the exception and return 0 as the default.
     */
    public static double parseDouble(String str)
    {
        double num = 0;
        try{num = Double.parseDouble(str.trim().replace("D","E"));}catch(Exception ex){};
        return num;
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

    /**
     * parses line into key value pair with key is string and value is double
     */
    protected static Object[] parseLine(String line, Object oa[])
    {
        if(oa == null)oa = new Object[2];
        String label = null;
        Double value = null;

        if(line.length() > 3)
        {
            char c = line.charAt(0);
            if(!(c=='!' || c=='#' || c=='/'))
            {
                String sa[] = line.split("=");
                label = sa[0];
                if(sa.length > 1)
                {
                    double val = parseDouble(sa[1].trim());
                    value = new Double(val);
                }
            }
        }
        
        oa[0] = label;
        oa[1] = value;

        return oa;
    }

    public static void readParameterFile(Parameters params, String file)
    {
        FileReader fr = null;
        try
        {
            fr = new FileReader(file);
            LineNumberReader lnr = new LineNumberReader(fr);
            String line = null;
            line = lnr.readLine();
            Object oa[] = null;
            String label = null;
            Double value = null;
            double val = 0;
            while(line != null)
            {
                oa = parseLine(line, oa);
                label = (String)oa[0];
                value = (Double)oa[1];

                if(label != null && value != null)
                {
                    val = value.doubleValue();
                    if(label.equals("potential_type")) params.potential_type = (int)val;
                    else if(label.equals("mass1")) params.mass1 = val;
                    else if(label.equals("mass2")) params.mass2 = val;
                    else if(label.equals("epsilon1")) params.epsilon1 = val;
                    else if(label.equals("epsilon2")) params.epsilon2 = val;
                    else if(label.equals("rin1")) params.rin1 = val;
                    else if(label.equals("rin2")) params.rin2 = val;
                    else if(label.equals("rout1")) params.rout1 = val;
                    else if(label.equals("rout2")) params.rout2 = val;
                    else if(label.equals("theta1")) params.theta1 = val;
                    else if(label.equals("theta2")) params.theta2 = val;
                    else if(label.equals("phi1")) params.phi1 = val;
                    else if(label.equals("phi2")) params.phi2 = val;
                    else if(label.equals("opt1")) params.opt1 = (int)val;
                    else if(label.equals("opt2")) params.opt2 = (int)val;
                    else if(label.equals("heat1")) params.heat1 = val;
                    else if(label.equals("heat2")) params.heat2 = val;
                    else if(label.equals("n1")) params.n1 = (int)val;
                    else if(label.equals("n2")) params.n2 = (int)val;
                    else if(label.equals("inclination_degree")) params.inclination_degree = val;
                    else if(label.equals("omega_degree")) params.omega_degree = val;
                    else if(label.equals("rmin")) params.rmin = val;
                    else if(label.equals("velocity_factor")) params.velocity_factor = val;
                    else if(label.equals("tstart")){ params.tstart = val; params.time = val; params.tIsSet = true;}
                    else if(label.equals("tend")) params.tend = val;
                    else if(label.equals("h")) params.h = val;
                    else if(label.equals("rx")) {params.sec_vec[0] = val;params.use_sec_vec=true;}
                    else if(label.equals("ry")) {params.sec_vec[1] = val;params.use_sec_vec=true;}
                    else if(label.equals("rz")) {params.sec_vec[2] = val;params.use_sec_vec=true;}
                    else if(label.equals("vx")) {params.sec_vec[3] = val;params.use_sec_vec=true;}
                    else if(label.equals("vy")) {params.sec_vec[4] = val;params.use_sec_vec=true;}
                    else if(label.equals("vz")) {params.sec_vec[5] = val;params.use_sec_vec=true;}
                    else if(label.equals("rscale1")) params.rscale1 = new double[]{val,val,val};
                    else if(label.equals("rscale2")) params.rscale2 = new double[]{val,val,val};
                    else if(label.equals("rscale11")) params.rscale1[0] = val;
                    else if(label.equals("rscale12")) params.rscale1[1] = val;
                    else if(label.equals("rscale13")) params.rscale1[2] = val;
                    else if(label.equals("rscale21")) params.rscale2[0] = val;
                    else if(label.equals("rscale22")) params.rscale2[1] = val;
                    else if(label.equals("rscale23")) params.rscale2[2] = val;
                    else System.out.println("Skipping line " + line);
                }

                line = lnr.readLine();
            }
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
        finally
        {
            if(fr != null) try{fr.close();}catch(Exception ex){};
        }
    }

    public static void writeParameterFile(Parameters params, String file)
    {
        FileWriter fw = null;
        try
        {
            fw = new FileWriter(file);

            fw.write("potential_type="+params.potential_type+"\n");
            fw.write("mass1="+params.mass1+"\n");
            fw.write("mass2="+params.mass2+"\n");
            fw.write("epsilon1="+params.epsilon1+"\n");
            fw.write("epsilon2="+params.epsilon2+"\n");
            fw.write("rin1="+params.rin1+"\n");
            fw.write("rin2="+params.rin2+"\n");
            fw.write("rout1="+params.rout1+"\n");
            fw.write("rout2="+params.rout2+"\n");
            fw.write("theta1="+params.theta1+"\n");
            fw.write("theta2="+params.theta2+"\n");
            fw.write("phi1="+params.phi1+"\n");
            fw.write("phi2="+params.phi2+"\n");
            fw.write("opt1="+params.opt1+"\n");
            fw.write("opt2="+params.opt2+"\n");
            fw.write("heat1="+params.heat1+"\n");
            fw.write("heat2="+params.heat2+"\n");
            fw.write("n1="+params.n1+"\n");
            fw.write("n2="+params.n2+"\n");
            fw.write("inclination_degree="+params.inclination_degree+"\n");
            fw.write("omega_degree="+params.omega_degree+"\n");
            fw.write("rmin="+params.rmin+"\n");
            fw.write("velocity_factor="+params.velocity_factor+"\n");
            fw.write("tstart="+params.tstart+"\n"); 
            fw.write("tend="+params.tend+"\n");
            fw.write("h="+params.h+"\n");
            fw.write("rx="+params.sec_vec[0]+"\n");
            fw.write("ry="+params.sec_vec[1]+"\n");
            fw.write("rz="+params.sec_vec[2]+"\n");
            fw.write("vx="+params.sec_vec[3]+"\n");
            fw.write("vy="+params.sec_vec[4]+"\n");
            fw.write("vz="+params.sec_vec[5]+"\n");
            fw.write("rscale11="+params.rscale1[0]+"\n");
            fw.write("rscale12="+params.rscale1[1]+"\n");
            fw.write("rscale13="+params.rscale1[2]+"\n");
            fw.write("rscale21="+params.rscale2[0]+"\n");
            fw.write("rscale22="+params.rscale2[1]+"\n");
            fw.write("rscale23="+params.rscale2[2]+"\n");
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
        finally
        {
            if(fw != null) try{fw.flush();fw.close();}catch(Exception ex){};
        }
    }

    public static void parseStateInfoString(Parameters params, String str)
    {
        double vals[] = new double[22];
        String sa[] = str.split(",");
        int len = Math.min(sa.length,vals.length);
        for(int i=0; i<len; i++)
        {
            vals[i] = parseDouble(sa[i]);
        }

        params.sec_vec[0] = vals[0];
        params.sec_vec[1] = vals[1];
        params.sec_vec[2] = vals[2];
        params.sec_vec[3] = vals[3];
        params.sec_vec[4] = vals[4];
        params.sec_vec[5] = vals[5];

        params.mass1 = vals[6];
        params.mass2 = vals[7];
        params.rout1 = vals[8];
        params.rout2 = vals[9];
        params.phi1 = vals[10];
        params.phi2 = vals[11];
        params.theta1 = vals[12];
        params.theta2 = vals[13];
        params.epsilon1 = vals[14];
        params.epsilon2 = vals[15];
        params.rscale1[0] = vals[16];
        params.rscale1[1] = vals[17];
        params.rscale1[2] = vals[18];
        params.rscale2[0] = vals[19];
        params.rscale2[1] = vals[20];
        params.rscale2[2] = vals[21];

        params.use_sec_vec = true;
    }
}
