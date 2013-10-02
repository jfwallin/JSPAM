package edu.gmu.cds.sim;

public class Run
{
    public Parameters params;
    public ForceModel forceModel;
    public Integrator integrator;
    public double x0[][] = null;
    public double xout[][] = null;

    public Run()
    {
        params = new Parameters();
        forceModel = new SPMModel();
        integrator = new Integrator(forceModel);
    }

    public void initRun(String args[])
    {
//Prof.startCall("initRun(args[])");
        SetupUtil su = new SetupUtil();

        su.setHelpers(forceModel, integrator, params);
        su.customCollision(args);
        // update the forceModel based upon passed in args
        forceModel = su.createForceModel(params.potential_type,true);

        forceModel.setParameters(params);
        integrator.initRKVar(params.n);

        x0  = new double[params.n+1][6];
        xout = new double[params.n+1][6];

        su.createCollision();

        copyParticles(integrator.x,x0);

        forceModel.setParameters(params);
//Prof.endCall("initRun(args[])");
    }

    public void initRun(Parameters params)
    {
        SetupUtil su = new SetupUtil();

        su.setHelpers(forceModel, integrator, params);
        // update the forceModel based upon passed in args
        forceModel = su.createForceModel(params.potential_type,true);

        forceModel.setParameters(params);
        integrator.initRKVar(params.n);

        x0  = new double[params.n+1][6];
        xout = new double[params.n+1][6];

        su.createCollision();

        copyParticles(integrator.x,x0);

        forceModel.setParameters(params);
    }

    public double[] getMinValues(Parameters params)
    {
        SetupUtil su = new SetupUtil();

        su.setHelpers(forceModel, integrator, params);
        // update the forceModel based upon passed in args
        forceModel = su.createForceModel(params.potential_type,true);

        forceModel.setParameters(params);
        integrator.initRKVar(params.n);

        x0  = new double[params.n+1][6];
        xout = new double[params.n+1][6];

        double mins[] = null;
        mins = su.perturberPositionVec(params.sec_vec, params.mass1, params.mass2, 
                                       params.eps1, params.eps2,
                                       params.h, params.n, params.n1, params.time, integrator.x);

        return mins;
    }

    /**
     * Copies the particles from x1 to x2.
     */
    public void copyParticles(double x1[][], double x2[][])
    {
//Prof.startCall("copyParticles");
        int n = x1.length;
        // assuming 6 members

        for(int i=0; i<n; i++)
        {
            x2[i][0] = x1[i][0];
            x2[i][1] = x1[i][1];
            x2[i][2] = x1[i][2];
            x2[i][3] = x1[i][3];
            x2[i][4] = x1[i][4];
            x2[i][5] = x1[i][5];
        }
//Prof.endCall("copyParticles");
    }

    public void takeAStep()
    {
//Prof.startCall("takeAStep");
        integrator.rk4(x0,xout,params.h);
        copyParticles(xout,x0);
        params.time = params.time + params.h;
//Prof.endCall("takeAStep");
    }

    public String getFilename(int i)
    {
        String str = String.valueOf(i);
        while(str.length() < 3)
        {
            str = "0" + str;
        }
        return "a." + str;
    }

    /**
     * Run the simulation with current parameters from tstart to tend
     */
    public void calculate(double tstart, double tend)
    {
        double t0 = 0;
        double time_interval = 0;

        int nstep_local = 7500;
        params.tstart = tstart;
        params.tend = tend;
        t0 = params.tstart;
        params.nstep = (int)((params.tend-t0)/params.h)+2;
        nstep_local = params.nstep;
        time_interval = (params.tend-t0)*2;
        for(int i=0; i< nstep_local; i++)
        {
            takeAStep();
        }
    }

    public static void main(String args[])
    {
//Prof.startCall("main");
        Run run = new Run();
        run.initRun(args);
        Parameters params = run.params;

        double t0 = 0;
        double time_interval = 0;

        int nstep_local = 7500;

        t0 = params.tstart;
        params.nstep = (int)((params.tend-t0)/params.h)+2;
        nstep_local = params.nstep;
        time_interval = (params.tend-t0)*2;
IOUtil.writeParameterFile(params,"tmp.p");
        for(int i=1; i<= nstep_local; i++)
        {
            run.takeAStep();
            if(i % 10 == 5)
            { 
                run.params.iout++;
                System.out.println(run.params.iout);
                //IOUtil.outputParticles(run.getFilename(run.params.iout), run.integrator.x);
            }
        } 

        IOUtil.outputParticles(run.getFilename(run.params.iout), run.integrator.x);
//Prof.endCall("main");
//Prof.report();
        System.exit(0);
    }
}
