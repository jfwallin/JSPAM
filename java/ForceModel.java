package edu.gmu.cds.sim;

/**
 * This interface defines the methods required for defining a force model.
 *
 */
public interface ForceModel
{
    /**
     * Sets the simulation parameters to this model.
     */
    public void setParameters(Parameters p);

    /**
     * Initialize temporary storage.
     */
    public void initVars(int n);

    /**
     * Cleanup temporary storage.
     */
    public void deallocateVars();

    /**
     * Compute the circular velocity for a particle at a distance r from the specified mass.
     * The rout scale of the disk and softening length, eps, are provided.
     */
    public double circularVelocity(double mass, double r, double rout, double eps);

    /**
     * For the given particle positions and velocities, calculate the accelerations.
     */
    public void diffeq(double x[][], double f[][]);    

    /**
     * Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
     */
    public void diffq1(double x[], double f[]);
}
