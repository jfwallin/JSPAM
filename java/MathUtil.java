package edu.gmu.cds.util;

public class MathUtil
{
    public static final int X_AXIS = 1;
    public static final int Y_AXIS = 2;
    public static final int Z_AXIS = 3;
    
    /**
     * Builds a rotation matrix for rotating about the
     * specified axis by the specified angle in radians.
     * 
     * @param angle
     * @param axis
     * @return
     */
    public static double[][] buildRot(double angle, int axis)
    {
        double mat[][] = new double[3][3];
        double cosa = Math.cos(angle);
        double sina = Math.sin(angle);
        
        if(axis == X_AXIS)
        {
            mat[0][0]=1.0d;
            mat[1][1]=cosa;
            mat[1][2]=sina;
            mat[2][1]=-sina;
            mat[2][2]=cosa;
        }
        else if(axis == Y_AXIS)
        {
            mat[0][0]=cosa;
            mat[0][2]=-sina;
            mat[1][1]=1.0d;
            mat[2][0]=sina;
            mat[2][2]=cosa;
        }
        else if(axis == Z_AXIS)
        {
            mat[0][0]=cosa;
            mat[0][1]=sina;
            mat[1][0]=-sina;
            mat[1][1]=cosa;
            mat[2][2]=1.0d;
        }
        
        return mat;
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

    /**
     * Determine the angle of rotation between two vectors.
     * 
     * @param v1
     * @param v2
     * @return
     */
    public static double angleBetween(double[] v1, double[] v2)
    {
        double m1 = mag(v1);
        double m2 = mag(v2);
        double ang = dot(v1,v2)/(m1*m2);
        ang = Math.acos(ang);
        return ang;
    }
    
    /**
     * Multiplies a 3x3 matrix by a 3 vector.
     * It is mildly optimized by unrolling loops.
     * 
     * @param m
     * @param v
     * @return
     */
    public static double[] mult3(double m[][], double v[])
    {
        double b[] = new double[3];
        if(v.length < 3)
        {
            return v;
        }
        b[0] = m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];
        b[1] = m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2];
        b[2] = m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2];
        
        return b;
    }
    
    /**
     * Multiplies an N vector by an MxN matrix.
     * 
     * @param m
     * @param v
     * @return
     */
    public static double[] mult(double m[][], double v[])
    {
        int col = m[0].length;
        int row = m.length;
        
        double b[] = new double[row];
        
        for(int i=0; i<row; i++)
        {
            for(int j=0; j<col; j++)
            {
                b[i]+=m[i][j]*v[j];
            }
        }
        
        return b;
    }
    
    /**
     * Multiply the two quaternions.
     * 
     * @param q1
     * @param q2
     * @return
     */
    public static double[] quatMult(double q1[], double q2[])
    {
        double q3[] = new double[4];
    
        q3[0] = q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
        q3[1] = q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2];
        q3[2] = q1[0]*q2[2]+q1[2]*q2[0]+q1[3]*q2[1]-q1[1]*q2[3];
        q3[3] = q1[0]*q2[3]+q1[3]*q2[0]+q1[1]*q2[2]-q1[2]*q2[1];
        
        return q3;
    }
    
    /**
     * Invert the quaternion.
     * 
     * @param q
     * @return
     */
    public static double[] quatInv(double q[])
    {
        double qInv[] = new double[4];
        
        double dot = dot(q,q);
        double scale = -1.0/dot;
        
        qInv[0] = -q[0]*scale;
        qInv[1] = q[1]*scale;
        qInv[2] = q[2]*scale;
        qInv[3] = q[3]*scale;
        
        return qInv;
    }
    
    /**
     * Build a rotation quaternion for the given angle
     * and axis.  This method will normalize the vector.
     * 
     * @param angle in radians
     * @param v
     * @return
     */
    public static double[] quatRot(double angle, double v[])
    {
        return quatRot(angle,v[0],v[1],v[2]);
    }
    
    /**
     * Build a rotation quaternion for the given angle
     * and axis.  This method will normalize the vector.
     * 
     * @param angle in radians
     * @param x
     * @param y
     * @param z
     * @return
     */
    public static double[] quatRot(double angle, double x, double y, double z)
    {
        double q[] = new double[4];
        
        double mag = Math.sqrt(x*x+y*y+z*z);
        double sinT = Math.sin(angle*0.5);
        double cosT = Math.cos(angle*0.5);
        
        double scale = sinT/mag;
        
        q[0] = cosT;
        q[1] = x*scale;
        q[2] = y*scale;
        q[3] = z*scale;
        
        return q;
    }
    
    /**
     * Convert the quaternion to a rotation matrix.
     * 
     * @param quat
     * @return
     */
    public static double[][] quatToMatrix(double quat[])
    {
        double mat[][] = new double[3][3];
        if(quat == null)
        {
            mat[0][0] = 1.0d;
            mat[1][1] = 1.0d;
            mat[2][2] = 1.0d;
            return mat;
        }
        double q0=quat[0];
        double q1=quat[1];
        double q2=quat[2];
        double q3=quat[3];
        
        mat[0][0] = 1.0-2.0*(q2*q2+q3*q3);
        mat[0][1] = 2.0*(q1*q2+q0*q3);
        mat[0][2] = 2.0*(q1*q3-q0*q2);
        mat[1][0] = 2.0*(q1*q2-q0*q3);
        mat[1][1] = 1.0-2.0*(q1*q1+q3*q3);
        mat[1][2] = 2.0*(q2*q3+q0*q1);
        mat[2][0] = 2.0*(q1*q3+q0*q2);
        mat[2][1] = 2.0*(q2*q3-q0*q1);
        mat[2][2] = 1.0-2.0*(q1*q1+q2*q2);

        return mat;
    }
}
