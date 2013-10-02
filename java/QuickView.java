package edu.gmu.cds.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import java.util.ArrayList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.gmu.cds.sim.IOUtil;

public class QuickView extends JPanel implements ChangeListener
{
    protected ScatterPanel sp = null;
    protected RotationController rc = null;
    protected JTextField txtAlpha = new JTextField(6);
    protected JTextField txtBeta = new JTextField(6);
    protected JSlider alphaSlider = new JSlider(JSlider.HORIZONTAL,0,360,0);
    protected JSlider betaSlider = new JSlider(JSlider.HORIZONTAL,0,360,0);

    protected double alpha = 0;
    protected double beta = 0;

    protected double[][] originalParticles = null;
    protected List frames = null;

    public QuickView()
    {
        super(new BorderLayout());

        sp = new ScatterPanel();
        rc = new RotationController();
        rc.setController(sp);

        add(sp,BorderLayout.CENTER);

        setMinSize(alphaSlider);
        setMinSize(betaSlider);
        setMinSize(txtAlpha);
        setMinSize(txtBeta);

        JPanel txtPanel = new JPanel(new GridBagLayout());

        Insets left = new Insets(2,5,2,2);
        Insets right = new Insets(2,2,2,5);
        GridBagConstraints gbc = null;
        int row = 0;

        // angles
        row++;
        gbc = getGBC(0,row,GridBagConstraints.EAST,left);
        txtPanel.add(new JLabel("Alpha"),gbc);
        
        gbc = getGBC(1,row,GridBagConstraints.WEST,right);
        txtPanel.add(txtAlpha,gbc);
        
        gbc = getGBC(2,row,GridBagConstraints.CENTER,right);
        gbc.gridwidth=5;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        txtPanel.add(alphaSlider,gbc);

        row++;
        gbc = getGBC(0,row,GridBagConstraints.EAST,left);
        txtPanel.add(new JLabel("Beta"),gbc);
        
        gbc = getGBC(1,row,GridBagConstraints.WEST,right);
        txtPanel.add(txtBeta,gbc);
        
        gbc = getGBC(2,row,GridBagConstraints.CENTER,right);
        gbc.gridwidth=5;
        gbc.fill = GridBagConstraints.HORIZONTAL;
        txtPanel.add(betaSlider,gbc);
        
        alphaSlider.addChangeListener(this);
        betaSlider.addChangeListener(this);

        add(txtPanel,BorderLayout.SOUTH);
    }

    public void setAlpha(double alpha)
    {
        this.alpha = alpha;
        alphaSlider.setValue((int)alpha);
        updateText();
    }

    public void setBeta(double beta)
    {
        this.beta = beta;
        betaSlider.setValue((int)beta);
        updateText();
    }

    /**
     * Force the minimum size to equal the preferred size.
     * 
     * @param c
     */
    protected void setMinSize(JComponent c)
    {
        c.setMinimumSize(c.getPreferredSize());    
    }
    
    public GridBagConstraints getGBC(int x, int y, int anchor, Insets insets)
    {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = x;
        gbc.gridy = y;
        gbc.anchor = anchor;
        gbc.fill = GridBagConstraints.NONE;
        gbc.gridheight = 1;
        gbc.gridwidth = 1;
        gbc.weightx = 1;
        gbc.weighty = 1;
        gbc.insets = insets;    
        
        return gbc;
    }

    public void updateText()
    {
        txtAlpha.setText(String.valueOf(alpha));
        txtBeta.setText(String.valueOf(beta));            
    }

    public double[][] applyViewingAngles(double alphaDeg, double betaDeg, double xin[][])
    { 
        double x,y,z,xv,yv,zv;
        double x2,y2,z2,xv2,yv2,zv2;
        double x3,y3,z3,xv3,yv3,zv3;
        double ctheta,stheta,cphi,sphi;

        stheta = Math.sin(Math.toRadians(alphaDeg));
        ctheta = Math.cos(Math.toRadians(alphaDeg));
        sphi = Math.sin(Math.toRadians(betaDeg));
        cphi = Math.cos(Math.toRadians(betaDeg));
        double ca = ctheta;
        double sa = stheta;
        double cb = cphi;
        double sb = sphi;
        
        int size = xin.length;
        double xout[][] = new double[size][6];

        for(int i=0; i<size; i++)
        {
            x=xin[i][0];
            y=xin[i][1];
            z=xin[i][2];

            xv=xin[i][3];
            yv=xin[i][4];
            zv=xin[i][5];

            x2 = x*cb - y*sb;
            y2 = x*sb + y*cb;
            z2 = z;
            x3 = x2*ca - z2*sa;
            y3 = -y2;
            z3 = -x2*sa-z2*ca;

/*
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
*/

            xout[i][0] = x3;
            xout[i][1] = y3;
            xout[i][2] = z3;
            //xout[i][3] = xv3;
            //xout[i][4] = yv3;
            //xout[i][5] = zv3;
        }

System.out.println(xout[size-1][0]+"\t"+xout[size-1][1]+"\t"+xout[size-1][2]+"\t"+xout[size-1][3]+"\t"+xout[size-1][4]+"\t"+xout[size-1][5]);

        return xout;
    }

    public double[][] applyOldStyleViewingAngles(double alphaDeg, double betaDeg, double x[][])
    { 
        double sa = Math.sin(Math.toRadians(alphaDeg));
        double ca = Math.cos(Math.toRadians(alphaDeg));
        double sb = Math.sin(Math.toRadians(betaDeg));
        double cb = Math.cos(Math.toRadians(betaDeg));
        double cg = 1.0;
        double sg = 0.0;

        double a11 = cg*ca-cb*sa*sg;
        double a12 = cg*sa+cb*ca*sg;
        double a13 = sg*sb;
        double a21 = -sg*ca-cb*sa*cg;
        double a22 = -sg*sa+cb*ca*cg;
        double a23 = cg*sb;

        double xp,yp,zt;

        int size = x.length;
        double xout[][] = new double[size][3];

        for(int i=0; i<size; i++)
        {
            xp=x[i][0];
            yp=x[i][1];
            zt=x[i][2];
            xout[i][0] = a11*xp+a12*yp+a13*zt;
            xout[i][1] = a21*xp+a22*yp+a23*zt;
        }

        return xout;
    }

    public void addSeries(double x[][])
    {
        sp.addSeries(x,Color.RED);
        sp.refreshDataRange();
    }
    /**
     * Handle slider events
     */
    public void stateChanged(ChangeEvent e)
    {
        Object obj = e.getSource();

        if(obj == alphaSlider)
        {
            alpha = alphaSlider.getValue();
        }
        else if(obj == betaSlider)
        {
            beta = betaSlider.getValue();
        }
        updateText();

        double particles[][] = applyViewingAngles(alpha,beta,originalParticles);
        sp.clear();
        addSeries(particles);

        repaint();
    }


    public static void usage()
    {
        System.out.println("\n\tSpecify a particle file.\n");
    }

    public static void main(String args[])
    {
        if(args.length < 1)
        {
            usage();
            System.exit(-1);
        }

        QuickView qv = new QuickView();

        double particles[][] = IOUtil.inputParticles(args[0]);
        qv.originalParticles = particles;
        if(args.length > 2)
        {
	    qv.setAlpha(Integer.parseInt(args[1]));
            qv.setBeta(Integer.parseInt(args[2]));
            particles = qv.applyViewingAngles(qv.alpha,qv.beta,particles);
        }
        qv.addSeries(particles);

        JFrame frame = new JFrame("Quick View");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(qv);
        frame.setSize(800,600);
        frame.setVisible(true);
    }
}
