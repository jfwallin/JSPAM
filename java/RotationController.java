package edu.gmu.cds.gui;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import edu.gmu.cds.util.MathUtil;

public class RotationController implements MouseListener, MouseMotionListener
{
    protected QuaternionRotator qr = null;
    
    protected boolean leftDown = false;
    protected double initV[] = null;
    
    public RotationController()
    {
        
    }
    
    public void clearRotation()
    {
        if(qr != null)
        {
            qr.clearQRot();
        }
    }
    
    public void setController(QuaternionRotator qr)
    {
        this.qr = qr;
        if(qr != null)
        {
            qr.addMouseListener(this);
            qr.addMouseMotionListener(this);
        }
    }
    
    /**
     * Take mouse x and y and get vector for
     * trackball location.
     * 
     * @param x
     * @param y
     * @return
     */
    public double[] getV(double x, double y)
    {
        int w = qr.getWidth();
        int h = qr.getHeight();

        x = x/(w*0.5);
        y = y/(h*0.5);
        
        x = x-1;
        y = 1-y;
        
        double z2 = 1-x*x-y*y;
        double z = 0;
        double scale = 1;

        if(z2 > 0)
        {
            z = Math.sqrt(z2);
        }
        scale = 1.0/Math.sqrt(x*x+y*y+z*z);

        double v[] = {x*scale,y*scale,z*scale};
        
        return v;
    }
    
    public void mouseDragged(MouseEvent e)
    {
        if(!leftDown)
        {
            return;
        }
        
        int x = e.getX();
        int y = e.getY();

        updateRotation(x,y);
    }
    
    public synchronized void updateRotation(int x, int y)
    {
        if(qr == null)
        {
            return;
        }

        double newV[] = getV(x,y);
        double v[] = MathUtil.cross(initV,newV);
        double angle = MathUtil.angleBetween(initV,newV);

        if(angle > 0.0d && !Double.isNaN(angle))
        {
            double qrot[] = MathUtil.quatRot(angle,v);
        
            qr.addQRot(qrot);
        
            initV = newV;
        
            qr.repaint();
        }
    }

    public void mouseMoved(MouseEvent e)
    {
        // do nothing
    }

    public void mouseClicked(MouseEvent e)
    {
        // do nothing
    }

    public void mousePressed(MouseEvent e)
    {
        if(e.getButton() == MouseEvent.BUTTON1)
        {
            leftDown = true;
            initV = getV(e.getX(),e.getY());
        }
    }

    public void mouseReleased(MouseEvent e)
    {
        if(e.getButton() == MouseEvent.BUTTON1)
        {
            leftDown = false;
        }
    }

    public void mouseEntered(MouseEvent e)
    {
        // do nothing
    }

    public void mouseExited(MouseEvent e)
    {
        // do nothing
    }

}
