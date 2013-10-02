package edu.gmu.cds.gui;

import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

public interface QuaternionRotator
{
    public void addQRot(double q[]);
    
    public void clearQRot();
 
    public double[] getQRot();
    
    public int getHeight();
    
    public int getWidth();
    
    public void addMouseListener(MouseListener ml);
    
    public void addMouseMotionListener(MouseMotionListener ml);
    
    public void repaint();
}
