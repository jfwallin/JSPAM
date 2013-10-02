/*
  Copyright 2008 by Anthony Holincheck
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/
package edu.gmu.cds.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import edu.gmu.cds.gui.QuaternionRotator;
import edu.gmu.cds.util.MathUtil;

/**
 * Used to display multiple series of x,y points.
 * 
 * @author Anthony Holincheck
 *
 */
public class ScatterPanel extends JPanel implements QuaternionRotator
{
    private static final long serialVersionUID = 1L;
    
    public static final Color MY_RED = new Color(255,0,0,100);
    public static final Color MY_BLUE = new Color(0,0,255,100);
    
    public List series = null;
    public List colors = null;
    public List maxInd = null;
    public List symbolSize = null;

    public Color background = null;

    public double minx = 0;
    public double miny = 0;
    public double maxx = 0;
    public double maxy = 0;

    public boolean bDoingPaint = false;
    
    public Color timeColor = Color.YELLOW;
    public boolean bDrawTime = false;
    public double timeVal = 0.0d;
    
    public String messageText;
    public boolean bDisplayMessage = false;
    public Color messageColor = Color.BLUE;
    
    protected double qrot[] = null;
    
    protected int defaultSymbolSize = 1;
    
    protected boolean showStereo = false;
    protected double eyeSep = 0.5;
    protected double z0 = 10;
    
    protected double zoomScale = 1.0;
    
    protected PaintHandler paintHandler = null;

    protected boolean paintAxes = true;
    
    public ScatterPanel()
    {
        super();
        series = new ArrayList();
        maxInd = new ArrayList();
        colors = new ArrayList();
        symbolSize = new ArrayList();
        background = Color.BLACK;
    }
    
    public void setPaintHandler(PaintHandler ph)
    {
        paintHandler = ph;
    }
    
    public void setDefaultSymbolSize(int size)
    {
        defaultSymbolSize = size;
    }
    
    /**
     * Add the new rotation quaternion to the previous rotation.
     * 
     * @param q
     */
    public void addQRot(double q[])
    {
        if(qrot == null)
        {
            qrot = q;
        }
        else
        {
            qrot = MathUtil.quatMult(qrot,q);
        }
    }
    
    /**
     * Unset the rotation.
     */
    public void clearQRot()
    {
        qrot = null;
    }
    
    /**
     * Get the rotation quaternion.
     * 
     * @return
     */
    public double[] getQRot()
    {
        return qrot;
    }
    
    /**
     * Point to the same series and colors as the
     * designated panel.
     */
    public void shareSeries(ScatterPanel sp)
    {
        clear();
        int size = sp.series.size();
        for(int i=0; i<size; i++)
        {
            series.add(sp.series.get(i));
            addColor(sp.getColor(i));
        }
    }
    
    public Object clone()
    {
        ScatterPanel sp = new ScatterPanel();
        sp.setDataRange(this.minx,this.maxx,this.miny,this.maxy);
        sp.series.addAll(this.series);
        sp.maxInd.addAll(this.maxInd);
        sp.colors.addAll(this.colors);
        sp.symbolSize.addAll(this.symbolSize);
        sp.background = this.background;
        sp.bDoingPaint = this.bDoingPaint;
        sp.timeColor = this.timeColor;
        sp.bDrawTime = this.bDrawTime;
        sp.timeVal = this.timeVal;
        sp.messageText = this.messageText;
        sp.bDisplayMessage = this.bDisplayMessage;
        sp.messageColor = this.messageColor;
        
        return sp;
    }
   
    public void refreshDataRange()
    {
        if(series.size() > 0)
        {
            double data[][] = (double[][])series.get(0);

            int size = data.length;
            double min1 = Double.MAX_VALUE;
            double min2 = Double.MAX_VALUE;
            double min3 = Double.MAX_VALUE;
            double max1 = -(min1-1);
            double max2 = -(min2-1);
            double max3 = -(min3-1);
            for(int i=0; i<size; i++)
            {
                if(data[i][0] <  min1)
                {
                    min1 = data[i][0];
                }
                if(data[i][0] >  max1)
                {
                    max1 = data[i][0];
                }
                if(data[i][1] <  min2)
                {
                    min2 = data[i][1];
                }
                if(data[i][1] >  max2)
                {
                    max2 = data[i][1];
                }
                if(data[i].length > 2)
                {
                    if(data[i][2] <  min3)
                    {
                        min3 = data[i][2];
                    }
                    if(data[i][2] >  max3)
                    {
                        max3 = data[i][2];
                    }
                }
            }

            //setDataRange(min1,max1,min2,max2);

            min1 = Math.min(min1,min2); 
            min1 = Math.min(min1,min3);
            max1 = Math.max(max1,max2); 
            max1 = Math.max(max1,max3);

            if(Math.abs(min1) > Math.abs(max1))
            {
                max1 = -min1;
            }
            else
            {
                min1 = -max1;
            }

            setDataRange(min1,max1,min1,max1);
        }
    }
 
    public void setDataRange(double minx, double maxx, double miny, double maxy)
    {
        this.minx = minx;
        this.maxx = maxx;
        this.miny = miny;
        this.maxy = maxy;
    }
    
    public void setMessage(String txt)
    {
        messageText = txt;
    }
    
    public void setDisplayMessage(boolean flag)
    {
        bDisplayMessage = flag;
    }
    
    public void setTimeValue(double t)
    {
        timeVal = t;
    }
    
    public void setDrawTime(boolean flag)
    {
        bDrawTime = flag;
    }
    
    public void clear()
    {
        series.clear();
        colors.clear();
        maxInd.clear();
        symbolSize.clear();
    }
    
    public void setSymbolSize(int series, int size)
    {
        symbolSize.set(series,new Integer(size));
    }
    
    public void addSeries(double data[][], Color color)
    {
        addSeries(data, color, false);
    }
    
    public void addSeries(double data[][], Color color, boolean doRepaint)
    {
        series.add(data);
        symbolSize.add(new Integer(defaultSymbolSize));
        maxInd.add(new Integer(data.length));
        addColor(color);
        
        if(doRepaint)
        {
            repaint();
        }
    }
    
    protected void addColor(Color color)
    {
        colors.add(color);
    }
    
    protected Color getColor(int ind)
    {
        return (Color)colors.get(ind);
    }
    
    public void setMaxInd(int seriesNum, int maxIndex)
    {
        if(seriesNum >=0 && seriesNum < maxInd.size())
        {
            maxInd.set(seriesNum,new Integer(maxIndex));
        }
    }
    
    public void setZoomScale(double zoom)
    {
        zoomScale = zoom;
    }
    
    public double getZoomScale()
    {
        return zoomScale;
    }
    
    public void paintComponent(Graphics gr)
    //public void paint(Graphics g)
    {
        if(bDoingPaint)
        {
            return;
        }
        
        bDoingPaint = true;
        
        Dimension dim = this.getSize();
        int width = (int)dim.getWidth();
        int height = (int)dim.getHeight();
        
        if(showStereo)
        {
            paintStereo(gr,width,height);
            bDoingPaint = false;
            return;
        }
        
        Graphics g = gr.create();
        
        double v[] = null;
        double rot[][] = null;
        rot = MathUtil.quatToMatrix(qrot);
        
        g.setColor(background);
        g.fillRect(0,0,width,height);
        
        // If we need to do a message, do it and return
        if(bDisplayMessage && messageText != null)
        {
            int y = height/2+10;
            int x = width/2-messageText.length()*4;
            g.setColor(messageColor);
            g.drawString(messageText,x,y);
            g.dispose();
            bDoingPaint = false;
            return;
        }
    
        double rngx = maxx-minx;
        double rngy = maxy-miny;
        
        rngx *= zoomScale;
        rngy *= zoomScale;

        if(paintAxes)
        {
            int sz = 101;
            double num = 0;
            double d[][] = new double[3*sz][3];
            for(int i=0; i<sz; i++)
            {
                num = i * (maxx/((double)(sz-1)) );
                d[i][0] = num;
                d[i][1] = 0;
                d[i][2] = 0;
                d[i+sz][0] = 0;
                d[i+sz][1] = num;
                d[i+sz][2] = 0;
                d[i+2*sz][0] = 0;
                d[i+2*sz][1] = 0;
                d[i+2*sz][2] = num;
            }
            series.add(d);
            colors.add(Color.yellow);
        }

        int size = series.size();
        
        for(int i=0; i<size; i++)
        {
            double data[][] = (double[][])series.get(i);
            g.setColor((Color)colors.get(i));
            
            int np = data.length;
            try
            {
                Integer mi = (Integer)maxInd.get(i);
                if(mi != null)
                {
                    np = mi.intValue();
                }
            }
            catch(Exception ex)
            {
                
            }
            
            int ss = 1;
            int soff = 0;
            
            try
            {
                int symsize = ((Integer)symbolSize.get(i)).intValue();
                if(symsize > 1)
                {
                    soff = symsize/2;
                    ss = symsize;
                }
            }
            catch(Exception ex)
            {
                
            }
            
            double x = 0;
            double y = 0;
            int px = 0;
            int py = 0;
            
            for(int j=0; j<np; j++)
            {
                //x = data[j][0];
               // y = data[j][1];
                   v = MathUtil.mult3(rot,data[j]);
                x = v[0];
                y = v[1];
                
                x = (x - minx*zoomScale)/rngx;
                y = (y - miny*zoomScale)/rngy;
            
                px = (int)(width*x);
                py = (int)(height*(1-y));
                
                if(rngy == 0)
                {
                    py = height/2;
                }
                
                if(rngx == 0)
                {
                    px = width/2;
                }
                
                if(px>=0 && py >=0)
                {
                    g.fillRect(px-soff,py-soff,ss,ss);
                }
            }
        }
        
        if(bDrawTime)
        {
            String timeStr = "t = "+String.valueOf(timeVal);
            if(timeStr.length() > 10)
            {
                timeStr = timeStr.substring(0,10);
            }
            
            g.setColor(timeColor);
            g.drawString(timeStr,10,20);
        }

        if(paintAxes)
        {
            series.remove(series.size()-1);
            colors.remove(colors.size()-1);
        }
    
        g.dispose();
        
        bDoingPaint = false;
    }
    
    /**
     * Inverts the data colors.
     */
    public void invertColors()
    {
        int size = colors.size();
    
        Color color = null;
        int r,g,b;
        for(int i=0; i<size; i++)
        {
            color = (Color)colors.get(i);
            r = 255-color.getRed();
            g = 255-color.getGreen();
            b = 255-color.getBlue();
            color = new Color(r,g,b);
            colors.set(i,color);
        }
    }
    
    public void paintStereo(Graphics g, int width, int height)
    {
        g.setColor(background);
        g.fillRect(0,0,width,height);
        int size = series.size();
        
        double rngx = maxx-minx;
        double rngy = maxy-miny;
        
        double rngxinv = 1.0d/rngx;
        double rngyinv = 1.0d/rngy;
        double zscale = Math.min(width,height);
        
        double v[] = null;
        double rot[][] = null;
        
        rot = MathUtil.quatToMatrix(qrot);
        
        double x = 0;
        double y = 0;
        double z = 0;
        double xproj = 0;
        double yproj = 0;
        int px = 0;
        int py = 0;
        double dx = 0;
        double dy = 0;
        double dz = 0;
        double den = 0;
        
        eyeSep = 0.5*(maxx-minx)/3.0;
        z0 = 10*(maxx-minx)/3.0;
        
        for(int i=0; i<size; i++)
        {
            double data[][] = (double[][])series.get(i);
            
            int np = data.length;
            int ss = ((Integer)symbolSize.get(i)).intValue();
            double d = 0;
            
            
            d = eyeSep;
            g.setColor(MY_RED);
            for(int j=1; j<np; j++)
            {
                x = data[j][0];
                y = data[j][1];
                
                v = MathUtil.mult3(rot,data[j]);
                x = v[0];
                y = v[1];
                z = v[2];
                
                // project coordinates
                den = 1.0d/(z0-z);
                xproj = (z0 * x - z * d)*den;
                yproj = (z0 * y)*den;
                
                x = (xproj - minx)*rngxinv;
                y = (yproj - miny)*rngyinv;
                
                dx = width*x;
                dy = height*(1.0d-y);
                                
                px = (int)dx;
                py = (int)dy;
                g.fillRect(px,py,ss,ss);   
            }
            
            d = -eyeSep;
            g.setColor(MY_BLUE);
            for(int j=1; j<np; j++)
            {
                x = data[j][0];
                y = data[j][1];
                
                v = MathUtil.mult3(rot,data[j]);
                x = v[0];
                y = v[1];
                z = v[2];
                
                // project coordinates
                den = 1.0d/(z0-z);
                xproj = (z0 * x - z * d)*den;
                yproj = (z0 * y)*den;
                
                x = (xproj - minx)*rngxinv;
                y = (yproj - miny)*rngyinv;
                
                dx = width*x;
                dy = height*(1.0d-y);
                                
                px = (int)dx;
                py = (int)dy;
                g.fillRect(px,py,ss,ss);   
            }            
        }
    }
    
    public boolean getShowStereo()
    {
        return showStereo;
    }
    
    public void setShowStereo(boolean flag)
    {
        showStereo = flag;
    }
    
    public void repaint()
    {
        if(paintHandler == null)
        {
            super.repaint();
        }
        else
        {
            paintHandler.repaint(this);
        }
    }

    public static interface PaintHandler
    {
        public void repaint(ScatterPanel sp);
    }
    
    public static void main(String args[])
    {
        ScatterPanel sp = new ScatterPanel();
        
        sp.minx = -3;
        sp.maxx = 3;
        sp.miny = -3;
        sp.maxy = 3;
        
        double data[][] = {{-3,-2.7,0},
        {-2.8,-2.1952,0},
        {-2.6,-1.7576,0},
        {-2.4,-1.3824,0},
        {-2.2,-1.0648,0},
        {-2,-0.8,0},
        {-1.8,-0.5832,0},
        {-1.6,-0.4096,0},
        {-1.4,-0.2744,0},
        {-1.2,-0.1728,0},
        {-1,-0.1,0},
        {-0.8,-0.0512,0},
        {-0.6,-0.0216,0},
        {-0.4,-0.0064,0},
        {-0.2,-0.0008,0},
        {0,0,0},
        {0.2,0.0008,0},
        {0.4,0.0064,0},
        {0.6,0.0216,0},
        {0.8,0.0512,0},
        {1,0.1,0},
        {1.2,0.1728,0},
        {1.4,0.2744,0},
        {1.6,0.4096,0},
        {1.8,0.5832,0},
        {2,0.8,0},
        {2.2,1.0648,0},
        {2.4,1.3824,0},
        {2.6,1.7576,0},
        {2.8,2.1952,0},
        {3,2.7,0}};

        sp.addSeries(data,Color.GREEN,true);

        RotationController rc = new RotationController();
        rc.setController(sp);
        JFrame frame = new JFrame("Test Plot");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(sp);
        frame.setSize(800,600);
        frame.setVisible(true);

    }
}
