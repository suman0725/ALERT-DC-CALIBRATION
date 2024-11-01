/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.t2d;

import java.math.RoundingMode;
import java.text.DecimalFormat;
//import static org.clas.detector.clas12calibration.dc.t2d.TableLoader.BfieldValues;
import static org.clas.detector.clas12calibration.dc.t2d.TableLoader.maxTBin;

public class TimeToDistanceEstimator {

    public TimeToDistanceEstimator() {
            // TODO Auto-generated constructor stub
    }
    
    /**
     * 
     * @param x value on grid
     * @param xa lower x bound on grid
     * @param xb upper x bound on grid
     * @param ya lower y bound on grid
     * @param yb upper y bound on grid
     * @return y value on grid from linear interpolation between a and b evaluated at x
     */
    private double interpolateLinear(double x0, double xa, double xb, double ya, double yb) {
        double x = x0;
        if(x>xb)
            x=xb;
        if(x<xa)
            x=xa;
        double y = (ya + yb)*0.5;
        if(xb - xa == 0) 
            return y;
        
        y = ya*(xb - x)/(xb - xa) + yb*(x - xa)/(xb - xa);
        
        return y;
    }

    /**
    * 
    * @param B B field in T
    * @param alpha is the local angle in degrees
    * @param t time in ns
    * @param SlyrIdx slyr index (0...5)
    * @return the distance to the wire in cm
    */
    public double interpolateOnGrid( double t,  int SecIdx, int SlyrIdx) {
        
        double f_B_alpha_t =0;
       /* interpolateLinear(Math.cos(Math.toRadians(30.))-Math.cos(Math.toRadians(alpha)), 
                    Math.cos(Math.toRadians(30.))-Math.cos(Math.toRadians(alpha1)), 
                    Math.cos(Math.toRadians(30.))-Math.cos(Math.toRadians(alpha2)), f_B_alpha1_t, f_B_alpha2_t);*/
        
        return f_B_alpha_t;

    }

    /**
     * 
     * @param binAlpha alpha parameter bin
     * @return value of alpha from alpha bin
     */
/*  We don't have alpha dependancy remove it -HA
 *
 *  private double getAlphaFromAlphaIdx(int binAlpha) {
        double cos30minusalpha = Math.cos(Math.toRadians(30.)) + (double) (binAlpha)*(1. - Math.cos(Math.toRadians(30.)))/5.;
        double alpha =  -(Math.toDegrees(Math.acos(cos30minusalpha)) - 30);
        double alpha1 = 0;
        double alpha2 = 30.;
        if(alpha<alpha1) {
            alpha=alpha1;
        }
         if(alpha>alpha2) {
             alpha=alpha2;
        }	
        return alpha;
    }*/

    /**
     * 
     * @param t1 time value in ns
     * @param is sector index (0...5)
     * @param ir superlayer index (0...5)
     * @param ibfield bfield bin (0...7)
     * @param icosalpha cosalpha bin (0...5)
     * @return time bin
     */
    //public int getTimeIdx(double t1, int is, int ir, int ibfield, int icosalpha) {
    public int getTimeIdx(double t1, int is, int ir) {
        
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);
       
        int binIdx =0;
        try{
            binIdx = Integer.parseInt(df.format(t1/2.) ) -1; 
        } catch (NumberFormatException e) {
            System.out.println(" time bin error "+t1+" ");
        }
        if(binIdx<0) {
            binIdx = TableLoader.minBinIdxT;
        }
        if(binIdx>maxTBin) {
            binIdx = maxTBin ;
        }

        return binIdx;
    }
    /**
     * 
     * @param b1 bfield value in T
     * @return B field bin
     */

/* No B dependency needed -HA

    public int getBIdx(double b1) {
        
//        int binIdx = (int) ((1+b1)*2) -2;
//        if(binIdx<0) {
//            binIdx = TableLoader.minBinIdxB;
//        }
//        if(binIdx>TableLoader.maxBinIdxB) {
//            binIdx = TableLoader.maxBinIdxB;
//        }
        int maxBinIdxB = TableLoader.BfieldValues.length-1;
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);
       
        int binIdx =0;
        try{
            binIdx = Integer.parseInt(df.format(b1*b1) ) -1; 
        } catch (NumberFormatException e) {
            System.out.println(" field bin error "+b1+" ");
        }
        if(binIdx<0) {
            binIdx = 0;
        }
        if(binIdx>maxBinIdxB)
            binIdx = maxBinIdxB;
        return binIdx;
    }

*/ 


    /**
     * 
     * @param alpha alpha parameter in deg
     * @return alpha bin
     */

 /* No ALpha dependency needed  - HA   
    private int getAlphaIdx(double alpha) {
        double Ccos30minusalpha = Math.cos(Math.toRadians(30.-alpha) ) ; 
        double Cicosalpha = (Ccos30minusalpha - Math.cos(Math.toRadians(30.)))/((1. - Math.cos(Math.toRadians(30.)))/5.);
        int binIdx = (int)  Cicosalpha; 
        if(binIdx<0) {
            binIdx = TableLoader.minBinIdxAlpha;
        }
        if(binIdx>TableLoader.maxBinIdxAlpha) {
            binIdx = TableLoader.maxBinIdxAlpha;
        } 
        return binIdx;
    }

    private int getTimeNextIdx(double t, int SecIdx, int SlyrIdx, int binlowB, int binlowAlpha) {
        int binlowT = this.getTimeIdx(t, SecIdx, SlyrIdx, binlowB, binlowAlpha);  
        int binhighT = binlowT + 1; 

        if(binhighT>TableLoader.maxBinIdxT[SecIdx][SlyrIdx][binlowB][binlowAlpha]) {
            binhighT=TableLoader.maxBinIdxT[SecIdx][SlyrIdx][binlowB][binlowAlpha];
        }
        return binhighT;
    }


    */
    
    /**
     * @param slyIdx superlayer index
     * @param time
     * @return test doca corr
     */
    public double addDOCACorr(double time, int slyIdx) {
        double dDoca = 0;
        if(slyIdx+1 == 5 || slyIdx+1 ==6) {
            if(time>600) {
                dDoca = 0.15;
            } else {
                dDoca = (7.6e-3 - 2.4e-4*time +9.8e-3*time*time - 3.8e-6*time*time*time)*5.5410595e-05;
            }
            //System.out.println("time "+time +" added doca "+(float)dDoca);
        }
        return dDoca;
    }
}
        


