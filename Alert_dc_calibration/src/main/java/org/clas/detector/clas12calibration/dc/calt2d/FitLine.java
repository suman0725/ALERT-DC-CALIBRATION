/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.math.Func1D;
import org.jlab.rec.dc.Constants;
/**
 *
 * @author ziegler
 * @author Hamza
 */
public class FitLine extends Func1D{
    public int i;
    private FitFunction fc ;
    public FitLine() {
        super("fcn", 0.0, 2.0);
        fc = new FitFunction();
    }
    public static final int nPars = 4;
    private double[] par = new double[nPars];
    public FitLine(String name, int i, MnUserParameters pars) {
        super(name, 0.0, 2.0);
        this.i = i;
        fc = new FitFunction();
        this.initParameters(pars);
    }

    private void initParameters(MnUserParameters pars) {
        for(int p = 0; p< nPars; p++) {
            par[p] = pars.value(p);
        }
    }
    @Override
    public double evaluate(double x) { 
        double calcTime = 0;
        double p1 = par[0]; 
        double p2 = par[1]; 
        double p3 = par[2]; 
        double p4 = par[3]; 

        calcTime = fc.polyFcnMac(x, p1,  p2,  p3,  p4, i+1);
        
        //System.out.println("ijk "+i+""+j+""+k+" b "+(float)T2DCalib.BfieldValues[k]+" ralpha "+(float)ralpha+" x "+x+" time "+(float)calcTime);
        return calcTime;
    }

    
}
