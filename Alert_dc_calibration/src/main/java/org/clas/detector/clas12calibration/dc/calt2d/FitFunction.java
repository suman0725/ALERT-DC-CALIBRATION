/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.freehep.math.minuit.FCNBase;
import org.jlab.groot.data.GraphErrors;
import org.jlab.rec.dc.Constants;

/**
 *
 * @author ziegler
 * @author Hamza
 */
public class FitFunction implements FCNBase{

  //  public double _beta = 1.0;
    
    private Map<Coordinate, GraphErrors> _tvstrkdocasProf;
    private int i;
    
    public FitFunction() {
        
    }
    public FitFunction(int i, Map<Coordinate, GraphErrors> tvstrkdocasProf) {
        this.i = i;
        _tvstrkdocasProf = tvstrkdocasProf;
    }
         
    public double eval(double x, double[] par) {
        double p1 = par[0]; 
        double p2 = par[1]; 
        double p3 = par[2]; 
        double p4 = par[3]; 

        double calcTime = this.polyFcnMac(x,p1,  p2,  p3,  p4, i+1);
        
        return calcTime;
    }

    
    public static double polyFcnMac(double x, double p1, double p2, double p3, double p4, int superlayer) {
      
        double time = 0;
        time = p1*x*x*x*x + p2*x*x*x + p3*x*x + p4*x ;
        
        return time;
    }


    @Override
    public double valueOf(double[] par) {
        double chisq = 0;
        double delta = 0;
                    if(_tvstrkdocasProf.get(new Coordinate(this.i)).getVectorX().size()>0){ 
                        GraphErrors gr = _tvstrkdocasProf.get(new Coordinate(this.i));
                            
                        for (int ix =0; ix< gr.getDataSize(0); ix++) {
                            double x = gr.getDataX(ix);
                            double time = gr.getDataY(ix);
                            double err = gr.getDataEY(ix);
                            if(err>0) {
                                double calcTime = this.eval(x, par);
                                delta = (time - calcTime) / err; 
                                chisq += delta * delta;
                            }
                        }
                    }
        return chisq;
        
    }
    
    
}
