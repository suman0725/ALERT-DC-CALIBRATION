/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.analysis.FitPanel;
import org.clas.detector.clas12calibration.dc.t2d.TimeToDistanceEstimator;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import static org.clas.detector.clas12calibration.viewer.T2DViewer.ccdb;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
import java.util.Arrays;

/**
 *
 * @author ziegler
 * @author Hamza 
 */
public class T2DCalib extends AnalysisMonitor{
    public HipoDataSync calwriter = null;
    public HipoDataSync writer = null;
    private HipoDataEvent calhipoEvent = null;
    private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    FitPanel fp;
    PrintWriter pw = null;
    File outfile = null;
    private int runNumber;



    ///////////////////////////////////////////////////////////////////////
    
    public T2DCalib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames("TrackDoca vs T","TrackDoca vs T Graphs","CalcDoca vs T","Time Residuals","Parameters");
        this.init(false, "p1:p2:p3:p4");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
       
        if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
            System.out.println(" BANK FOUND........");
        } else {
            System.out.println(" BANK NOT FOUND........");
        }
        calwriter = new HipoDataSync(schemaFactory);
        calwriter.setCompressionType(2);
        writer = new HipoDataSync(schemaFactory);
        writer.setCompressionType(2);
        calhipoEvent = (HipoDataEvent) calwriter.createEvent();
        hipoEvent = (HipoDataEvent) writer.createEvent();
        calwriter.open("TestCalOutPut.hipo");
        calwriter.writeEvent(calhipoEvent);
        writer.open("TestOutPut.hipo");
        writer.writeEvent(hipoEvent);
 
    }



    private Map<Coordinate, H2F> Tvstrkdocas                = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tvscalcdocas               = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, GraphErrors> TvstrkdocasProf    = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> TvstrkdocasFit             = new HashMap<Coordinate, FitFunction>();
    private Map<Coordinate, MnUserParameters> TvstrkdocasFitPars    = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> TvstrkdocasFits                = new HashMap<Coordinate, FitLine>();
    private Map<Coordinate, H1F> timeResi               = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiFromFile       = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiNew            = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> fitResi        = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> ParsVsIter    = new HashMap<Coordinate, H1F>();
    int nsl = 6;






    ///////////////////////////////////////////////////////////////////////
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup td = new DataGroup(7,2);
        DataGroup tdp = new DataGroup(14,8);
        DataGroup cd = new DataGroup(7,2);
        DataGroup tr = new DataGroup(6,1);
        DataGroup fr = new DataGroup(6,1);
        
        for (int i = 0; i < nsl; i++) {
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            timeResi.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.5, 0.5)); 
            timeResiFromFile.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.5, 0.5)); 
            timeResiNew.put(new Coordinate(i), new H1F("time residual for sly " + (i+1), 100, -0.5, 0.5)); 
            fitResi.put(new Coordinate(i), new H1F("fit residual for sly " + (i+1), 100, -0.5, 0.5));
            
            tr.addDataSet(timeResi.get(new Coordinate(i)), i);
            tr.addDataSet(timeResiFromFile.get(new Coordinate(i)), i);
            tr.addDataSet(timeResiNew.get(new Coordinate(i)), i);
            fr.addDataSet(fitResi.get(new Coordinate(i)), i);
            
                DataGroup trkdvst = new DataGroup(1,1);
                DataGroup dvst = new DataGroup(1,1);
                
                    DataGroup prfdvst = new DataGroup(1,1);
                    //Track Doca vs T histogram
                    Tvstrkdocas.put(new Coordinate(i), new H2F("trkDocavsT" + (i + 1)*1000+26, "superlayer" + (i + 1)
                                , 200, 0, 2.0, 200, 0, 500.0+(int)(i/2)*450.0));
                    
                    TvstrkdocasProf.put(new Coordinate(i), new GraphErrors());
                    TvstrkdocasProf.get(new Coordinate(i)).setMarkerColor(i+1);
                    
                    //Calc Doca vs T histogram
                    Tvscalcdocas.put(new Coordinate(i), new H2F("calcDocavsT" + (i + 1)*1000+26, "superlayer" + (i + 1)
                            , 200, 0, 2.0, 200, 0, 500.0+(int)(i/2)*450.0));
                    tdp.addDataSet(TvstrkdocasProf.get(new Coordinate(i)), i);
                    prfdvst.addDataSet(TvstrkdocasProf.get(new Coordinate(i)), 0);
                    TvstrkdocasFits.put(new Coordinate(i), new FitLine());
                    prfdvst.addDataSet(TvstrkdocasFits.get(new Coordinate(i)), 0);
                    this.getDataGroup().add(prfdvst, 0, i+1, 0);
                    
        }
        
        this.getDataGroup().add(td, 0,0,0);
        this.getDataGroup().add(tdp,1,0,0);
        this.getDataGroup().add(cd, 2,0,0);
        this.getDataGroup().add(tr, 3,0,0);
        this.getDataGroup().add(fr, 4,0,0);
        
        for (int i = 0; i < nsl; i++) {
                this.getCalib().addEntry(0,i+1, 0);
                //blank out
                this.getCalib().setDoubleValue((double)999, "p1", 0, i+1, 0);
                this.getCalib().setDoubleValue((double)999, "p2", 0, i+1, 0);
                this.getCalib().setDoubleValue((double)999, "p3", 0, i+1, 0);
                this.getCalib().setDoubleValue((double)999, "p4", 0, i+1, 0);
        }
        
        this.getCalib().fireTableDataChanged();
        System.out.println("Pass createHistos function");
    }




    ///////////////////////////////////////////////////////////////////////
    private void updateTable(int i) {
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(0), "p1", 0, i+1, 0);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(1), "p2", 0, i+1, 0);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(2), "p3", 0, i+1, 0);
        this.getCalib().setDoubleValue(TvstrkdocasFitPars.get(new Coordinate(i)).value(3), "p4", 0, i+1, 0);
        System.out.println("Pass updateTable function");
    }    



    ///////////////////////////////////////////////////////////////////////
    @Override
    public void plotHistos() {
        String[] Names = {"TrackDoca vs T","TrackDoca vs T Graphs","CalcDoca vs T","Time Residuals","Parameters"};
        for(int s = 0; s<3; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
            int NumPads = 
            this.getAnalysisCanvas().getCanvas(Names[s]).getCanvasPads().size();
            for (int n = 0; n < NumPads; n++) {
                this.getAnalysisCanvas().getCanvas(Names[s]).getPad(n).getAxisZ().setLog(true);
            }
        }
        
        for(int s = 3; s<5; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
        }
        this.getAnalysisCanvas().getCanvas(Names[3]).divide(this.nsl, 2);
        this.getAnalysisCanvas().getCanvas(Names[4]).divide(6, 6);
        
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").update();
        this.getAnalysisCanvas().getCanvas("CalcDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("Time Residuals").update();
        this.getAnalysisCanvas().getCanvas("Parameters").update();
        System.out.println("Pass plotHistos function");
    }




    @Override
    public void timerUpdate() {
    }
    

    ///////////////////////////////////////////////////////////////////////
    @Override
    public void analysis() {
        calwriter.close();
        writer.close();
       // this.UpdateBBinCenters();
        for (int i = 0; i < this.nsl; i++) {
            
            this.filltrkDocavsTGraphs(i);
            //runFit(i); 
        }
        this.getAnalysisCanvas().getCanvas("Time Residuals").divide(nsl, 3);
        //
        for(int i = 0; i<this.nsl; i++) {
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i);
            this.fitTimeResPlot(timeResiFromFile.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
        }
        fp.updateFitButton();
        
        System.out.println("Pass analysis function");
    }



    public void plotFits(boolean fitted) throws FileNotFoundException {
        if(fitted==true) {
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "Files/ccdb_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
           
            pw = new PrintWriter(fileName);
            pw.printf("#& sector superlayer component p1 p2 p3 p4\n");
        
            int ij =0;
            int ip =0;
            NbRunFit++;
            for (int i = 0; i < this.nsl; i++) {
               
                //for(int p = 0; p<6; p++) {
                    for(int p = 0; p<4; p++) {
                    ParsVsIter.get(new Coordinate(i,p)).setBinContent(NbRunFit, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
                    ParsVsIter.get(new Coordinate(i,p)).setBinError(NbRunFit, TvstrkdocasFitPars.get(new Coordinate(i)).error(p));
                    ParsVsIter.get(new Coordinate(i,p)).setOptStat(0);
                    this.getAnalysisCanvas().getCanvas("Parameters").cd(ip);
                    GraphErrors gr = new GraphErrors();
                    double min = ParsVsIter.get(new Coordinate(i,p)).getMin()-2*ParsVsIter.get(new Coordinate(i,p)).getBinError(1);
                    double max = ParsVsIter.get(new Coordinate(i,p)).getMax()+2*ParsVsIter.get(new Coordinate(i,p)).getBinError(1);
                    if(Math.abs(min)<1.e-06 && Math.abs(max)<1.e-06) {
                        min = -0.1;
                        max = 0.1;
                    }
                    
                    this.getAnalysisCanvas().getCanvas("Parameters").
                            draw(ParsVsIter.get(new Coordinate(i,p)));

                    ip++;
                }
            }

            for (int i = 0; i < this.nsl; i++) {


                                this.updateTable(i);
                                TvstrkdocasFits.put(new Coordinate(i), new FitLine("f"+""+i, i, 
                                TvstrkdocasFitPars.get(new Coordinate(i))));
                                TvstrkdocasFits.get(new Coordinate(i)).setLineStyle(4);
                                TvstrkdocasFits.get(new Coordinate(i)).setLineWidth(5);
                                TvstrkdocasFits.get(new Coordinate(i)).setLineColor(i+1);
            }
            this.getCalib().fireTableDataChanged();    
            for(int isec = 0; isec < 6; isec++) {
                for(int i = 0; i<6; i++) {
                    pw.printf("%d\t %d\t %d\t %.6f\t %d\t %.6f\t %.6f\t %.6f\t %d\n",
                    (isec+1), (i+1), 0,
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(0),
                    0,
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(1),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(2),
                    TvstrkdocasFitPars.get(new Coordinate(i)).value(3),
                    0);
                    
                }
            }
            pw.close();
            //this.reProcess();
        }
        System.out.println("Pass plotfits function");
    }
 
    



    private int maxIter = 10;
    double[][] fixPars = new double[6][4];
    
    private MnScan  scanner = null;
    private MnMigrad migrad = null;
    
    public int NbRunFit = 0;



    public void runFit(int i, boolean fixFit[][]) {
        // i = superlayer - 1;
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** ");
        
        TvstrkdocasFit.put(new Coordinate(i), 
                new FitFunction(i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
        
        scanner = new MnScan((FCNBase) TvstrkdocasFit.get(new Coordinate(i)), 
                TvstrkdocasFitPars.get(new Coordinate(i)),2);
	
        //scanner.fix(4);
        for (int p = 0; p < 4; p++) {
            if(fixFit[p][i]==true) {
                scanner.fix(p); 
            }
        }
        FunctionMinimum scanmin = scanner.minimize();
        if(scanmin.isValid())
            TvstrkdocasFitPars.put(new Coordinate(i),scanmin.userParameters());
        
        migrad = new MnMigrad((FCNBase) TvstrkdocasFit.get(new Coordinate(i)),
                TvstrkdocasFitPars.get(new Coordinate(i)),1);
        migrad.setCheckAnalyticalDerivatives(true);
        
        FunctionMinimum min ;
        
        
        for(int it = 0; it<maxIter; it++) {
            
            min = migrad.minimize();
            System.err.println("****************************************************");
            System.err.println("*   FIT RESULTS  FOR SUPERLAYER  "+(i+1)+" at iteration "+(it+1)+"  *");
            System.err.println("****************************************************");  
            for(int pi = 0; pi<4; pi++) 
                System.out.println("par["+pi+"]="+(TvstrkdocasFitPars.get(new Coordinate(i)).value(pi)-min.userParameters().value(pi)));
            
            if(min.isValid()) {
                TvstrkdocasFitPars.put(new Coordinate(i),min.userParameters());  
            }
            
            
            System.err.println(min);
            
        }
        
        
        //release all so they can be fixed again
       // TvstrkdocasFitPars.get(new Coordinate(i)).release(4);
        for (int p = 0; p < 4; p++) {
            if(fixFit[p][i]==true) {
                TvstrkdocasFitPars.get(new Coordinate(i)).release(p);
            }
        }
        if(i==5) {
            System.err.println("***************************************************");
            System.err.println("*                  FITS ARE DONE!                 *");
            System.err.println("***************************************************");  
            
        }
    }
    






    int counter = 0;
    private int iterationNum = 1;
    public  HipoDataSource calreader = new HipoDataSource();
    public  HipoDataSource reader = new HipoDataSource();




    
    public void reCook() {
        iterationNum++;
        for (int i = 0; i < nsl; i++) {
                    Tvstrkdocas.get(new Coordinate(i)).reset();
                    Tvscalcdocas.get(new Coordinate(i)).reset();
                    //TvstrkdocasProf.get(new Coordinate(i,j,k)).reset();
        }
        reLoadFitPars();
        
        reader = new HipoDataSource();
        reader.open("TestOutPut.hipo");
        while (reader.hasEvent()) { 
            hits.clear();
            DataEvent event = reader.getNextEvent();
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        hits.add(this.getCalHit(bnkHits, i));
                }
                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit);
                }
                //refit with new constants
                Refit rf = new Refit(hits);
                rf.reFit();    
            }
        }
        reader.close();
        
        calreader = new HipoDataSource();
        calreader.open("TestCalOutPut.hipo");
        while (calreader.hasEvent()) { 
            calhits.clear();
            DataEvent event = calreader.getNextEvent();
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        calhits.add(this.getCalHit(bnkHits, i));
                }
                for(FittedHit hit : calhits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit);
                }
                //refit with new constants
                Refit rf = new Refit(calhits);
                rf.reFit();    
                for(FittedHit hit : calhits) {
                         
                        Tvstrkdocas.get(new Coordinate(hit.get_Superlayer() - 1)).fill(hit.get_ClusFitDoca(), hit.get_Time());
                        Tvscalcdocas.get(new Coordinate(hit.get_Superlayer() - 1)).fill(hit.get_Doca(), hit.get_Time());

                    }
                }
            }
        //}
        for (int i = 0; i < this.nsl; i++) {
                this.filltrkDocavsTGraphs(i);
        }
        System.out.println("RECOOKING DONE!");
        calreader.close();
        fp.updateFitButton();
    }





    public void reProcess() {
        //reset histos to refill
        for (int i = 0; i < this.nsl; i++) {
               timeResi.get(new Coordinate(i)).reset();
               timeResiNew.get(new Coordinate(i)).reset();
                
        }
        reader = new HipoDataSource();
        reader.open("TestOutPut.hipo");
        while (reader.hasEvent()) { 
            hits.clear();
            DataEvent event = reader.getNextEvent();
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        hits.add(this.getCalHit(bnkHits, i));
                }

                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit);
                }
                //refit with new constants
                Refit rf = new Refit(hits);
                rf.reFit();
                // fill calibrated plot
                for(FittedHit hit : hits) {
                    timeResi.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                }

            } 
        }
        


        reader = new HipoDataSource();
        reader.open("TestOutPut.hipo");
        while (reader.hasEvent()) { 
            hits.clear();
            DataEvent event = reader.getNextEvent();
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
                
                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        hits.add(this.getCalHit(bnkHits, i));
                }

                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit);
                }
                //refit with new constants
                Refit rf = new Refit(hits);
                rf.reFit();
                // fill calibrated plot
                for(FittedHit hit : hits) {
                    timeResi.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                }

            } 
        }
        
        //--------------------------------------------
        System.out.println("************ reloading Fit Parameters");
        reLoadFitPars();
        //--------------------------------------------
        reader.gotoEvent(0);
        
        hits.clear();
        while (reader.hasEvent()) {
            DataEvent event = reader.getNextEvent();
            if(event.hasBank("TimeBasedTrkg::TBHits")) {
                DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");

                for (int i = 0; i < bnkHits.rows(); i++) {
                    if(this.getCalHit(bnkHits, i)!=null)
                        hits.add(this.getCalHit(bnkHits, i));
                }
                
                for(FittedHit hit : hits) {
                    hit.set_TimeResidual(-999);
                    updateHit(hit);
                }
                //refit with new constants
                Refit rf = new Refit(hits);
                rf.reFit();
                // fill calibrated plot
                for(FittedHit hit : rf.hits) {
                    timeResiNew.get(new Coordinate(hit.get_Superlayer()-1)).fill(hit.get_TimeResidual());
                }
                //DataBank bnkHitsUpd = this.fillTBHitsBank(event, rf.calhits);
                //this.processEventIterate(bnkHitsUpd);
                hits.clear();
            }
        }
        this.getAnalysisCanvas().getCanvas("Time Residuals").clear();
        this.getAnalysisCanvas().getCanvas("Time Residuals").divide(nsl, 3);
        //
        for(int i = 0; i<this.nsl; i++) {
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i);
            this.fitTimeResPlot(timeResiFromFile.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i+6);
            this.fitTimeResPlot(timeResi.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
            this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i+12);
            //this.getAnalysisCanvas().getCanvas("Time Residuals").draw(timeResiNew.get(new Coordinate(i)));
            this.fitTimeResPlot(timeResiNew.get(new Coordinate(i)), 
                    this.getAnalysisCanvas().getCanvas("Time Residuals"));
        }
        //for(int i = 0; i<this.nsl; i++) {
        //    this.getAnalysisCanvas().getCanvas("Fit Residuals").cd(i);
        //    this.getAnalysisCanvas().getCanvas("Fit Residuals").draw(fitResi.get(new Coordinate(i)));
        //}
        System.out.println("REPROCESSING DONE!");
        reader.close();
        this.getCalib().fireTableDataChanged();  
    }
         




    private void fitTimeResPlot(H1F h1, EmbeddedCanvas canvasRes) {
        if (h1==null) return;
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+[amp2]*gaus(x,[mean],[sigma2])", -0.5, 0.5); 
        gausFunc.setLineColor(4);
        gausFunc.setLineStyle(1);
        gausFunc.setLineWidth(2);
        gausFunc.setParameter(0, h1.getMax());
        gausFunc.setParameter(1, -0.0);
        gausFunc.setParameter(2, 0.05);
        gausFunc.setParameter(3, h1.getMax()/2.);
        gausFunc.setParameter(4, 0.5);
        gausFunc.setOptStat(1110);
        h1.setOptStat(0);
        //canvasRes.clear();
        
        DataFitter.fit(gausFunc, h1, "Q");
        //gausFunc.setOptStat(101100);
        gausFunc.setOptStat(0);
        
        int effSig = (int) (10000*this.getError(gausFunc.getParameter(0), gausFunc.getParameter(1), gausFunc.getParameter(2), 
                gausFunc.getParameter(3), gausFunc.getParameter(1), gausFunc.getParameter(4)));
        String t = "eff. Sig. "+effSig + "microns";
        h1.setTitle(t);
        canvasRes.draw(h1, "E1");
        //canvasRes.draw(gausFunc, "same");
        
        
    }







    private int MINENTRIES = 10;
    F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])+[p0]", 0, 1.8);
    F1D f2 = new F1D("f2","[amp1]*gaus(x,[mean1],[sigma1])+[amp2]*gaus(x,[mean2],[sigma2])+[p02]", 0, 1.8);
    
    private void filltrkDocavsTGraphs(int i) {
        
        if(TvstrkdocasProf.get(new Coordinate(i))!=null) {
            
            TvstrkdocasProf.get(new Coordinate(i)).reset();
            H2F h2 = Tvstrkdocas.get(new Coordinate(i));
            ArrayList<H1F> hslice = h2.getSlicesX();
            
            for(int si=0; si<hslice.size(); si++) {
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                double mean = hslice.get(si).getDataX(hslice.get(si).getMaximumBin());
                if(amp<this.MINENTRIES) {
                    
                } else {
                    
                    double x = h2.getXAxis().getBinCenter(si);
                    double y = hslice.get(si).getMean();
                    double sigma = hslice.get(si).getRMS();
                   
                    if(x/(2.*Constants.wpdist[i])>0.9) { // for large docas (90% of dmax take fit with a double gauss and take the peak of the one at smaller distance as y
                        f2.setParameter(0, amp);
                        f2.setParameter(1, mean);
                        f2.setParameter(4, y+10);
                        f2.setParameter(2, sigma);
                        f2.setParameter(5, sigma+10);
                        f2.setParameter(6, 0);
                    
                        DataFitter.fit(f2, hslice.get(si), "Q"); //No options uses error for sigma 
                        if(f2.getChiSquare()<200 && f2.getParameter(1)>0 && f2.parameter(1).error()<50) {
                            //TvstrkdocasProf.get(new Coordinate(i, j, k)).
                            //        addPoint(x, f2.getParameter(1), 0, this.getError(f2.getParameter(0),
                            //                f2.getParameter(1),f2.getParameter(2),f2.getParameter(3),
                            //                f2.getParameter(4), f2.getParameter(5) ));
                            TvstrkdocasProf.get(new Coordinate(i)).
                                    addPoint(x, f2.getParameter(1), 0, sigma);
                        }
                    } else {
                        f1.setParameter(0, amp);
                        f1.setParameter(1, mean);
                        f1.setParameter(2, sigma);
                        f1.setParameter(3, 0);
                    
                        DataFitter.fit(f1, hslice.get(si), "Q"); //No options uses error for sigma 
                        if(f1.getChiSquare()<200 && f1.getParameter(1)>0 && f1.parameter(1).error()<50) {
                            TvstrkdocasProf.get(new Coordinate(i)).
                                    addPoint(x, f1.getParameter(1), 0, f1.getParameter(2));
                        }
                    }
                }
            }
        }
    }



    int count = 0;
    public static int polarity =-1;
    public List<FittedHit> calhits = new ArrayList<>();
    public List<FittedHit> hits = new ArrayList<>();
    Map<Integer, ArrayList<Integer>> segMapTBHits = new HashMap<Integer, ArrayList<Integer>>();
    Map<Integer, SegmentProperty> segPropMap = new HashMap<Integer, SegmentProperty>();
    List<FittedHit> calhitlist = new ArrayList<>();
    List<FittedHit> hitlist = new ArrayList<>();
    @Override


    
    public void processEvent(DataEvent event) {
        
        if (!event.hasBank("RUN::config")) {
            return ;
        }
        
        DataBank bank = event.getBank("RUN::config");
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) {
           return ;
        } else {
           count++;
        }
       
        if(count==1) {
            Constants.Load();
            String newVar = String.valueOf(T2DViewer.calVariation.getSelectedItem());
            System.out.println("* VARIATION *"+newVar);
            ccdb.setVariation(newVar);
            TableLoader.FillT0Tables(newRun, newVar);
            TableLoader.Fill(T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/time2dist")); 
            this.loadFitPars(); 
            polarity = (int)Math.signum(event.getBank("RUN::config").getFloat("torus",0));
            runNumber = newRun;
            
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
        // get segment property
        
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        this.getSegProperty(bnkHits);
        
        for (int i = 0; i < bnkHits.rows(); i++) {
            int superlayer = bnkHits.getInt("superlayer", i);

            if(this.passResiCuts(event, bnkHits, i))
                hitlist.add(this.getHit(bnkHits, i)); 
            if(this.passCalibCuts(event,bnkHits, i)){
                calhitlist.add(this.getHit(bnkHits, i));
                double calibTime = (double) (bnkHits.getInt("TDC", i) - bnkHits.getFloat("TProp", i)
                                        - bnkHits.getFloat("TFlight", i) - bnkHits.getFloat("TStart", i) 
                                        - bnkHits.getFloat("T0", i));

                Tvstrkdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1)).fill(bnkHits.getFloat("trkDoca", i), calibTime);
                Tvscalcdocas.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1)).fill(bnkHits.getFloat("doca", i), calibTime);
                
                // fill uncalibrated plot
               timeResiFromFile.get(new Coordinate(bnkHits.getInt("superlayer", i) - 1))
                                .fill(bnkHits.getFloat("timeResidual", i)); 
            //                    System.out.println("====================== filling the residuals");
            //                }
                
            }
        }
        calhipoEvent = (HipoDataEvent) calwriter.createEvent();
        hipoEvent = (HipoDataEvent) writer.createEvent();
        //selected.show();
        calhipoEvent.appendBank(this.fillTBHitsBank(event, calhitlist));
        hipoEvent.appendBank(this.fillTBHitsBank(event, hitlist));
        calwriter.writeEvent(calhipoEvent);
        writer.writeEvent(hipoEvent);
        calhitlist.clear();
        hitlist.clear();
      // System.out.println("Pass processevent function");
    }
    
   // public double[][] resetPars = new double[6][11];
    //private String[] parNames = {"v0", "vmid", "R", "tmax", "distbeta", "delBf","b1", "b2", "b3", "b4", "dmax"};
    //private double[] errs = {0.001,0.001,0.01,1.0,0.01,0.001,0.001,0.001,0.001,0.001,0.00001};
    public double[][] resetPars = new double[6][4];
    private String[] parNames = {"p1", "p2", "p3", "p4"};
    private double[] errs = {0.01,0.001,0.001,0.001};
    
    public void resetPars() {
        for (int i = 0; i < this.nsl; i++) {
            double[] pars = resetPars[i];
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            for(int p = 0; p < 4; p++) {
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                //create graphs of parameters for various iterations
               ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p+": " +parNames[p]+", superlayer "+(i+1),this.maxIter+1, 0.5,this.maxIter+1.5));
            }
            //TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[3], pars[3], errs[3]);
        }
        System.out.println("======================it is working now!!");
        fp.openFitPanel("fit panel", TvstrkdocasFitPars);
    }
    public void loadFitPars() {
        for (int i = 0; i < this.nsl; i++) {
            double[] pars = new double[4];
        //    T2DFunctions.polyFcnMac(x, p1, p2, p3, p4, superlayer) ;
 
            pars[0] = TableLoader.p1[0][i];
            pars[1] = TableLoader.p2[0][i];
            pars[2] = TableLoader.p3[0][i];
            pars[3] = TableLoader.p4[0][i];

            resetPars[i] = pars;
            TvstrkdocasFitPars.put(new Coordinate(i), new MnUserParameters());
            for(int p = 0; p < 4; p++) {
                TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[p], pars[p], errs[p]);
                //create graphs of parameters for various iterations
                ParsVsIter.put(new Coordinate(i,p), new H1F("h"+p, "superlayer "+(i+1)+" par "+p,this.maxIter+1, 0.5,this.maxIter+1.5));
            }
           //TvstrkdocasFitPars.get(new Coordinate(i)).add(parNames[3], pars[3], errs[3]);
        }   
        // Fit panel
        fp = new FitPanel(this);
        fp.openFitPanel("fit panel", TvstrkdocasFitPars);
    } 

    
    



    private void reLoadFitPars() {
        for (int s =0; s < 6; s++) {
            for (int i = 0; i < this.nsl; i++) {
              
                TableLoader.p1[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(0);
                TableLoader.p2[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(1);
                TableLoader.p3[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(2);
                TableLoader.p4[s][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(3);
            }
        }
        for (int i = 0; i < this.nsl; i++) {
            for(int p = 0; p<4; p++) {
                ParsVsIter.get(new Coordinate(i,p)).setBinContent(0, TvstrkdocasFitPars.get(new Coordinate(i)).value(p));
                ParsVsIter.get(new Coordinate(i,p)).setBinError(0, TvstrkdocasFitPars.get(new Coordinate(i)).error(p));
            }
        }
        TableLoader.ReFill();
    }
    
    public void Plot(int i,int j) {
        
            if(TvstrkdocasProf.get(new Coordinate(i)).getVectorX().size()>0) {
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T").cd(0);
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T").draw(Tvstrkdocas.get(new Coordinate(i)));
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").cd(0);
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").draw(Tvstrkdocas.get(new Coordinate(i)));
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").cd(0);
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").draw(Tvscalcdocas.get(new Coordinate(i)));
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").draw(TvstrkdocasProf.get(new Coordinate(i)), "same");
                this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").draw(TvstrkdocasFits.get(new Coordinate(i)), "same");
                this.getAnalysisCanvas().getCanvas("CalcDoca vs T").draw(TvstrkdocasFits.get(new Coordinate(i)), "same");
            }
    }





    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_component = (String) cc.getValueAt(row, 2);
        System.out.println(str_sector + " " + str_layer + " " + str_component);
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);
       int component = Integer.parseInt(str_component);

       if(group.hasItem(sector,layer,component)==true){
           this.Plot(layer-1, component-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       //System.out.println("Pass");
   
    }





    private FittedHit getCalHit(DataBank bnkHits, int i) {
        FittedHit hit = null;
        int id = bnkHits.getShort("id", i);;
        int sector = bnkHits.getByte("sector", i);
        int superlayer = bnkHits.getByte("superlayer", i);
        int layer = bnkHits.getByte("layer", i);
        int wire = bnkHits.getShort("wire", i);
        int TDC = bnkHits.getInt("TDC", i);
        double doca = bnkHits.getFloat("doca", i);
        double docaError = bnkHits.getFloat("docaError", i);
        double trkDoca = bnkHits.getFloat("trkDoca", i);
        int LR = bnkHits.getByte("LR", i);
        double X = bnkHits.getFloat("X", i);
        double Z = bnkHits.getFloat("Z", i);
        //double B = bnkHits.getFloat("B", i);
        //double Alpha = bnkHits.getFloat("Alpha", i);
        double TProp = bnkHits.getFloat("TProp", i);
        double TFlight = bnkHits.getFloat("TFlight", i);
        double T0 = bnkHits.getFloat("T0", i);
        double TStart = bnkHits.getFloat("TStart", i);
        int clusterID = bnkHits.getShort("clusterID", i);
        int trkID = bnkHits.getByte("trkID", i);
        double time = bnkHits.getFloat("time", i);
        double beta = bnkHits.getFloat("beta", i);
        double tBeta = bnkHits.getFloat("tBeta", i);
        double resiTime = bnkHits.getFloat("timeResidual", i);
        double resiFit = bnkHits.getFloat("fitResidual", i);
        

        hit = new FittedHit(sector, superlayer, layer, wire, TDC, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        double T0Sub = (TDC - TProp - TFlight - T0);
        hit.set_Time(T0Sub-TStart);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(T2DViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(T2DViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
            //        hit.setAlpha(alphaRadUncor);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        this.getSegProperty(bnkHits);
        
        return hit;
    }

    private FittedHit getHit(DataBank bnkHits, int i) {
        FittedHit hit = null;
        int id = bnkHits.getShort("id", i);;
        int sector = bnkHits.getByte("sector", i);
        int superlayer = bnkHits.getByte("superlayer", i);
        int layer = bnkHits.getByte("layer", i);
        int wire = bnkHits.getShort("wire", i);
        int TDC = bnkHits.getInt("TDC", i);
        double doca = bnkHits.getFloat("doca", i);
        double docaError = bnkHits.getFloat("docaError", i);
        double trkDoca = bnkHits.getFloat("trkDoca", i);
        int LR = bnkHits.getByte("LR", i);
        double X = bnkHits.getFloat("X", i);
        double Z = bnkHits.getFloat("Z", i);
      //  double B = bnkHits.getFloat("B", i);
      //  double Alpha = bnkHits.getFloat("Alpha", i);
        double TProp = bnkHits.getFloat("TProp", i);
        double TFlight = bnkHits.getFloat("TFlight", i);
        double T0 = bnkHits.getFloat("T0", i);
        double TStart = bnkHits.getFloat("TStart", i);
        int clusterID = bnkHits.getShort("clusterID", i);
        int trkID = bnkHits.getByte("trkID", i);
        double time = bnkHits.getFloat("time", i);
        double beta = bnkHits.getFloat("beta", i);
        double tBeta = bnkHits.getFloat("tBeta", i);
        double resiTime = bnkHits.getFloat("timeResidual", i);
        double resiFit = bnkHits.getFloat("fitResidual", i);
        
      //  double bFieldVal = (double) bnkHits.getFloat("B", i);
        //int region = (int) (superlayer + 1) / 2;
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
     //   double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
     //   double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;

        hit = new FittedHit(sector, superlayer, layer, wire, TDC, id);
        hit.set_Id(id); // use event number as id to recompose the clusters
       // hit.setB(B);
        hit.setT0(T0);
        hit.setTStart(TStart);
        hit.setTProp(TProp);
        hit.set_Beta(beta);
        hit.setTFlight(TFlight);
        double T0Sub = (TDC - TProp - TFlight - T0);
        hit.set_Time(T0Sub-TStart);
        hit.set_LeftRightAmb(LR);
        hit.calc_CellSize(T2DViewer.dcDetector);
        hit.set_X(X);
        hit.set_Z(Z);
        hit.calc_GeomCorr(T2DViewer.dcDetector, 0);
        hit.set_ClusFitDoca(trkDoca);
        hit.set_DeltaTimeBeta(tBeta);
        hit.set_Doca(doca);
        hit.set_Time(time);
      //  hit.setAlpha(alphaRadUncor);
        hit.set_DocaErr(docaError);
        hit.set_AssociatedClusterID(clusterID);
        hit.set_AssociatedHBTrackID(trkID); 
        hit.set_TimeResidual(resiTime);
        hit.set_Residual(resiFit);
        this.getSegProperty(bnkHits);
        
        return hit;
    }


    private DataBank fillTBHitsBank(DataEvent event, List<FittedHit> hitlist) {
        //if(event.hasBank("TimeBasedTrkg::TBHits")) { // for second pass tracking
        //     event.removeBank("TimeBasedTrkg::TBHits");
        //}
        DataBank bank = event.createBank("TimeBasedTrkg::TBHits", hitlist.size());

        for (int i = 0; i < hitlist.size(); i++) {
            
            bank.setShort("id", i, (short) hitlist.get(i).get_Id());
            bank.setShort("status", i, (short) hitlist.get(i).get_QualityFac());
            bank.setByte("superlayer", i, (byte) hitlist.get(i).get_Superlayer());
            bank.setByte("layer", i, (byte) hitlist.get(i).get_Layer());
            bank.setByte("sector", i, (byte) hitlist.get(i).get_Sector());
            bank.setShort("wire", i, (short) hitlist.get(i).get_Wire());

            bank.setFloat("X", i, (float) hitlist.get(i).get_X());
            bank.setFloat("Z", i, (float) hitlist.get(i).get_Z());
            bank.setByte("LR", i, (byte) hitlist.get(i).get_LeftRightAmb());
            bank.setFloat("time", i, (float) (hitlist.get(i).get_Time() + 0*hitlist.get(i).get_DeltaTimeBeta()));
            bank.setFloat("tBeta", i, (float) hitlist.get(i).get_DeltaTimeBeta());
            bank.setFloat("fitResidual", i, (float) hitlist.get(i).get_Residual());
           // bank.setFloat("Alpha", i, (float) hitlist.get(i).getAlpha());
            
            bank.setFloat("doca", i, (float) hitlist.get(i).get_Doca());
            bank.setFloat("docaError", i, (float) hitlist.get(i).get_DocaErr());
            bank.setFloat("trkDoca", i, (float) hitlist.get(i).get_ClusFitDoca());

            bank.setShort("clusterID", i, (short) hitlist.get(i).get_AssociatedClusterID());
            bank.setByte("trkID", i, (byte) hitlist.get(i).get_AssociatedHBTrackID());
            bank.setFloat("timeResidual", i, (float) hitlist.get(i).get_TimeResidual());
            
            bank.setInt("TDC",i,hitlist.get(i).get_TDC());
            //bank.setFloat("B", i, (float) hitlist.get(i).getB());
            bank.setFloat("TProp", i, (float) hitlist.get(i).getTProp());
            bank.setFloat("TFlight", i, (float) hitlist.get(i).getTFlight());
            bank.setFloat("T0", i, (float) hitlist.get(i).getT0());
            bank.setFloat("TStart", i, (float) hitlist.get(i).getTStart());
            
            bank.setFloat("beta", i, (float) hitlist.get(i).get_Beta());
            
            

        }
        //System.out.println(" Created Bank "); bank.show();
        return bank;

    }
    private int readPID(DataEvent event, int trkId) {
        int pid = 0;
        //fetch the track associated pid from the REC tracking bank
        if (!event.hasBank("REC::Particle") || !event.hasBank("REC::Track"))
            return pid;
        DataBank bank = event.getBank("REC::Track");
        //match the index and return the pid
        int rows = bank.rows();
        for (int i = 0; i < rows; i++) {
            if (bank.getByte("detector", i) == 6 &&
                    bank.getShort("index", i) == trkId - 1) {
                DataBank bank2 = event.getBank("REC::Particle");
                if(bank2.getByte("charge", bank.getShort("pindex", i))!=0) {
                    pid = bank2.getInt("pid", bank.getShort("pindex", i));
                }
            }
        }
        
        return pid;
    } 
    
    private void updateHit(FittedHit hit) { 
        double distbeta = TvstrkdocasFitPars.get(new Coordinate(hit.get_Superlayer()-1)).value(4);
        double deltatime_beta = (Math.sqrt(hit.get_ClusFitDoca() * hit.get_ClusFitDoca() 
                + (distbeta * hit.get_Beta() * hit.get_Beta()) 
                * (distbeta* hit.get_Beta() * hit.get_Beta())) - hit.get_ClusFitDoca()) / Constants.V0AVERAGED;
        hit.set_DeltaTimeBeta(deltatime_beta);
        hit.set_Doca(this.timeToDistance(hit.get_Sector(), hit.get_Superlayer(), hit.get_Beta(), hit.get_ClusFitDoca(), hit.get_Time(), distbeta));
        double x = hit.get_XWire(); 
       // double alphaRadUncor = Math.toRadians(hit.getAlpha());
       // double trkAngle = Math.tan(alphaRadUncor);
       // double cosTkAng = 1./Math.sqrt(trkAngle*trkAngle + 1.);
       // hit.set_X(x + hit.get_LeftRightAmb() * (hit.get_Doca() / cosTkAng) );
        
    }
    
    private TimeToDistanceEstimator tde = new TimeToDistanceEstimator();
    
    private double timeToDistance(int sector, int superlayer, double beta, double x, double time,
                                double distbeta) {
        double distance = 0;
        
         
        if(beta>1.0)
            beta = 1.0;
        double deltatime_beta = (Math.sqrt(x * x + (distbeta * beta * beta) 
                * (distbeta* beta * beta)) - x) / Constants.V0AVERAGED;
        double correctedTime = time -deltatime_beta;
        if(correctedTime<=0)
            correctedTime=0.01; // fixes edge effects ... to be improved
        distance = tde.interpolateOnGrid(correctedTime,sector-1, superlayer-1) ; 
        return distance;
    }


    private void getSegProperty(DataBank bnkHits) {
        
        segMapTBHits.clear();
        segPropMap.clear();
               
        for (int j = 0; j < bnkHits.rows(); j++){
            Integer cID = bnkHits.getInt("clusterID", j);
            if(segMapTBHits.containsKey(cID)==false) {
                segMapTBHits.put(cID, new ArrayList<Integer>());
            }
            
            segMapTBHits.get(cID).add((int)bnkHits.getShort("wire", j));
            
        }
        
        Iterator<Map.Entry<Integer, ArrayList<Integer>>> itr = segMapTBHits.entrySet().iterator(); 
          
        while(itr.hasNext()) { 
            Map.Entry<Integer, ArrayList<Integer>> entry = itr.next(); 
            segPropMap.put(entry.getKey() , 
                    new SegmentProperty(entry.getKey(),entry.getValue(),Integer.parseInt(T2DViewer.deltaWire.getText())));
        } 

    }

    private double getError(double a1, double mu1, double sig_1, double a2, double mu2, double sig_2) {
        
        double sig1 = Math.abs(sig_1);
        double sig2 = Math.abs(sig_2);
        
        double A1 = Math.sqrt(2*Math.PI)*sig1*a1;
        double A2 = Math.sqrt(2*Math.PI)*sig2*a2;
        
        double Momt1 = A1/(A1+A2) * mu1 + A2/(A1+A2) * mu2;
        double Momt2 = A1/(A1+A2) * (mu1*mu1 + sig1*sig1) + A2/(A1+A2) * (mu2*mu2 + sig2*sig2);
        
        double var = Momt2 - Momt1*Momt1;
        
        return Math.sqrt(var);
    }

    private boolean passPID(DataEvent event, DataBank bnkHits, int rowIdxinTrkBank) {
        boolean pass = false;
        int trkID = bnkHits.getByte("trkID", rowIdxinTrkBank);
        int pid = this.readPID(event, trkID);
        //pass if the track is identified as an electron or as a hadron
        if(pid==11 || Math.abs(pid)==211 || Math.abs(pid)==2212 || Math.abs(pid)==321) {
            pass = true;
        }
        return pass;
    }
    
    private boolean passResiCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        
        double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
        double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
        double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;

        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("TFlight", i)>0 
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && this.passPID(event, bnkHits, i)==true)
            {
                pass = true;
            }
        return pass;
    }

    private boolean passCalibCuts(DataEvent event, DataBank bnkHits, int i) {
        boolean pass = false;
        
        int superlayer = bnkHits.getInt("superlayer", i);

        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("beta", i)> Double.parseDouble(T2DViewer.betaCut.getText()) 
                    && bnkHits.getFloat("TFlight", i)>0 
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getNumWireWithinDW()<=Integer.parseInt(T2DViewer.npassWires.getText())
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()>Integer.parseInt(T2DViewer.nWires.getText())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && this.passPID(event, bnkHits, i)==true)
            {
                pass = true;
            }
        return pass;
    } 
    
    private boolean passCalibCuts(DataBank bnkHits, int i) {
        boolean pass = false;
        
        //double bFieldVal = (double) bnkHits.getFloat("B", i);
        int superlayer = bnkHits.getInt("superlayer", i);
        // alpha in the bank is corrected for B field.  To fill the alpha bin use the uncorrected value
       // double theta0 = Math.toDegrees(Math.acos(1-0.02*bFieldVal));
       // double alphaRadUncor = bnkHits.getFloat("Alpha", i)+(double)T2DCalib.polarity*theta0;

        if (bnkHits.getByte("trkID", i) >0 
                    && bnkHits.getFloat("beta", i)> Double.parseDouble(T2DViewer.betaCut.getText()) 
                    && bnkHits.getFloat("TFlight", i)>0 
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getNumWireWithinDW()<=Integer.parseInt(T2DViewer.npassWires.getText())
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.0001*Double.parseDouble(T2DViewer.fitresiCut.getText()) 
                    && segPropMap.get(bnkHits.getInt("clusterID", i)).getSize()>Integer.parseInt(T2DViewer.nWires.getText()))
            {
                pass = true;
            }
        return pass;
    } 
}


