/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.analysis;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.plaf.metal.MetalButtonUI;
import org.clas.detector.clas12calibration.dc.calt2d.T2DCalib;
import org.freehep.math.minuit.MnUserParameters;


public class FitPanel {
    
    
    private Map<Integer, ArrayList<Double>> pars     = new HashMap<Integer, ArrayList<Double>>();
    private double[]          range    = new double[2];
    private JFrame            frame    = new JFrame();
    private CustomPanel2      panel    = null;
    private T2DCalib _pM;

    public FitPanel(T2DCalib pM) {
        this._pM = pM;
        //init pars container
        for(int j = 0; j<6; j++) {
            pars.put(j, new ArrayList<Double>());
        }
    }
    
    public void openFitPanel(String title, Map<Coordinate, MnUserParameters> TvstrkdocasFitPars){
        
        panel = new CustomPanel2(TvstrkdocasFitPars);
        frame.setSize(400, 400); 
        frame.setTitle(title);
        frame.add(panel);
        frame.pack();
        frame.setVisible(true);
            //        frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
            
    }
    public void updateFitButton() {
        panel.fitButton.setBackground(Color.GREEN);
    }
    public boolean fitted = false;
    public void refit(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) throws FileNotFoundException{
        System.out.println("Refit method called.");
        boolean[][] fixedPars = new boolean[4][6];
        for(int j = 0; j<6; j++) {
            pars.get(j).clear();
            System.out.println("Cleared pars array.");
        }
        int npar = 4;
        for(int j = 0; j<6; j++) {
            for(int i=0; i<npar; i++){   
                if(panel.params[i][j].getText().isEmpty()){
                    this.pars.get(j).add(TvstrkdocasFitPars.get(new Coordinate(j)).value(i));
                }
                else { 
                    this.pars.get(j).add(Double.parseDouble(panel.params[i][j].getText()));
                }
            }
        }
        if(!panel.minRange.getText().isEmpty())this.range[0] = Double.parseDouble(panel.minRange.getText());
        else this.range[0] = 0.0;
        if(!panel.maxRange.getText().isEmpty())this.range[1] = Double.parseDouble(panel.maxRange.getText());
        else this.range[1] = 2.0;
        for(int j = 0; j<6; j++) {    
            for(int i=0; i<npar; i++){
                TvstrkdocasFitPars.get(new Coordinate(j)).setValue(i,this.pars.get(j).get(i));
                //TvstrkdocasFitPars.get(new Coordinate(j)).setLimits(i, 
                //        this.pars.get(i)*0.8-0.1,this.pars.get(i)*1.2+0.1);
            }
        }
        
        for(int j = 0; j<6; j++) {
            for(int i=0; i<npar; i++){ 
                System.out.println("j "+j+" par "+this.pars.get(j).get(i));
                if(panel.fixFit[i][j].isSelected()==true)
                    fixedPars[i][j] = true;
            }
            
            this._pM.runFit(j, fixedPars);
            for(int p = 0; p<4; p++) {
                panel.pars[p][j] = TvstrkdocasFitPars.get(new Coordinate(j)).value(p);
                //if(p!=3) {
                    panel.params[p][j].setText(String.format("%.5f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                //} else {
               //     panel.params[p][j].setText(String.format("%.3f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
              //  }
            }
        }
        fitted = true;
        this._pM.plotFits(fitted);
        System.out.println("Refit method finished.");
    }
    
    
    public void plotResiduals() {
        this._pM.reProcess();
    }
    public void reCook() {
        this._pM.reCook();
    }
    public void reset() {
        this._pM.resetPars();
    }
    private final class CustomPanel2 extends JPanel {
        JLabel label;
        JPanel panel;
    	JTextField minRange = new JTextField(5);
	JTextField maxRange = new JTextField(5);
	JTextField[][] params = new JTextField[4][6];
        JCheckBox[][]    fixFit ;
        
        JButton   fitButton = null;
        JButton   resetButton = null;
        JButton   resButton = null;
        JButton   reCookButton = null;
        
        String[] parNames = new String[] {"p1", "p2", "p3", "p4"};
        double[][] pars = new double[4][6];
        
        public CustomPanel2(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) {        
            super(new BorderLayout());
            for(int i = 0; i < 6; i++) {
                for(int p = 0; p<4; p++) {
                    pars[p][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(p);
                }
            }
            int npar = 4;
            panel = new JPanel(new GridLayout(npar+1, 13));            
            fixFit = new JCheckBox[npar][6];

            for (int i = 0; i < npar; i++) {  
                JLabel l = new JLabel("      "+parNames[i], JLabel.LEADING);
                panel.add(l);
                
                for (int j = 0; j < 6; j++) {
                    fixFit[i][j] = new JCheckBox("Fix");
                   /* if(i==2 || i>4) {
                    fixFit[i][j].setSelected(true);
                    } else {
                        fixFit[i][j].setSelected(false);
                   } */
                    params[i][j] = new JTextField(3);

                    if(i!=3) {
                        params[i][j].setText(String.format("%.5f", pars[i][j]));
                    } else {
                        params[i][j].setText(String.format("%.3f", pars[i][j]));
                   }
                    panel.add(params[i][j]);
                    panel.add(fixFit[i][j]);
                }
            }



            panel.add(new JLabel("    Fit range min"));
            minRange.setText(Double.toString(0));
            panel.add(minRange);
            panel.add(new JLabel("    Fit range max"));
            maxRange.setText(Double.toString(2.0));
            panel.add(maxRange);
            
            resetButton = new JButton("RESET");
            resetButton.setUI(new MetalButtonUI());
            resetButton.setBackground(Color.CYAN);
            resetButton.setContentAreaFilled(false);
            resetButton.setOpaque(true);
            resetButton.setFont(bBold);
            resetButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    reset();
                    return;
                }
            });
            panel.add(resetButton);
            
            resButton = new JButton("Residuals");
            resButton.setUI(new MetalButtonUI());
            resButton.setBackground(Color.YELLOW);
            resButton.setContentAreaFilled(false);
            resButton.setOpaque(true);
            resButton.setFont(bBold);
            resButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    plotResiduals();
                    return;
                }
            });
            panel.add(resButton);
            
            reCookButton = new JButton("Reprocess");
            reCookButton.setUI(new MetalButtonUI());
            reCookButton.setBackground(Color.ORANGE);
            reCookButton.setContentAreaFilled(false);
            reCookButton.setOpaque(true);
            reCookButton.setFont(bBold);
            reCookButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    fitButton.setBackground(Color.RED);
                    fitButton.setContentAreaFilled(false);
                    fitButton.setOpaque(true);
                    System.out.println("******************************************");
                    System.out.println("* REFITTING SEGMENTS WITH NEW PARAMETERS *");
                    System.out.println("******************************************");
                    reCook();
                    System.out.println("******************************************");
                    System.out.println("*      READY TO REFIT THE HISTOGRAMS     *");
                    System.out.println("******************************************");
                    return;
                }
            });
            panel.add(reCookButton);
            
            fitButton = new JButton("FIT");
            fitButton.setUI(new MetalButtonUI());
            fitButton.setBackground(Color.RED);
            fitButton.setContentAreaFilled(false);
            fitButton.setOpaque(true);
            fitButton.setFont(bBold);
            fitButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        refit(TvstrkdocasFitPars);
                    } catch (FileNotFoundException ex) {
                        Logger.getLogger(FitPanel.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    return;
                }
            });
            
            this.add(panel, BorderLayout.CENTER);
            this.add(fitButton, BorderLayout.PAGE_END);
            
            
            label = new JLabel("Click the \"Show it!\" button"
                           + " to bring up the selected dialog.",
                           JLabel.CENTER);
        }
        private Font bBold = new Font("Arial", Font.BOLD, 16);
        void setLabel(String newText) {
            label.setText(newText);
        }

    }
}
