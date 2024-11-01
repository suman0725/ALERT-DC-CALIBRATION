/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.viewer;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dialog;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.plaf.metal.MetalButtonUI;

/**
 *
 * @author ziegler
 */
public class Driver {
    
      
    public static void main(String[] args) throws FileNotFoundException {

        System.out.println("Driver class started");

        
        JFrame    frame    = new JFrame();
        JButton   T2DButton = null;
        JButton   T0Button = null;
        JButton   TDCButton = null;
        JPanel panel = new JPanel(new GridLayout(4, 1)); 
        frame.setSize(1400, 800); 
        frame.setTitle("DC CALIBRATIONS");
        ImageIcon imageIcon = new ImageIcon("ALERT.png");
        imageIcon.getImage().getScaledInstance(800, 400, java.awt.Image.SCALE_SMOOTH);
        JLabel imgLabel = new JLabel(imageIcon);
        frame.add(imgLabel, BorderLayout.PAGE_START);
        frame.add(panel);
        frame.pack();
        frame.setVisible(true);
        
        T2DButton = new JButton("T2D");
        T2DButton.setUI(new MetalButtonUI());
        T2DButton.setBackground(Color.YELLOW);
        T2DButton.setContentAreaFilled(false);
        T2DButton.setOpaque(true);
        T2DButton.setFont(new Font("Arial", Font.BOLD, 18));
        T2DButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JFrame frame = new JFrame("DC Calibration");
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                T2DViewer viewer = null;
                try {
                    viewer = new T2DViewer();
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(Driver.class.getName()).log(Level.SEVERE, null, ex);
                }
                frame.add(viewer.mainPanel);
                frame.setJMenuBar(viewer.menuBar);
                frame.setSize(300, 300);
                frame.setVisible(true);
                viewer.configFrame.setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
                viewer.configure();
                return;
            }
        });
        panel.add(T2DButton);
        
        //Now only for T2D button 


        /* T0Button = new JButton("T0");
        T0Button.setUI(new MetalButtonUI());
        T0Button.setBackground(Color.RED);
        T0Button.setContentAreaFilled(false);
        T0Button.setOpaque(true);
        T0Button.setFont(new Font("Arial", Font.BOLD, 18));
        T0Button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JFrame frame = new JFrame("DC Calibration");
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                T0Viewer viewer = null;
                try {
                    viewer = new T0Viewer();
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(Driver.class.getName()).log(Level.SEVERE, null, ex);
                }
                frame.add(viewer.mainPanel);
                frame.setJMenuBar(viewer.menuBar);
                frame.setSize(300, 300);
                frame.setVisible(true);
                viewer.configFrame.setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
                viewer.configure();
                return;
            }
        });
        panel.add(T0Button);
        
        
        TDCButton = new JButton("TDC CUTS");
        TDCButton.setUI(new MetalButtonUI());
        TDCButton.setBackground(Color.YELLOW);
        TDCButton.setContentAreaFilled(false);
        TDCButton.setOpaque(true);
        TDCButton.setFont(new Font("Arial", Font.BOLD, 18));
        TDCButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JFrame frame = new JFrame("DC Calibration");
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                TDCViewer viewer = null;
                try {
                    viewer = new TDCViewer();
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(Driver.class.getName()).log(Level.SEVERE, null, ex);
                }
                frame.add(viewer.mainPanel);
                frame.setJMenuBar(viewer.menuBar);
                frame.setSize(300, 300);
                frame.setVisible(true);
                viewer.configFrame.setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
                viewer.configure();
                return;
            }
        }); */
      //  panel.add(TDCButton); // remove the TDC cuts Hamza
        
        
        frame.add(panel, BorderLayout.PAGE_END);
        //frame.add(T2DButton, BorderLayout.PAGE_END);
        

    }
}
