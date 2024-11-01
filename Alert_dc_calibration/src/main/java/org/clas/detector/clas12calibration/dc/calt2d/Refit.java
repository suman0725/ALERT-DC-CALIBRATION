/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import org.jlab.rec.dc.cluster.Cluster;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.hit.FittedHit;

/**
 *
 * @author ziegler
 */
public class Refit {
    
    public List<FittedHit> hits;
    public Refit(List<FittedHit> fhits) {
        hits = new ArrayList<FittedHit>();
        hits.addAll(fhits);
    }
    public List<FittedCluster> recomposeClusters(List<FittedHit> fhits) {
        Map<Integer, ArrayList<FittedHit>> grpHits = new HashMap<Integer, ArrayList<FittedHit>>();
        List<FittedCluster> clusters = new ArrayList<FittedCluster>();
        for (FittedHit hit : fhits) { 
            if (hit.get_AssociatedClusterID() == -1 || hit.get_AssociatedHBTrackID() == -1) {
                continue;
            }
            if (hit.get_AssociatedClusterID() != -1 &&
                    hit.get_AssociatedHBTrackID() != -1) {
                int index = hit.get_AssociatedHBTrackID()*10000+hit.get_AssociatedClusterID();
               
                if(grpHits.get(index)==null) { // if the list not yet created make it
                    grpHits.put(index, new ArrayList<FittedHit>()); 
                    grpHits.get(index).add(hit); // append hit
                    //System.out.println("appended first hit "+hit.get_Sector()+","+hit.get_Superlayer()+", "+hit.get_Layer()+","+hit.get_TDC()+
                    //        " to cluster "+index);
                } else {
                    grpHits.get(index).add(hit); // append hit
                    //System.out.println("appended subs hit "+hit.get_Sector()+","+hit.get_Superlayer()+", "+hit.get_Layer()+","+hit.get_TDC()+
                    //        " to cluster "+index);
                }
            }
        }
        Iterator<Map.Entry<Integer, ArrayList<FittedHit>>> itr = grpHits.entrySet().iterator(); 
          
        while(itr.hasNext()) {
            Map.Entry<Integer, ArrayList<FittedHit>> entry = itr.next(); 
             
            if(entry.getValue().size()>3) {
                Cluster cluster = new Cluster(entry.getValue().get(0).get_Sector(), 
                        entry.getValue().get(0).get_Superlayer(), entry.getValue().get(0).get_AssociatedClusterID());
                FittedCluster fcluster = new FittedCluster(cluster);
                
                fcluster.addAll(entry.getValue());
                clusters.add(fcluster);
            }
        }
        return clusters;
    }
    private ClusterFitter cf = new ClusterFitter();
    public void reFit() {
        List<FittedCluster> clusters = this.recomposeClusters(hits);
        for(FittedCluster clus : clusters) {
            cf.SetFitArray(clus, "TSC");
            cf.Fit(clus, true);
            cf.SetResidualDerivedParams(clus, true, false, T2DViewer.dcDetector); //calcTimeResidual=false, resetLRAmbig=false 
            
         /*   double trkAngle = clus.get_clusterLineFitSlope();
            //local angle updated
            double alpha = Math.toDegrees(Math.atan(trkAngle));
            for(FittedHit h : clus) {
                h.setAlpha(alpha);
            }*/
        }
        clusters.clear();
    }
    
    
    
}
