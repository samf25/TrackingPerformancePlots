#include "SimplePlotTemplate.h"

// Load EDM4hep
R__LOAD_LIBRARY(libedm4hep)
#include "edm4hep/TrackerHitPlaneData.h"
#include "edm4hep/Vector3d.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

void HitDensity_MuCollStyle(const char* filename) {
    // Initialize Muon Collider style
    MuCollStyle::InitializeStyle();
    
    if (!filename) {
        printf("Usage: HitDensity_MuCollStyle(\"<root_file>\")\n");
        return;
    }
    
    // Open the ROOT file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        printf("Error: Could not open file %s\n", filename);
        return;
    }
    
    // Get the event TTree
    TTree* tree = (TTree*)file->Get("events");
    if (!tree) {
        printf("Error: Could not find 'events' TTree in file\n");
        file->Close();
        return;
    }

    // Set up branches for different hit collections
    std::vector<std::vector<edm4hep::TrackerHitPlaneData>*> hitCollections(6, nullptr);
    const char* branchNames[] = {"VXDBarrelHits", "VXDEndcapHits", "ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits"};
    const int layerMapping[] = {0, 0, 1, 1, 2, 2}; // Map to layers (VXD=0, IT=1, OT=2)
    
    for (int i = 0; i < 6; i++) {
        tree->SetBranchAddress(branchNames[i], &hitCollections[i]);
    }
    
    // Define theta bins
    const int nThetaBins = 20;
    const double thetaMin = 0.0;
    const double thetaMax = 3.14159;
    const double thetaBinWidth = (thetaMax - thetaMin) / nThetaBins;
    
    // Create histograms for each layer
    const char* layerNames[] = {"VXD", "IT", "OT"};
    const char* layerLabels[] = {"Vertex Detector", "Inner Tracker", "Outer Tracker"};
    TH1D* histograms[3];
    
    for (int layer = 0; layer < 3; layer++) {
        histograms[layer] = new TH1D(layerNames[layer], "", nThetaBins, thetaMin, thetaMax);
    }
    
    // Arrays to store hit counts per theta bin for each layer
    int hitCounts[3][nThetaBins] = {0};
    
    Long64_t nEntries = tree->GetEntries();
    printf("Processing %lld events...\n", nEntries);
    
    // Process events
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Process all hit collections
        for (int collIdx = 0; collIdx < 6; collIdx++) {
            if (hitCollections[collIdx]) {
                int layer = layerMapping[collIdx];
                for (const auto &hit : *hitCollections[collIdx]) {
                    edm4hep::Vector3d position = hit.position;
                    double theta = atan2(sqrt(position.x*position.x + position.y*position.y), position.z);
                    
                    int thetaBin = (int)((theta - thetaMin) / thetaBinWidth);
                    if (thetaBin >= 0 && thetaBin < nThetaBins) {
                        hitCounts[layer][thetaBin]++;
                    }
                }
            }
        }
    }
    
    // Fill histograms with hit density values
    for (int layer = 0; layer < 3; layer++) {
        for (int bin = 0; bin < nThetaBins; bin++) {
            double thetaLow = thetaMin + bin * thetaBinWidth;
            double thetaHigh = thetaMin + (bin + 1) * thetaBinWidth;
            double angularArea = 2 * 3.14159 * (cos(thetaLow) - cos(thetaHigh));
            
            double density = hitCounts[layer][bin] / angularArea;
            histograms[layer]->SetBinContent(bin + 1, density);
        }
    }
    
    // Create output ROOT file
    TFile* outputFile = new TFile("hit_density_results_styled.root", "RECREATE");
    
    // Create individual plots for each layer
    for (int layer = 0; layer < 3; layer++) {
        TCanvas* canvas = MuCollStyle::CreateCanvas(Form("c_%s", layerNames[layer]), 
                                                  Form("%s Hit Density", layerLabels[layer]));
        
        // Style the histogram
        MuCollStyle::StyleHist(histograms[layer], MuCollStyle::GetColor(layer));
        histograms[layer]->GetXaxis()->SetTitle("#theta [rad]");
        histograms[layer]->GetYaxis()->SetTitle("Hit Density [hits/cm^{2}]");
        
        histograms[layer]->Draw("PE");
        
        // Add standard labels
        MuCollStyle::AddStandardLabels(canvas, "10 TeV");
        
        // Create legend for individual plot
        std::vector<TH1*> histVector = {histograms[layer]};
        TLegend* legend = MuCollStyle::CreateSmartLegend(histVector);
        legend->AddEntry(histograms[layer], layerLabels[layer], "pe");
        legend->Draw();
        // Write canvas to ROOT file
        canvas->Write();
        
        // Write histogram to ROOT file
        histograms[layer]->Write();
    }
    
    // Create combined plot
    TCanvas* combinedCanvas = MuCollStyle::CreateCanvas("c_combined", "Hit Density Comparison");
    
    // Find maximum for scaling
    double maxVal = 0;
    for (int layer = 0; layer < 3; layer++) {
        double layerMax = histograms[layer]->GetMaximum();
        if (layerMax > maxVal) maxVal = layerMax;
    }
    
    // Draw all histograms
    for (int layer = 0; layer < 3; layer++) {
        MuCollStyle::StyleHist(histograms[layer], MuCollStyle::GetColor(layer));
        histograms[layer]->GetXaxis()->SetTitle("#theta [rad]");
        histograms[layer]->GetYaxis()->SetTitle("Hit Density [hits/cm^{2}]");
        histograms[layer]->SetMaximum(maxVal * 1.1);
        
        if (layer == 0) {
            histograms[layer]->Draw("PE");
        } else {
            histograms[layer]->Draw("PE SAME");
        }
    }
    
    // Create legend
    std::vector<TH1*> histVector(histograms, histograms + 3);
    TLegend* legend = MuCollStyle::CreateSmartLegend(histVector);
    for (int layer = 0; layer < 3; layer++) {
        legend->AddEntry(histograms[layer], layerLabels[layer], "pe");
    }
    legend->Draw();
    
    // Add standard labels
    MuCollStyle::AddStandardLabels(combinedCanvas, "10 TeV");
    
    // Write combined canvas to ROOT file
    combinedCanvas->Write();
    
    outputFile->Close();
    delete outputFile;
    
    printf("Muon Collider styled hit density plots created!\n");
    printf("All plots and histograms saved to: hit_density_results_styled.root\n");

    file->Close();
}
