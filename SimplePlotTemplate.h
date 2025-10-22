// Simple Muon Collider Plot Template - Header Only
// Include this at the top of your ROOT macros for consistent styling

#ifndef SIMPLE_PLOT_TEMPLATE_H
#define SIMPLE_PLOT_TEMPLATE_H

#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TROOT.h"

namespace MuCollStyle {

// Initialize Muon Collider style (call once at beginning of macro)
void InitializeStyle() {
    gROOT->SetStyle("Plain");
    
    // Canvas settings
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetCanvasDefH(600);
    gStyle->SetCanvasDefW(800);
    
    // Pad settings
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(kWhite);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    
    // Frame settings
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(kBlack);
    gStyle->SetFrameLineWidth(1);
    
    // Axis settings
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    // Label settings
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    
    // Title settings
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.2);
    
    // Statistics and title
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // Margins
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.05);
    
    gROOT->ForceStyle();
}

// Color palette
int GetColor(int index) {
    int colors[] = {kBlue+1, kRed+1, kGreen+2, kOrange+1, kMagenta+1, 
                   kCyan+1, kYellow+2, kViolet+1, kSpring-1, kAzure+2};
    return colors[index % 10];
}

// Style a histogram
void StyleHist(TH1* hist, int color = -1, int markerStyle = 20, double markerSize = 0.8) {
    if (!hist) return;
    
    if (color == -1) color = GetColor(0);
    
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
    hist->SetLineWidth(2);
    hist->SetTitle("");
    
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(1.3);
}

// Style a TEfficiency
void StyleEff(TEfficiency* eff, int color = -1, int markerStyle = 20, double markerSize = 0.8) {
    if (!eff) return;
    
    if (color == -1) color = GetColor(0);
    
    eff->SetLineColor(color);
    eff->SetMarkerColor(color);
    eff->SetMarkerStyle(markerStyle);
    eff->SetMarkerSize(markerSize);
    eff->SetLineWidth(2);
}

// Create a styled canvas
TCanvas* CreateCanvas(const char* name, const char* title, int width = 800, int height = 600, 
                     bool logY = false, bool logX = false) {
    TCanvas* canvas = new TCanvas(name, title, width, height);
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameBorderMode(0);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    
    if (logY) canvas->SetLogy();
    if (logX) canvas->SetLogx();
    
    return canvas;
}

// Find optimal legend position to avoid data overlap
std::pair<double, double> FindOptimalLegendPosition(std::vector<TH1*> hists, double legendWidth = 0.25, double legendHeight = 0.25) {
    if (hists.empty()) return std::make_pair(0.65, 0.65);
    
    // Define possible legend positions (x1, y1) - avoiding top left where label is
    std::vector<std::pair<double, double>> positions = {
        {0.7, 0.7},   // Top right (default)
        {0.7, 0.2},   // Bottom right
        {0.2, 0.2},   // Bottom left
        {0.45, 0.7},  // Top center
        {0.7, 0.45},  // Middle right
        {0.45, 0.2}   // Bottom center
    };
    
    double minOverlap = 1e9;
    std::pair<double, double> bestPosition = {0.7, 0.7};
    
    for (const auto& pos : positions) {
        double x1 = pos.first;
        double y1 = pos.second;
        double x2 = x1 + legendWidth;
        double y2 = y1 + legendHeight;
        
        // Check if legend goes outside canvas
        if (x2 > 0.95 || y2 > 0.95 || x1 < 0.05 || y1 < 0.05) continue;
        
        // Check if legend overlaps with label area (top left region)
        // Label is at x=0.2, y=0.87, extends roughly to x=0.6, y=0.75
        if ((x1 < 0.6 && x2 > 0.15) && (y1 < 0.9 && y2 > 0.72)) continue;
        
        // Calculate overlap with all histogram data
        double overlap = 0;
        
        for (TH1* hist : hists) {
            if (!hist) continue;
            
            int nBins = hist->GetNbinsX();
            double histMax = hist->GetMaximum();
            double histMin = hist->GetMinimum();
            
            for (int i = 1; i <= nBins; i++) {
                double binCenter = hist->GetBinCenter(i);
                double binContent = hist->GetBinContent(i);
                
                if (binContent <= histMin) continue;
                
                // Convert bin position to NDC coordinates
                double xNDC = (binCenter - hist->GetXaxis()->GetXmin()) / 
                             (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin());
                double yNDC = (binContent - histMin) / (histMax - histMin);
                
                // Apply pad margins
                xNDC = 0.16 + xNDC * (0.95 - 0.16);
                yNDC = 0.13 + yNDC * (0.95 - 0.13);
                
                // Check if data point is inside legend area
                if (xNDC >= x1 && xNDC <= x2 && yNDC >= y1 && yNDC <= y2) {
                    overlap += binContent;
                }
            }
        }
        
        if (overlap < minOverlap) {
            minOverlap = overlap;
            bestPosition = pos;
        }
    }
    
    return bestPosition;
}

// Create an intelligently positioned legend
TLegend* CreateSmartLegend(std::vector<TH1*> hists, double width = 0.25, double height = 0.25) {
    std::pair<double, double> position;
    
    if (!hists.empty()) {
        position = FindOptimalLegendPosition(hists, width, height);
    } else {
        position = std::make_pair(0.65, 0.65);
    }
    
    TLegend* legend = new TLegend(position.first, position.second, 
                                 position.first + width, position.second + height);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    return legend;
}

// Create a styled legend
TLegend* CreateLegend(double x1 = 0.65, double y1 = 0.65, double x2 = 0.9, double y2 = 0.9) {
    TLegend* legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    return legend;
}

// Add Muon Collider label
void AddLabel(TCanvas* canvas, const char* energy = "10 TeV", 
              double x = 0.2, double y = 0.87, double textSize = 0.04) {
    canvas->cd();
    
    // Combined "Muon Collider" and "Simulation" labels with different styles
    TLatex* latex1 = new TLatex();
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(textSize);
    latex1->SetTextColor(kBlack);
    latex1->DrawLatex(x, y, "#bf{Muon Collider} Simulation");
    
    // Energy label with geometry on same line
    TLatex* latex2 = new TLatex();
    latex2->SetNDC();
    latex2->SetTextFont(42);
    latex2->SetTextSize(textSize * 0.8);
    latex2->SetTextColor(kBlack);
    latex2->DrawLatex(x, y - 0.06, Form("#sqrt{s} = %s, MAIA Geometry", energy));
}

// Convenient function to add all standard labels
void AddStandardLabels(TCanvas* canvas, const char* energy = "10 TeV", const char* lumiText = "Simulation") {
    AddLabel(canvas, energy);
}

// Save plot in multiple formats
void SavePlot(TCanvas* canvas, const char* filename, const char* formats = "pdf,png") {
    if (!canvas) return;
    
    std::string formatStr(formats);
    std::string format;
    size_t pos = 0;
    
    while ((pos = formatStr.find(",")) != std::string::npos) {
        format = formatStr.substr(0, pos);
        canvas->SaveAs(Form("%s.%s", filename, format.c_str()));
        formatStr.erase(0, pos + 1);
    }
    // Save the last format
    if (!formatStr.empty()) {
        canvas->SaveAs(Form("%s.%s", filename, formatStr.c_str()));
    }
}

} // namespace MuCollStyle

#endif
