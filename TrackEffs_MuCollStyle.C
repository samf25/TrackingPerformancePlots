#include "SimplePlotTemplate.h"

void run_collection(const char* filename, const char* branchName, bool isSubsetCollection = false) {
    // Initialize Muon Collider style
    MuCollStyle::InitializeStyle();
    
    if (!filename) {
        printf("Usage: TrackEffs_MuCollStyle(\"<root_file>\")\n");
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
    
    // Set up branches for MC particles, tracks, track states, and relations
    std::vector<edm4hep::TrackState>* trackStates = nullptr;
    std::vector<edm4hep::MCParticleData>* mcParticles = nullptr;
    std::vector<edm4hep::TrackData>* tracks = nullptr;
    std::vector<podio::ObjectID>* toRelations = nullptr;
    std::vector<podio::ObjectID>* fromRelations = nullptr;
    std::vector<podio::ObjectID>* subsetTrackIndices = nullptr;
    
    tree->SetBranchAddress("MCParticles", &mcParticles);
    if (isSubsetCollection) {
        tree->SetBranchAddress("AllTracks", &tracks);
        tree->SetBranchAddress("_AllTracks_trackStates", &trackStates);
        tree->SetBranchAddress((std::string(branchName)+"s_objIdx").c_str(), &subsetTrackIndices);
    } else {
        tree->SetBranchAddress(branchName, &tracks);
    }
    tree->SetBranchAddress(("_" + std::string(branchName) + "Relations_to").c_str(), &toRelations);
    tree->SetBranchAddress(("_" + std::string(branchName) + "Relations_from").c_str(), &fromRelations);
    
    // Define binning for efficiency plots
    const int nPtBins = 20;
    const double ptMin = 0.5;
    const double ptMax = 110.0;
    
    const int nThetaBins = 20;
    const double thetaMin = 0.0;
    const double thetaMax = 3.14159;
    
    const double matchProb = 0.5; // Minimum matching probability
    
    // Create efficiency histograms
    TH1D* h_allTruths_pt = new TH1D("allTruths_pt", "", nPtBins, ptMin, ptMax);
    TH1D* h_realTruths_pt = new TH1D("realTruths_pt", "", nPtBins, ptMin, ptMax);
    
    TH1D* h_allTruths_theta = new TH1D("allTruths_theta", "", nThetaBins, thetaMin, thetaMax);
    TH1D* h_realTruths_theta = new TH1D("realTruths_theta", "", nThetaBins, thetaMin, thetaMax);
    
    TH1D* h_allTracks = new TH1D("allTracks", "", nPtBins, ptMin, ptMax);
    TH1D* h_realTracks = new TH1D("realTracks", "", nPtBins, ptMin, ptMax);
    TH1D* h_fakeTracks = new TH1D("fakeTracks", "", nPtBins, ptMin, ptMax);
    
    TH1D* h_numberOfTracks = new TH1D("numberOfTracks", "", 100, 0, 1000);

    // Create resolution histograms
    TH1D* h_resolutions_q_over_pt = new TH1D("resolutions_q/pt", "", 100, -0.5, 0.5);
    TH1D* h_resolutions_d0 = new TH1D("resolutions_d0", "", 100, -1., 1.);
    TH1D* h_resolutions_z0 = new TH1D("resolutions_z0", "", 100, -1., 1.);

    // Create Track Data Histograms
    TH1D* h_allTracks_nHits = new TH1D("allTracks_nHits", "", 20, 0, 20);
    TH1D* h_allTracks_nHoles = new TH1D("allTracks_nHoles", "", 10, 0, 10);
    TH1D* h_allTracks_chi2ndof = new TH1D("allTracks_chi2ndof", "", 50, 0, 10);

    TH1D* h_realTracks_nHits = new TH1D("realTracks_nHits", "", 20, 0, 20);
    TH1D* h_realTracks_nHoles = new TH1D("realTracks_nHoles", "", 10, 0, 10);
    TH1D* h_realTracks_chi2ndof = new TH1D("realTracks_chi2ndof", "", 50, 0, 10);

    Long64_t nEntries = tree->GetEntries();
    printf("Processing %lld events for %s...\n", nEntries, branchName);
    
    // Process events
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Select acceptable MC particles
        std::vector<edm4hep::MCParticleData> mcpSet;
        if (mcParticles) {
            for (const auto& mcp : *mcParticles) {
                if (mcp.generatorStatus != 1) { continue; }
                if (mcp.charge == 0) { continue; }
              //if (mcp.isDecayedInTracker) { continue; }
                
                // Tracker Acceptance
                const edm4hep::Vector3d& mom = mcp.momentum;
                double pt = std::sqrt(std::pow(mom.x, 2) + std::pow(mom.y, 2));
                std::cout << pt << std::endl;
                double lambda = std::atan2(mom.z, pt);
                
                //if (fabs(lambda) > 75. / 180 * M_PI) { continue; }
                
                mcpSet.push_back(mcp);
                
                // Fill all truth histograms
                h_allTruths_pt->Fill(pt);
                double theta = std::atan2(std::sqrt(mom.x*mom.x + mom.y*mom.y), mom.z);
                h_allTruths_theta->Fill(theta);
            }
        }
        
        // Process tracks
        std::vector<edm4hep::TrackData> trkSet;
        if (tracks) {
            if (isSubsetCollection) {
            printf("Track Collection Size: %lu\n", tracks->size());
            // Subset collection: build track set from indices
            for (const auto& idxObj : *subsetTrackIndices) {
                const auto& trk = (*tracks)[idxObj.index];
                trkSet.push_back(trk);
                
                // Calculate track pT (approximate from track parameters)
                const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
                double trackPt = fabs(0.3 * 3.57 / firstTrackState.omega / 1000);
                h_allTracks->Fill(trackPt);
                
                // Fill track quality histograms
                h_allTracks_nHits->Fill(trk.trackerHits_end - trk.trackerHits_begin);
                h_allTracks_nHoles->Fill(trk.Nholes);
                if (trk.ndf > 0) {
                h_allTracks_chi2ndof->Fill(trk.chi2 / trk.ndf);
                }
            }
            h_numberOfTracks->Fill(trkSet.size());
            } else {
            printf("Track Collection Size: %lu\n", tracks->size());
            for (const auto& trk : *tracks) {
                trkSet.push_back(trk);
                
                // Calculate track pT (approximate from track parameters)
                const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
                double trackPt = fabs(0.3 * 3.57 / firstTrackState.omega / 1000);
                h_allTracks->Fill(trackPt);
                
                // Fill track quality histograms
                h_allTracks_nHits->Fill(trk.trackerHits_end - trk.trackerHits_begin);
                h_allTracks_nHoles->Fill(trk.Nholes);
                if (trk.ndf > 0) {
                h_allTracks_chi2ndof->Fill(trk.chi2 / trk.ndf);
                }
            }
            h_numberOfTracks->Fill(trkSet.size());
            }
        }
        
        // Process MC Relations and save matched objects
        if (toRelations) {
            for (size_t i = 0; i < toRelations->size(); ++i) {
                const edm4hep::MCParticleData& mcpObj = mcParticles->at((*toRelations)[i].index);
                const edm4hep::TrackData& trkObj = tracks->at((*fromRelations)[i].index);
                
                // Look for mcpObj in mcpSet
                auto itMC = std::find_if(mcpSet.begin(), mcpSet.end(), 
                    [&mcpObj](const edm4hep::MCParticleData& obj) { 
                        // Compare by momentum and charge as proxy for equality
                        return (obj.momentum.x == mcpObj.momentum.x && 
                               obj.momentum.y == mcpObj.momentum.y && 
                               obj.momentum.z == mcpObj.momentum.z &&
                               obj.charge == mcpObj.charge);
                    });
                
                if (itMC != mcpSet.end()) { // Truth particle accepted
                    auto itTRK = std::find_if(trkSet.begin(), trkSet.end(), 
                        [&trkObj, &trackStates](const edm4hep::TrackData& obj) {
                            // Compare by track parameters as proxy for equality
                            if (trackStates && obj.trackStates_begin < trackStates->size() && 
                                trkObj.trackStates_begin < trackStates->size()) {
                                const auto& objState = (*trackStates)[obj.trackStates_begin];
                                const auto& trkState = (*trackStates)[trkObj.trackStates_begin];
                                return (objState.omega == trkState.omega &&
                                       objState.phi == trkState.phi &&
                                       objState.tanLambda == trkState.tanLambda);
                            }
                            return false;
                    });
                        
                    if (itTRK != trkSet.end()) {
                        // Calculate quantities for matched particles
                        const edm4hep::Vector3d& mom = mcpObj.momentum;
                        double pt = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                        double theta = std::atan2(std::sqrt(mom.x*mom.x + mom.y*mom.y), mom.z);
                        const auto& firstTrackState = (*trackStates)[trkObj.trackStates_begin];
                        double trackPt = fabs(0.3 * 3.57 / firstTrackState.omega / 1000);
                        double trackD0 = firstTrackState.D0;
                        double trackZ0 = firstTrackState.Z0;
                        
                        // Fill matched histograms
                        h_realTruths_pt->Fill(pt);
                        h_realTruths_theta->Fill(theta);
                        h_realTracks->Fill(trackPt);
                        
                        // Fill resolution histograms
                        double true_q_over_pt = (mcpObj.charge) / pt;
                        double reco_q_over_pt = (mcpObj.charge) / trackPt;
                        h_resolutions_q_over_pt->Fill((reco_q_over_pt - true_q_over_pt) / true_q_over_pt);
                        // Calculate MC particle d0 and z0 from vertex
                        double mcpD0 = std::sqrt(mcpObj.vertex.x * mcpObj.vertex.x + mcpObj.vertex.y * mcpObj.vertex.y);
                        double mcpZ0 = mcpObj.vertex.z;
                        h_resolutions_d0->Fill(trackD0 - mcpD0);
                        h_resolutions_z0->Fill(trackZ0 - mcpZ0);

                        // Fill track quality histograms for real tracks
                        h_realTracks_nHits->Fill(trkObj.trackerHits_end - trkObj.trackerHits_begin);
                        h_realTracks_nHoles->Fill(trkObj.Nholes);
                        if (trkObj.ndf > 0) {
                            h_realTracks_chi2ndof->Fill(trkObj.chi2 / trkObj.ndf);
                        }

                        trkSet.erase(itTRK);
                    }
                }
            }
        }
        
        // Save unmatched tracks as fake tracks
        for (const auto& trk : trkSet) {
            const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
            double trackPt = fabs(0.3 * 3.57 / firstTrackState.omega / 1000);
            h_fakeTracks->Fill(trackPt);
        }
    }
    
    // Create efficiency plots
    TEfficiency* eff_pt = new TEfficiency(*h_realTruths_pt, *h_allTruths_pt);
    TEfficiency* eff_theta = new TEfficiency(*h_realTruths_theta, *h_allTruths_theta);
    TEfficiency* fake_rate = new TEfficiency(*h_fakeTracks, *h_allTracks);
    
    // Create and style efficiency vs pT plot
    TCanvas* c_eff_pt = MuCollStyle::CreateCanvas("c_eff_pt", "Tracking Efficiency vs pT");
    TH1D background("background", "", 1, 0, 110);
    background.SetMaximum(1.3);
    background.SetMinimum(0.0);
    background.GetXaxis()->SetTitle("p_{T} [GeV]");
    background.GetYaxis()->SetTitle("Tracking Efficiency");
    background.Draw();
    
    MuCollStyle::StyleEff(eff_pt, MuCollStyle::GetColor(0));
    eff_pt->Draw("Same");
    
    MuCollStyle::AddStandardLabels(c_eff_pt, "10 TeV");
    
    // Create and style efficiency vs theta plot
    TCanvas* c_eff_theta = MuCollStyle::CreateCanvas("c_eff_theta", "Tracking Efficiency vs #theta");
    
    TH1D background_theta("background_theta", "", 1, 0, 3.14159);
    background_theta.SetMaximum(1.3);
    background_theta.SetMinimum(0.0);
    background_theta.GetXaxis()->SetTitle("#theta [rad]");
    background_theta.GetYaxis()->SetTitle("Tracking Efficiency");
    background_theta.Draw();
    
    MuCollStyle::StyleEff(eff_theta, MuCollStyle::GetColor(1));
    eff_theta->Draw("Same");
    
    MuCollStyle::AddStandardLabels(c_eff_theta, "10 TeV");
    
    // Create fake rate plot
    TCanvas* c_fake = MuCollStyle::CreateCanvas("c_fake", "Fake Rate vs pT");
    
    TH1D background_fake("background_fake", "", 1, 0, 110);
    background_fake.SetMaximum(1.3);
    background_fake.SetMinimum(0.0);
    background_fake.GetXaxis()->SetTitle("p_{T} [GeV]");
    background_fake.GetYaxis()->SetTitle("Fake Rate");
    background_fake.Draw();
    
    MuCollStyle::StyleEff(fake_rate, MuCollStyle::GetColor(2));
    fake_rate->Draw("Same");
    
    // Create number of tracks histogram
    TCanvas* c_ntracks = MuCollStyle::CreateCanvas("c_ntracks", "Number of Tracks");
    MuCollStyle::StyleHist(h_numberOfTracks, MuCollStyle::GetColor(3));
    h_numberOfTracks->GetXaxis()->SetTitle("Number of Tracks");
    h_numberOfTracks->GetYaxis()->SetTitle("Events");
    h_numberOfTracks->Draw("HIST");
    
    MuCollStyle::AddStandardLabels(c_ntracks, "10 TeV");

    // Create resolution canvases
    TCanvas* c_res_pt = MuCollStyle::CreateCanvas("c_res_pt", "Resolution of q/p_{T}");
    MuCollStyle::StyleHist(h_resolutions_q_over_pt, MuCollStyle::GetColor(4));
    h_resolutions_q_over_pt->GetXaxis()->SetTitle("Reconstructed q/p_{T} - True q/p_{T}");
    h_resolutions_q_over_pt->GetYaxis()->SetTitle("Entries");
    h_resolutions_q_over_pt->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_pt, "10 TeV"); 
    TCanvas* c_res_d0 = MuCollStyle::CreateCanvas("c_res_d0", "Resolution of d_{0}");
    MuCollStyle::StyleHist(h_resolutions_d0, MuCollStyle::GetColor(5));
    h_resolutions_d0->GetXaxis()->SetTitle("Reconstructed d - True d_{0} [mm]");
    h_resolutions_d0->GetYaxis()->SetTitle("Entries");
    h_resolutions_d0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_d0, "10 TeV");
    TCanvas* c_res_z0 = MuCollStyle::CreateCanvas("c_res_z0", "Resolution of z_{0}");
    MuCollStyle::StyleHist(h_resolutions_z0, MuCollStyle::GetColor(6));
    h_resolutions_z0->GetXaxis()->SetTitle("Reconstructed z - True z_{0} [mm]");
    h_resolutions_z0->GetYaxis()->SetTitle("Entries");
    h_resolutions_z0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_z0, "10 TeV");

    // Create Track Data quality canvases
    TCanvas* c_nHits = MuCollStyle::CreateCanvas("c_nHits", "Number of Hits");
    MuCollStyle::StyleHist(h_allTracks_nHits, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_realTracks_nHits, MuCollStyle::GetColor(1));
    h_allTracks_nHits->GetXaxis()->SetTitle("Number of Hits");
    h_allTracks_nHits->GetYaxis()->SetTitle("Entries");
    h_allTracks_nHits->Draw("HIST");
    h_realTracks_nHits->Draw("HIST SAME");
    MuCollStyle::AddStandardLabels(c_nHits, "10 TeV");

    TCanvas* c_nHoles = MuCollStyle::CreateCanvas("c_nHoles", "Number of Holes");
    MuCollStyle::StyleHist(h_allTracks_nHoles, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_realTracks_nHoles, MuCollStyle::GetColor(1));
    h_allTracks_nHoles->GetXaxis()->SetTitle("Number of Holes");
    h_allTracks_nHoles->GetYaxis()->SetTitle("Entries");
    h_allTracks_nHoles->Draw("HIST");
    h_realTracks_nHoles->Draw("HIST SAME");
    MuCollStyle::AddStandardLabels(c_nHoles, "10 TeV");

    TCanvas* c_chi2ndof = MuCollStyle::CreateCanvas("c_chi2ndof", "#chi^{2}/ndof");
    MuCollStyle::StyleHist(h_allTracks_chi2ndof, MuCollStyle::GetColor(0));
    MuCollStyle::StyleHist(h_realTracks_chi2ndof, MuCollStyle::GetColor(1));
    h_allTracks_chi2ndof->GetXaxis()->SetTitle("#chi^{2}/ndof");
    h_allTracks_chi2ndof->GetYaxis()->SetTitle("Entries");
    h_allTracks_chi2ndof->Draw("HIST");
    h_realTracks_chi2ndof->Draw("HIST SAME");
    MuCollStyle::AddStandardLabels(c_chi2ndof, "10 TeV");
    
    // Create output ROOT file and save histograms
    TFile* outputFile = new TFile("track_efficiency_stylized.root", "RECREATE");
    
    
    // Save canvases
    c_eff_pt->Write();
    c_eff_theta->Write();
    c_fake->Write();
    c_ntracks->Write();
    c_res_pt->Write();
    c_res_d0->Write();
    c_res_z0->Write();
    c_nHits->Write();
    c_nHoles->Write();
    c_chi2ndof->Write();
    
    outputFile->Close();
    delete outputFile;

    file->Close();
}

void TrackEffs_MuCollStyle(const char* filename) {
    run_collection(filename, "SiTrack", true);
}
