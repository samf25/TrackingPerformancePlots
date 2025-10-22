#include "SimplePlotTemplate.h"
#include "DDSegmentation/BitFieldCoder.h"
#include <unordered_map>


void run_collection(const char* filename, const char* branchName, bool isSubsetCollection = false) {
    // Initialize Muon Collider style
    MuCollStyle::InitializeStyle();
    dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder("system:5,side:-2,layer:6,module:11,sensor:8");
    
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

    std::vector<std::string> trackerHitCollections = {
        "ITBarrelHits", "ITEndcapHits", "VXDBarrelHits", 
        "VXDEndcapHits", "OTBarrelHits", "OTEndcapHits"
    };    
    std::vector<std::string> simTrackerHitCollections = {
        "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", 
        "VertexBarrelCollection", "VertexEndcapCollection", 
        "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"
    };

    // Use podio metadata to get the true numeric collectionIDs for each collection
    TTree* meta = static_cast<TTree*>(file->Get("podio_metadata"));
    std::vector<podio::root_utils::CollectionWriteInfo>* infos = nullptr;
    std::unordered_map<std::string, uint32_t> name2id;
    if (meta) {
        meta->SetBranchAddress("events___CollectionTypeInfo", &infos);
        meta->GetEntry(0);
        if (infos) {
            for (const auto& ci : *infos) {
                // ci.name is the collection name; ci.id is the numeric collectionID used in ObjectID
                name2id[ci.name] = ci.collectionID;
            }
        }
    }
    // Build list of hit collection IDs and a fast id->slot mapping aligned with trackerHitCollections order
    std::vector<uint32_t> trackerHitIDs;
    std::unordered_map<uint32_t, size_t> idToSlot;
    trackerHitIDs.reserve(trackerHitCollections.size());
    for (size_t slot = 0; slot < trackerHitCollections.size(); ++slot) {
        const auto& cname = trackerHitCollections[slot];
        auto it = name2id.find(cname);
        if (it != name2id.end()) {
            trackerHitIDs.push_back(it->second);
            idToSlot[it->second] = slot;
        } else {
            // collection absent in this file; mark as -1 to keep alignment
            trackerHitIDs.push_back(-1);
        }
    }

    // Set up branches for MC particles, tracks, track states, and relations
    std::vector<edm4hep::MCParticleData>* mcParticles = nullptr;
    std::vector<edm4hep::TrackData>* tracks = nullptr;
    std::vector<edm4hep::TrackState>* trackStates = nullptr;
    std::vector<podio::ObjectID>* seedTrackerHits = nullptr;
    std::vector<std::vector<edm4hep::TrackerHitPlaneData>*> trackerHits(trackerHitCollections.size(), nullptr);
    std::vector<std::vector<podio::ObjectID>*> toRelations(trackerHitCollections.size(), nullptr);
    std::vector<std::vector<podio::ObjectID>*> fromRelations(trackerHitCollections.size(), nullptr);
    std::vector<std::vector<podio::ObjectID>*> simParticles(trackerHitCollections.size(), nullptr);

    tree->SetBranchAddress("MCParticles", &mcParticles);
    tree->SetBranchAddress(branchName, &tracks);
    tree->SetBranchAddress(("_" + std::string(branchName) + "_trackStates").c_str(), &trackStates);
    tree->SetBranchAddress(("_" + std::string(branchName) + "_trackerHits").c_str(), &seedTrackerHits);
    for (size_t i = 0; i < trackerHitCollections.size(); ++i) {
        const auto& collectionName = trackerHitCollections[i];
        const auto& simCollectionName = simTrackerHitCollections[i];
        tree->SetBranchAddress((collectionName).c_str(), &trackerHits[i]);
        tree->SetBranchAddress(("_" + std::string(collectionName) + "Relations_to").c_str(), &toRelations[i]);
        tree->SetBranchAddress(("_" + std::string(collectionName) + "Relations_from").c_str(), &fromRelations[i]);
        tree->SetBranchAddress(("_" + std::string(simCollectionName)+"_particle").c_str(), &simParticles[i]);
    }

    // Define binning for efficiency plots
    const int nPtBins = 20;
    const double ptMin = 0.5;
    const double ptMax = 110.0;
    
    const int nThetaBins = 20;
    const double thetaMin = 0.0;
    const double thetaMax = 3.14159;

    // Create seed information histograms
    TH1D* h_seed_theta = new TH1D("seed_theta", "", nThetaBins, thetaMin, thetaMax);
    TH1D* h_seed_layer_barrel = new TH1D("seed_layer_barrel", "", 10, 0, 10);
    TH1D* h_seed_layer_endcap = new TH1D("seed_layer_endcap", "", 15, 0, 15);
    TH1D* h_seed_layer_barrel_matched = new TH1D("seed_layer_barrel_matched", "", 10, 0, 10);
    TH1D* h_seed_layer_endcap_matched = new TH1D("seed_layer_endcap_matched", "", 15, 0, 15);
    TH1D* h_seed_layer_barrel_unmatched = new TH1D("seed_layer_barrel_unmatched", "", 10, 0, 10);
    TH1D* h_seed_layer_endcap_unmatched = new TH1D("seed_layer_endcap_unmatched", "", 15, 0, 15);
    TH1D* h_seeds_per_MCP = new TH1D("seed_per_MCP", "", 10, 0, 10);
    TH1D* h_seed_number = new TH1D("seed_number", "", 100, 0, 1000);
    TH1D* h_seed_unmatched = new TH1D("seed_unmatched", "", 100, 0, 1000);
    TH1D* h_seed_matched = new TH1D("seed_matched", "", 100, 0, 1000);

    // Create resolution histograms
    TH1D* h_seed_resolutions_q_over_pt = new TH1D("seed_resolutions_q/pt", "", 100, -0.5, 0.5);
    TH1D* h_seed_resolutions_d0 = new TH1D("seed_resolutions_d0", "", 100, -1., 1.);
    TH1D* h_seed_resolutions_z0 = new TH1D("seed_resolutions_z0", "", 100, -1., 1.);

    // Loop over events
    Long64_t nEntries = tree->GetEntries();
    printf("Processing %lld events for %s...\n", nEntries, branchName);
    
    // Process events
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // counters
        int matchedSeeds = 0;
        int unmatchedSeeds = 0;
        int totalSeeds = 0;

        // Process Seeds
        // Track MCP matching counts
        std::unordered_map<unsigned int, int> mcpMatchCount;
        for (const auto& trk : *tracks) {
            totalSeeds++;
            const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
            // Fill seed info histogram
            double tanLambda = firstTrackState.tanLambda;
            double theta = atan2(1.0, tanLambda);
            h_seed_theta->Fill(theta);

            std::vector<int> matchedMCPIndices;
            // Loop over tracker hits for this track
            for (unsigned int hitIdx = trk.trackerHits_begin; hitIdx < trk.trackerHits_end; ++hitIdx) {
                if (!seedTrackerHits || hitIdx >= seedTrackerHits->size()) continue;
                const auto& hitID = (*seedTrackerHits)[hitIdx];
                // Resolve collection slot via metadata-provided ID mapping
                auto itSlot = idToSlot.find(hitID.collectionID);
                if (itSlot == idToSlot.end()) {
                    std::cout << "Warning: Could not find collection ID " << hitID.collectionID 
                              << " in idToSlot mapping" << std::endl;
                    continue;
                }
                const size_t collIdx = itSlot->second;

                // get relations for this hit collection
                const auto* fromVec = (collIdx < fromRelations.size() ? fromRelations[collIdx] : nullptr);
                const auto* toVec   = (collIdx < toRelations.size() ? toRelations[collIdx] : nullptr);
                const auto* mcVec  = (collIdx < simParticles.size() ? simParticles[collIdx] : nullptr);
                if (!fromVec || !toVec || !mcVec) {
                    std::cout << "Warning: Relations vectors are null for collection index " << collIdx 
                              << " (collection ID " << hitID.collectionID << ")" << std::endl;
                    continue;
                }
                
                for (size_t relIdx = 0; relIdx < fromVec->size(); ++relIdx) {
                    const auto& fromOID = fromVec->at(relIdx);
                    if (fromOID.index == hitID.index && fromOID.collectionID == hitID.collectionID) {
                        matchedMCPIndices.push_back(mcVec->at(toVec->at(relIdx).index).index);
                        break;
                    }
                }

            }

            // Check if seed is matched to an MC particle
            std::cout << "Seed has " << matchedMCPIndices.size() << " matched MCP indices." << std::endl;
            if (matchedMCPIndices.size() >= 3) {
                // Count occurrences of each MCP index
                std::unordered_map<unsigned int, int> mcpCount;
                for (unsigned int mcpIdx : matchedMCPIndices) {
                    mcpCount[mcpIdx]++;
                }
                
                // Check if any MCP index appears 2 or more times
                bool isMatched = false;
                unsigned int matchedMCPIndex = 0;
                for (const auto& pair : mcpCount) {
                    if (pair.second >= 2) {
                        isMatched = true;
                        matchedMCPIndex = pair.first;
                        mcpMatchCount[matchedMCPIndex]++;
                        break;
                    }
                }

                // Fill layer histograms based on matching
                for (unsigned int hitIdx = trk.trackerHits_begin; hitIdx < trk.trackerHits_end; ++hitIdx) {
                    if (!seedTrackerHits || hitIdx >= seedTrackerHits->size()) continue;
                    const auto& hitID = (*seedTrackerHits)[hitIdx];
                    auto itSlot = idToSlot.find(hitID.collectionID);
                    if (itSlot == idToSlot.end()) continue;
                    const size_t collIdx = itSlot->second;
                    const auto* hitVec  = (collIdx < trackerHits.size() ? trackerHits[collIdx] : nullptr);
                    if (!hitVec) continue;
                    if (static_cast<size_t>(hitID.index) >= hitVec->size()) continue;
                    // Fill layer histograms
                    int layer = bitFieldCoder.get(hitVec->at(hitID.index).cellID, "layer");
                    h_seed_layer_barrel->Fill(layer);
                    if (isMatched) {
                        h_seed_layer_barrel_matched->Fill(layer);
                    } else {
                        h_seed_layer_barrel_unmatched->Fill(layer);
                    }
                }
                
                if (isMatched) {
                    matchedSeeds++;
                    // Access the matched MCP
                    const auto& matchedMCP = (*mcParticles)[matchedMCPIndex];
                    // Calculate quantities for matched particles
                    const edm4hep::Vector3d& mom = matchedMCP.momentum;
                    double pt = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                    double theta = std::atan2(std::sqrt(mom.x*mom.x + mom.y*mom.y), mom.z);
                    const auto& firstTrackState = (*trackStates)[trk.trackStates_begin];
                    double trackPt = fabs(0.3 * 3.57 / firstTrackState.omega / 1000);
                    double trackD0 = firstTrackState.D0;
                    double trackZ0 = firstTrackState.Z0;

                    // Fill response histograms
                    double true_q_over_pt = (matchedMCP.charge) / pt;
                    double reco_q_over_pt = (matchedMCP.charge) / trackPt;
                    // Debug print statements for resolution histograms
                    // std::cout << "True pT: " << pt << ", Reco pT: " << trackPt << std::endl;
                    // std::cout << "True q/pT: " << true_q_over_pt << ", Reco q/pT: " << reco_q_over_pt << std::endl;
                    // std::cout << "True d0: " << std::sqrt(matchedMCP.vertex.x * matchedMCP.vertex.x + matchedMCP.vertex.y * matchedMCP.vertex.y) 
                    //           << ", Reco d0: " << trackD0 << std::endl;
                    // std::cout << "True z0: " << matchedMCP.vertex.z << ", Reco z0: " << trackZ0 << std::endl;
                    // std::cout << "q/pT resolution: " << (reco_q_over_pt - true_q_over_pt) / true_q_over_pt << std::endl;
                    // std::cout << "d0 resolution: " << trackD0 - std::sqrt(matchedMCP.vertex.x * matchedMCP.vertex.x + matchedMCP.vertex.y * matchedMCP.vertex.y) << std::endl;
                    // std::cout << "z0 resolution: " << trackZ0 - matchedMCP.vertex.z << std::endl;
                    h_seed_resolutions_q_over_pt->Fill((reco_q_over_pt - true_q_over_pt) / true_q_over_pt);
                    h_seed_resolutions_d0->Fill(trackD0 - std::sqrt(matchedMCP.vertex.x * matchedMCP.vertex.x + matchedMCP.vertex.y * matchedMCP.vertex.y));
                    h_seed_resolutions_z0->Fill(trackZ0 - matchedMCP.vertex.z);
                } else {
                    unmatchedSeeds++;
                }
            } else {
                unmatchedSeeds++;
            }
        }
        // Fill matched and unmatched seed histograms
        h_seed_matched->Fill(matchedSeeds);
        h_seed_unmatched->Fill(unmatchedSeeds);
        h_seed_number->Fill(totalSeeds);

        // Fill seeds per MCP histogram with average value
        if (!mcpMatchCount.empty()) {
            double totalCount = 0;
            for (const auto& pair : mcpMatchCount) {
            totalCount += pair.second;
            }
            double average = totalCount / mcpMatchCount.size();
            h_seeds_per_MCP->Fill(average);
        }
    }
    // Create canvas for seed theta distribution
    TCanvas* c_seed_theta = MuCollStyle::CreateCanvas("c_seed_theta", "Seed Theta Distribution");
    MuCollStyle::StyleHist(h_seed_theta, MuCollStyle::GetColor(1));
    h_seed_theta->GetXaxis()->SetTitle("Theta [rad]");
    h_seed_theta->GetYaxis()->SetTitle("Entries");
    h_seed_theta->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_seed_theta, "10 TeV");

    // Create canvas for seed layer barrel distribution
    TCanvas* c_seed_layer_barrel = MuCollStyle::CreateCanvas("c_seed_layer_barrel", "Seed Layer Distribution (Barrel)");
    MuCollStyle::StyleHist(h_seed_layer_barrel, MuCollStyle::GetColor(2));
    h_seed_layer_barrel->GetXaxis()->SetTitle("Layer");
    h_seed_layer_barrel->GetYaxis()->SetTitle("Entries");
    h_seed_layer_barrel->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_seed_layer_barrel, "10 TeV");

    // Create canvas for seed layer endcap distribution
    TCanvas* c_seed_layer_endcap = MuCollStyle::CreateCanvas("c_seed_layer_endcap", "Seed Layer Distribution (Endcap)");
    MuCollStyle::StyleHist(h_seed_layer_endcap, MuCollStyle::GetColor(3));
    h_seed_layer_endcap->GetXaxis()->SetTitle("Layer");
    h_seed_layer_endcap->GetYaxis()->SetTitle("Entries");
    h_seed_layer_endcap->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_seed_layer_endcap, "10 TeV");

    // Create canvas for matched vs unmatched barrel layers
    TCanvas* c_seed_layer_barrel_comparison = MuCollStyle::CreateCanvas("c_seed_layer_barrel_comparison", "Barrel Layer Distribution: Matched vs Unmatched");
    MuCollStyle::StyleHist(h_seed_layer_barrel_matched, MuCollStyle::GetColor(2));
    MuCollStyle::StyleHist(h_seed_layer_barrel_unmatched, MuCollStyle::GetColor(4));
    h_seed_layer_barrel_matched->GetXaxis()->SetTitle("Layer");
    h_seed_layer_barrel_matched->GetYaxis()->SetTitle("Entries");
    h_seed_layer_barrel_matched->Draw("HIST");
    h_seed_layer_barrel_unmatched->Draw("HIST SAME");
    
    std::vector<TH1*> barrelHistVector = {h_seed_layer_barrel_matched, h_seed_layer_barrel_unmatched};
    TLegend* barrelLegend = MuCollStyle::CreateSmartLegend(barrelHistVector);
    barrelLegend->AddEntry(h_seed_layer_barrel_matched, "Matched", "l");
    barrelLegend->AddEntry(h_seed_layer_barrel_unmatched, "Unmatched", "l");
    barrelLegend->Draw();
    
    MuCollStyle::AddStandardLabels(c_seed_layer_barrel_comparison, "10 TeV");

    // Create canvas for matched vs unmatched endcap layers
    TCanvas* c_seed_layer_endcap_comparison = MuCollStyle::CreateCanvas("c_seed_layer_endcap_comparison", "Endcap Layer Distribution: Matched vs Unmatched");
    MuCollStyle::StyleHist(h_seed_layer_endcap_matched, MuCollStyle::GetColor(2));
    MuCollStyle::StyleHist(h_seed_layer_endcap_unmatched, MuCollStyle::GetColor(4));
    h_seed_layer_endcap_matched->GetXaxis()->SetTitle("Layer");
    h_seed_layer_endcap_matched->GetYaxis()->SetTitle("Entries");
    h_seed_layer_endcap_matched->Draw("HIST");
    h_seed_layer_endcap_unmatched->Draw("HIST SAME");
    
    std::vector<TH1*> endcapHistVector = {h_seed_layer_endcap_matched, h_seed_layer_endcap_unmatched};
    TLegend* endcapLegend = MuCollStyle::CreateSmartLegend(endcapHistVector);
    endcapLegend->AddEntry(h_seed_layer_endcap_matched, "Matched", "l");
    endcapLegend->AddEntry(h_seed_layer_endcap_unmatched, "Unmatched", "l");
    endcapLegend->Draw();
    
    MuCollStyle::AddStandardLabels(c_seed_layer_endcap_comparison, "10 TeV");

    // Create canvas for seeds per MCP
    TCanvas* c_seeds_per_MCP = MuCollStyle::CreateCanvas("c_seeds_per_MCP", "Seeds per MC Particle");
    MuCollStyle::StyleHist(h_seeds_per_MCP, MuCollStyle::GetColor(5));
    h_seeds_per_MCP->GetXaxis()->SetTitle("Average Seeds per MCP");
    h_seeds_per_MCP->GetYaxis()->SetTitle("Entries");
    h_seeds_per_MCP->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_seeds_per_MCP, "10 TeV");

    // Create canvas for seed number distribution
    TCanvas* c_seed_number = MuCollStyle::CreateCanvas("c_seed_number", "Number of Seeds per Event");
    MuCollStyle::StyleHist(h_seed_number, MuCollStyle::GetColor(6));
    h_seed_number->GetXaxis()->SetTitle("Number of Seeds");
    h_seed_number->GetYaxis()->SetTitle("Entries");
    h_seed_number->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_seed_number, "10 TeV");

    // Create canvas for matched vs unmatched seeds
    TCanvas* c_seed_matching = MuCollStyle::CreateCanvas("c_seed_matching", "Matched vs Unmatched Seeds");
    MuCollStyle::StyleHist(h_seed_matched, MuCollStyle::GetColor(2));
    MuCollStyle::StyleHist(h_seed_unmatched, MuCollStyle::GetColor(4));
    h_seed_matched->GetXaxis()->SetTitle("Number of Seeds");
    h_seed_matched->GetYaxis()->SetTitle("Entries");
    h_seed_matched->Draw("HIST");
    h_seed_unmatched->Draw("HIST SAME");
    
    std::vector<TH1*> matchingHistVector = {h_seed_matched, h_seed_unmatched};
    TLegend* matchingLegend = MuCollStyle::CreateSmartLegend(matchingHistVector);
    matchingLegend->AddEntry(h_seed_matched, "Matched", "l");
    matchingLegend->AddEntry(h_seed_unmatched, "Unmatched", "l");
    matchingLegend->Draw();
    
    MuCollStyle::AddStandardLabels(c_seed_matching, "10 TeV");


    // Create resolution canvases
    TCanvas* c_res_q_over_pt = MuCollStyle::CreateCanvas("c_res_q_over_pt", "Resolution of q/p_{T}");
    MuCollStyle::StyleHist(h_seed_resolutions_q_over_pt, MuCollStyle::GetColor(4));
    h_seed_resolutions_q_over_pt->GetXaxis()->SetTitle("(Reco q/p_{T} - True q/p_{T}) / True q/p_{T}");
    h_seed_resolutions_q_over_pt->GetYaxis()->SetTitle("Entries");
    h_seed_resolutions_q_over_pt->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_q_over_pt, "10 TeV");

    TCanvas* c_res_d0 = MuCollStyle::CreateCanvas("c_res_d0", "Resolution of d_{0}");
    MuCollStyle::StyleHist(h_seed_resolutions_d0, MuCollStyle::GetColor(3));
    h_seed_resolutions_d0->GetXaxis()->SetTitle("Reco d_{0} - True d_{0} [mm]");
    h_seed_resolutions_d0->GetYaxis()->SetTitle("Entries");
    h_seed_resolutions_d0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_d0, "10 TeV");

    TCanvas* c_res_z0 = MuCollStyle::CreateCanvas("c_res_z0", "Resolution of z_{0}");
    MuCollStyle::StyleHist(h_seed_resolutions_z0, MuCollStyle::GetColor(5));
    h_seed_resolutions_z0->GetXaxis()->SetTitle("Reco z_{0} - True z_{0} [mm]");
    h_seed_resolutions_z0->GetYaxis()->SetTitle("Entries");
    h_seed_resolutions_z0->Draw("HIST");
    MuCollStyle::AddStandardLabels(c_res_z0, "10 TeV");

    // Create output ROOT file and save histograms
    TFile* outputFile = new TFile("seed_analysis_stylized.root", "RECREATE");

    // Save canvases
    c_seed_theta->Write();
    c_seed_layer_barrel->Write();
    c_seed_layer_endcap->Write();
    c_seed_layer_barrel_comparison->Write();
    c_seed_layer_endcap_comparison->Write();
    c_seeds_per_MCP->Write();
    c_seed_number->Write();
    c_seed_matching->Write();
    c_res_q_over_pt->Write();
    c_res_d0->Write();
    c_res_z0->Write();

    outputFile->Close();
    delete outputFile;

    file->Close();
}

void Seeds_MuCollStyle(const char* filename) {
    run_collection(filename, "SeedTracks");
}
