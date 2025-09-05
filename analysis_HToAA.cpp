// Safiye Eker //sep4/2025
//Execute the code by using the command:
//TFile *file = TFile::Open("HToAATo2Mu2B.root"); TTree *tree = (TTree*)file->Get("Events");
//Then run the code by using:
// ./run_analyses
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <iostream>
#include <TF1.h>
#include <RooRealVar.h>
#include <RooVoigtian.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <TStyle.h>

int main() {
    TChain *chain = new TChain("Events");
    chain->Add("HToAATo2Mu2B.root");

    TTreeReader reader(chain);
    
    // --- Muon branches ---
    TTreeReaderValue<UInt_t> nMuon(reader, "nMuon");
    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<int> Muon_charge(reader, "Muon_charge");
    
    // --- Jet branches ---
    TTreeReaderValue<UInt_t> nJet(reader, "nJet");
    TTreeReaderArray<float> Jet_pt(reader, "Jet_pt");
    TTreeReaderArray<float> Jet_eta(reader, "Jet_eta");
    TTreeReaderArray<float> Jet_phi(reader, "Jet_phi");
    TTreeReaderArray<float> Jet_mass(reader, "Jet_mass");

    // --- Histograms ---
    // Single object histograms
    TH1F* h_muon_pt = new TH1F("h_muon_pt", ";Muon p_{T} [GeV];Events", 50, 0, 200);
    TH1F* h_muon_eta = new TH1F("h_muon_eta", ";Muon #eta;Events", 50, -3, 3);
    TH1F* h_muon_phi = new TH1F("h_muon_phi", ";Muon #phi;Events", 50, -TMath::Pi(), TMath::Pi());
    
    TH1F* h_jet_pt = new TH1F("h_jet_pt", ";Jet p_{T} [GeV];Events", 25, 10, 160);
    TH1F* h_jet_eta = new TH1F("h_jet_eta", ";Jet #eta;Events", 50, -3, 3);
    TH1F* h_jet_phi = new TH1F("h_jet_phi", ";Jet #phi;Events", 50, -TMath::Pi(), TMath::Pi());
    
    // Tight selection histograms
    TH1F* h_muon_pt_tight = new TH1F("h_muon_pt_tight", ";Muon p_{T} [GeV] (pT > 20 GeV);Events", 50, 0, 200);
    TH1F* h_jet_pt_tight = new TH1F("h_jet_pt_tight", ";Jet p_{T} [GeV] (pT > 30 GeV);Events", 25, 10, 160);

    // Di-object mass histograms
    TH1F* h_dimu_mass = new TH1F("h_dimu_mass", ";Di-muon Mass [GeV];Events", 80, 0, 200);
    TH1F* h_dijet_mass = new TH1F("h_dijet_mass", ";Di-jet Mass [GeV];Events", 80, 0, 200);
    TH1F* h_higgs_mass = new TH1F("h_higgs_mass", ";Higgs Candidate Mass [GeV];Events", 80, 0, 300);

    // 2D correlation histograms
    TH2F* h2_dimu_dijet = new TH2F("h2_dimu_dijet", "Dimuon Mass vs Dijet Mass;M_{μμ} [GeV];M_{bb} [GeV]", 
                                  80, 0, 200, 80, 0, 200);
    TH2F* h2_jet_pt_eta = new TH2F("h2_jet_pt_eta", "Jet p_{T} vs #eta;Jet #eta;Jet p_{T} [GeV]", 
                                  50, -3, 3, 25, 10, 160);
    TH2F* h2_muon_pt_eta = new TH2F("h2_muon_pt_eta", "Muon p_{T} vs #eta;Muon #eta;Muon p_{T} [GeV]", 
                                  50, -3, 3, 50, 0, 200);
    
    // Phi correlation histograms
    TH2F* h2_dimu_phi = new TH2F("h2_dimu_phi", "Di-muon #phi correlation;#phi_{μ1};#phi_{μ2}", 
                               50, -TMath::Pi(), TMath::Pi(), 50, -TMath::Pi(), TMath::Pi());
    TH2F* h2_dijet_phi = new TH2F("h2_dijet_phi", "Di-jet #phi correlation;#phi_{j1};#phi_{j2}", 
                               50, -TMath::Pi(), TMath::Pi(), 50, -TMath::Pi(), TMath::Pi());
    TH2F* h2_dimu_dijet_phi = new TH2F("h2_dimu_dijet_phi", "#Delta#phi(μμ) vs #Delta#phi(bb);#Delta#phi_{μμ};#Delta#phi_{bb}", 
                                     50, 0, TMath::Pi(), 50, 0, TMath::Pi());

    // --- Event loop ---
    int nEvents = 0;
    int nPassed = 0;
    
    while (reader.Next()) {
        nEvents++;
        
        // Basic selection
        if (*nMuon < 2 || *nJet < 2) continue;
        
        // Muon selection
        if (Muon_pt[0] < 5.0 || fabs(Muon_eta[0]) > 2.5) continue;
        if (Muon_pt[1] < 5.0 || fabs(Muon_eta[1]) > 2.5) continue;
        if (Muon_charge[0] * Muon_charge[1] > 0) continue;
        
        // Jet selection
        if (Jet_pt[0] < 10.0 || fabs(Jet_eta[0]) > 2.5) continue;
        if (Jet_pt[1] < 10.0 || fabs(Jet_eta[1]) > 2.5) continue;
        
        nPassed++;

        // Fill single object histograms
        for (size_t i = 0; i < Muon_pt.GetSize(); i++) {
            h_muon_pt->Fill(Muon_pt[i]);
            h_muon_eta->Fill(Muon_eta[i]);
            h_muon_phi->Fill(Muon_phi[i]);
            h2_muon_pt_eta->Fill(Muon_eta[i], Muon_pt[i]);
            
            if (Muon_pt[i] > 20.0) {
                h_muon_pt_tight->Fill(Muon_pt[i]);
            }
        }

        for (size_t i = 0; i < Jet_pt.GetSize(); i++) {
            h_jet_pt->Fill(Jet_pt[i]);
            h_jet_eta->Fill(Jet_eta[i]);
            h_jet_phi->Fill(Jet_phi[i]);
            h2_jet_pt_eta->Fill(Jet_eta[i], Jet_pt[i]);
            
            if (Jet_pt[i] > 30.0) {
                h_jet_pt_tight->Fill(Jet_pt[i]);
            }
        }

        // Dimuon reconstruction
        TLorentzVector mu1, mu2;
        mu1.SetPtEtaPhiM(Muon_pt[0], Muon_eta[0], Muon_phi[0], 0.105);
        mu2.SetPtEtaPhiM(Muon_pt[1], Muon_eta[1], Muon_phi[1], 0.105);
        float dimu_mass = (mu1 + mu2).M();
        float dphi_dimu = mu1.DeltaPhi(mu2);
        
        h_dimu_mass->Fill(dimu_mass);
        h2_dimu_phi->Fill(Muon_phi[0], Muon_phi[1]);

        // Dijet reconstruction
        TLorentzVector j1, j2;
        j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
        j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
        float dijet_mass = (j1 + j2).M();
        float dphi_dijet = j1.DeltaPhi(j2);
        
        h_dijet_mass->Fill(dijet_mass);
        h2_dijet_phi->Fill(Jet_phi[0], Jet_phi[1]);
        h2_dimu_dijet_phi->Fill(fabs(dphi_dimu), fabs(dphi_dijet));

        // Dimuon vs Dijet correlation
        h2_dimu_dijet->Fill(dimu_mass, dijet_mass);

        // Higgs candidate: di-muon + di-jet
        TLorentzVector higgs = mu1 + mu2 + j1 + j2;
        h_higgs_mass->Fill(higgs.M());
    }

    // --- Save output ---
    TFile fout("effnow_output.root", "RECREATE");
    
    // Single object histograms
    h_muon_pt->Write();
    h_muon_eta->Write();
    h_muon_phi->Write();
    h_jet_pt->Write();
    h_jet_eta->Write();
    h_jet_phi->Write();
    
    // Tight selection histograms
    h_muon_pt_tight->Write();
    h_jet_pt_tight->Write();
    
    // Di-object mass histograms
    h_dimu_mass->Write();
    h_dijet_mass->Write();
    h_higgs_mass->Write();
    
    // 2D correlation histograms
    h2_dimu_dijet->Write();
    h2_jet_pt_eta->Write();
    h2_muon_pt_eta->Write();
    h2_dimu_phi->Write();
    h2_dijet_phi->Write();
    h2_dimu_dijet_phi->Write();
    
    fout.Close();

    std::cout << "Processed " << nEvents << " events" << std::endl;
    std::cout << "Passed selection: " << nPassed << " events" << std::endl;
    std::cout << "Efficiency: " << (nPassed * 100.0 / nEvents) << "%" << std::endl;
    std::cout << "All histograms written to effnow_output.root" << std::endl;

    return 0;
}
