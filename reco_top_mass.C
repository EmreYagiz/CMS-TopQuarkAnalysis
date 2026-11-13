#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

const double LUMI = 8605.0;
const double BTAG_CUT = 0.2489;

struct Lepton {
    double pt, eta, phi, mass;
    int charge;
    int type; // 0 muon, 1 electron
};

double DeltaR(double eta1, double phi1, double eta2, double phi2) {
    double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    double deta = eta1 - eta2;
    return sqrt(deta*deta + dphi*dphi);
}

double SumGenWeight(string filename) {
    TFile *f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) return 1.0;

    TTree *t = (TTree*)f->Get("Events");
    if (!t) {
        f->Close();
        return 1.0;
    }

    float genWeight = 1.0;
    t->SetBranchAddress("genWeight", &genWeight);

    double sumw = 0.0;
    Long64_t n = t->GetEntries();

    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);
        sumw += genWeight;
    }

    f->Close();
    if (sumw == 0) return 1.0;
    return sumw;
}

TLorentzVector BuildNeutrino(TLorentzVector lep, double met_px, double met_py) {
    const double MW = 80.379;

    double px = met_px;
    double py = met_py;

    double mu = (MW * MW) / 2.0 + lep.Px() * px + lep.Py() * py;
    double a = lep.E() * lep.E() - lep.Pz() * lep.Pz();
    double b = -2.0 * mu * lep.Pz();
    double c = lep.E() * lep.E() * (px * px + py * py) - mu * mu;

    double disc = b * b - 4.0 * a * c;
    double pz = 0.0;

    if (fabs(a) < 1e-6) {
        pz = 0.0;
    } else if (disc >= 0) {
        double pz1 = (-b + sqrt(disc)) / (2.0 * a);
        double pz2 = (-b - sqrt(disc)) / (2.0 * a);
        pz = (fabs(pz1) < fabs(pz2)) ? pz1 : pz2;
    } else {
        pz = -b / (2.0 * a);
    }

    double e = sqrt(px * px + py * py + pz * pz);

    TLorentzVector nu;
    nu.SetPxPyPzE(px, py, pz, e);
    return nu;
}

void ProcessSample(
    string filename,
    TH1D *h,
    bool isData,
    double xsec,
    int requiredChannel
) {
    cout << "Processing " << filename << endl;

    double norm = 1.0;

    if (!isData) {
        double sumw = SumGenWeight(filename);
        norm = LUMI * xsec / sumw;
        cout << "  sum genWeight = " << sumw << endl;
        cout << "  norm = " << norm << endl;
    }

    TFile *f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) {
        cout << "Cannot open " << filename << endl;
        return;
    }

    TTreeReader r("Events", f);

    TTreeReaderValue<Float_t> MET_pt(r, "MET_pt");
    TTreeReaderValue<Float_t> MET_phi(r, "MET_phi");

    TTreeReaderValue<Bool_t> HLT_IsoMu24(r, "HLT_IsoMu24");
    TTreeReaderValue<Bool_t> HLT_Ele27_WPTight_Gsf(r, "HLT_Ele27_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(
        r, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"
    );

    TTreeReaderArray<Float_t> Muon_pt(r, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(r, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(r, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_mass(r, "Muon_mass");
    TTreeReaderArray<Int_t>   Muon_charge(r, "Muon_charge");
    TTreeReaderArray<Bool_t>  Muon_tightId(r, "Muon_tightId");
    TTreeReaderArray<Float_t> Muon_pfRelIso04_all(r, "Muon_pfRelIso04_all");

    TTreeReaderArray<Float_t> Electron_pt(r, "Electron_pt");
    TTreeReaderArray<Float_t> Electron_eta(r, "Electron_eta");
    TTreeReaderArray<Float_t> Electron_phi(r, "Electron_phi");
    TTreeReaderArray<Float_t> Electron_mass(r, "Electron_mass");
    TTreeReaderArray<Int_t>   Electron_charge(r, "Electron_charge");
    TTreeReaderArray<Int_t>   Electron_cutBased(r, "Electron_cutBased");
    TTreeReaderArray<Float_t> Electron_pfRelIso03_all(r, "Electron_pfRelIso03_all");

    TTreeReaderArray<Float_t> Jet_pt(r, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_eta(r, "Jet_eta");
    TTreeReaderArray<Float_t> Jet_phi(r, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_btagDeepFlavB(r, "Jet_btagDeepFlavB");

    TTreeReaderValue<Float_t>* genWeight = 0;
    if (!isData) {
        genWeight = new TTreeReaderValue<Float_t>(r, "genWeight");
    }

    while (r.Next()) {
        double weight = 1.0;
        if (!isData) weight = (**genWeight) * norm;

        vector<Lepton> leptons;

        for (unsigned int i = 0; i < Muon_pt.GetSize(); i++) {
            if (Muon_pt[i] < 20) continue;
            if (fabs(Muon_eta[i]) > 2.4) continue;
            if (!Muon_tightId[i]) continue;
            if (Muon_pfRelIso04_all[i] > 0.15) continue;

            Lepton l;
            l.pt = Muon_pt[i];
            l.eta = Muon_eta[i];
            l.phi = Muon_phi[i];
            l.mass = Muon_mass[i];
            l.charge = Muon_charge[i];
            l.type = 0;
            leptons.push_back(l);
        }

        for (unsigned int i = 0; i < Electron_pt.GetSize(); i++) {
            if (Electron_pt[i] < 20) continue;
            if (fabs(Electron_eta[i]) > 2.5) continue;
            if (Electron_cutBased[i] < 4) continue;
            if (Electron_pfRelIso03_all[i] > 0.15) continue;

            Lepton l;
            l.pt = Electron_pt[i];
            l.eta = Electron_eta[i];
            l.phi = Electron_phi[i];
            l.mass = Electron_mass[i];
            l.charge = Electron_charge[i];
            l.type = 1;
            leptons.push_back(l);
        }

        if (leptons.size() < 2) continue;

        int best_i = -1;
        int best_j = -1;
        int channel = -1;
        double best_ptsum = -1.0;

        for (unsigned int i = 0; i < leptons.size(); i++) {
            for (unsigned int j = i + 1; j < leptons.size(); j++) {
                if (leptons[i].charge * leptons[j].charge >= 0) continue;

                int ch = -1;
                if (leptons[i].type == 0 && leptons[j].type == 0) ch = 0;
                if (leptons[i].type == 1 && leptons[j].type == 1) ch = 1;
                if (leptons[i].type != leptons[j].type) ch = 2;

                if (requiredChannel != -1 && ch != requiredChannel) continue;

                double ptsum = leptons[i].pt + leptons[j].pt;
                if (ptsum > best_ptsum) {
                    best_ptsum = ptsum;
                    best_i = i;
                    best_j = j;
                    channel = ch;
                }
            }
        }

        if (best_i < 0 || best_j < 0) continue;

        bool passTrigger = false;

        if (channel == 0 && *HLT_IsoMu24) passTrigger = true;
        if (channel == 1 && *HLT_Ele27_WPTight_Gsf) passTrigger = true;
        if (channel == 2 && *HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) passTrigger = true;

        if (!passTrigger) continue;

        Lepton l1 = leptons[best_i];
        Lepton l2 = leptons[best_j];

        TLorentzVector lep1, lep2;
        lep1.SetPtEtaPhiM(l1.pt, l1.eta, l1.phi, l1.mass);
        lep2.SetPtEtaPhiM(l2.pt, l2.eta, l2.phi, l2.mass);

        double mll = (lep1 + lep2).M();

        if (mll < 10) continue;

        if (channel == 0 || channel == 1) {
            if (fabs(mll - 91.1876) < 15.0) continue;
            if (*MET_pt < 30) continue;
        }

        if (channel == 2) {
            if (*MET_pt < 20) continue;
        }

        vector<TLorentzVector> selectedJets;
        vector<double> selectedBtags;

        int njet = 0;
        int nbtag = 0;

        for (unsigned int j = 0; j < Jet_pt.GetSize(); j++) {
            if (Jet_pt[j] < 30) continue;
            if (fabs(Jet_eta[j]) > 2.5) continue;

            if (DeltaR(Jet_eta[j], Jet_phi[j], l1.eta, l1.phi) < 0.4) continue;
            if (DeltaR(Jet_eta[j], Jet_phi[j], l2.eta, l2.phi) < 0.4) continue;

            TLorentzVector jet;
            jet.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], 0.0);

            selectedJets.push_back(jet);
            selectedBtags.push_back(Jet_btagDeepFlavB[j]);

            njet++;
            if (Jet_btagDeepFlavB[j] > BTAG_CUT) nbtag++;
        }

        if (njet < 2) continue;
        if (nbtag < 1) continue;

        int j1 = -1;
        int j2 = -1;
        double bestB1 = -999.0;
        double bestB2 = -999.0;

        for (unsigned int j = 0; j < selectedJets.size(); j++) {
            if (selectedBtags[j] > bestB1) {
                bestB2 = bestB1;
                j2 = j1;
                bestB1 = selectedBtags[j];
                j1 = j;
            } else if (selectedBtags[j] > bestB2) {
                bestB2 = selectedBtags[j];
                j2 = j;
            }
        }

        if (j1 < 0 || j2 < 0) continue;

        double met_px = (*MET_pt) * cos(*MET_phi);
        double met_py = (*MET_pt) * sin(*MET_phi);

        double nu1_px = 0.5 * met_px;
        double nu1_py = 0.5 * met_py;
        double nu2_px = 0.5 * met_px;
        double nu2_py = 0.5 * met_py;

        TLorentzVector nu1 = BuildNeutrino(lep1, nu1_px, nu1_py);
        TLorentzVector nu2 = BuildNeutrino(lep2, nu2_px, nu2_py);

        TLorentzVector b1 = selectedJets[j1];
        TLorentzVector b2 = selectedJets[j2];

        double mt_a = ((b1 + lep1 + nu1).M() + (b2 + lep2 + nu2).M()) / 2.0;
        double mt_b = ((b1 + lep2 + nu2).M() + (b2 + lep1 + nu1).M()) / 2.0;

        double mt_reco = fabs(mt_a - 172.5) < fabs(mt_b - 172.5) ? mt_a : mt_b;

        if (mt_reco > 100 && mt_reco < 400) {
            h->Fill(mt_reco, weight);
        }
    }

    if (genWeight) delete genWeight;
    f->Close();
}

void DrawTopMass(
    TH1D *h_data,
    TH1D *h_tt,
    TH1D *h_dy,
    TH1D *h_vv,
    TH1D *h_st
) {
    double dataIntegral = h_data->Integral();
    double mcIntegral =
        h_tt->Integral() +
        h_dy->Integral() +
        h_vv->Integral() +
        h_st->Integral();

    double scale = 1.0;
    if (mcIntegral > 0) scale = dataIntegral / mcIntegral;

    h_tt->Scale(scale);
    h_dy->Scale(scale);
    h_vv->Scale(scale);
    h_st->Scale(scale);

    TCanvas *c = new TCanvas("c_topmass", "c_topmass", 800, 800);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.30, 1, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.30);

    pad1->SetBottomMargin(0.03);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.30);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();

    h_tt->SetFillColor(kRed - 7);
    h_dy->SetFillColor(kAzure - 9);
    h_vv->SetFillColor(kGreen - 6);
    h_st->SetFillColor(kMagenta - 7);

    h_tt->SetLineColor(kBlack);
    h_dy->SetLineColor(kBlack);
    h_vv->SetLineColor(kBlack);
    h_st->SetLineColor(kBlack);

    THStack *st = new THStack("st", ";Reconstructed top mass [GeV];Events");
    st->Add(h_dy);
    st->Add(h_vv);
    st->Add(h_st);
    st->Add(h_tt);

    st->Draw("hist");

    double ymax = st->GetMaximum();
    if (h_data->GetMaximum() > ymax) ymax = h_data->GetMaximum();
    st->SetMaximum(1.5 * ymax);

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.0);
    h_data->SetLineColor(kBlack);
    h_data->Draw("E SAME");

    TLegend *leg = new TLegend(0.62, 0.62, 0.88, 0.88);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_tt, "t#bar{t}", "f");
    leg->AddEntry(h_dy, "Drell-Yan", "f");
    leg->AddEntry(h_vv, "VV", "f");
    leg->AddEntry(h_st, "Single top", "f");
    leg->Draw();

    pad2->cd();

    TH1D *h_mc = (TH1D*)h_tt->Clone("h_mc");
    h_mc->Add(h_dy);
    h_mc->Add(h_vv);
    h_mc->Add(h_st);

    TH1D *ratio = (TH1D*)h_data->Clone("ratio");
    ratio->Divide(h_mc);

    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->GetYaxis()->SetRangeUser(0, 2);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.09);
    ratio->GetYaxis()->SetLabelSize(0.08);
    ratio->GetXaxis()->SetTitle("Reconstructed top mass [GeV]");
    ratio->GetXaxis()->SetTitleSize(0.11);
    ratio->GetXaxis()->SetLabelSize(0.09);
    ratio->Draw("E");

    TLine *line = new TLine(100, 1.0, 400, 1.0);
    line->SetLineStyle(2);
    line->Draw("SAME");

    gSystem->mkdir("plots", true);
    c->SaveAs("plots/reconstructed_top_mass.png");
    c->SaveAs("plots/reconstructed_top_mass.pdf");
}

void reco_top_mass() {
    gStyle->SetOptStat(0);

    TH1D *h_data = new TH1D("h_data", ";Reconstructed top mass [GeV];Events", 30, 100, 400);
    TH1D *h_tt   = new TH1D("h_tt",   ";Reconstructed top mass [GeV];Events", 30, 100, 400);
    TH1D *h_dy   = new TH1D("h_dy",   ";Reconstructed top mass [GeV];Events", 30, 100, 400);
    TH1D *h_vv   = new TH1D("h_vv",   ";Reconstructed top mass [GeV];Events", 30, 100, 400);
    TH1D *h_st   = new TH1D("h_st",   ";Reconstructed top mass [GeV];Events", 30, 100, 400);

    h_data->Sumw2();
    h_tt->Sumw2();
    h_dy->Sumw2();
    h_vv->Sumw2();
    h_st->Sumw2();

    ProcessSample("muonmuon.root", h_data, true, 0.0, 0);
    ProcessSample("electronelectron.root", h_data, true, 0.0, 1);
    ProcessSample("muonelectron.root", h_data, true, 0.0, 2);

    ProcessSample("TTTo2L2Nu.root", h_tt, false, 88.29, -1);

    ProcessSample("dy.root", h_dy, false, 6077.22, -1);

    ProcessSample("WWTo2L2Nu.root", h_vv, false, 12.178, -1);
    ProcessSample("WZ.root", h_vv, false, 47.13, -1);
    ProcessSample("ZZ.root", h_vv, false, 16.523, -1);

    ProcessSample("ST_tW_top.root", h_st, false, 35.85, -1);
    ProcessSample("ST_tW_antitop.root", h_st, false, 35.85, -1);

    cout << endl;
    cout << "Top mass selected yields:" << endl;
    cout << "Data      = " << h_data->Integral() << endl;
    cout << "ttbar     = " << h_tt->Integral() << endl;
    cout << "DY        = " << h_dy->Integral() << endl;
    cout << "VV        = " << h_vv->Integral() << endl;
    cout << "SingleTop = " << h_st->Integral() << endl;

    DrawTopMass(h_data, h_tt, h_dy, h_vv, h_st);

    TFile *out = new TFile("reconstructed_top_mass.root", "RECREATE");
    h_data->Write();
    h_tt->Write();
    h_dy->Write();
    h_vv->Write();
    h_st->Write();
    out->Close();

    cout << endl;
    cout << "Saved: plots/reconstructed_top_mass.png" << endl;
    cout << "Saved: plots/reconstructed_top_mass.pdf" << endl;
    cout << "Saved: reconstructed_top_mass.root" << endl;
}
