#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TVector2.h>

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
    int type; // 0 = muon, 1 = electron
};

struct HistSet {
    TH1D* h_njet;
    TH1D* h_btag;
    TH1D* h_met;
    TH1D* h_mll;
    TH1D* h_ht;
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

HistSet MakeHistSet(string name) {
    HistSet h;

    h.h_njet = new TH1D((name+"_njet").c_str(),
        ";Number of jets;Events", 5, -0.5, 4.5);

    h.h_btag = new TH1D((name+"_btag").c_str(),
        ";Number of b-tagged jets;Events", 3, -0.5, 2.5);

    h.h_met = new TH1D((name+"_met").c_str(),
        ";MET [GeV];Events", 20, 0, 200);

    h.h_mll = new TH1D((name+"_mll").c_str(),
        ";m_{ll} [GeV];Events", 30, 0, 300);

    h.h_ht = new TH1D((name+"_ht").c_str(),
        ";H_{T} [GeV];Events", 25, 0, 500);

    h.h_njet->Sumw2();
    h.h_btag->Sumw2();
    h.h_met->Sumw2();
    h.h_mll->Sumw2();
    h.h_ht->Sumw2();

    return h;
}

void ProcessSample(
    string filename,
    HistSet &h,
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

    TTreeReaderValue<Bool_t> HLT_IsoMu24(r, "HLT_IsoMu24");
    TTreeReaderValue<Bool_t> HLT_Ele27_WPTight_Gsf(r, "HLT_Ele27_WPTight_Gsf");
    TTreeReaderValue<Bool_t> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL(
        r, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL");

    TTreeReaderArray<Float_t> Muon_pt(r, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(r, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(r, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_mass(r, "Muon_mass");
    TTreeReaderArray<Int_t> Muon_charge(r, "Muon_charge");
    TTreeReaderArray<Bool_t> Muon_tightId(r, "Muon_tightId");
    TTreeReaderArray<Float_t> Muon_pfRelIso04_all(r, "Muon_pfRelIso04_all");

    TTreeReaderArray<Float_t> Electron_pt(r, "Electron_pt");
    TTreeReaderArray<Float_t> Electron_eta(r, "Electron_eta");
    TTreeReaderArray<Float_t> Electron_phi(r, "Electron_phi");
    TTreeReaderArray<Float_t> Electron_mass(r, "Electron_mass");
    TTreeReaderArray<Int_t> Electron_charge(r, "Electron_charge");
    TTreeReaderArray<Int_t> Electron_cutBased(r, "Electron_cutBased");
    TTreeReaderArray<Float_t> Electron_pfRelIso03_all(r, "Electron_pfRelIso03_all");

    TTreeReaderArray<Float_t> Jet_pt(r, "Jet_pt");
    TTreeReaderArray<Float_t> Jet_eta(r, "Jet_eta");
    TTreeReaderArray<Float_t> Jet_phi(r, "Jet_phi");
    TTreeReaderArray<Float_t> Jet_btagDeepFlavB(r, "Jet_btagDeepFlavB");

    TTreeReaderValue<Float_t>* genWeight = nullptr;

    if (!isData) {
        genWeight = new TTreeReaderValue<Float_t>(r, "genWeight");
    }

    while (r.Next()) {

        double weight = 1.0;

        if (!isData) {
            weight = (**genWeight) * norm;
        }

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
        double best_ptsum = -1.0;
        int channel = -1;

        for (unsigned int i = 0; i < leptons.size(); i++) {
            for (unsigned int j = i + 1; j < leptons.size(); j++) {

                if (leptons[i].charge * leptons[j].charge >= 0) continue;

                int ch = -1;

                if (leptons[i].type == 0 && leptons[j].type == 0)
                    ch = 0; // mumu

                if (leptons[i].type == 1 && leptons[j].type == 1)
                    ch = 1; // ee

                if (leptons[i].type != leptons[j].type)
                    ch = 2; // emu

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

        if (channel == 0 && *HLT_IsoMu24)
            passTrigger = true;

        if (channel == 1 && *HLT_Ele27_WPTight_Gsf)
            passTrigger = true;

        if (channel == 2 && *HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL)
            passTrigger = true;

        if (!passTrigger) continue;

        Lepton l1 = leptons[best_i];
        Lepton l2 = leptons[best_j];

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(l1.pt, l1.eta, l1.phi, l1.mass);
        v2.SetPtEtaPhiM(l2.pt, l2.eta, l2.phi, l2.mass);

        double mll = (v1 + v2).M();

        if (mll < 10.0) continue;

        if (channel == 0 || channel == 1) {
            if (fabs(mll - 91.1876) < 15.0) continue;
            if (*MET_pt < 30.0) continue;
        }

        if (channel == 2) {
            if (*MET_pt < 20.0) continue;
        }

        int njet = 0;
        int nbtag = 0;
        double ht = 0.0;

        for (unsigned int j = 0; j < Jet_pt.GetSize(); j++) {
            if (Jet_pt[j] < 30.0) continue;
            if (fabs(Jet_eta[j]) > 2.5) continue;

            if (DeltaR(Jet_eta[j], Jet_phi[j], l1.eta, l1.phi) < 0.4) continue;
            if (DeltaR(Jet_eta[j], Jet_phi[j], l2.eta, l2.phi) < 0.4) continue;

            njet++;
            ht += Jet_pt[j];

            if (Jet_btagDeepFlavB[j] > BTAG_CUT)
                nbtag++;
        }

        int njet_plot = njet;
        if (njet_plot >= 4) njet_plot = 4;

        int nbtag_plot = nbtag;
        if (nbtag_plot >= 2) nbtag_plot = 2;

        // Figure 1 paper logic:
        // fill jet multiplicity before requiring >=2 jets
        h.h_njet->Fill(njet_plot, weight);

        // Final paper event selection ends here
        if (njet < 2) continue;

        // Figure 3 paper logic:
        // b-tag multiplicity is only a validation plot, not a cut
        h.h_btag->Fill(nbtag_plot, weight);

        // Final selected distributions
        h.h_met->Fill(*MET_pt, weight);
        h.h_mll->Fill(mll, weight);
        h.h_ht->Fill(ht, weight);
    }

    if (genWeight) delete genWeight;
    f->Close();
}

void DrawStack(
    TH1D* h_data,
    TH1D* h_tt,
    TH1D* h_dy,
    TH1D* h_vv,
    TH1D* h_st,
    string cname
) {
    TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 800, 800);

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

    double dataIntegral = h_data->Integral();

    double mcIntegral =
        h_tt->Integral()
        + h_dy->Integral()
        + h_vv->Integral()
        + h_st->Integral();

    if (mcIntegral > 0) {
        double scale = dataIntegral / mcIntegral;

        h_tt->Scale(scale);
        h_dy->Scale(scale);
        h_vv->Scale(scale);
        h_st->Scale(scale);
    }

    THStack *st = new THStack("st", "");

    st->Add(h_dy);
    st->Add(h_vv);
    st->Add(h_st);
    st->Add(h_tt);

    st->Draw("hist");
    st->GetYaxis()->SetTitle("Events");

    double max1 = st->GetMaximum();
    double max2 = h_data->GetMaximum();

    if (max2 > max1) st->SetMaximum(max2 * 1.5);
    else st->SetMaximum(max1 * 1.5);

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.0);
    h_data->SetLineColor(kBlack);
    h_data->Draw("E SAME");

    TLegend *leg = new TLegend(0.65, 0.62, 0.88, 0.88);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_tt, "t#bar{t}", "f");
    leg->AddEntry(h_dy, "Drell-Yan", "f");
    leg->AddEntry(h_vv, "VV", "f");
    leg->AddEntry(h_st, "Single top", "f");
    leg->Draw();

    pad2->cd();

    TH1D *h_mc = (TH1D*)h_tt->Clone((cname+"_mc").c_str());
    h_mc->Add(h_dy);
    h_mc->Add(h_vv);
    h_mc->Add(h_st);

    TH1D *ratio = (TH1D*)h_data->Clone((cname+"_ratio").c_str());
    ratio->Divide(h_mc);

    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->GetYaxis()->SetRangeUser(0, 2);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.09);
    ratio->GetYaxis()->SetLabelSize(0.08);
    ratio->GetXaxis()->SetTitleSize(0.11);
    ratio->GetXaxis()->SetLabelSize(0.09);

    ratio->Draw("E");

    TLine *line = new TLine(
        ratio->GetXaxis()->GetXmin(),
        1.0,
        ratio->GetXaxis()->GetXmax(),
        1.0
    );

    line->SetLineStyle(2);
    line->Draw("SAME");

    gSystem->mkdir("plots", true);

    c->SaveAs(("plots/" + cname + ".png").c_str());
    c->SaveAs(("plots/" + cname + ".pdf").c_str());
}

void ttbar_dilepton_analysis1() {
    gStyle->SetOptStat(0);

    HistSet h_data = MakeHistSet("data");
    HistSet h_tt   = MakeHistSet("tt");
    HistSet h_dy   = MakeHistSet("dy");
    HistSet h_vv   = MakeHistSet("vv");
    HistSet h_st   = MakeHistSet("st");

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
    cout << "========== FINAL YIELDS AFTER PAPER-LIKE SELECTION ==========" << endl;
    cout << "Selection ends at >=2 jets. No b-tag cut is applied." << endl;
    cout << "Data      = " << h_data.h_met->Integral() << endl;
    cout << "ttbar     = " << h_tt.h_met->Integral() << endl;
    cout << "DY        = " << h_dy.h_met->Integral() << endl;
    cout << "VV        = " << h_vv.h_met->Integral() << endl;
    cout << "SingleTop = " << h_st.h_met->Integral() << endl;
    cout << "Total MC  = "
         << h_tt.h_met->Integral()
          + h_dy.h_met->Integral()
          + h_vv.h_met->Integral()
          + h_st.h_met->Integral()
         << endl;

    DrawStack(h_data.h_njet, h_tt.h_njet, h_dy.h_njet, h_vv.h_njet, h_st.h_njet,
        "jet_multiplicity_before_2jet_cut");

    DrawStack(h_data.h_btag, h_tt.h_btag, h_dy.h_btag, h_vv.h_btag, h_st.h_btag,
        "btag_multiplicity_validation");

    DrawStack(h_data.h_met, h_tt.h_met, h_dy.h_met, h_vv.h_met, h_st.h_met,
        "met_after_paper_selection");

    DrawStack(h_data.h_mll, h_tt.h_mll, h_dy.h_mll, h_vv.h_mll, h_st.h_mll,
        "dilepton_mass_after_paper_selection");

    DrawStack(h_data.h_ht, h_tt.h_ht, h_dy.h_ht, h_vv.h_ht, h_st.h_ht,
        "ht_after_paper_selection");

    TFile *out = new TFile("ttbar_output_paper_like.root", "RECREATE");

    h_data.h_njet->Write();
    h_data.h_btag->Write();
    h_data.h_met->Write();
    h_data.h_mll->Write();
    h_data.h_ht->Write();

    h_tt.h_njet->Write();
    h_tt.h_btag->Write();
    h_tt.h_met->Write();
    h_tt.h_mll->Write();
    h_tt.h_ht->Write();

    h_dy.h_njet->Write();
    h_dy.h_btag->Write();
    h_dy.h_met->Write();
    h_dy.h_mll->Write();
    h_dy.h_ht->Write();

    h_vv.h_njet->Write();
    h_vv.h_btag->Write();
    h_vv.h_met->Write();
    h_vv.h_mll->Write();
    h_vv.h_ht->Write();

    h_st.h_njet->Write();
    h_st.h_btag->Write();
    h_st.h_met->Write();
    h_st.h_mll->Write();
    h_st.h_ht->Write();

    out->Close();

    cout << endl;
    cout << "Saved plots in plots/ folder." << endl;
    cout << "Saved histograms to ttbar_output_paper_like.root" << endl;
}
