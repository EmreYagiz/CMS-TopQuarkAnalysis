#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TVector2.h>

#include <iostream>
#include <fstream>
#include <iomanip>
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

struct CutFlow {
    double all = 0;
    double trigger = 0;
    double twoLeptons = 0;
    double oppositeSign = 0;
    double mll = 0;
    double zveto = 0;
    double met = 0;
    double twoJets = 0;
    double oneBtag = 0;
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

void AddCutFlow(CutFlow &a, const CutFlow &b) {
    a.all += b.all;
    a.trigger += b.trigger;
    a.twoLeptons += b.twoLeptons;
    a.oppositeSign += b.oppositeSign;
    a.mll += b.mll;
    a.zveto += b.zveto;
    a.met += b.met;
    a.twoJets += b.twoJets;
    a.oneBtag += b.oneBtag;
}

CutFlow ProcessSample(
    string filename,
    bool isData,
    double xsec,
    int requiredChannel
) {
    cout << "Processing " << filename << endl;

    CutFlow cf;

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
        return cf;
    }

    TTreeReader r("Events", f);

    TTreeReaderValue<Float_t> MET_pt(r, "MET_pt");

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

        cf.all += weight;

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
        cf.twoLeptons += weight;

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
        cf.oppositeSign += weight;

        bool passTrigger = false;

        if (channel == 0 && *HLT_IsoMu24) passTrigger = true;
        if (channel == 1 && *HLT_Ele27_WPTight_Gsf) passTrigger = true;
        if (channel == 2 && *HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL) passTrigger = true;

        if (!passTrigger) continue;
        cf.trigger += weight;

        Lepton l1 = leptons[best_i];
        Lepton l2 = leptons[best_j];

        TLorentzVector v1, v2;
        v1.SetPtEtaPhiM(l1.pt, l1.eta, l1.phi, l1.mass);
        v2.SetPtEtaPhiM(l2.pt, l2.eta, l2.phi, l2.mass);

        double mll = (v1 + v2).M();

        if (mll < 10) continue;
        cf.mll += weight;

        if (channel == 0 || channel == 1) {
            if (fabs(mll - 91.1876) < 15.0) continue;
        }
        cf.zveto += weight;

        if (channel == 0 || channel == 1) {
            if (*MET_pt < 30) continue;
        }

        if (channel == 2) {
            if (*MET_pt < 20) continue;
        }
        cf.met += weight;

        int njet = 0;
        int nbtag = 0;

        for (unsigned int j = 0; j < Jet_pt.GetSize(); j++) {
            if (Jet_pt[j] < 30) continue;
            if (fabs(Jet_eta[j]) > 2.5) continue;

            if (DeltaR(Jet_eta[j], Jet_phi[j], l1.eta, l1.phi) < 0.4) continue;
            if (DeltaR(Jet_eta[j], Jet_phi[j], l2.eta, l2.phi) < 0.4) continue;

            njet++;

            if (Jet_btagDeepFlavB[j] > BTAG_CUT) {
                nbtag++;
            }
        }

        if (njet < 2) continue;
        cf.twoJets += weight;

        if (nbtag < 1) continue;
        cf.oneBtag += weight;
    }

    if (genWeight) delete genWeight;
    f->Close();

    return cf;
}

void PrintRow(ofstream &csv, string name, CutFlow cf) {
    cout << left << setw(22) << name
         << setw(14) << fixed << setprecision(1) << cf.all
         << setw(14) << cf.twoLeptons
         << setw(14) << cf.oppositeSign
         << setw(14) << cf.trigger
         << setw(14) << cf.mll
         << setw(14) << cf.zveto
         << setw(14) << cf.met
         << setw(14) << cf.twoJets
         << setw(14) << cf.oneBtag
         << endl;

    csv << name << ","
        << cf.all << ","
        << cf.twoLeptons << ","
        << cf.oppositeSign << ","
        << cf.trigger << ","
        << cf.mll << ","
        << cf.zveto << ","
        << cf.met << ","
        << cf.twoJets << ","
        << cf.oneBtag << endl;
}

void cutflow() {
    CutFlow data;
    CutFlow ttbar;
    CutFlow dy;
    CutFlow vv;
    CutFlow st;

    AddCutFlow(data, ProcessSample("muonmuon.root", true, 0.0, 0));
    AddCutFlow(data, ProcessSample("electronelectron.root", true, 0.0, 1));
    AddCutFlow(data, ProcessSample("muonelectron.root", true, 0.0, 2));

    AddCutFlow(ttbar, ProcessSample("TTTo2L2Nu.root", false, 88.29, -1));

    AddCutFlow(dy, ProcessSample("dy.root", false, 6077.22, -1));

    AddCutFlow(vv, ProcessSample("WWTo2L2Nu.root", false, 12.178, -1));
    AddCutFlow(vv, ProcessSample("WZ.root", false, 47.13, -1));
    AddCutFlow(vv, ProcessSample("ZZ.root", false, 16.523, -1));

    AddCutFlow(st, ProcessSample("ST_tW_top.root", false, 35.85, -1));
    AddCutFlow(st, ProcessSample("ST_tW_antitop.root", false, 35.85, -1));

    CutFlow totalMC;
    AddCutFlow(totalMC, ttbar);
    AddCutFlow(totalMC, dy);
    AddCutFlow(totalMC, vv);
    AddCutFlow(totalMC, st);

    ofstream csv("cutflow.csv");

    csv << "Sample,All,2 leptons,Opposite sign,Trigger,mll > 10,Z veto,MET,>=2 jets,>=1 btag\n";

    cout << endl;
    cout << "============================== CUTFLOW TABLE ==============================" << endl;

    cout << left << setw(22) << "Sample"
         << setw(14) << "All"
         << setw(14) << "2 lep"
         << setw(14) << "OS"
         << setw(14) << "Trigger"
         << setw(14) << "mll>10"
         << setw(14) << "Z veto"
         << setw(14) << "MET"
         << setw(14) << ">=2 jets"
         << setw(14) << ">=1 btag"
         << endl;

    cout << string(148, '-') << endl;

    PrintRow(csv, "Data", data);
    PrintRow(csv, "ttbar", ttbar);
    PrintRow(csv, "DY", dy);
    PrintRow(csv, "VV", vv);
    PrintRow(csv, "SingleTop", st);
    PrintRow(csv, "Total MC", totalMC);

    csv.close();

    cout << endl;
    cout << "Saved cutflow.csv" << endl;
}
