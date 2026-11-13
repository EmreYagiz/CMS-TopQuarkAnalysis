#include "cutflow.C"

void cross_section() {

    double LUMI = 8605.0; // pb^-1, Run2016H yaklaşık değer

    double sigma_ttbar_theory = 831.76; // pb, 13 TeV inclusive ttbar
    double sigma_TTTo2L2Nu = 88.29;     // pb, TTTo2L2Nu sample cross section

    CutFlow data;
    CutFlow ttbar;
    CutFlow dy;
    CutFlow vv;
    CutFlow st;

    AddCutFlow(data, ProcessSample("muonmuon.root", true, 0.0, 0));
    AddCutFlow(data, ProcessSample("electronelectron.root", true, 0.0, 1));
    AddCutFlow(data, ProcessSample("muonelectron.root", true, 0.0, 2));

    AddCutFlow(ttbar, ProcessSample("TTTo2L2Nu.root", false, sigma_TTTo2L2Nu, -1));

    AddCutFlow(dy, ProcessSample("dy.root", false, 6077.22, -1));

    AddCutFlow(vv, ProcessSample("WWTo2L2Nu.root", false, 12.178, -1));
    AddCutFlow(vv, ProcessSample("WZ.root", false, 47.13, -1));
    AddCutFlow(vv, ProcessSample("ZZ.root", false, 16.523, -1));

    AddCutFlow(st, ProcessSample("ST_tW_top.root", false, 35.85, -1));
    AddCutFlow(st, ProcessSample("ST_tW_antitop.root", false, 35.85, -1));

    double Ndata = data.oneBtag;

    double Nbkg = dy.oneBtag + vv.oneBtag + st.oneBtag;

    double Nsignal = Ndata - Nbkg;

    double Nttbar_pred = ttbar.oneBtag;

    double scale_factor = Nsignal / Nttbar_pred;

    double sigma_TTTo2L2Nu_measured = scale_factor * sigma_TTTo2L2Nu;

    double sigma_ttbar_inclusive_measured = scale_factor * sigma_ttbar_theory;

    double stat_unc = sqrt(Ndata) / Nttbar_pred * sigma_ttbar_theory;

    cout << endl;
    cout << "========== CROSS SECTION ESTIMATE ==========" << endl;

    cout << "Luminosity = " << LUMI << " pb^-1" << endl;
    cout << endl;

    cout << "N_data              = " << Ndata << endl;
    cout << "N_background        = " << Nbkg << endl;
    cout << "N_signal = data-bkg = " << Nsignal << endl;
    cout << "N_ttbar_MC_expected = " << Nttbar_pred << endl;
    cout << endl;

    cout << "Scale factor = (Data - Background) / ttbar_MC" << endl;
    cout << "Scale factor = " << scale_factor << endl;
    cout << endl;

    cout << "Measured TTTo2L2Nu cross section:" << endl;
    cout << "sigma(TTTo2L2Nu) = "
         << sigma_TTTo2L2Nu_measured
         << " pb" << endl;

    cout << endl;

    cout << "Approximate inclusive ttbar cross section:" << endl;
    cout << "sigma(pp -> ttbar) = "
         << sigma_ttbar_inclusive_measured
         << " +/- "
         << stat_unc
         << " pb (stat only, approximate)" << endl;

    cout << endl;
    cout << "Theory reference used: sigma_ttbar = "
         << sigma_ttbar_theory << " pb" << endl;

    cout << endl;
    cout << "WARNING:" << endl;
    cout << "This is not a publication-level measurement." << endl;
    cout << "No pileup reweighting, lepton SF, trigger SF, b-tag SF, JES/JER, fake lepton estimate, or systematic uncertainty is included." << endl;
}
