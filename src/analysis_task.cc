// Created by Aleksandra on 11/02/21

#include "analysis_task.h"
#include <cmath>
#include <math.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCutG.h>

TASK_IMPL(AnalysisTask)
boost::program_options::options_description AnalysisTask::GetBoostOptions() {
    return UserTask::GetBoostOptions();
}

void AnalysisTask::UserInit(std::map<std::string, void *> &branch_map) {
    // linking pointers with branch fields
    event_header_ = GetInBranch("RecEventHeader");
    vtx_tracks_ = GetInBranch("VtxTracks");
    VtxTracks_ = static_cast<AnalysisTree::TrackDetector *>(branch_map.at("VtxTracks"));
    vtx_tracks_->UseFields({
                                   {"dcax", dcax},
                                   {"dcay", dcay},
                                   {"dcaz", dcaz},
                                   {"nhits_total", nhits_total},
                                   {"nhits_vtpc1", nhits_vtpc1},
                                   {"nhits_vtpc2", nhits_vtpc2},
                                   {"nhits_mtpc", nhits_mtpc},
                                   {"nhits_pot_total", nhits_pot_total},
                                   {"nhits_pot_vtpc1", nhits_pot_vtpc1},
                                   {"nhits_pot_vtpc2", nhits_pot_vtpc2},
                                   {"nhits_pot_mtpc", nhits_pot_mtpc},
                                   {"p", p},
                                   { "dedx_total", dedx},
                                   { "q", q}
                           });
    psd_modules_ = GetInBranch("PsdModules");
    psd_modules_->UseFields({
                                    {"signal", psd_signal},
                                    {"number", psd_number}
                            });
    std::tie(multiplicity_var_, energy_psd_var_, vtxX_var_, vtxY_var_, vtxZ_var_, wfa_t4_, wfa_s1_) =
            event_header_->GetVars( "M", "Epsd", "vtx_x", "vtx_y", "vtx_z", "wfa_t4", "wfa_s1");
    event_header_ ->UseFields({{"t1", evt_t1}, {"t2", evt_t2}, {"t4", evt_t4}});

    // Histo
    multiplicity_distribution_ = new TH1F( "multiplicity_distribution", ";M ;entries", 500, 0, 500 );
    psd_energy_distribution_ = new TH1F( "psd_energy_distribution", ";E_{PSD} [GeV];entries", 500, 0., 6000.0 );
    psd_energy_vs_mult_ = new TH2F ("psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 250, 0, 500, 600, 0., 6000.0);
    pT_distribution_ = new TH1F( "pT_distribution", ";p_{T} [GeV/c]; entries", 500, 0., 5 );
    phi_distribution_ = new TH1F( "phi_distribution", ";#varphi ; entries", 350, -3.5, 3.5 );
    eta_distribution_ = new TH1F ("eta_distribution", "; #eta; entries", 350, 0, 7.0);
    pT_vs_eta_ = new TH2F ("pT_vs_eta", "; #eta; p_{T} [GeV/c]",350, 0, 7.0, 300, 0., 3  );
    pT_vs_phi_ = new TH2F ("pT_vs_phi", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 300, 0., 3  );
    phi_vs_eta_ =  new TH2F ("phi_vs_eta", "; #varphi; #eta", 350, -3.5, 3.5, 350, 0, 7.0 );

    vtxX_distribution_ = new TH1F ( "vtxX_distribution", ";Vtx_X [cm]; entries", 300, -3.0, 3.0 );
    vtxY_distribution_ = new TH1F ( "vtxY_distribution", ";Vtx_Y [cm]; entries", 300, -4.0, 2.0 );
    vtxZ_distribution_ = new TH1F ( "vtxZ_distribution", ";Vtx_Z [cm]; entries", 400, -700.0, -500.0 );
    vtxX_vs_vtxY_ = new TH2F( "vtxX_vtxY_distribution", ";Vtx_X [cm];Vtx_Y [cm];", 100, -3.0, 3.0, 100, -4.0, 2.0);
    vtxZ_vs_mult_ = new TH2F ("vtxZ_vs_mult", ";Vtx_Z [cm]; M", 400, -700.0, -500.0 , 250, 0, 500 );
    cut1_mult_dist_ = new TH1F( "cut1_multiplicity_distribution", ";M ;entries", 500, 0, 500);
    cut2_mult_dist_ = new TH1F( "cut2_multiplicity_distribution", ";M ;entries", 500, 0, 500 );
    cut12_mult_dist_ = new TH1F( "cut12_multiplicity_distribution", ";M ;entries", 500, 0, 500);
    cut2_vtxX_distribution_ = new TH1F ( "cut2_vtxX_distribution", ";Vtx_X [cm]; entries", 150, -1.5, 1.5 );
    cut2_vtxY_distribution_= new TH1F ( "cut2_vtxY_distribution", ";Vtx_Y [cm]; entries", 150, -1.5, 1.5 );
    cut1_vtxZ_distribution_ = new TH1F ("cut1_vtxZ_distribution", ";Vtx_Z [cm]; entries", 20, -596.0, -586.0 );
    cut2_vtxX_vs_vtxY_ = new TH2F( "cut1_vtxX_vs_vtxY", ";Vtx_X [cm];Vtx_Y [cm];", 150, -1.5, 1.5, 150, -1.5, 1.5);
    cut1_vtxZ_vs_mult_ = new TH2F ("cut1_vtxZ_vs_mult", ";Vtx_Z [cm]; M", 100, -596.0, -586.0, 250, 0, 500 );
    cut1_psd_energy_vs_mult_ = new TH2F ("cut1_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 250, 0, 500, 600, 0., 6000.0);
    cut2_psd_energy_vs_mult_ = new TH2F ("cut2_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 250, 0, 500, 600, 0., 6000.0);
    cut12_psd_energy_vs_mult_ = new TH2F ("cut12_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 250, 0, 500, 600, 0., 6000.0);

    wfa_t4_distribution_ = new TH1F ("WFA_t4", ";WFA_{t4}; entries", 50, 0, 11000000);
    wfa_s1_distribution_ = new TH1F ("WFA_s1", ";WFA_{s1}; entries", 50, 0, 11000000);
    cut3_wfa_t4_distribution_ = new TH1F ("cut3_WFA_t4", ";WFA_{t4}; entries", 200, 10000000, 10020000);
    cut3_wfa_s1_distribution_ = new TH1F ("cut3_WFA_s1", ";WFA_{s1}; entries", 55, 0, 11000000);
    cut3_mult_dist_ = new TH1F( "cut3_multiplicity_distribution", ";M ;entries", 500, 0, 500 );
    cut3_psd_energy_vs_mult_ = new TH2F ("cut3_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 250, 0, 500, 600, 0., 6000.0);

    cut4_mult_dist_ = new TH1F( "cut4_multiplicity_distribution", ";M ;entries", 500, 0, 500 );
    cut4_psd_energy_vs_mult_ = new TH2F ("cut4_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 250, 0, 500, 600, 0., 6000.0);

    good_events_mult_dist_ = new TH1F( "good_events_multiplicity_distribution", ";M ;entries", 450, 0, 450 );
    good_events_psd_energy_vs_mult_ = new TH2F ("good_events_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 225, 0, 450, 600, 0., 3000.0);

    dcaX_distribution_ = new TH1F ( "dcaX", ";DCA_{x} [cm]; entries", 200, -5.0, 5.0 );
    dcaY_distribution_ = new TH1F ( "dcaY", ";DCA_{y} [cm]; entries", 200, -5.0, 5.0 );
    dcaX_vs_dcaY_ = new TH2F ( "dcaX_vs_dcaY", ";DCA_{x} [cm]; DCA_{y} [cm]", 200, -5.0, 5.0, 200, -5.0, 5.0 );
    cut5_dcaX_distribution_ = new TH1F ( "cut5_dcaX", ";DCA_{x} [cm]; entries", 80, -2.0, 2.0 );
    cut5_dcaY_distribution_ = new TH1F ( "cut5_dcaY", ";DCA_{y} [cm]; entries", 50, -1.0, 1.0 );
    cut5_dcaX_vs_dcaY_ = new TH2F ( "cut5_dcaX_vs_dcaY", ";DCA_{x} [cm]; DCA_{y} [cm]", 80, -2.0, 2.0, 50, -1.0, 1.0 );
    cut5_mult_dist_ = new TH1F( "cut5_multiplicity_distribution", ";M ;entries", 350, 0, 350  );
    cut5_mult_dist_scaled_ = new TH1F( "cut5_scaled_multiplicity_distribution", ";M ;entries", 100, 0, 500);
    cut5_psd_energy_vs_mult_ = new TH2F ("cut5_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 200, 0, 400 , 600, 0., 3000.0);

    nhits_total_distribution_ = new TH1F( "nhits_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    nhits_vtpc1_distribution_ = new TH1F( "nhits_vtpc1_distribution", ";N_{hits} ;entries", 80, 0, 80 );
    nhits_vtpc2_distribution_ = new TH1F( "nhits_vtpc2_distribution", ";N_{hits} ;entries", 80, 0, 80 );
    nhits_mtpc_distribution_ = new TH1F( "nhits_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //nhits_pot_total_distribution_ = new TH1F( "nhits_pot_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //nhits_pot_vtpc1_distribution_ = new TH1F( "nhits_pot_vtpc1_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //nhits_pot_vtpc2_distribution_ = new TH1F( "nhits_pot_vtpc2_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //nhits_pot_mtpc_distribution_ = new TH1F( "nhits_pot_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    nhits_divided_nhits_pot_ = new TH1F ("Nhits/Nhits_pot", "; N_{hits}/N_{hits_pot}; entries", 150, 0, 3.);
    cut6_nhits_divided_nhits_pot_ = new TH1F ("cut6_Nhits/Nhits_pot", "; N_{hits}/N_{hits_pot}; entries", 150, 0, 3.);
    nhits_total_vs_ratio_ = new TH2F ("Nhits_vs_ratio", "; N_{hits}; N_{hits}/N_{hits_pot}; entries", 250, 0, 250, 60, 0, 3.);
    cut6_nhits_total_vs_ratio_ = new TH2F ("cut6_Nhits_vs_ratio", "; N_{hits}; N_{hits}/N_{hits_pot}; entries", 250, 0, 250, 20, 0.5, 1.5);
    cut6_nhits_total_distribution_ = new TH1F( "cut6_nhits_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    cut6_nhits_vtpc1_distribution_ = new TH1F( "cut6_nhits_vtpc1_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    cut6_nhits_vtpc2_distribution_ = new TH1F( "cut6_nhits_vtpc2_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    cut6_nhits_mtpc_distribution_ = new TH1F( "cut6_nhits_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //cut6_nhits_pot_total_distribution_ = new TH1F( "cut6_nhits_pot_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //cut6_nhits_pot_vtpc1_distribution_ = new TH1F( "cut6_nhits_pot_vtpc1_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //cut6_nhits_pot_vtpc2_distribution_ = new TH1F( "cut6_nhits_pot_vtpc2_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //cut6_nhits_pot_mtpc_distribution_ = new TH1F( "cut6_nhits_pot_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    cut6_mult_dist_ = new TH1F( "cut6_multiplicity_distribution", ";M ;entries", 250, 0, 250 );
    cut6_mult_dist_scaled_ = new TH1F( "cut6_scaled_multiplicity_distribution", ";M ;entries", 100, 0, 500);
    cut6_psd_energy_vs_mult_ = new TH2F ("cut6_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 125, 0, 250 , 600, 0., 3000.0);

    dedx_distribution_ = new TH1F ("dEdx", "; dE/dx; entries", 150, 0.5, 3.5);
    p_distribution_ = new TH1F( "p_distribution", ";p [GeV/c]; entries", 100, 0, 30 );
    q_distribution_ = new TH1F( "q_distribution", ";q ; entries", 100, -5, 5 );
    dedx_vs_p_dist_ = new TH2F ("dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, -8., 8., 150, 0.5, 3.5);

    all_cut_mult_dist_ = new TH1F( "all_cut_multiplicity_distribution", ";M ;entries", 250, 0, 250 );
    all_cut_mult_dist_scaled_ = new TH1F( "all_cut_scaled_multiplicity_distribution", ";M ;entries", 100, 0, 500);
    all_cut_psd_energy_vs_mult_ = new TH2F ("all_cut_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 125, 0, 250, 600, 0., 3000.0);
    all_cut_p_distribution_ = new TH1F( "all_cut_p_distribution", ";p [GeV/c]; entries", 300, 0., 15 );
    all_cut_dedx_vs_p_dist_ = new TH2F ("all_cut_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, -8., 8., 150, 0.5, 3.5);

    pions_multiplicity_dist_ = new TH1F( "pions_multiplicity_distribution", ";M ;entries", 100, 0, 100 );
    pions_multiplicity_dist_scaled_ = new TH1F( "pions_scaled_multiplicity_distribution", ";M ;entries", 50, 0, 250 );
    protons_multiplicity_dist_ = new TH1F( "protons_multiplicity_distribution", ";M ;entries", 100, 0, 100 );
    protons_multiplicity_dist_scaled_ = new TH1F( "protons_scaled_multiplicity_distribution", ";M ;entries", 50, 0, 250 );
    protons_all_multiplicity_dist_ = new TH1F( "protons_all_multiplicity_distribution", ";M ;entries", 100, 0, 100 );
    protons_all_multiplicity_dist_scaled_ = new TH1F( "protons_all_scaled_multiplicity_distribution", ";M ;entries", 50, 0, 250 );

    pions_dedx_vs_p_dist_ = new TH2F ("pions_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, -8., 0.,150, 0.5, 3.5);
    pions_pT_distribution_ = new TH1F( "pions_pT_distribution", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    pions_phi_distribution_ = new TH1F( "pions_phi_distribution", ";#varphi ; entries", 350, -3.5, 3.5 );
    pions_eta_distribution_ = new TH1F ("pions_eta_distribution", "; #eta; entries", 300, 1.0, 7.0);
    pions_pT_vs_eta_ = new TH2F ("pions_pT_vs_eta", "; #eta; p_{T} [GeV/c]",300, 1.0, 7.0, 150, 0., 3  );
    pions_pT_vs_phi_ = new TH2F ("pions_pT_vs_phi", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_phi_vs_eta_ =  new TH2F ("pions_phi_vs_eta", " ; #varphi; #eta", 350, -3.5, 3.5, 300, 1.0, 7.0 );
    pions_rapidity_ = new TH1F ("pions_rapidity", "; #y; entries", 250, 0, 5.0);
    pions_rapidity_vs_pT_ = new TH2F ("pions_pT_vs_y", "; y; p_{T} [GeV/c]",100, 0, 5.0, 150, 0., 3 );
    pions_rapidity_vs_phi = new TH2F ("pions_phi_vs_y", " ; #varphi; y", 350, -3.5, 3.5, 100, 0, 5.0);
    protons_dedx_vs_p_dist_ = new TH2F ("protons_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, 0., 8.,150, 0.5, 3.5);
    protons_pT_distribution_ = new TH1F( "protons_pT_distribution", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    protons_phi_distribution_ = new TH1F( "protons_phi_distribution", ";#varphi ; entries", 350, -3.5, 3.5 );
    protons_eta_distribution_ = new TH1F ("protons_eta_distribution", "; #eta; entries", 300, 1.0, 7.0);
    protons_pT_vs_eta_ = new TH2F ("protons_pT_vs_eta", "; #eta; p_{T} [GeV/c]",300, 1.0, 7.0, 150, 0., 3  );
    protons_pT_vs_phi_ = new TH2F ("protons_pT_vs_phi", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    protons_phi_vs_eta_ =  new TH2F ("protons_phi_vs_eta", " ; #varphi; #eta", 350, -3.5, 3.5, 300, 1.0, 7.0 );
    protons_rapidity_ = new TH1F ("protons_rapidity", "; #y; entries", 250, 0, 5.0);
    protons_rapidity_vs_pT_ = new TH2F ("protons_pT_vs_y", "; y; p_{T} [GeV/c]",100, 0, 5.0, 150, 0., 3 );
    protons_rapidity_vs_phi = new TH2F ("protons_phi_vs_y", " ; #varphi; y", 350, -3.5, 3.5, 100, 0, 5.0);
    protons_all_dedx_vs_p_dist_ = new TH2F ("protons_all_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, 0., 8.,150, 0.5, 3.5);
    protons_all_pT_distribution_ = new TH1F( "protons_all_pT_distribution", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    protons_all_phi_distribution_ = new TH1F( "protons_all_phi_distribution", ";#varphi ; entries", 350, -3.5, 3.5 );
    protons_all_eta_distribution_ = new TH1F ("protons_all_eta_distribution", "; #eta; entries", 300, 1.0, 7.0);
    protons_all_pT_vs_eta_ = new TH2F ("protons_all_pT_vs_eta", "; #eta; p_{T} [GeV/c]",300, 1.0, 7.0, 150, 0., 3  );
    protons_all_pT_vs_phi_ = new TH2F ("protons_all_pT_vs_phi", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    protons_all_phi_vs_eta_ =  new TH2F ("protons_all_phi_vs_eta", " ; #varphi; #eta", 350, -3.5, 3.5, 300, 1.0, 7.0 );
    protons_all_rapidity_ = new TH1F ("protons_all_rapidity", "; #y; entries", 250, 0, 5.0);
    protons_all_rapidity_vs_pT_ = new TH2F ("protons_all_pT_vs_y", "; y; p_{T} [GeV/c]",100, 0, 5.0, 150, 0., 3 );
    protons_all_rapidity_vs_phi = new TH2F ("protons_all_phi_vs_y", " ; #varphi; y", 350, -3.5, 3.5, 100, 0, 5.0);

}

void AnalysisTask::UserExec() {
    using AnalysisTree::Track;
    using AnalysisTree::Module;

    TCutG *cut_pi = new TCutG("cut_pi",6);
    cut_pi->SetPoint(0, -5.26214, 1.41964);
    cut_pi->SetPoint(1, -5.26214, 1.19643 );
    cut_pi->SetPoint(2, -2.25243, 0.883929);
    cut_pi->SetPoint(3, -2.25243, 1.10714 );
    cut_pi->SetPoint(4, -5.24272, 1.41964);
    cut_pi->SetPoint(5, -5.26214, 1.41964);

    TCutG *cut_all_p = new TCutG("cut_all_p", 15);
    cut_all_p->SetPoint(0, 2.27184, 3.48065);
    cut_all_p->SetPoint(1, 2.54369, 2.20833);
    cut_all_p->SetPoint(2, 3.08738, 1.27083);
    cut_all_p->SetPoint(3, 4.07767, 0.883929);
    cut_all_p->SetPoint(4, 5.06796, 0.958333);
    cut_all_p->SetPoint(5, 5.78641, 1.08482);
    cut_all_p->SetPoint(6, 5.80583, 1.20387);
    cut_all_p->SetPoint(7, 4.58252, 1.08482);
    cut_all_p->SetPoint(8, 3.98058, 1.06994);
    cut_all_p->SetPoint(9, 3.30097, 1.29315);
    cut_all_p->SetPoint(10, 2.73786, 2.14137);
    cut_all_p->SetPoint(11, 2.4466, 3.47321);
    cut_all_p->SetPoint(12, 2.31068, 3.46577);
    cut_all_p->SetPoint(13, 2.29126, 3.46577);
    cut_all_p->SetPoint(14, 2.27184, 3.48065);

    TCutG *cut_p = new TCutG("cut_protons",6);
    cut_p->SetPoint(0, 4.10163, 1.03274);
    cut_p->SetPoint(1, 5.80762, 1.20387 );
    cut_p->SetPoint(2, 5.80762, 1.0625);
    cut_p->SetPoint(3, 4.11978, 0.89881 );
    cut_p->SetPoint(4, 4.11978, 1.0253);
    cut_p->SetPoint(5, 4.10163, 1.03274);

    // work with minbias events (T4); commit(change) trigger condition, if you want to analyse all(other trigger) events
    if ((*event_header_)[evt_t4].GetBool()){
        int n_tracks_ = 0, n_tracks_cut5 = 0, n_tracks_cut6 = 0, n_tracks_all_cut = 0, n_pions_minus=0, n_protons_all = 0, n_protons =0;    //to calculate multiplicity of different events
        bool cut1_ = false, cut2_ = false, cut12_ = false, cut3_ = false, cut4_= false;     //cut flag
        //calculate multiplicity of event without cuts and fill kinematic histo
        for( auto& track : vtx_tracks_->Loop() ) {
            n_tracks_++;
            auto pT = track.DataT<Track>()->GetPt();
            auto phi = track.DataT<Track>()->GetPhi();
            auto eta = track.DataT<Track>()->GetEta();
            pT_distribution_->Fill(pT);
            eta_distribution_->Fill(eta);
            phi_distribution_->Fill(phi);
            pT_vs_eta_->Fill(eta, pT);
            pT_vs_phi_->Fill(phi, pT);
            phi_vs_eta_->Fill(phi, eta);
        }
        // Variables
        auto energy_psd = (*event_header_)[energy_psd_var_].GetVal();
        auto vtxZ_ = ((*event_header_)[vtxZ_var_].GetVal());
        auto vtxX_ = (*event_header_)[vtxX_var_].GetVal();
        auto vtxY_ = (*event_header_)[vtxY_var_].GetVal();
        auto wfaT4_ = (*event_header_)[wfa_t4_].GetVal();
        auto wfaS1_ = (*event_header_)[wfa_s1_].GetVal();
        //histo before event cuts:
        multiplicity_distribution_->Fill(n_tracks_);
        psd_energy_distribution_ ->Fill(energy_psd);
        psd_energy_vs_mult_ -> Fill(n_tracks_, energy_psd);
        vtxZ_distribution_->Fill(vtxZ_);
        vtxX_distribution_->Fill(vtxX_);
        vtxY_distribution_->Fill(vtxY_);
        vtxZ_vs_mult_->Fill(vtxZ_,n_tracks_);
        vtxX_vs_vtxY_->Fill(vtxX_,vtxY_);
        wfa_t4_distribution_->Fill(wfaT4_);
        wfa_s1_distribution_->Fill(wfaS1_);
        //Events cuts:
        //cut1 vtxZ
        if ( abs(vtxZ_ + 591.8) <= 2. ){
            cut1_ = true;
            cut1_vtxZ_distribution_->Fill(vtxZ_);
            cut1_vtxZ_vs_mult_->Fill(vtxZ_ ,n_tracks_);
            cut1_mult_dist_->Fill(n_tracks_);
            cut1_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
        }
        //cut2 vtx_X & vtx_Y
        double x_min = -0.3746 - 0.2157, x_max = -0.3746 + 0.2157;
        double y_min = -0.2835 - 0.1531, y_max = -0.2835 + 0.1531;
        if ( (vtxX_>=x_min && vtxX_<=x_max) && (vtxY_>=y_min && vtxY_<=y_max) ){
            cut2_ = true;
            cut2_vtxX_vs_vtxY_-> Fill(vtxX_, vtxY_);
            cut2_vtxX_distribution_-> Fill(vtxX_);
            cut2_vtxY_distribution_-> Fill(vtxY_);
            cut2_mult_dist_-> Fill(n_tracks_);
            cut2_psd_energy_vs_mult_ -> Fill(n_tracks_, energy_psd);
        }
        //cut 1 && 2
        if (cut1_ && cut2_){
            cut12_ = true;
            cut12_mult_dist_->Fill(n_tracks_);
            cut12_psd_energy_vs_mult_ ->Fill(n_tracks_, energy_psd);
        }
        //cut3 wfa_t4 & wfa_s1
        if ( ((*event_header_)[wfa_t4_].GetVal()>=25000) && (wfaS1_>=4000)){
            cut3_ = true;
            cut3_mult_dist_->Fill(n_tracks_);
            cut3_psd_energy_vs_mult_ ->Fill(n_tracks_, energy_psd);
            cut3_wfa_t4_distribution_->Fill(wfaT4_);
            cut3_wfa_s1_distribution_->Fill(wfaS1_);
        }
        //cut 4: graphical cut
        if ( energy_psd <= (-310./37)*n_tracks_+3100){
            cut4_ = true;
            cut4_mult_dist_->Fill(n_tracks_);
            cut4_psd_energy_vs_mult_ ->Fill(n_tracks_, energy_psd);
        }
        //For good events:
        if (cut12_ && cut3_ && cut4_) {
            //Track cuts
            for (auto &track: vtx_tracks_->Loop()) {
                //initialise variables
                bool cut5_ = false, cut6_ = false;
                auto p_ = track[p].GetVal();
                auto dEdx_ = track[dedx].GetVal();
                auto q_ = track[q].GetVal();
                auto dcaX = track[dcax].GetVal();
                auto dcaY = track[dcay].GetVal();
                int nhits[4];
                nhits[0] = track[nhits_vtpc1].GetVal();
                nhits[1] = track[nhits_vtpc2].GetVal();
                nhits[2] = track[nhits_mtpc].GetVal();
                nhits[3] = track[nhits_total].GetVal();
                int nhits_pot[4];
                nhits_pot[0] = track[nhits_pot_vtpc1].GetVal();
                nhits_pot[1] = track[nhits_pot_vtpc2].GetVal();
                nhits_pot[2] = track[nhits_pot_mtpc].GetVal();
                nhits_pot[3] = track[nhits_pot_total].GetVal();
                float ratio = 1. * nhits[3] / nhits_pot[3];

                q_distribution_->Fill(q_);
                p_distribution_->Fill(p_);
                dedx_distribution_->Fill(dEdx_);

                dcaX_distribution_->Fill(dcaX);
                dcaY_distribution_->Fill(dcaY);
                dcaX_vs_dcaY_->Fill(dcaX, dcaY);
                dedx_vs_p_dist_->Fill(q_ * log(20 * p_), dEdx_);

                nhits_total_distribution_->Fill(nhits[3]);
                nhits_vtpc1_distribution_->Fill(nhits[0]);
                nhits_vtpc2_distribution_->Fill(nhits[1]);
                nhits_mtpc_distribution_->Fill(nhits[2]);
                nhits_divided_nhits_pot_->Fill(ratio);
                nhits_total_vs_ratio_->Fill(nhits[3], ratio);
                //cut 5: dcaX vs dcaY
                if (abs(dcaX) <= 2. && abs(dcaY) <= 1.) {
                    cut5_ = true;
                    n_tracks_cut5++;
                    cut5_dcaX_distribution_->Fill(dcaX);
                    cut5_dcaY_distribution_->Fill(dcaY);
                    cut5_dcaX_vs_dcaY_->Fill(dcaX, dcaY);
                }
                //cut 6: Nhits
                if (((nhits[0] + nhits[1]) > 15) && (nhits[3] > 30) && (ratio > 0.55) && (ratio < 1.1)) {
                    cut6_ = true;
                    n_tracks_cut6++;
                    cut6_nhits_total_distribution_->Fill(nhits[3]);
                    cut6_nhits_vtpc1_distribution_->Fill(nhits[0]);
                    cut6_nhits_vtpc2_distribution_->Fill(nhits[1]);
                    cut6_nhits_mtpc_distribution_->Fill(nhits[2]);
                    cut6_nhits_divided_nhits_pot_->Fill(ratio);
                    cut6_nhits_total_vs_ratio_->Fill(nhits[3], ratio);
                }
                auto pT = track.DataT<Track>()->GetPt();
                auto phi = track.DataT<Track>()->GetPhi();
                auto eta = track.DataT<Track>()->GetEta();
                //all event and track cuts + pions
                if (cut5_ && cut6_) {
                    n_tracks_all_cut++;
                    all_cut_p_distribution_->Fill(p_);
                    all_cut_dedx_vs_p_dist_->Fill(q_ * log(20 * p_), dEdx_);
                    auto p_beam = 13;
                    double beta = p_beam / sqrt(p_beam * p_beam + 0.9315 * 0.9315);
                    double y_cm = atanh(beta / 2);
                    double log_p = q_ * log(20 * p_);
                    if (cut_pi->IsInside(log_p, dEdx_)) {
                        n_pions_minus++;
                        auto y = track.DataT<Track>()->GetRapidity(-211);
                        pions_dedx_vs_p_dist_->Fill(log_p, dEdx_);
                        pions_rapidity_->Fill(y - y_cm);
                        pions_rapidity_vs_phi ->Fill (phi, y);
                        pions_rapidity_vs_pT_ ->Fill(y, pT);
                        pions_pT_distribution_->Fill(pT);
                        pions_phi_distribution_->Fill(phi);
                        pions_eta_distribution_->Fill(eta);
                        pions_pT_vs_eta_->Fill(eta, pT);
                        pions_pT_vs_phi_->Fill(phi, pT);
                        pions_phi_vs_eta_->Fill(phi, eta);
                    }
                    if (cut_all_p->IsInside(log_p, dEdx_)) {
                        n_protons_all++;
                        auto y = track.DataT<Track>()->GetRapidity(2212);
                        protons_all_dedx_vs_p_dist_->Fill(log_p, dEdx_);
                        protons_all_rapidity_->Fill(y - y_cm);
                        protons_all_eta_distribution_->Fill(eta);
                        protons_all_rapidity_vs_pT_ ->Fill(y, pT);
                        protons_all_rapidity_vs_phi ->Fill(phi, y);
                        protons_all_pT_distribution_->Fill(pT);
                        protons_all_phi_distribution_->Fill(phi);
                        protons_all_pT_vs_eta_->Fill(eta, pT);
                        protons_all_pT_vs_phi_->Fill(phi, pT);
                        protons_all_phi_vs_eta_->Fill(phi, eta);
                    }
                    if (cut_p->IsInside(log_p, dEdx_)) {
                        n_protons++;
                        auto y = track.DataT<Track>()->GetRapidity(2212);
                        protons_dedx_vs_p_dist_->Fill(log_p, dEdx_);
                        protons_eta_distribution_->Fill(eta);
                        protons_rapidity_->Fill(y - y_cm);
                        protons_rapidity_vs_pT_ ->Fill(y, pT);
                        protons_rapidity_vs_phi ->Fill(phi, y);
                        protons_pT_distribution_->Fill(pT);
                        protons_phi_distribution_->Fill(phi);
                        protons_pT_vs_eta_->Fill(eta, pT);
                        protons_pT_vs_phi_->Fill(phi, pT);
                        protons_phi_vs_eta_->Fill(phi, eta);
                    }

                }

            }

            good_events_mult_dist_->Fill(n_tracks_);
            good_events_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
            pions_multiplicity_dist_->Fill(n_pions_minus);
            pions_multiplicity_dist_scaled_->Fill(1.0*n_pions_minus*140/45);
            protons_all_multiplicity_dist_->Fill (n_protons_all);
            protons_all_multiplicity_dist_scaled_->Fill(1.0*n_protons_all*140/40);
            protons_multiplicity_dist_->Fill (n_protons);
            protons_multiplicity_dist_scaled_->Fill(1.0*n_protons*140/40);
            good_events_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
            //fill multiplicity histo for track cut
            cut5_mult_dist_->Fill(n_tracks_cut5);
            cut5_mult_dist_scaled_->Fill(1.0*n_tracks_cut5*320/230);
            cut5_psd_energy_vs_mult_->Fill(n_tracks_cut5, energy_psd);
            cut6_mult_dist_->Fill(n_tracks_cut6);
            cut6_mult_dist_scaled_->Fill(1.0*n_tracks_cut6*320/160);
            cut6_psd_energy_vs_mult_->Fill(n_tracks_cut6, energy_psd);
            all_cut_mult_dist_->Fill(n_tracks_all_cut);
            all_cut_mult_dist_scaled_->Fill(1.0*n_tracks_all_cut*320/140);
            all_cut_psd_energy_vs_mult_->Fill(n_tracks_all_cut, energy_psd);
        }
    }


}

void AnalysisTask::UserFinish() {
    // Writing histograms to file
    out_file_->cd();

    TCutG *cut_pi = new TCutG("cut_pi",6);
    cut_pi->SetPoint(0, -5.26214, 1.41964);
    cut_pi->SetPoint(1, -5.26214, 1.19643 );
    cut_pi->SetPoint(2, -2.25243, 0.883929);
    cut_pi->SetPoint(3, -2.25243, 1.10714 );
    cut_pi->SetPoint(4, -5.24272, 1.41964);
    cut_pi->SetPoint(5, -5.26214, 1.41964);

    TCutG *cut_all_p = new TCutG("cut_all_p", 15);
    cut_all_p->SetPoint(0, 2.27184, 3.48065);
    cut_all_p->SetPoint(1, 2.54369, 2.20833);
    cut_all_p->SetPoint(2, 3.08738, 1.27083);
    cut_all_p->SetPoint(3, 4.07767, 0.883929);
    cut_all_p->SetPoint(4, 5.06796, 0.958333);
    cut_all_p->SetPoint(5, 5.78641, 1.08482);
    cut_all_p->SetPoint(6, 5.80583, 1.20387);
    cut_all_p->SetPoint(7, 4.58252, 1.08482);
    cut_all_p->SetPoint(8, 3.98058, 1.06994);
    cut_all_p->SetPoint(9, 3.30097, 1.29315);
    cut_all_p->SetPoint(10, 2.73786, 2.14137);
    cut_all_p->SetPoint(11, 2.4466, 3.47321);
    cut_all_p->SetPoint(12, 2.31068, 3.46577);
    cut_all_p->SetPoint(13, 2.29126, 3.46577);
    cut_all_p->SetPoint(14, 2.27184, 3.48065);

    TCutG *cut_p = new TCutG("cut_protons",6);
    cut_p->SetPoint(0, 4.10163, 1.03274);
    cut_p->SetPoint(1, 5.80762, 1.20387 );
    cut_p->SetPoint(2, 5.80762, 1.0625);
    cut_p->SetPoint(3, 4.11978, 0.89881 );
    cut_p->SetPoint(4, 4.11978, 1.0253);
    cut_p->SetPoint(5, 4.10163, 1.03274);

    cut_pi->Write();
    cut_p->Write();
    cut_all_p->Write();

    multiplicity_distribution_->Write();
    psd_energy_distribution_->Write();
    psd_energy_vs_mult_->Write();
    pT_distribution_->Write();
    eta_distribution_->Write();
    phi_distribution_->Write();
    pT_vs_eta_->Write();
    pT_vs_phi_->Write();
    phi_vs_eta_->Write();
    vtxX_distribution_->Write();
    vtxY_distribution_->Write();
    vtxZ_distribution_->Write();
    vtxX_vs_vtxY_->Write();
    vtxZ_vs_mult_->Write();

    cut1_mult_dist_->Write();
    cut1_psd_energy_vs_mult_ ->Write();
    cut1_vtxZ_distribution_->Write();
    cut1_vtxZ_vs_mult_->Write();

    cut2_mult_dist_->Write();
    cut2_psd_energy_vs_mult_ ->Write();

    cut12_mult_dist_->Write();
    cut12_psd_energy_vs_mult_->Write();

    wfa_t4_distribution_->Write();
    wfa_s1_distribution_->Write();
    cut3_wfa_t4_distribution_->Write();
    cut3_wfa_s1_distribution_->Write();
    cut3_mult_dist_->Write();
    cut3_psd_energy_vs_mult_ ->Write();

    cut4_mult_dist_ ->Write();
    cut4_psd_energy_vs_mult_ ->Write();
    good_events_psd_energy_vs_mult_->Write();
    good_events_mult_dist_->Write();

    dcaX_distribution_ ->Write();
    dcaY_distribution_ ->Write();
    dcaX_vs_dcaY_ ->Write();
    cut5_dcaX_distribution_ ->Write();
    cut5_dcaY_distribution_ ->Write();
    cut5_dcaX_vs_dcaY_ ->Write();
    cut5_mult_dist_->Write();
    cut5_mult_dist_scaled_->Write();
    cut5_psd_energy_vs_mult_->Write();

    nhits_total_distribution_ ->Write();
    nhits_vtpc1_distribution_ ->Write();
    nhits_vtpc2_distribution_ ->Write();
    nhits_mtpc_distribution_ ->Write();
    nhits_divided_nhits_pot_ ->Write();
    nhits_total_vs_ratio_ ->Write();
    cut6_nhits_total_distribution_ ->Write();
    cut6_nhits_vtpc1_distribution_ ->Write();
    cut6_nhits_vtpc2_distribution_ ->Write();
    cut6_nhits_mtpc_distribution_ ->Write();
    cut6_nhits_divided_nhits_pot_ ->Write();
    cut6_nhits_total_vs_ratio_ ->Write();

    cut6_mult_dist_ ->Write();
    cut6_mult_dist_scaled_->Write();
    cut6_psd_energy_vs_mult_ ->Write();
    all_cut_mult_dist_->Write();
    all_cut_mult_dist_scaled_->Write();
    all_cut_psd_energy_vs_mult_->Write();
    all_cut_dedx_vs_p_dist_->Write();
    q_distribution_->Write();
    p_distribution_->Write();
    dedx_distribution_->Write();
    dedx_vs_p_dist_->Write();
    all_cut_p_distribution_->Write();

    pions_multiplicity_dist_->Write();
    pions_multiplicity_dist_scaled_->Write();
    pions_pT_distribution_->Write();
    pions_eta_distribution_->Write();
    pions_phi_distribution_->Write();
    pions_pT_vs_eta_->Write();
    pions_pT_vs_phi_->Write();
    pions_phi_vs_eta_->Write();
    pions_dedx_vs_p_dist_->Write();
    pions_rapidity_->Write();
    pions_rapidity_vs_phi ->Write();
    pions_rapidity_vs_pT_ ->Write();

    protons_multiplicity_dist_->Write();
    protons_multiplicity_dist_scaled_->Write();
    protons_pT_distribution_->Write();
    protons_eta_distribution_->Write();
    protons_phi_distribution_->Write();
    protons_pT_vs_eta_->Write();
    protons_pT_vs_phi_->Write();
    protons_phi_vs_eta_->Write();
    protons_dedx_vs_p_dist_->Write();
    protons_rapidity_->Write();
    protons_rapidity_vs_phi ->Write();
    protons_rapidity_vs_pT_ ->Write();

    protons_all_multiplicity_dist_->Write();
    protons_all_multiplicity_dist_scaled_->Write();
    protons_all_pT_distribution_->Write();
    protons_all_eta_distribution_->Write();
    protons_all_phi_distribution_->Write();
    protons_all_pT_vs_eta_->Write();
    protons_all_pT_vs_phi_->Write();
    protons_all_phi_vs_eta_->Write();
    protons_all_dedx_vs_p_dist_->Write();
    protons_all_rapidity_->Write();
    protons_all_rapidity_vs_phi ->Write();
    protons_all_rapidity_vs_pT_ ->Write();
}
