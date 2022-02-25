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
    //psd_modules_ = GetInBranch("PsdModules");

    auto branch_Mch_good_events = tree->Branch("M_good_events",  &M_good, "M_good/I");
    auto branch_Mch_cut_dca = tree->Branch("M_cut_dca",  &M_cut_dca, "M_cut_dca/I");
    auto branch_Mch_all_cut = tree->Branch("Mch",  &M_ch, "M_ch/I");
    auto branch_Mpi = tree->Branch("Mpions",  &M_pi, "M_pi/I");
    auto branch_Mp= tree->Branch("Mprotons",  &M_p, "M_p/I");

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
    psd_energy_distribution_central_ = new TH1F( "central_psd_energy_distribution", ";E_{PSD_central} [GeV];entries", 500, 0., 5000.0 );
    psd_energy_distribution_pherical_ = new TH1F( "pherical_psd_energy_distribution", ";E_{PSD_pherical} [GeV];entries", 100, 0., 600.0 );
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

    good_events_mult_dist_scaled_ = new TH1F( "good_events_multiplicity_distribution", ";M/300 ;entries", 200, 0., 2. );
    good_events_psd_energy_vs_mult_ = new TH2F ("good_events_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 225, 0, 450, 720, 0., 3600.0);

    dcaX_distribution_ = new TH1F ( "dcaX", ";DCA_{x} [cm]; entries", 200, -5.0, 5.0 );
    dcaY_distribution_ = new TH1F ( "dcaY", ";DCA_{y} [cm]; entries", 200, -5.0, 5.0 );
    dcaX_vs_dcaY_ = new TH2F ( "dcaX_vs_dcaY", ";DCA_{x} [cm]; DCA_{y} [cm]", 200, -5.0, 5.0, 200, -5.0, 5.0 );
    cut5_dcaX_distribution_ = new TH1F ( "cut5_dcaX", ";DCA_{x} [cm]; entries", 80, -2.0, 2.0 );
    cut5_dcaY_distribution_ = new TH1F ( "cut5_dcaY", ";DCA_{y} [cm]; entries", 50, -1.0, 1.0 );
    cut5_dcaX_vs_dcaY_ = new TH2F ( "cut5_dcaX_vs_dcaY", ";DCA_{x} [cm]; DCA_{y} [cm]", 80, -2.0, 2.0, 50, -1.0, 1.0 );
    cut5_mult_dist_ = new TH1F( "cut5_multiplicity_distribution", ";M ;entries", 350, 0, 350  );
    cut5_mult_dist_scaled_ = new TH1F( "cut5_scaled_multiplicity_distribution", ";M/230 ;entries", 200, 0., 2.);
    cut5_psd_energy_vs_mult_ = new TH2F ("cut5_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 200, 0, 400 , 720, 0., 3600.0);

    nhits_total_distribution_ = new TH1F( "nhits_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    nhits_vtpc1_distribution_ = new TH1F( "nhits_vtpc1_distribution", ";N_{hits} ;entries", 80, 0, 80 );
    nhits_vtpc2_distribution_ = new TH1F( "nhits_vtpc2_distribution", ";N_{hits} ;entries", 80, 0, 80 );
    nhits_mtpc_distribution_ = new TH1F( "nhits_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //nhits_pot_total_distribution_ = new TH1F( "nhits_pot_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //nhits_pot_vtpc1_distribution_ = new TH1F( "nhits_pot_vtpc1_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //nhits_pot_vtpc2_distribution_ = new TH1F( "nhits_pot_vtpc2_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //nhits_pot_mtpc_distribution_ = new TH1F( "nhits_pot_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    nhits_divided_nhits_pot_ = new TH1F ("Nhits/Nhits_pot", "; N_{hits}/N_{hits_pot}; entries", 150, 0, 3.);
    nhits_total_vs_potential_ = new TH2F ("Nhits_vs_Nhits_pot", "; N_{hits}; N_{hits_pot}", 250, 0, 250, 250, 0, 250);
    cut6_nhits_divided_nhits_pot_ = new TH1F ("cut6_Nhits/Nhits_pot", "; N_{hits}/N_{hits_pot}; entries", 150, 0, 3.);
    cut6_nhits_total_distribution_ = new TH1F( "cut6_nhits_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    cut6_nhits_vtpc1_distribution_ = new TH1F( "cut6_nhits_vtpc1_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    cut6_nhits_vtpc2_distribution_ = new TH1F( "cut6_nhits_vtpc2_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    cut6_nhits_mtpc_distribution_ = new TH1F( "cut6_nhits_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    cut6_nhits_total_vs_potential_ = new TH2F ("cut6_Nhits_vs_Nhits_pot", "; N_{hits}; N_{hits_pot}", 250, 0, 250, 250, 0, 250);
    //cut6_nhits_pot_total_distribution_ = new TH1F( "cut6_nhits_pot_total_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //cut6_nhits_pot_vtpc1_distribution_ = new TH1F( "cut6_nhits_pot_vtpc1_distribution", ";N_{hits} ;entries", 250, 0, 250 );
    //cut6_nhits_pot_vtpc2_distribution_ = new TH1F( "cut6_nhits_pot_vtpc2_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    //cut6_nhits_pot_mtpc_distribution_ = new TH1F( "cut6_nhits_pot_mptc_distribution", ";N_{hits} ;entries", 100, 0, 100 );
    cut6_mult_dist_ = new TH1F( "cut6_multiplicity_distribution", ";M ;entries", 250, 0, 250 );
    cut6_mult_dist_scaled_ = new TH1F( "cut6_scaled_multiplicity_distribution", ";M/160 ;entries", 200, 0., 2.);
    cut6_psd_energy_vs_mult_ = new TH2F ("cut6_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 125, 0, 250 , 720, 0., 3600.0);

    dedx_distribution_ = new TH1F ("dEdx", "; dE/dx; entries", 150, 0.5, 3.5);
    p_distribution_ = new TH1F( "p_distribution", ";p [GeV/c]; entries", 600, 0, 30 );
    q_distribution_ = new TH1F( "q_distribution", ";q ; entries", 100, -5, 5 );
    dedx_vs_p_dist_ = new TH2F ("dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, -8., 8., 150, 0.5, 3.5);

    all_cut_mult_dist_ = new TH1F( "all_cut_multiplicity_distribution", ";M ;entries", 250, 0, 250 );
    all_cut_mult_dist_scaled_ = new TH1F( "all_cut_scaled_multiplicity_distribution", ";M/140 ;entries", 200, 0., 2.);
    all_cut_psd_energy_vs_mult_ = new TH2F ("all_cut_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 125, 0, 250, 720, 0., 3600.0);
    all_cut_p_distribution_ = new TH1F( "all_cut_p_distribution", ";p [GeV/c]; entries", 600, 0., 30 );
    all_cut_dedx_vs_p_dist_ = new TH2F ("all_cut_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, -8., 8., 150, 0.5, 3.5);
    all_cut_psd_energy_central_vs_mult_ = new TH2F ("all_cut_central_psd_energy_vs_mult", "; M;E_{PSD_central} [GeV]", 125, 0, 250, 720, 0., 3600.0);
    all_cut_psd_energy_pherical_vs_mult_ = new TH2F ("all_cut_pherical_psd_energy_vs_mult", "; M;E_{PSD_pherical} [GeV]", 125, 0, 250, 100, 0., 600.0 );

    pions_multiplicity_dist_ = new TH1F( "pions_multiplicity_distribution", ";M ;entries", 100, 0, 100 );
    pions_multiplicity_dist_scaled_ = new TH1F( "pions_scaled_multiplicity_distribution", ";M/50 ;entries", 200, 0., 2. );
    pions_mult_vs_psd_energy_ = new TH2F ("pions_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 100, 0, 100, 720, 0., 3600.0);
    pions_mult_vs_psd_energy_central_ = new TH2F ("pions_central_psd_energy_vs_mult", "; M;E_{PSD_central} [GeV]", 100, 0, 100, 720, 0., 3600.0);
    pions_mult_vs_psd_energy_pherical_ = new TH2F ("pions_pherical_psd_energy_vs_mult", "; M;E_{PSD_pherical} [GeV]", 100, 0, 100, 100, 0., 600.0 );
    protons_multiplicity_dist_ = new TH1F( "protons_multiplicity_distribution", ";M ;entries", 100, 0, 100 );
    protons_multiplicity_dist_scaled_ = new TH1F( "protons_scaled_multiplicity_distribution", ";M/45 ;entries", 200, 0., 2. );
    protons_mult_vs_psd_energy_ = new TH2F ("protons_psd_energy_vs_mult", "; M;E_{PSD} [GeV]", 100, 0, 100, 720, 0., 3600.0);
    protons_mult_vs_psd_energy_central_ = new TH2F ("protons_central_psd_energy_vs_mult", "; M;E_{PSD_central} [GeV]", 100, 0, 100, 720, 0., 3600.0);
    protons_mult_vs_psd_energy_pherical_ = new TH2F ("protons_pherical_psd_energy_vs_mult", "; M;E_{PSD_pherical} [GeV]", 100, 0, 100, 100, 0., 600.0 );

    pions_dedx_vs_p_dist_ = new TH2F ("pions_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, -8., 0.,150, 0.5, 3.5);
    pions_p_distribution_ = new TH1F( "pions_p_distribution", ";p [GeV/c]; entries", 300, 0., 15 );
    pions_pT_distribution_ = new TH1F( "pions_pT_distribution", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    pions_pT_distribution_cut_y_pt1_ = new TH1F( "pions_pT_cut_y_pt1_ ", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    pions_pT_distribution_cut_y_pt2_ = new TH1F( "pions_pT_cut_y_pt2_ ", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    pions_phi_distribution_ = new TH1F( "pions_phi_distribution", ";#varphi ; entries", 350, -3.5, 3.5 );
    pions_eta_distribution_ = new TH1F ("pions_eta_distribution", "; #eta; entries", 300, 1.0, 7.0);
    pions_pT_vs_eta_ = new TH2F ("pions_pT_vs_eta", "; #eta; p_{T} [GeV/c]",300, 1.0, 7.0, 150, 0., 3  );
    pions_pT_vs_phi_ = new TH2F ("pions_pT_vs_phi", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_pT_vs_phi_cut_y1_ = new TH2F ("pions_pT_vs_phi_cut_y1_", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_pT_vs_phi_cut_y2_ = new TH2F ("pions_pT_vs_phi_cut_y2_", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_pT_vs_phi_cut_y3_ = new TH2F ("pions_pT_vs_phi_cut_y3_", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_pT_vs_phi_cut_y_pt1_ = new TH2F ("pions_pT_vs_phi_cut_y_pt1_", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_pT_vs_phi_cut_y_pt2_ = new TH2F ("pions_pT_vs_phi_cut_y_pt2_", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    pions_phi_vs_eta_ =  new TH2F ("pions_phi_vs_eta", " ; #varphi; #eta", 350, -3.5, 3.5, 300, 1.0, 7.0 );
    pions_rapidity_ = new TH1F ("pions_rapidity", "; y; entries", 250, -1.50, 3.50);
    pions_rapidity_cut_y1_ = new TH1F ("pions_rapidity_cut_y1_", "; y; entries", 75, -0.2, 1.3);
    pions_rapidity_cut_y2_ = new TH1F ("pions_rapidity_cut_y2_", "; y; entries", 75, -0.2, 1.3);
    pions_rapidity_cut_y3_ = new TH1F ("pions_rapidity_cut_y3_", "; y; entries", 75, -0.2, 1.3);
    pions_rapidity_vs_pT_ = new TH2F ("pions_pT_vs_y", "; y; p_{T} [GeV/c]",250, -1.50, 3.50, 150, 0., 3 );
    pions_rapidity_vs_pT_cut_y1_ = new TH2F ("pions_pT_vs_y_cut_y1_", "; y; p_{T} [GeV/c]",75, -0.2, 1.3, 150, 0., 3 );
    pions_rapidity_vs_pT_cut_y2_ = new TH2F ("pions_pT_vs_y_cut_y2_", "; y; p_{T} [GeV/c]",75, -0.2, 1.3, 150, 0., 3 );
    pions_rapidity_vs_pT_cut_y3_ = new TH2F ("pions_pT_vs_y_cut_y3_", "; y; p_{T} [GeV/c]",75, -0.2, 1.3, 150, 0., 3 );
    pions_rapidity_vs_pT_cut_y_pt1_ = new TH2F ("pions_pT_vs_y_cut_y_pt1_", "; y; p_{T} [GeV/c]",75, -0.2, 1.3, 150, 0., 3 );
    pions_rapidity_vs_pT_cut_y_pt2_ = new TH2F ("pions_pT_vs_y_cut_y_pt2_", "; y; p_{T} [GeV/c]",75, -0.2, 1.3, 150, 0., 3 );
    pions_rapidity_vs_phi_ = new TH2F ("pions_phi_vs_y_", " ; #varphi; y", 350, -3.5, 3.5, 75, -0.2, 1.3);
    pions_rapidity_vs_phi_cut_y1_ = new TH2F ("pions_phi_vs_y_cut_y1_", " ; #varphi; y", 350, -3.5, 3.5, 75, -0.2, 1.3);
    pions_rapidity_vs_phi_cut_y2_ = new TH2F ("pions_phi_vs_y_cut_y2_", " ; #varphi; y", 350, -3.5, 3.5, 75, -0.2, 1.3);
    pions_rapidity_vs_phi_cut_y3_ = new TH2F ("pions_phi_vs_y_cut_y3_", " ; #varphi; y", 350, -3.5, 3.5, 75, -0.2, 1.3);
    pions_rapidity_vs_phi_cut_y_pt1_ = new TH2F ("pions_phi_vs_y_cut_y_pt1_", " ; #varphi; y", 350, -3.5, 3.5, 75, -0.2, 1.3);
    pions_rapidity_vs_phi_cut_y_pt2_ = new TH2F ("pions_phi_vs_y_cut_y_pt2_", " ; #varphi; y", 350, -3.5, 3.5, 75, -0.2, 1.3);
    protons_dedx_vs_p_dist_ = new TH2F ("protons_dEdx_vs_p", "; q*log(20p/1 GeV/c); dE/dx ", 320, 0., 8.,150, 0.5, 3.5);
    protons_p_distribution_ = new TH1F( "protons_p_distribution", ";p [GeV/c]; entries", 600, 0., 30 );
    protons_pT_distribution_ = new TH1F( "protons_pT_distribution", ";p_{T} [GeV/c]; entries", 300, 0., 3 );
    protons_phi_distribution_ = new TH1F( "protons_phi_distribution", ";#varphi ; entries", 350, -3.5, 3.5 );
    protons_eta_distribution_ = new TH1F ("protons_eta_distribution", "; #eta; entries", 300, 1.0, 7.0);
    protons_pT_vs_eta_ = new TH2F ("protons_pT_vs_eta", "; #eta; p_{T} [GeV/c]",300, 1.0, 7.0, 150, 0., 3  );
    protons_pT_vs_phi_ = new TH2F ("protons_pT_vs_phi", "; #varphi; p_{T} [GeV/c]",350, -3.5, 3.5, 150, 0., 3   );
    protons_phi_vs_eta_ =  new TH2F ("protons_phi_vs_eta", " ; #varphi; #eta", 350, -3.5, 3.5, 300, 1.0, 7.0 );
    protons_rapidity_ = new TH1F ("protons_rapidity", "; y; entries", 250, 0, 5.0);
    protons_rapidity_vs_pT_ = new TH2F ("protons_pT_vs_y", "; y; p_{T} [GeV/c]",100, 0, 5.0, 150, 0., 3 );
    protons_rapidity_vs_phi = new TH2F ("protons_phi_vs_y", " ; #varphi; y", 350, -3.5, 3.5, 100, 0, 5.0);

    pions_mult_cut_y1 = new TH1F( "pions_mult_cut_y1", ";M ;entries", 60, 0, 60 );
    pions_mult_cut_y2 = new TH1F( "pions_mult_cut_y2", ";M ;entries", 60, 0, 60 );
    pions_mult_cut_y3 = new TH1F( "pions_mult_cut_y3", ";M ;entries", 60, 0, 60 );
    pions_mult_cut_y_pt1 = new TH1F( "pions_mult_cut_y_pt1", ";M ;entries", 60, 0, 60 );
    pions_mult_cut_y_pt2 = new TH1F( "pions_mult_cut_y_pt2", ";M ;entries", 60, 0, 60 );
    pions_mult_vs_psd_energy_cut_y1 = new TH2F ("pions_psd_energy_vs_mult_cut_y1", "; M;E_{PSD} [GeV]", 60, 0, 60, 720, 0., 3600.0);
    pions_mult_vs_psd_energy_cut_y2 = new TH2F ("pions_psd_energy_vs_mult_cut_y2", "; M;E_{PSD} [GeV]", 60, 0, 60, 720, 0., 3600.0);
    pions_mult_vs_psd_energy_cut_y3 = new TH2F ("pions_psd_energy_vs_mult_cut_y3", "; M;E_{PSD} [GeV]", 60, 0, 60, 720, 0., 3600.0);
    pions_mult_vs_psd_energy_cut_y_pt1 = new TH2F ("pions_psd_energy_vs_mult_cut_y_pt1", "; M;E_{PSD} [GeV]", 60, 0, 60, 720, 0., 3600.0);
    pions_mult_vs_psd_energy_cut_y_pt2 = new TH2F ("pions_psd_energy_vs_mult_cut_y_pt2", "; M;E_{PSD} [GeV]", 60, 0, 60, 720, 0., 3600.0);
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
        int n_pions_minus_cut_y1=0, n_pions_minus_cut_y2=0, n_pions_minus_cut_y3=0, n_pions_minus_cut_y_pt1=0, n_pions_minus_cut_y_pt2=0, n_pions_minus_cut_y_pt3=0;
        //auto psd_modules_posotions = data_header_->GetModulePositions(0);
        std::array<float, 3> psd_subevents_energy{0, 0, 0};
        for (auto psd_ch : psd_modules_->Loop()) {
            auto mod_no = psd_ch[psd_number].GetInt();
            auto mod_signal = psd_ch[psd_signal].GetVal();
            int subevent_id;
            if ((0 <= mod_no && mod_no < 16) || mod_no == 44)
                subevent_id = 0;
            else if (16 <= mod_no && mod_no < 29)
                subevent_id = 1;
            else if (29 <= mod_no && mod_no < 44)
                subevent_id = 2;
            else
                 assert(false);
            psd_subevents_energy[subevent_id] += mod_signal;
        }
        psd_energy_distribution_central_->Fill(psd_subevents_energy[0]);
        psd_energy_distribution_pherical_->Fill(psd_subevents_energy[2]);
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
        auto vtxZ_ = (*event_header_)[vtxZ_var_].GetVal();
        auto vtxX_ = (*event_header_)[vtxX_var_].GetVal();
        auto vtxY_ = (*event_header_)[vtxY_var_].GetVal();
        auto wfaT4_ = (*event_header_)[wfa_t4_].GetVal();
        auto wfaS1_ = (*event_header_)[wfa_s1_].GetVal();
        //histo before event cuts:
        multiplicity_distribution_->Fill(n_tracks_);
        psd_energy_distribution_ ->Fill(energy_psd);
        psd_energy_vs_mult_ ->Fill(n_tracks_, energy_psd);
        vtxZ_distribution_->Fill(vtxZ_);
        vtxZ_vs_mult_->Fill(vtxZ_,n_tracks_);
        wfa_t4_distribution_->Fill(wfaT4_);
        wfa_s1_distribution_->Fill(wfaS1_);

        //Events cuts:
        bool cut1_ = false, cut2_ = false, cut3_ = false, cut4_ = false, event_cuts_ = false;     //cut flag
        //cut1 vtxZ
        if ( abs(vtxZ_ + 591.8) <= 2.) cut1_ = true;
        //cut2 vtx_X & vtx_Y
        double x_min = -0.38 - 2*0.21, x_max = -0.38+ 2*0.21;
        double y_min = -0.29 - 2*0.15, y_max = -0.29 + 2*0.15;
        if ( (vtxX_ >= x_min && vtxX_ <= x_max) && (vtxY_ >= y_min && vtxY_ <= y_max)) cut2_ = true;
        //cut3 wfa_t4 & wfa_s1
        if (((*event_header_)[wfa_t4_].GetVal() >= 25000) && (wfaS1_ >= 4000)) cut3_ = true;
        //graphical cut
        if (energy_psd <= (-360. / 43) * n_tracks_ + 3600) cut4_ = true;

        if ( cut1_ ) {
            cut1_vtxZ_distribution_->Fill(vtxZ_);
            cut1_vtxZ_vs_mult_->Fill(vtxZ_, n_tracks_);
            cut1_mult_dist_->Fill(n_tracks_);
            cut1_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
            vtxX_distribution_->Fill(vtxX_);
            vtxY_distribution_->Fill(vtxY_);
            vtxX_vs_vtxY_->Fill(vtxX_,vtxY_);
        }

        if (cut1_ && cut2_) {
            cut2_vtxX_vs_vtxY_->Fill(vtxX_, vtxY_);
            cut2_vtxX_distribution_->Fill(vtxX_);
            cut2_vtxY_distribution_->Fill(vtxY_);
            cut2_mult_dist_->Fill(n_tracks_);
            cut2_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
        }
        if (cut1_ && cut2_ && cut3_) {
            cut3_mult_dist_->Fill(n_tracks_);
            cut3_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
            cut3_wfa_t4_distribution_->Fill(wfaT4_);
            cut3_wfa_s1_distribution_->Fill(wfaS1_);
        }
        if (cut1_ && cut2_ && cut3_ && cut4_){
            event_cuts_ = true;
            cut4_mult_dist_->Fill(n_tracks_);
            cut4_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
        }

        //For good events:
        if (event_cuts_) {
            //Track cuts
            for (auto &track: vtx_tracks_->Loop()) {
                //initialise variables
                bool cut5_ = false, cut56_ = false;
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
                nhits_total_vs_potential_->Fill(nhits[3], nhits_pot[3]);
                //cut 5: dcaX vs dcaY
                if (abs(dcaX) <= 2. && abs(dcaY) <= 1.) {
                    cut5_ = true;
                    n_tracks_cut5++;
                    cut5_dcaX_distribution_->Fill(dcaX);
                    cut5_dcaY_distribution_->Fill(dcaY);
                    cut5_dcaX_vs_dcaY_->Fill(dcaX, dcaY);
                }
                //cut 6: Nhits
                if (cut5_ && ((nhits[0] + nhits[1]) > 15) && (nhits[3] > 30) && (ratio > 0.55) && (ratio < 1.1)){
                    cut56_ = true;
                    n_tracks_cut6++;
                    cut6_nhits_total_distribution_->Fill(nhits[3]);
                    cut6_nhits_vtpc1_distribution_->Fill(nhits[0]);
                    cut6_nhits_vtpc2_distribution_->Fill(nhits[1]);
                    cut6_nhits_mtpc_distribution_->Fill(nhits[2]);
                    cut6_nhits_divided_nhits_pot_->Fill(ratio);
                    cut6_nhits_total_vs_potential_->Fill(nhits[3], nhits_pot[3]);
                }

                auto p = track.DataT<Track>()->GetP();
                auto pT = track.DataT<Track>()->GetPt();
                auto pz = track.DataT<Track>()->GetPz();
                auto phi = track.DataT<Track>()->GetPhi();
                auto eta = track.DataT<Track>()->GetEta();
                //all event and track cuts + pions
                if (cut56_) {
                    n_tracks_all_cut++;
                    all_cut_p_distribution_->Fill(p_);
                    all_cut_dedx_vs_p_dist_->Fill(q_ * log(20 * p_), dEdx_);
                    auto p_beam = 13;
                    double beta_beam = p_beam / sqrt(p_beam * p_beam + 0.9315 * 0.9315);
                    double y_beam = atanh(0.5 * beta_beam);
                    double log_p = q_ * log(20 * p_);
                    if (cut_pi->IsInside(log_p, dEdx_)) {
                        n_pions_minus++;
                        //TDatabasePDG* DatabasePDG=TDatabasePDG::Instance();
                        //double fMass = DatabasePDG->GetParticle(-211)->Mass();
                        auto y = track.DataT<Track>()->GetRapidity(-211);
                        double y_cm = y - y_beam;
                        //Double_t fRapidity = 0.5*log((sqrt(fMass*fMass+p*p)+pz)/(sqrt(fMass*fMass+p*p)-pz));
                        pions_dedx_vs_p_dist_->Fill(log_p, dEdx_);
                        pions_rapidity_->Fill(y_cm);
                        pions_rapidity_vs_phi_ ->Fill (phi, y_cm);
                        pions_rapidity_vs_pT_ ->Fill(y_cm, pT);
                        pions_pT_distribution_->Fill(pT);
                        pions_p_distribution_->Fill(p);
                        pions_phi_distribution_->Fill(phi);
                        pions_eta_distribution_->Fill(eta);
                        pions_pT_vs_eta_->Fill(eta, pT);
                        pions_pT_vs_phi_->Fill(phi, pT);
                        pions_phi_vs_eta_->Fill(phi, eta);

                        if ( y_cm >= 0 ) {
                            if ( y_cm <= 0.8){
                                n_pions_minus_cut_y1++;
                                pions_rapidity_cut_y1_->Fill(y_cm);
                                pions_rapidity_vs_pT_cut_y1_ ->Fill(y_cm, pT);
                                pions_rapidity_vs_phi_cut_y1_ ->Fill (phi, y_cm);
                                pions_pT_vs_phi_cut_y1_->Fill(phi, pT);
                            }
                            if ( y_cm <= 1.0){
                                n_pions_minus_cut_y2++;
                                pions_rapidity_cut_y2_->Fill(y_cm);
                                pions_rapidity_vs_pT_cut_y2_ ->Fill(y_cm, pT);
                                pions_rapidity_vs_phi_cut_y2_ ->Fill (phi, y_cm);
                                pions_pT_vs_phi_cut_y2_->Fill(phi, pT);
                                if ( pT>0.2){
                                    n_pions_minus_cut_y_pt1++;
                                    pions_pT_distribution_cut_y_pt1_ ->Fill (pT);
                                    pions_rapidity_vs_pT_cut_y_pt1_ ->Fill(y_cm, pT);
                                    pions_rapidity_vs_phi_cut_y_pt1_ ->Fill (phi, y_cm);
                                    pions_pT_vs_phi_cut_y_pt1_->Fill(phi, pT);
                                }
                                if ( pT>0.3){
                                    n_pions_minus_cut_y_pt2++;
                                    pions_pT_distribution_cut_y_pt2_ ->Fill (pT);
                                    pions_rapidity_vs_pT_cut_y_pt2_ ->Fill(y_cm, pT);
                                    pions_rapidity_vs_phi_cut_y_pt2_ ->Fill (phi, y_cm);
                                    pions_pT_vs_phi_cut_y_pt2_->Fill(phi, pT);
                                }
                            }
                            if ( y_cm <= 1.2){
                                n_pions_minus_cut_y3++;
                                pions_rapidity_cut_y3_->Fill(y_cm);
                                pions_rapidity_vs_pT_cut_y3_ ->Fill(y_cm, pT);
                                pions_rapidity_vs_phi_cut_y3_ ->Fill (phi, y_cm);
                                pions_pT_vs_phi_cut_y3_->Fill(phi, pT);
                            }
                        }
                    }
                    if (cut_p->IsInside(log_p, dEdx_)) {
                        n_protons++;
                        auto y = track.DataT<Track>()->GetRapidity(2212);
                        double y_cm = y - y_beam;
                        protons_dedx_vs_p_dist_->Fill(log_p, dEdx_);
                        protons_eta_distribution_->Fill(eta);
                        protons_rapidity_->Fill(y_cm);
                        protons_rapidity_vs_pT_ ->Fill(y_cm, pT);
                        protons_rapidity_vs_phi ->Fill(phi, y_cm);
                        protons_pT_distribution_->Fill(pT);
                        protons_p_distribution_->Fill(p);
                        protons_phi_distribution_->Fill(phi);
                        protons_pT_vs_eta_->Fill(eta, pT);
                        protons_pT_vs_phi_->Fill(phi, pT);
                        protons_phi_vs_eta_->Fill(phi, eta);
                    }
                }

            }

            good_events_mult_dist_scaled_->Fill(1.0*n_tracks_/300);
            good_events_psd_energy_vs_mult_->Fill(n_tracks_, energy_psd);
            pions_multiplicity_dist_->Fill(n_pions_minus);
            pions_multiplicity_dist_scaled_->Fill(1.0*n_pions_minus/50);
            pions_mult_vs_psd_energy_ ->Fill (n_pions_minus, energy_psd);
            pions_mult_vs_psd_energy_central_ ->Fill (n_pions_minus, psd_subevents_energy[0]);
            M_good = n_tracks_;
            M_cut_dca = n_tracks_cut5;
            M_ch = n_tracks_all_cut;
            M_pi = n_pions_minus;
            M_p = n_protons;
            tree->Fill();
            pions_mult_vs_psd_energy_pherical_ ->Fill(n_pions_minus, psd_subevents_energy[2]);
            protons_multiplicity_dist_->Fill (n_protons);
            protons_multiplicity_dist_scaled_->Fill(1.0*n_protons/40);
            protons_mult_vs_psd_energy_ ->Fill(n_protons, energy_psd);
            protons_mult_vs_psd_energy_central_ ->Fill(n_protons, psd_subevents_energy[0]);
            protons_mult_vs_psd_energy_pherical_ ->Fill(n_protons, psd_subevents_energy[2]);
            //fill multiplicity histo for track cut
            cut5_mult_dist_->Fill(n_tracks_cut5);
            cut5_mult_dist_scaled_->Fill(1.0*n_tracks_cut5/230);
            cut5_psd_energy_vs_mult_->Fill(n_tracks_cut5, energy_psd);
            cut6_mult_dist_->Fill(n_tracks_cut6);
            cut6_mult_dist_scaled_->Fill(1.0*n_tracks_cut6/160);
            cut6_psd_energy_vs_mult_->Fill(n_tracks_cut6, energy_psd);
            all_cut_mult_dist_->Fill(n_tracks_all_cut);
            all_cut_mult_dist_scaled_->Fill(1.0*n_tracks_all_cut/160);
            all_cut_psd_energy_vs_mult_->Fill(n_tracks_all_cut, energy_psd);
            all_cut_psd_energy_central_vs_mult_->Fill(n_tracks_all_cut, psd_subevents_energy[0]);
            all_cut_psd_energy_pherical_vs_mult_->Fill(n_tracks_all_cut, psd_subevents_energy[2]);

            pions_mult_cut_y1 ->Fill(n_pions_minus_cut_y1);
            pions_mult_cut_y2 ->Fill(n_pions_minus_cut_y2);
            pions_mult_cut_y3 ->Fill(n_pions_minus_cut_y3);
            pions_mult_cut_y_pt1 ->Fill(n_pions_minus_cut_y_pt1);
            pions_mult_cut_y_pt2 ->Fill(n_pions_minus_cut_y_pt2);
            pions_mult_vs_psd_energy_cut_y1->Fill(n_pions_minus_cut_y1, energy_psd);
            pions_mult_vs_psd_energy_cut_y2->Fill(n_pions_minus_cut_y2, energy_psd);
            pions_mult_vs_psd_energy_cut_y3->Fill(n_pions_minus_cut_y3, energy_psd);
            pions_mult_vs_psd_energy_cut_y_pt1->Fill(n_pions_minus_cut_y_pt1, energy_psd);
            pions_mult_vs_psd_energy_cut_y_pt2->Fill(n_pions_minus_cut_y_pt2, energy_psd);
        }
    }
}

void AnalysisTask::UserFinish() {
    // Writing histograms to file
    out_file_->cd();

    tree->Write();

    TCutG *cut_pi = new TCutG("cut_pi",6);
    cut_pi->SetPoint(0, -5.26214, 1.41964);
    cut_pi->SetPoint(1, -5.26214, 1.19643 );
    cut_pi->SetPoint(2, -2.25243, 0.883929);
    cut_pi->SetPoint(3, -2.25243, 1.10714 );
    cut_pi->SetPoint(4, -5.24272, 1.41964);
    cut_pi->SetPoint(5, -5.26214, 1.41964);

    TCutG *cut_p = new TCutG("cut_protons",6);
    cut_p->SetPoint(0, 4.10163, 1.03274);
    cut_p->SetPoint(1, 5.80762, 1.20387 );
    cut_p->SetPoint(2, 5.80762, 1.0625);
    cut_p->SetPoint(3, 4.11978, 0.89881 );
    cut_p->SetPoint(4, 4.11978, 1.0253);
    cut_p->SetPoint(5, 4.10163, 1.03274);

    cut_pi->Write();
    cut_p->Write();

    multiplicity_distribution_->Write();
    psd_energy_distribution_->Write();
    psd_energy_distribution_central_->Write();
    psd_energy_distribution_pherical_->Write();
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
    good_events_mult_dist_scaled_->Write();
    good_events_psd_energy_vs_mult_->Write();

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
    nhits_total_vs_potential_ ->Write();
    cut6_nhits_total_distribution_ ->Write();
    cut6_nhits_vtpc1_distribution_ ->Write();
    cut6_nhits_vtpc2_distribution_ ->Write();
    cut6_nhits_mtpc_distribution_ ->Write();
    cut6_nhits_divided_nhits_pot_ ->Write();
    cut6_nhits_total_vs_potential_ ->Write();

    cut6_mult_dist_ ->Write();
    cut6_mult_dist_scaled_->Write();
    cut6_psd_energy_vs_mult_ ->Write();
    all_cut_mult_dist_->Write();
    all_cut_mult_dist_scaled_->Write();
    all_cut_psd_energy_vs_mult_->Write();
    all_cut_psd_energy_central_vs_mult_->Write();
    all_cut_psd_energy_pherical_vs_mult_->Write();
    all_cut_dedx_vs_p_dist_->Write();
    q_distribution_->Write();
    p_distribution_->Write();
    dedx_distribution_->Write();
    dedx_vs_p_dist_->Write();
    all_cut_p_distribution_->Write();

    TDirectory *dir_2 = out_file_->mkdir("pi_minus");
    dir_2->cd();
    pions_multiplicity_dist_->Write();
    pions_multiplicity_dist_scaled_->Write();
    pions_mult_vs_psd_energy_ ->Write();
    pions_mult_vs_psd_energy_central_ ->Write();
    pions_mult_vs_psd_energy_pherical_->Write();
    pions_p_distribution_->Write();
    pions_pT_distribution_->Write();
    pions_eta_distribution_->Write();
    pions_phi_distribution_->Write();
    pions_rapidity_->Write();
    pions_rapidity_cut_y2_->Write();
    pions_rapidity_vs_phi_ ->Write();
    pions_rapidity_vs_phi_cut_y_pt1_ ->Write();
    pions_rapidity_vs_pT_ ->Write();
    pions_pT_vs_eta_->Write();
    pions_pT_vs_phi_->Write();
    pions_phi_vs_eta_->Write();
    pions_dedx_vs_p_dist_->Write();

    TDirectory *dir_1 = out_file_->mkdir("protons");
    dir_1->cd();
    protons_multiplicity_dist_->Write();
    protons_multiplicity_dist_scaled_->Write();
    protons_mult_vs_psd_energy_ ->Write();
    protons_mult_vs_psd_energy_central_ ->Write();
    protons_mult_vs_psd_energy_pherical_->Write();
    protons_p_distribution_->Write();
    protons_pT_distribution_->Write();
    protons_eta_distribution_->Write();
    protons_phi_distribution_->Write();
    protons_rapidity_->Write();
    protons_rapidity_vs_phi ->Write();
    protons_rapidity_vs_pT_ ->Write();
    protons_pT_vs_eta_->Write();
    protons_pT_vs_phi_->Write();
    protons_phi_vs_eta_->Write();
    protons_dedx_vs_p_dist_->Write();

    TDirectory *dir1 = out_file_->mkdir("pions_cut_y1");
    dir1->cd();
    pions_rapidity_cut_y1_->Write();
    pions_rapidity_vs_phi_cut_y1_ ->Write();
    pions_pT_vs_phi_cut_y1_->Write();
    pions_rapidity_vs_pT_cut_y1_->Write();
    pions_mult_cut_y1->Write();
    pions_mult_vs_psd_energy_cut_y1->Write();
    TDirectory *dir2 = out_file_->mkdir("pions_cut_y2");
    dir2->cd();
    pions_rapidity_cut_y2_->Write();
    pions_rapidity_vs_phi_cut_y2_ ->Write();
    pions_pT_vs_phi_cut_y2_->Write();
    pions_rapidity_vs_pT_cut_y2_->Write();
    pions_mult_cut_y2->Write();
    pions_mult_vs_psd_energy_cut_y2->Write();
    TDirectory *dir3 = out_file_->mkdir("pions_cut_y3");
    dir3->cd();
    pions_rapidity_cut_y3_->Write();
    pions_rapidity_vs_phi_cut_y3_ ->Write();
    pions_pT_vs_phi_cut_y3_->Write();
    pions_rapidity_vs_pT_cut_y3_->Write();
    pions_mult_cut_y3->Write();
    pions_mult_vs_psd_energy_cut_y3->Write();
    TDirectory *dir4 = out_file_->mkdir("pions_cut_y_pt1");
    dir4->cd();
    pions_rapidity_vs_phi_cut_y_pt1_ ->Write();
    pions_pT_vs_phi_cut_y_pt1_->Write();
    pions_rapidity_vs_pT_cut_y_pt1_->Write();
    pions_mult_cut_y_pt1->Write();
    pions_mult_vs_psd_energy_cut_y_pt1->Write();
    TDirectory *dir5 = out_file_->mkdir("pions_cut_y_pt2");
    dir5->cd();
    pions_rapidity_vs_phi_cut_y_pt2_ ->Write();
    pions_pT_vs_phi_cut_y_pt2_->Write();
    pions_rapidity_vs_pT_cut_y_pt2_->Write();
    pions_mult_cut_y_pt2->Write();
    pions_mult_vs_psd_energy_cut_y_pt2->Write();
}
