// Created by Aleksandra 11/02/2021

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <TFile.h>
#include <TTree.h>


#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/DataHeader.hpp>
#include "TStyle.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <at_task/Task.h>
#include <memory>
#include <string>
#include <TCutG.h>

class AnalysisTask : public UserFillTask{
public:
    AnalysisTask() = default;
    ~AnalysisTask() override = default;
    void UserInit( std::map<std::string, void*>& branch_map ) override;
    void UserExec() override;
    void UserFinish() override;
    boost::program_options::options_description GetBoostOptions() override;
    void PreInit() override {};
    void PostFinish() override {
        UserTask::PostFinish();
    }
private:
    //TFile *file_M{TFile::Open("out_tree", "create")};
    TTree *tree{new TTree("tree", "tree")};
    Int_t M_good, M_cut_dca, M_cut_tpc, M_ch, M_pi, M_p;
    //ATI2::Branch* branch_Mch_cut_dca, *branch_Mch_all_cuts, *branch_Mpi, *branch_Mp;

    /* pointers to link tree's branches with */
    ATI2::Branch* event_header_{nullptr}; 		// event info
    ATI2::Branch* vtx_tracks_{nullptr}; 		        // reconstructed tracks
    ATI2::Branch* psd_modules_{nullptr}; 		// modules of FhCal branch
    AnalysisTree::TrackDetector* VtxTracks_{nullptr};

    ATI2::Variable psd_signal;
    ATI2::Variable psd_number;
    ATI2::Variable dcax, dcay, dcaz, p, dedx, q;
    ATI2::Variable nhits_total, nhits_vtpc1, nhits_vtpc2, nhits_mtpc;
    ATI2::Variable nhits_pot_total, nhits_pot_vtpc1, nhits_pot_vtpc2, nhits_pot_mtpc;

    ATI2::Variable evt_t1;
    ATI2::Variable evt_t2;
    ATI2::Variable evt_t4;
    ATI2::Variable multiplicity_var_;
    ATI2::Variable vtxX_var_;
    ATI2::Variable vtxY_var_;
    ATI2::Variable vtxZ_var_;
    ATI2::Variable energy_psd_var_;
    ATI2::Variable pT_var_;
    ATI2::Variable wfa_t4_;
    ATI2::Variable wfa_s1_;

    //TCutG *cutg;

    //before cut
    TH1F* multiplicity_distribution_, * psd_energy_distribution_, *psd_energy_distribution_central_, *psd_energy_distribution_pherical_;
    TH2F* psd_energy_vs_mult_;
    TH1F* pT_distribution_, * eta_distribution_, * phi_distribution_;
    TH2F* pT_vs_eta_, * pT_vs_phi_, * phi_vs_eta_;

    // cut1-2: vtxZ && vtxX_Y
    TH1F* vtxX_distribution_, * vtxY_distribution_, * vtxZ_distribution_;
    TH2F* vtxZ_vs_mult_, * vtxX_vs_vtxY_;
    TH1F* cut2_vtxX_distribution_, * cut2_vtxY_distribution_, * cut1_vtxZ_distribution_;
    TH2F* cut1_vtxZ_vs_mult_, * cut2_vtxX_vs_vtxY_;
    TH1F* cut1_mult_dist_, * cut2_mult_dist_, * cut12_mult_dist_;
    TH2F* cut1_psd_energy_vs_mult_, * cut2_psd_energy_vs_mult_, * cut12_psd_energy_vs_mult_;

    // cut3: wfa_t4 & wfa_s1
    TH1F* wfa_t4_distribution_,  * wfa_s1_distribution_;
    TH1F* cut3_wfa_t4_distribution_, * cut3_wfa_s1_distribution_;
    TH1F* cut3_mult_dist_;
    TH2F* cut3_psd_energy_vs_mult_;

    // cut4: graphical cut
    TH1F* cut4_mult_dist_;
    TH2F* cut4_psd_energy_vs_mult_;

    // histo for good event without track cuts
    TH1F* good_events_mult_dist_scaled_;
    TH2F* good_events_psd_energy_vs_mult_, *all_cut_psd_energy_central_vs_mult_, *all_cut_psd_energy_pherical_vs_mult_ ;

    // cut5: dcaX && dcaY
    TH1F* dcaX_distribution_, * dcaY_distribution_;
    TH2F* dcaX_vs_dcaY_;
    TH1F* cut5_dcaX_distribution_,* cut5_dcaY_distribution_;
    TH2F* cut5_dcaX_vs_dcaY_;
    TH1F* cut5_mult_dist_, *cut5_mult_dist_scaled_;
    TH2F* cut5_psd_energy_vs_mult_;

    // cut6: Nhits
    TH1F* nhits_total_distribution_, * nhits_vtpc1_distribution_, * nhits_vtpc2_distribution_, * nhits_mtpc_distribution_;
    //TH1F* nhits_pot_total_distribution_, * nhits_pot_vtpc1_distribution_, * nhits_pot_vtpc2_distribution_, * nhits_pot_mtpc_distribution_;
    TH1F* cut6_nhits_total_distribution_, * cut6_nhits_vtpc1_distribution_, * cut6_nhits_vtpc2_distribution_, * cut6_nhits_mtpc_distribution_;
    //TH1F* cut6_nhits_pot_total_distribution_, * cut6_nhits_pot_vtpc1_distribution_, * cut6_nhits_pot_vtpc2_distribution_, * cut6_nhits_pot_mtpc_distribution_;
    TH1F* nhits_divided_nhits_pot_;
    TH2F* nhits_total_vs_potential_, *cut6_nhits_total_vs_potential_;
    TH2F* nhits_vtpc1_vs_vtpc2, * nhits_vtpc_vs_mtpc, * nhits_mtpc_vs_total;
    TH2F* cut6_nhits_vtpc1_vs_vtpc2, * cut6_nhits_vtpc_vs_mtpc, * cut6_nhits_mtpc_vs_total;
    TH1F* cut6_nhits_divided_nhits_pot_;
    TH1F* cut6_mult_dist_, *cut6_mult_dist_scaled_ ;
    TH2F* cut6_psd_energy_vs_mult_;

    // dEdx_vs_q
    TH1F* dedx_distribution_, * p_distribution_, * q_distribution_;
    TH2F* dedx_vs_p_dist_, *pions_dedx_vs_p_dist_, *protons_dedx_vs_p_dist_, *protons_all_dedx_vs_p_dist_;

    // all cuts
    TH1F* all_cut_mult_dist_, *all_cut_mult_dist_scaled_;
    TH1F* pions_multiplicity_dist_, *pions_multiplicity_dist_scaled_;
    TH1F* protons_multiplicity_dist_, *protons_multiplicity_dist_scaled_;
    TH1F* protons_all_multiplicity_dist_, *protons_all_multiplicity_dist_scaled_;
    TH2F* all_cut_psd_energy_vs_mult_;
    TH1F* all_cut_p_distribution_;
    TH2F* all_cut_dedx_vs_p_dist_;

    // kinematic for pions
    TH1F* pions_pT_distribution_, *pions_p_distribution_, *pions_eta_distribution_, *pions_phi_distribution_, *pions_rapidity_;
    TH2F* pions_mult_vs_psd_energy_, *pions_mult_vs_psd_energy_central_, *pions_mult_vs_psd_energy_pherical_;
    TH2F* pions_pT_vs_eta_, * pions_pT_vs_phi_, * pions_phi_vs_eta_, * pions_rapidity_vs_pT_, * pions_rapidity_vs_phi_;
    // pions with cut (y,pT)
    TH1F* pions_rapidity_cut_y1_, *pions_rapidity_cut_y2_, *pions_rapidity_cut_y3_, *pions_pT_cut_y_pt1_, *pions_pT_cut_y_pt2_;
    TH1F* pions_pT_distribution_cut_y_pt1_, *pions_pT_distribution_cut_y_pt2_;
    TH2F* pions_pT_vs_phi_cut_y1_, *pions_rapidity_vs_pT_cut_y1_, *pions_rapidity_vs_phi_cut_y1_;
    TH2F* pions_pT_vs_phi_cut_y2_, *pions_rapidity_vs_pT_cut_y2_, *pions_rapidity_vs_phi_cut_y2_;
    TH2F* pions_pT_vs_phi_cut_y_pt1_, *pions_rapidity_vs_pT_cut_y_pt1_, *pions_rapidity_vs_phi_cut_y_pt1_;
    TH2F* pions_pT_vs_phi_cut_y_pt2_, *pions_rapidity_vs_pT_cut_y_pt2_, *pions_rapidity_vs_phi_cut_y_pt2_;
    TH2F* pions_pT_vs_phi_cut_y3_, *pions_rapidity_vs_pT_cut_y3_, *pions_rapidity_vs_phi_cut_y3_;
    TH1F * pions_mult_cut_y1, *pions_mult_cut_y2, *pions_mult_cut_y3, *pions_mult_cut_y_pt1, *pions_mult_cut_y_pt2;
    TH2F *pions_mult_vs_psd_energy_cut_y1, *pions_mult_vs_psd_energy_cut_y2, *pions_mult_vs_psd_energy_cut_y3, *pions_mult_vs_psd_energy_cut_y_pt1, *pions_mult_vs_psd_energy_cut_y_pt2;

    TH1F* protons_pT_distribution_, *protons_p_distribution_, *protons_eta_distribution_, *protons_phi_distribution_, * protons_rapidity_;
    TH2F* protons_mult_vs_psd_energy_, *protons_mult_vs_psd_energy_central_, *protons_mult_vs_psd_energy_pherical_;
    TH2F* protons_pT_vs_eta_, *protons_pT_vs_phi_, *protons_phi_vs_eta_, *protons_rapidity_vs_pT_, *protons_rapidity_vs_phi;
    TASK_DEF(AnalysisTask, 0)
};
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
