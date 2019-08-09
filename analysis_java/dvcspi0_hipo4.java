import java.io.*;
import java.util.*;
import org.jlab.groot.math.*;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;
import java.util.Random;

public class dvcspi0_hipo4 {
	public boolean isMC, isDVCS, isPI0, usePCAL, useFTCAL;
	public Random myMC;
	public double EB, HWP;
	public int HELI, NPHI, haz_elec, haz_prot, haz_g1, haz_g2, runNum, Ndvcs, Npi0;
	public int e_sect, p_sect;
	public double elec_mom, elec_theta, elec_phi, elec_vz, elec_phi_sect;
	public double prot_mom, prot_theta, prot_phi, prot_vz, prot_beta, prot_phi_sect;
	public double g1_E, g1_theta, g1_phi, g1_beta;
	public double g2_E, g2_theta, g2_phi, g2_beta;
	public double elast_W, elast_EB, elast_bb;

	public LorentzVector VB, VT, VE, VGS, VPROT, Vmand, VG1, VG2, VPI0, VMISS, VmissP, VmissG;
	public Vector3 Vlept, Vhadr, Vhad2;

	public H2F H_elec_theta_mom, H_elec_theta_phi, H_elec_phi_mom;
	public H2F H_prot_theta_mom, H_prot_theta_phi, H_prot_phi_mom;
	public H2F H_g_theta_E, H_g_theta_phi, H_g_phi_E;
	public H2F H_g1_theta_E, H_g1_theta_phi, H_g1_phi_E;
	public H2F H_g2_theta_E, H_g2_theta_phi, H_g2_phi_E;
	public H2F H_prot_beta_p, H_prot_e_vz, H_phot_b_p;

	public H2F H_prot_miss_phi, H_phot_elec_angle, H_phot_virt_angle;
	public H2F H_prot_pi0_miss_phi, H_pi0_elec_angle, H_pi0_virt_angle;

	public H2F H_dvcs_Q2xB, H_dvcs_tphi;
	public H1F H_dvcs_pT, H_dvcs_cone, H_dvcs_copl, H_dvcs_ME, H_dvcs_MM_ep, H_dvcs_MM_eg, H_dvcs_MM_epg, H_dvcs_Phi;
	public H2F H_pi0_Q2xB, H_pi0_tphi, H_gg_open_E;
	public H1F H_pi0_pT, H_pi0_cone, H_pi0_copl, H_pi0_ME, H_pi0_MM_ep, H_pi0_MM_epi0, H_pi0_MM_eppi0, H_pi0_Phi, H_IM_pi0;
	GraphErrors g_dvcs_bsa;

	public H2F[] H_elast_e_th_mom, H_elast_e_th_phi, H_elast_e_phi_mom;
	public H2F[] H_elast_p_th_mom, H_elast_p_th_phi, H_elast_p_phi_mom;
        public H2F[] H_elast_vz_vz, H_elast_W_th, H_elast_W_phi, H_elast_th_th, H_elast_ph_ph, H_elast_EB_bb;

	public H1F H_dvcs_phi, H_dvcs_phi_plus, H_dvcs_phi_minus;
	public H1F H_pi0_phi, H_pi0_phi_plus, H_pi0_phi_minus;

	public dvcspi0_hipo4(int reqrunNum, float EBreq, int reqNphi, int reqCal) {
		isMC = false;usePCAL=true;useFTCAL=true;
		myMC = new Random();
		NPHI = 72;//12;
		if(reqNphi>0)NPHI=reqNphi;
		runNum = reqrunNum;
		EB = 6.42f;
		if(EBreq>8f)EB = 10.6;
		if(reqCal>0){
			if(reqCal==1)usePCAL=false;
			if(reqCal==2)useFTCAL=false;
		}
		HWP = 1;
		Ndvcs  = 0;
		Npi0 = 0;
		VB = new LorentzVector(0,0,EB,EB);
		VT = new LorentzVector(0,0,0,0.938);
		H_elec_theta_mom = new H2F("H_elec_theta_mom","H_elec_theta_mom",100,0,EB,100,0,40);
		H_elec_theta_mom.setTitle("elec #theta vs p");
		H_elec_theta_mom.setTitleX("p (GeV)");
		H_elec_theta_mom.setTitleY("#theta (^o)");
		H_elec_theta_phi = new H2F("H_elec_theta_phi","H_elec_theta_phi",100,-180,180,100,0,40);
		H_elec_theta_phi.setTitle("elec #theta vs #phi");
		H_elec_theta_phi.setTitleX("#phi (^o)");
		H_elec_theta_phi.setTitleY("#theta (^o)");
		H_elec_phi_mom = new H2F("H_elec_phi_mom","H_elec_phi_mom",100,0,EB,100,-180,180);
		H_elec_phi_mom.setTitle("elec #phi vs p");
		H_elec_phi_mom.setTitleX("p (GeV)");
		H_elec_phi_mom.setTitleY("#phi (^o)");
		H_prot_theta_mom = new H2F("H_prot_theta_mom","H_prot_theta_mom",100,0,11,100,0,90);
		H_prot_theta_mom.setTitle("prot #theta vs p");
		H_prot_theta_mom.setTitleX("p (GeV)");
		H_prot_theta_mom.setTitleY("#theta (^o)");
		H_prot_theta_phi = new H2F("H_prot_theta_phi","H_prot_theta_phi",100,-180,180,100,0,90);
		H_prot_theta_phi.setTitle("prot #theta vs #phi");
		H_prot_theta_phi.setTitleX("#phi (^o)");
		H_prot_theta_phi.setTitleY("#theta (^o)");
		H_prot_phi_mom = new H2F("H_prot_phi_mom","H_prot_phi_mom",100,0,EB,100,-180,180);
		H_prot_phi_mom.setTitle("prot #phi vs p");
		H_prot_phi_mom.setTitleX("p (GeV)");
		H_prot_phi_mom.setTitleY("#phi (^o)");
		H_g_theta_E = new H2F("H_g_theta_E","H_g_theta_E",100,0,EB,100,0,40);
		H_g_theta_E.setTitle("#gamma #theta vs p");
		H_g_theta_E.setTitleX("p (GeV)");
		H_g_theta_E.setTitleY("#theta (^o)");
		H_g_theta_phi = new H2F("H_g_theta_phi","H_g_theta_phi",100,-180,180,100,0,40);
		H_g_theta_phi.setTitle("#gamma #theta vs #phi");
		H_g_theta_phi.setTitleX("#phi (^o)");
		H_g_theta_phi.setTitleY("#theta (^o)");
		H_g_phi_E = new H2F("H_g_phi_E","H_g_phi_E",100,0,EB,100,-180,180);
		H_g_phi_E.setTitle("#gamma #phi vs p");
		H_g_phi_E.setTitleX("p (GeV)");
		H_g_phi_E.setTitleY("#phi (^o)");
		H_g1_theta_E = new H2F("H_g1_theta_E","H_g1_theta_E",100,0,EB,100,0,40);
		H_g1_theta_E.setTitle("#gamma1 #theta vs p");
		H_g1_theta_E.setTitleX("p (GeV)");
		H_g1_theta_E.setTitleY("#theta (^o)");
		H_g1_theta_phi = new H2F("H_g1_theta_phi","H_g1_theta_phi",100,-180,180,100,0,40);
		H_g1_theta_phi.setTitle("#gamma1 #theta vs #phi");
		H_g1_theta_phi.setTitleX("#phi (^o)");
		H_g1_theta_phi.setTitleY("#theta (^o)");
		H_g1_phi_E = new H2F("H_g1_phi_E","H_g1_phi_E",100,0,EB,100,-180,180);
		H_g1_phi_E.setTitle("#gamma1 #phi vs p");
		H_g1_phi_E.setTitleX("p (GeV)");
		H_g1_phi_E.setTitleY("#phi (^o)");
		H_g2_theta_E = new H2F("H_g2_theta_E","H_g2_theta_E",100,0,EB,100,0,40);
		H_g2_theta_E.setTitle("#gamma2 #theta vs p");
		H_g2_theta_E.setTitleX("p (GeV)");
		H_g2_theta_E.setTitleY("#theta (^o)");
		H_g2_theta_phi = new H2F("H_g2_theta_phi","H_g2_theta_phi",100,-180,180,100,0,40);
		H_g2_theta_phi.setTitle("#gamma2 #theta vs #phi");
		H_g2_theta_phi.setTitleX("#phi (^o)");
		H_g2_theta_phi.setTitleY("#theta (^o)");
		H_g2_phi_E = new H2F("H_g2_phi_E","H_g2_phi_E",100,0,EB,100,-180,180);
		H_g2_phi_E.setTitle("#gamma2 #phi vs p");
		H_g2_phi_E.setTitleX("p (GeV)");
		H_g2_phi_E.setTitleY("#phi (^o)");
		H_prot_beta_p = new H2F("H_prot_beta_p","H_prot_beta_p",100,0,4,100,0,1.2);
		H_prot_beta_p.setTitle("prot #beta vs p");
		H_prot_beta_p.setTitleX("p (GeV)");
		H_prot_beta_p.setTitleY("#beta");
		H_phot_b_p = new H2F("H_phot_b_p","H_phot_b_p",100,0.5,EB,100,0.75,1.25);
		H_phot_b_p.setTitle("#gamma #beta vs p");
		H_phot_b_p.setTitleX("p (GeV)");
		H_phot_b_p.setTitleY("#beta");
		H_prot_e_vz = new H2F("H_prot_e_vz","H_prot_e_vz",100,-25,25,100,-25,25);
		H_prot_e_vz.setTitle("prot vz vs elec vz");
		H_prot_e_vz.setTitleX("elec vz (cm)");
		H_prot_e_vz.setTitleY("prot vz (cm)");

		H_prot_miss_phi = new H2F("H_prot_miss_phi","H_prot_miss_phi",100,-180,180,100,-180,180);
		H_prot_miss_phi.setTitle("DVCS proton #phi vs miss #phi");
		H_prot_miss_phi.setTitleX("proton #phi (^o)");
		H_prot_miss_phi.setTitleY("miss #phi (^o)");
		H_phot_elec_angle = new H2F("H_phot_elec_angle","H_phot_elec_angle",100,-180,180,100,0,40);
	       	H_phot_elec_angle.setTitle("angle #gamma e vs #phi");
		H_phot_elec_angle.setTitleX("#phi Trento (^o)");
		H_phot_elec_angle.setTitleY("angle #gamma e (^o)");
		H_phot_virt_angle = new H2F("H_phot_virt_angle","H_phot_virt_angle",100,0,5,100,0,35);
		H_phot_virt_angle.setTitle("angle #gamma#gamma^* vs -t");
		H_phot_virt_angle.setTitleX("-t (GeV^2)");
		H_phot_virt_angle.setTitleY("angle #gamma #gamma^* (^o)");

		H_prot_pi0_miss_phi = new H2F("H_prot_pi0_miss_phi","H_prot_pi0_miss_phi",100,-180,180,100,-180,180);
		H_prot_pi0_miss_phi.setTitle("#pi^0 proton #phi vs miss #phi");
		H_prot_pi0_miss_phi.setTitleX("proton #phi (^o)");
		H_prot_pi0_miss_phi.setTitleY("miss #phi (^o)");
		H_pi0_elec_angle = new H2F("H_pi0_elec_angle","H_pi0_elec_angle",100,-180,180,100,0,40);
	       	H_pi0_elec_angle.setTitle("angle #pi^0 e vs #phi");
		H_pi0_elec_angle.setTitleX("#phi Trento (^o)");
		H_pi0_elec_angle.setTitleY("angle #pi^0 e (^o)");
		H_pi0_virt_angle = new H2F("H_pi0_virt_angle","H_pi0_virt_angle",100,0,5,100,0,35);
		H_pi0_virt_angle.setTitle("angle #pi^0#gamma^* vs -t");
		H_pi0_virt_angle.setTitleX("-t (GeV^2)");
		H_pi0_virt_angle.setTitleY("angle #pi^0 #gamma^* (^o)");

		H_dvcs_Q2xB = new H2F("H_dvcs_Q2xB","H_dvcs_Q2xB",100,0,1,100,0,12);
		H_dvcs_Q2xB.setTitle("DVCS Q2 vs xB");
		H_dvcs_Q2xB.setTitleX("xB");
		H_dvcs_Q2xB.setTitleY("Q2");
		H_dvcs_tphi = new H2F("H_dvcs_tphi","H_dvcs_tphi",100,-180,180,100,0,5);
		H_dvcs_tphi.setTitle("DVCS -t vs #phi");
		H_dvcs_tphi.setTitleX("#phi (^o)");
		H_dvcs_tphi.setTitleY("-t (GeV^2)");
		H_dvcs_pT = new H1F("H_dvcs_pT","H_dvcs_pT",100,0,1.0);
		H_dvcs_pT.setTitle("DVCS missing pT");
		H_dvcs_pT.setTitleX("pT (GeV)");
		H_dvcs_cone = new H1F("H_dvcs_cone","H_dvcs_cone",100,0,5);
		H_dvcs_cone.setTitle("DVCS cone angle");
		H_dvcs_cone.setTitleX("cone ang (^o)");
		H_dvcs_copl = new H1F("H_dvcs_copl","H_dvcs_copl",100,0,30);
		H_dvcs_copl.setTitle("DVCS copl. angle");
		H_dvcs_copl.setTitleX("copl angle (^o)");
		H_dvcs_ME = new H1F("H_dvcs_ME","H_dvcs_ME",100,-1,2);
		H_dvcs_ME.setTitle("DVCS miss E");
		H_dvcs_ME.setTitleX("miss E (GeV)");
		H_dvcs_MM_ep = new H1F("H_dvcs_MM_ep","H_dvcs_MM_ep",100,-1,1.2);
		H_dvcs_MM_ep.setTitle("DVCS ep MM^2");
		H_dvcs_MM_ep.setTitleX("ep MM^2 (GeV^2)");
		H_dvcs_MM_eg = new H1F("H_dvcs_MM_eg","H_dvcs_MM_eg",100,-0.5,4);
		H_dvcs_MM_eg.setTitle("DVCS e#gamma MM^2");
		H_dvcs_MM_eg.setTitleX("e#gamma MM^2 (GeV^2)");
		H_dvcs_MM_epg = new H1F("H_dvcs_MM_epg","H_dvcs_MM_epg",100,-0.25,0.25);
		H_dvcs_MM_epg.setTitle("DVCS ep#gamma MM^2");
		H_dvcs_MM_epg.setTitleX("ep#gamma MM^2 (GeV^2)");
		H_dvcs_Phi = new H1F("H_dvcs_Phi","H_dvcs_Phi",100,-180,180);
		H_dvcs_Phi.setTitle("DVCS phi");
		H_dvcs_Phi.setTitleX("#phi (^o)");

		H_pi0_Q2xB = new H2F("H_pi0_Q2xB","H_pi0_Q2xB",100,0,1,100,0,12);
		H_pi0_Q2xB.setTitle("DV#pi^0 Q2 vs xB");
		H_pi0_Q2xB.setTitleX("xB");
		H_pi0_Q2xB.setTitleY("Q2");
		H_pi0_tphi = new H2F("H_pi0_tphi","H_pi0_tphi",100,-180,180,100,0,5);
		H_pi0_tphi.setTitle("DV#pi^0 -t vs #phi");
		H_pi0_tphi.setTitleX("#phi (^o)");
		H_pi0_tphi.setTitleY("-t (GeV)");
		H_gg_open_E = new H2F("H_gg_open_E","H_gg_open_E",100,0,14,100,0,20);
		H_gg_open_E.setTitle("DV#pi^0 #gamma#gamma open angle vs E");
		H_gg_open_E.setTitleX("E #pi^0 (GeV)");
		H_gg_open_E.setTitleY("#gamma#gamma open angle (^o)");
		H_pi0_pT = new H1F("H_pi0_pT","H_pi0_pT",100,0,2.5);
		H_pi0_pT.setTitle("DV#pi^0 missing pT");
		H_pi0_pT.setTitleX("pT (GeV)");
		H_pi0_cone = new H1F("H_pi0_cone","H_pi0_cone",100,0,25);
		H_pi0_cone.setTitle("DV#pi^0 cone angle");
		H_pi0_cone.setTitleX("cone ang (^o)");
		H_pi0_copl = new H1F("H_pi0_copl","H_pi0_copl",100,0,60);
		H_pi0_copl.setTitle("DV#pi^0 copl angle");
		H_pi0_copl.setTitleX("copl ang (^o)");
		H_pi0_ME = new H1F("H_pi0_ME","H_pi0_ME",100,-5,5);
		H_pi0_ME.setTitle("DV#pi^0 miss E");
		H_pi0_ME.setTitleX("miss E (GeV)");
		H_pi0_MM_ep = new H1F("H_pi0_MM_ep","H_pi0_MM_ep",100,-3.5,2.5);
		H_pi0_MM_ep.setTitle("DV#pi^0 ep MM^2");
		H_pi0_MM_ep.setTitleX("ep MM^2 (GeV^2)");
		H_pi0_MM_epi0 = new H1F("H_pi0_MM_epi0","H_pi0_MM_epi0",100,-2.5,7.5);
		H_pi0_MM_epi0.setTitle("DV#pi^0 e#pi^0 MM^2");
		H_pi0_MM_epi0.setTitleX("e#pi^0 MM^2 (GeV^2)");
		H_pi0_MM_eppi0 = new H1F("H_pi0_MM_eppi0","H_pi0_MM_eppi0",100,-2,1);
		H_pi0_MM_eppi0.setTitle("DV#pi^0 ep#pi^0 MM");
		H_pi0_MM_eppi0.setTitleX("ep#pi^0 MM^2 (GeV^2)");
		H_pi0_Phi = new H1F("H_pi0_Phi","H_pi0_Phi",100,-180,180);
		H_pi0_Phi.setTitle("DV#pi^0 #phi");
		H_pi0_Phi.setTitleX("#phi (^0)");
		H_IM_pi0 = new H1F("H_IM_pi0","H_IM_pi0",100,0,0.25);
		H_IM_pi0.setTitle("DV#pi^0 #gamma#gamma IM");
		H_IM_pi0.setTitleX("#gamma#gamma IM (GeV)");

		H_elast_e_th_mom = new H2F[7];
		H_elast_e_th_phi = new H2F[7];
		H_elast_e_phi_mom = new H2F[7];
		H_elast_W_th = new H2F[7];
		H_elast_W_phi = new H2F[7];
		H_elast_p_th_mom = new H2F[7];
		H_elast_p_th_phi = new H2F[7];
		H_elast_p_phi_mom = new H2F[7];
		H_elast_th_th = new H2F[7];
		H_elast_ph_ph = new H2F[7];
		H_elast_EB_bb = new H2F[7];
		H_elast_vz_vz = new H2F[7];
		for(int s=0;s<7;s++){
			H_elast_e_th_mom[s] = new H2F(String.format("H_elast_e_th_mom_S%d",s),String.format("S%d e #theta vs mom",s),100,0,EB,100,0,40);
			H_elast_e_th_phi[s] = new H2F(String.format("H_elast_e_th_phi_S%d",s),String.format("S%d e #theta vs #phi",s),100,-60,60,100,0,40);
			H_elast_e_phi_mom[s] = new H2F(String.format("H_elast_e_phi_mom_S%d",s),String.format("S%d e #phi vs mom",s),100,0,EB,100,-60,60);
			H_elast_W_th[s] = new H2F(String.format("H_elast_W_th_S%d",s),String.format("S%d W vs #theta",s),100,0,40,100,0,4);
			H_elast_W_phi[s] = new H2F(String.format("H_elast_W_phi_S%d",s),String.format("S%d W vs #phi",s),100,-60,60,100,0,4);
			H_elast_p_th_mom[s] = new H2F(String.format("H_elast_p_th_mom_S%d",s),String.format("S%d p #theta vs mom",s),100,0,4,100,0,70);
			H_elast_p_th_phi[s] = new H2F(String.format("H_elast_p_th_phi_S%d",s),String.format("S%d p #theta vs mom",s),100,-90,30,100,0,70);
			H_elast_p_phi_mom[s] = new H2F(String.format("H_elast_p_phi_mom_S%d",s),String.format("S%d p #phi vs mom",s),100,0,4,100,-90,30);
			H_elast_th_th[s] = new H2F(String.format("H_elast_th_th_S%d",s),String.format("S%d #theta vs #theta",s),100,0,70,100,0,40);
			H_elast_ph_ph[s] = new H2F(String.format("H_elast_ph_ph_S%d",s),String.format("S%d #phi vs #phi",s),100,-180,180,100,-180,180);
			H_elast_EB_bb[s] = new H2F(String.format("H_elast_EB_bb_S%d",s),String.format("S%d Ebeam vs #Delta#phi",s),100,-180,180,100,0,50);
			H_elast_vz_vz[s] = new H2F(String.format("H_elast_vz_vz_S%d",s),String.format("S%d vz p vs vz e",s),100,-25,25,100,-25,25);
		}
		H_dvcs_phi = new H1F("H_dvcs_phi","H_dvcs_phi",NPHI,-180,180);
		H_dvcs_phi.setTitle("DVCS #phi counts");
		H_dvcs_phi.setTitleX("#phi (^o)");
		H_dvcs_phi_plus = new H1F("H_dvcs_phi_plus","H_dvcs_phi_plus",NPHI,-180,180);
		H_dvcs_phi_plus.setTitle("DVCS #phi counts h+");
		H_dvcs_phi_plus.setTitleX("#phi (^o)");
		H_dvcs_phi_minus = new H1F("H_dvcs_phi_minus","H_dvcs_phi_minus",NPHI,-180,180);
		H_dvcs_phi_minus.setTitle("DVCS #phi counts h-");
		H_dvcs_phi_minus.setTitleX("#phi (^o)");
		H_pi0_phi = new H1F("H_pi0_phi","H_pi0_phi",NPHI,-180,180);
		H_pi0_phi.setTitle("#pi^0 #phi counts");
		H_pi0_phi.setTitleX("#phi (^o)");
		H_pi0_phi_plus = new H1F("H_pi0_phi_plus","H_pi0_phi_plus",NPHI,-180,180);
		H_pi0_phi_plus.setTitle("#pi^0 #phi counts h+");
		H_pi0_phi_plus.setTitleX("#phi (^o)");
		H_pi0_phi_minus = new H1F("H_pi0_phi_minus","H_pi0_phi_minus",NPHI,-180,180);
		H_pi0_phi_minus.setTitle("#pi^0 #phi counts h-");
		H_pi0_phi_minus.setTitleX("#phi (^o)");
	}
	
	public void setHWP(double val){HWP=val;}

	public double Vangle(Vector3 v1, Vector3 v2){ 
                double res = 0; 
		double l1 = v1.mag();
                double l2 = v2.mag();
                double prod = v1.dot(v2);
                if( l1 * l2 !=0 && Math.abs(prod)<l1*l2 )res = Math.toDegrees( Math.acos(prod/(l1*l2) ) ); 
                return res; 
	}
	public void HasElectron(DataBank recPart){
		for(int p=0;p<recPart.rows() && haz_elec==-1 ;p++){
			//if(recPart.getInt("charge",p)<0){}
			if(recPart.getInt("pid",p)==11 || (isMC && recPart.getInt("charge",p)<0) ){
				float px = recPart.getFloat("px", p);
				float py = recPart.getFloat("py", p);
				float pz = recPart.getFloat("pz", p);
				float vz = recPart.getFloat("vz", p);
				float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
				float theta = (float)Math.toDegrees(Math.acos(pz/mom));
				float phi = (float)Math.toDegrees(Math.atan2(py,px));
				phi+=20;
				if(phi<0)phi+=180;
				//System.out.println("FOUND PID 11 mom="+mom+" , theta="+theta+" , vz="+vz);
				if(mom>1.75 && theta>7 && Math.abs(vz)<15 && theta>17*(1-mom/7) ){
					e_sect = (int)phi/60;
					haz_elec=p;
				}
			}
		}
	}
	public void MakeProton(DataBank recPart){
		int Npos=0, Nneg=0;
		for(int p=0;p<recPart.rows();p++){
			if(recPart.getInt("charge",p)>0)Npos++;
			if(recPart.getInt("charge",p)<0)Nneg++;
		}
		if(true /*Npos==1 && Nneg==1*/ ){
			for(int p=0;p<recPart.rows();p++){
				//if(recPart.getInt("charge",p)>0 ){}
				if(recPart.getInt("pid",p)==2212 || (isMC && recPart.getInt("charge",p)>0) ){
					float px = recPart.getFloat("px", p);
					float py = recPart.getFloat("py", p);
					float pz = recPart.getFloat("pz", p);
					float vz = recPart.getFloat("vz", p);
					prot_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
					prot_theta = (float)Math.toDegrees(Math.acos(pz/prot_mom));
					if(prot_mom>0.4 && prot_mom<3 && prot_theta>15 && prot_theta<75 && Math.abs(vz)<15){
						haz_prot=p;
					}
				}
			}
		}
	}
	public void MakePhotons(DataBank recPart){
		float select_g_E = 0;
		for(int p=0;p<recPart.rows();p++){
			//if(recPart.getInt("pid",p)==22 || recPart.getInt("charge",p)==0)}{
			if(recPart.getInt("charge",p)==0){
				float be = recPart.getFloat("beta", p);
				float px = recPart.getFloat("px", p);
				float py = recPart.getFloat("py", p);
				float pz = recPart.getFloat("pz", p);
				float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
				float theta = (float)Math.toDegrees(Math.acos(pz/mom));
				float phi = (float)Math.toDegrees(Math.atan2(py,px));
				phi+=30;
				if(phi<0)phi+=180;
				int sect = (int)phi/60;
				float phi_sect = phi - 60f*sect;
				phi_sect-=30;
				if( /*sect != e_sect &&*/ mom>1.5 /*&& (be>0.9||sect==2||theta<5) */ 
					&& (theta>8 ||theta<5) && select_g_E<mom && ( Math.abs(phi_sect)<23 || theta<5)
					&& (useFTCAL||theta>8)
					&& (usePCAL ||theta<5)
				){
					select_g_E = mom;haz_g1=p;
				}
			}
		}
		select_g_E = 0;
		if(haz_g1>-1)for(int p=0;p<recPart.rows();p++){
			//if( (recPart.getInt("pid",p)==22 || recPart.getInt("charge",p)==0) && p!=haz_g1){}
			if( recPart.getInt("charge",p)==0 && p!=haz_g1){
				float be = recPart.getFloat("beta", p);
				float px = recPart.getFloat("px", p);
				float py = recPart.getFloat("py", p);
				float pz = recPart.getFloat("pz", p);
				float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
				float theta = (float)Math.toDegrees(Math.acos(pz/mom));
				float phi = (float)Math.toDegrees(Math.atan2(py,px));
				phi+=30;
				if(phi<0)phi+=180;
				int sect = (int)phi/60;
				float phi_sect = phi - 60f*sect;
				phi_sect-=30;
				if( /*sect != e_sect &&*/ mom>0.75 /* && (be>0.9||sect==2||theta<5) */ 
					&& (theta>8||theta<5) && select_g_E<mom && ( Math.abs(phi_sect)<23 || theta<5)
					&& (useFTCAL||theta>8)
					&& (usePCAL ||theta<5)
				){
					select_g_E = mom;haz_g2=p;
				}
			}
		}
	}
        public void fillEBTrack(DataEvent event){
                DataBank bank = event.getBank("REC::Track");
		for(int k = 0; k < bank.rows(); k++){
                        short pind = bank.getShort("pindex",k);
                        if(pind==haz_elec){
                                e_sect = 1;//event.getBank("TimeBasedTrkg::TBTracks").getInt("sector",bank.getShort("index",k));
                        }    
                        if(pind==haz_prot){
                                p_sect = 1;//event.getBank("TimeBasedTrkg::TBTracks").getInt("sector",bank.getShort("index",k));
                        }    
                }    
        }    
	public void MakeParticles(DataBank recPart){
		int p = haz_elec;
		float px = recPart.getFloat("px", p);
		float py = recPart.getFloat("py", p);
		float pz = recPart.getFloat("pz", p);
		float vz = recPart.getFloat("vz", p);
		float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
		VE = new LorentzVector(px,py,pz,mom);
		VGS = new LorentzVector(0,0,0,0);
		VGS.add(VB);
		VGS.sub(VE);
		Vlept = (VB.vect()).cross(VE.vect());
		elec_mom = mom;
		elec_theta = Math.toDegrees(VE.theta());
		elec_vz = recPart.getFloat("vz", p);
		float solenoid_scale = -1.0f;
		float kick_1GeV = 30f;//35f
		elec_phi_sect = Math.toDegrees(VE.phi());
		if(e_sect>3 && elec_phi_sect<0)elec_phi_sect+=360;
		elec_phi_sect +=30f  + solenoid_scale * kick_1GeV/elec_mom ;
		elec_phi_sect -= 60f * (e_sect-1);
		while(elec_phi_sect>60)elec_phi_sect-=60;
		elec_phi_sect -= 30f + solenoid_scale * kick_1GeV/elec_mom;

		p = haz_prot;
		px = recPart.getFloat("px", p);
		py = recPart.getFloat("py", p);
		pz = recPart.getFloat("pz", p);
		float E = (float)Math.sqrt(px*px+py*py+pz*pz + 0.938*0.938);
		VPROT = new LorentzVector(px,py,pz,E);
		Vmand = new LorentzVector(0,0,0,0);
		Vmand.sub(VT);
		Vmand.add(VPROT);
		VmissG = new LorentzVector(0,0,0,0);
		VmissG.add(VB);
		VmissG.add(VT);
		VmissG.sub(VE);
		VmissG.sub(VPROT);
		Vhadr = (VPROT.vect()).cross(VGS.vect());
		prot_beta = recPart.getFloat("beta", p);
		prot_vz = recPart.getFloat("vz", p);
		prot_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
		prot_phi_sect = Math.toDegrees(VPROT.phi());
		if(p_sect>3 && prot_phi_sect<0)prot_phi_sect+=360;
		prot_phi_sect += 30f - solenoid_scale * kick_1GeV/prot_mom;
		prot_phi_sect -= 60f * (p_sect-1);
		while(prot_phi_sect>60)prot_phi_sect-=60;
		prot_phi_sect -= 30f - solenoid_scale * kick_1GeV/prot_mom;
		if(haz_g1>-1){
			p = haz_g1;
			px = recPart.getFloat("px", p);
			py = recPart.getFloat("py", p);
			pz = recPart.getFloat("pz", p);
			E = (float)Math.sqrt(px*px+py*py+pz*pz);
			VG1 = new LorentzVector(px,py,pz,E);
			Vhad2 = (VGS.vect()).cross(VG1.vect());
			VmissP = new LorentzVector(0,0,0,0);
			VmissP.add(VB);
			VmissP.add(VT);
			VmissP.sub(VE);
			VmissP.sub(VG1);
			VMISS = new LorentzVector(0,0,0,0);
			VMISS.add(VB);
			VMISS.add(VT);
			VMISS.sub(VE);
			VMISS.sub(VPROT);
			VMISS.sub(VG1);
			g1_beta = recPart.getFloat("beta", p);
		}
		if(haz_g2>-1){
			p = haz_g2;
			px = recPart.getFloat("px", p);
			py = recPart.getFloat("py", p);
			pz = recPart.getFloat("pz", p);
			E = (float)Math.sqrt(px*px+py*py+pz*pz);
			VG2 = new LorentzVector(px,py,pz,E);
			VmissP.sub(VG2);
			VMISS.sub(VG2);
			VPI0 = new LorentzVector(0,0,0,0);
			VPI0.add(VG1);
			VPI0.add(VG2);
			Vhad2 = (VGS.vect()).cross(VPI0.vect());
			g2_beta = recPart.getFloat("beta", p);
		}
	}
	// VB ; VT ; VE ; VPROT ; Vhadr ; VG1 ; VG2 ; VPI0 ; VmissG
	// VMISS ; VmissP ; vhad2
	public boolean KineCut(){
		boolean res = false;
		LorentzVector W = new LorentzVector(0,0,0,0);
		W.add(VB);
		W.add(VT);
		W.sub(VE);
		//double checkQ2 = 4f*EB*elec_mom*Math.pow(Math.sin(Math.toRadians(elec_theta)/2),2);
		//float checkW = (float)Math.sqrt(0.93827*0.93827 + 2*0.93827*(EB-elec_mom) - checkQ2);
		//System.out.println(checkW + " = " + W.mass() );
		if( haz_g2==-1
		   && VG1.e()>3
		   && Vangle(VG1.vect(),VE.vect())>4
		   && VMISS.e()<1.5 
		   && VMISS.mass2() <0.2 && VMISS.mass2() >-0.2 
		   && VmissP.mass2() < 3 && VmissP.mass2() > -0.25
		   && VmissG.mass2() < 1 && VmissG.mass2() > -1
		   && Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) < 0.3
		   && Vangle(VG1.vect(),VmissG.vect()) < 3
		   && Vangle(Vhad2,Vhadr) < 25
		   //&& Vangle(Vhad2,Vhadr) < 90
		  ){
			isDVCS=true;
			res=true;
			Ndvcs++;
		}
		if(haz_g2>-1 
		   && Vangle(VG1.vect(),VG2.vect())>1.75
		   && Vangle(VG1.vect(),VG2.vect())>5*(1-VPI0.e()/5) 
		   && VPI0.mass()>0.02 && VPI0.mass() < 0.25
		   && Vangle(Vhad2,Vhadr) < 60
		   && VmissP.mass2()>-2 && VmissP.mass2()<7
		   && Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) < 1.5
		   && Vangle(VPI0.vect(),VmissG.vect()) < 20
		   && VMISS.e()>-2 && VMISS.e()<5 
		   && VmissG.mass2()>-2.5 && VmissG.mass2()<2
		   && VMISS.mass2() > -2 && VMISS.mass2() < 1
		   && Vangle(VPI0.vect(),VE.vect())>12
		  ) {
			isPI0=true;
			res=true;
			Npi0++;
		}
		return res;
	}
	public void FillHists(){
		H_elec_theta_mom.fill(  VE.p(),Math.toDegrees(VE.theta()));
		H_prot_theta_mom.fill(VPROT.p(),Math.toDegrees(VPROT.theta()));
		if(isDVCS)H_g_theta_E.fill( VG1.p(),Math.toDegrees(VG1.theta()));
		if(isPI0)H_g1_theta_E.fill( VG1.p(),Math.toDegrees(VG1.theta()));
		H_elec_theta_phi.fill(Math.toDegrees(VE.phi())  ,Math.toDegrees(VE.theta()));
		H_prot_theta_phi.fill(Math.toDegrees(VPROT.phi()),Math.toDegrees(VPROT.theta()));
		if(isDVCS)H_g_theta_phi.fill(Math.toDegrees(VG1.phi()) ,Math.toDegrees(VG1.theta()));
		if(isPI0)H_g1_theta_phi.fill(Math.toDegrees(VG1.phi()) ,Math.toDegrees(VG1.theta()));
		H_elec_phi_mom.fill(  VE.p(),Math.toDegrees(VE.phi()));
		H_prot_phi_mom.fill(VPROT.p(),Math.toDegrees(VPROT.phi()));
		if(isDVCS)H_g_phi_E.fill( VG1.p(),Math.toDegrees(VG1.phi()));
		if(isPI0)H_g1_phi_E.fill( VG1.p(),Math.toDegrees(VG1.phi()));
		H_prot_beta_p.fill(VPROT.p(),prot_beta);
		H_phot_b_p.fill( VG1.p(),g1_beta);
		H_prot_e_vz.fill(elec_vz,prot_vz);
		float TrentoAng = (float)Vangle(Vlept,Vhadr);
		if((VPROT.vect()).dot(Vlept)<0)TrentoAng=-TrentoAng;
		if(isMC){
			float ThBSA = 0.8f*0.18f * (float)Math.sin(Math.toRadians(TrentoAng));
			float myRand = 2*myMC.nextFloat()-1;
			if(myRand<ThBSA)HELI=1;
			else HELI = -1;
		}
		if(haz_g2==-1){}
		if(isDVCS){
			H_dvcs_Q2xB.fill(-VGS.mass2()/(2*0.938*VGS.e()),-VGS.mass2());
			H_dvcs_tphi.fill(TrentoAng,-Vmand.mass2());
			H_dvcs_copl.fill(Vangle(Vhad2,Vhadr));
			H_dvcs_MM_eg.fill(VmissP.mass2());
			H_dvcs_pT.fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()));
			//H_dvcs_cone.fill(Vangle(VPROT.vect(),VmissP.vect()));
			H_dvcs_cone.fill(Vangle(VG1.vect(),VmissG.vect()));
			H_dvcs_ME.fill(VMISS.e());
			H_dvcs_MM_ep.fill(VmissG.mass2());
			H_dvcs_MM_epg.fill(VMISS.mass2());
			H_dvcs_Phi.fill(TrentoAng);
			H_prot_miss_phi.fill(Math.toDegrees(VPROT.phi()),Math.toDegrees(VmissP.phi()));
			H_phot_elec_angle.fill(TrentoAng,Vangle(VG1.vect(),VE.vect()));
			H_phot_virt_angle.fill(-Vmand.mass2(),Vangle(VG1.vect(),VGS.vect()));
			H_dvcs_phi.fill(TrentoAng);
			if(HWP*HELI==1)H_dvcs_phi_plus.fill(TrentoAng);
			if(HWP*HELI==-1)H_dvcs_phi_minus.fill(TrentoAng);

		}
		if(haz_g2>-1){}
		if(isPI0){
			H_g2_theta_E.fill( VG2.p(),Math.toDegrees(VG2.theta()));
			H_g2_theta_phi.fill(Math.toDegrees(VG2.phi()) ,Math.toDegrees(VG2.theta()));
			H_g2_phi_E.fill( VG2.p(),Math.toDegrees(VG2.phi()));
		        H_phot_b_p.fill( VG2.p(),g2_beta);
			H_gg_open_E.fill(VPI0.e() , Vangle(VG1.vect(),VG2.vect()));
			H_IM_pi0.fill(VPI0.mass());
			H_pi0_Q2xB.fill(-VGS.mass2()/(2*0.938*VGS.e()),-VGS.mass2());
			H_pi0_tphi.fill(TrentoAng,-Vmand.mass2());
			H_pi0_copl.fill(Vangle(Vhad2,Vhadr));
			H_pi0_MM_epi0.fill(VmissP.mass2());
			H_pi0_pT.fill(Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()));
			//H_pi0_cone.fill(Vangle(VPROT.vect(),VmissP.vect()));
			H_pi0_cone.fill(Vangle(VPI0.vect(),VmissG.vect()));
			H_pi0_ME.fill(VMISS.e());
			H_pi0_MM_ep.fill(VmissG.mass2());
			H_pi0_MM_eppi0.fill(VMISS.mass2());
			H_pi0_Phi.fill(TrentoAng);
			H_prot_pi0_miss_phi.fill(Math.toDegrees(VPROT.phi()),Math.toDegrees(VmissP.phi()));
			H_pi0_elec_angle.fill(TrentoAng,Vangle(VPI0.vect(),VE.vect()));
			H_pi0_virt_angle.fill(-Vmand.mass2(),Vangle(VPI0.vect(),VGS.vect()));
			H_pi0_phi.fill(TrentoAng);
			if(HWP*HELI==1)H_pi0_phi_plus.fill(TrentoAng);
			if(HWP*HELI==-1)H_pi0_phi_minus.fill(TrentoAng);
		}
	}
	public boolean IsElastic(){
                boolean res = false;
                LorentzVector W = new LorentzVector(0,0,0,0);
                W.add(VB);
                W.add(VT);
                W.sub(VE);
		elast_W = W.mass();
		double tantan = Math.tan(VPROT.theta()) * Math.tan(VE.theta()/2);
		elast_EB = 0.938*(1-tantan)/tantan;
		elast_bb = Math.toDegrees( VE.phi()-VPROT.phi() ) + 180;
		while(elast_bb>180)elast_bb-=360;
		while(elast_bb<-180)elast_bb+=360;
		return true;
	}
	public void FillElastic(int s){
		H_elast_e_th_mom[s].fill(VE.p(),Math.toDegrees(VE.theta()));
		H_elast_e_th_phi[s].fill(elec_phi_sect  ,Math.toDegrees(VE.theta()));
		H_elast_e_phi_mom[s].fill(VE.p(),elec_phi_sect);
		H_elast_p_th_mom[s].fill(VPROT.p(),Math.toDegrees(VPROT.theta()));
		H_elast_p_th_phi[s].fill(prot_phi_sect,Math.toDegrees(VPROT.theta()));
		H_elast_p_phi_mom[s].fill(VPROT.p(),prot_phi_sect);
		H_elast_th_th[s].fill(Math.toDegrees(VPROT.theta()),Math.toDegrees(VE.theta()));
		H_elast_ph_ph[s].fill(Math.toDegrees(VPROT.phi()),Math.toDegrees(VE.phi()));
		H_elast_W_th[s].fill(Math.toDegrees(VE.theta()),elast_W);
		H_elast_W_phi[s].fill(elec_phi_sect,elast_W);
		H_elast_EB_bb[s].fill(elast_bb,elast_EB);
		H_elast_vz_vz[s].fill(elec_vz,prot_vz);
	}

        public void processEvent(DataEvent event) {
		isDVCS=false;
		isPI0=false;
		HELI=0;haz_elec=-1;haz_prot=-1;haz_g1=-1;haz_g2=-1;
		if(event.hasBank("MC::Event"))isMC=true;
		if(event.hasBank("HEL::flip")){
			
			if(event.getBank("HEL::flip").getInt("helicity",0)==1)HELI=1;
			if(event.getBank("HEL::flip").getInt("helicity",0)==-1)HELI=-1;
		}
		if(event.hasBank("REC::Particle"))HasElectron(event.getBank("REC::Particle"));
		if(haz_elec>-1)MakeProton(event.getBank("REC::Particle"));
		if(haz_elec>-1 && haz_prot>-1 )MakePhotons(event.getBank("REC::Particle"));
		if(haz_elec>-1 && haz_prot>-1 && haz_g1>-1){
			//fillEBTrack(event);
			MakeParticles(event.getBank("REC::Particle"));
			//if(IsElastic()){FillElastic(e_sect-1);FillElastic(6);}
			if(haz_g1>-1)if(KineCut())FillHists();
		}
	}
		
        public void plot() {
                EmbeddedCanvas can_dvcspi0_parts  = new EmbeddedCanvas();
                can_dvcspi0_parts.setSize(3500,1500);
                can_dvcspi0_parts.divide(7,3);
                can_dvcspi0_parts.setAxisTitleSize(18);
                can_dvcspi0_parts.setAxisFontSize(18);
                can_dvcspi0_parts.setTitleSize(18);
		can_dvcspi0_parts.cd(0);can_dvcspi0_parts.draw(H_elec_theta_mom);
		can_dvcspi0_parts.cd(1);can_dvcspi0_parts.draw(H_prot_theta_mom);
		can_dvcspi0_parts.cd(2);can_dvcspi0_parts.draw(H_g_theta_E);
		can_dvcspi0_parts.cd(3);can_dvcspi0_parts.draw(H_g1_theta_E);
		can_dvcspi0_parts.cd(4);can_dvcspi0_parts.draw(H_g2_theta_E);
		can_dvcspi0_parts.cd(5);can_dvcspi0_parts.draw(H_prot_miss_phi);
		can_dvcspi0_parts.cd(6);can_dvcspi0_parts.draw(H_prot_pi0_miss_phi);
	
		can_dvcspi0_parts.cd(7);can_dvcspi0_parts.draw(H_elec_theta_phi);
		can_dvcspi0_parts.cd(8);can_dvcspi0_parts.draw(H_prot_theta_phi);
		can_dvcspi0_parts.cd(9);can_dvcspi0_parts.draw(H_g_theta_phi);
		can_dvcspi0_parts.cd(10);can_dvcspi0_parts.draw(H_g1_theta_phi);
		can_dvcspi0_parts.cd(11);can_dvcspi0_parts.draw(H_g2_theta_phi);
		can_dvcspi0_parts.cd(12);can_dvcspi0_parts.draw(H_phot_elec_angle);
		can_dvcspi0_parts.cd(13);can_dvcspi0_parts.draw(H_pi0_elec_angle);
		
		can_dvcspi0_parts.cd(14);can_dvcspi0_parts.draw(H_elec_phi_mom);
		can_dvcspi0_parts.cd(15);can_dvcspi0_parts.draw(H_prot_phi_mom);
		can_dvcspi0_parts.cd(16);can_dvcspi0_parts.draw(H_g_phi_E);
		can_dvcspi0_parts.cd(17);can_dvcspi0_parts.draw(H_g1_phi_E);
		can_dvcspi0_parts.cd(18);can_dvcspi0_parts.draw(H_g2_phi_E);
		can_dvcspi0_parts.cd(19);can_dvcspi0_parts.draw(H_phot_virt_angle);
		can_dvcspi0_parts.cd(20);can_dvcspi0_parts.draw(H_pi0_virt_angle);
                
		can_dvcspi0_parts.save(String.format("plots"+runNum+"/dvcspi0_parts.png"));
                System.out.println(String.format("saved plots"+runNum+"/dvcspi0_parts.png"));

		EmbeddedCanvas can_dvcspi0_id = new EmbeddedCanvas();
		can_dvcspi0_id.setSize(2500,2500);
		can_dvcspi0_id.divide(5,5);
		can_dvcspi0_id.setAxisTitleSize(18);
		can_dvcspi0_id.setAxisFontSize(18);
		can_dvcspi0_id.setTitleSize(18);
		can_dvcspi0_id.cd(0);can_dvcspi0_id.draw(H_prot_beta_p);
		F1D prot_beta_line = new F1D("prot_beta_line","x/sqrt(x*x + 0.93827*0.93827)",0.5,3);
		can_dvcspi0_id.draw(prot_beta_line,"same");
		can_dvcspi0_id.cd(1);can_dvcspi0_id.draw(H_prot_e_vz);
		can_dvcspi0_id.cd(2);can_dvcspi0_id.draw(H_gg_open_E);
		can_dvcspi0_id.cd(3);can_dvcspi0_id.draw(H_IM_pi0);
		can_dvcspi0_id.cd(4);can_dvcspi0_id.draw(H_phot_b_p);

		can_dvcspi0_id.cd(5);can_dvcspi0_id.draw(H_dvcs_Q2xB);
		can_dvcspi0_id.cd(6);can_dvcspi0_id.draw(H_dvcs_tphi);
		can_dvcspi0_id.cd(7);can_dvcspi0_id.draw(H_dvcs_copl);
		can_dvcspi0_id.cd(8);can_dvcspi0_id.draw(H_dvcs_MM_eg);
		can_dvcspi0_id.cd(9);can_dvcspi0_id.draw(H_dvcs_pT);
		can_dvcspi0_id.cd(10);can_dvcspi0_id.draw(H_dvcs_cone);
		can_dvcspi0_id.cd(11);can_dvcspi0_id.draw(H_dvcs_ME);
		can_dvcspi0_id.cd(12);can_dvcspi0_id.draw(H_dvcs_MM_ep);
		can_dvcspi0_id.cd(13);can_dvcspi0_id.draw(H_dvcs_MM_epg);
		can_dvcspi0_id.cd(14);can_dvcspi0_id.draw(H_dvcs_Phi);

		can_dvcspi0_id.cd(15);can_dvcspi0_id.draw(H_pi0_Q2xB);
		can_dvcspi0_id.cd(16);can_dvcspi0_id.draw(H_pi0_tphi);
		can_dvcspi0_id.cd(17);can_dvcspi0_id.draw(H_pi0_copl);
		can_dvcspi0_id.cd(18);can_dvcspi0_id.draw(H_pi0_MM_epi0);
		can_dvcspi0_id.cd(19);can_dvcspi0_id.draw(H_pi0_pT);
		can_dvcspi0_id.cd(20);can_dvcspi0_id.draw(H_pi0_cone);
		can_dvcspi0_id.cd(21);can_dvcspi0_id.draw(H_pi0_ME);
		can_dvcspi0_id.cd(22);can_dvcspi0_id.draw(H_pi0_MM_ep);
		can_dvcspi0_id.cd(23);can_dvcspi0_id.draw(H_pi0_MM_eppi0);
		can_dvcspi0_id.cd(24);can_dvcspi0_id.draw(H_pi0_Phi);
		can_dvcspi0_id.save(String.format("plots"+runNum+"/dvcspi0.png"));
		System.out.println(String.format("saved plots"+runNum+"/dvcspi0.png"));

		g_dvcs_bsa = new GraphErrors();
		g_dvcs_bsa.setName("g_dvcs_bsa");
		String stringBSA = "double BSA[NPHI] = {";
		String strindBSA = "double dBSA[NPHI] = {";
		String stringNpos = "double NPOS[NPHI] = {";
		String stringNneg = "double NNEG[NPHI] = {";
		for(int p=0;p<NPHI;p++){
			float PHI = (float)H_dvcs_phi.getDataX(p);
			float Np = (float)H_dvcs_phi_plus.getDataY(p);
			float Nm = (float)H_dvcs_phi_minus.getDataY(p);
			float BSA = (Np-Nm) / (Np+Nm) / 0.8f;
			float dBSA = 1/0.8f * (float)Math.sqrt( (1 - BSA*BSA*0.8f*0.8f) / (Np+Nm) );
			g_dvcs_bsa.addPoint(PHI,BSA,0,dBSA);
			System.out.println(p+ " : DVCS phi="+PHI+" , N+="+Np+" , N-="+Nm+" , BSA="+BSA+" , dBSA="+dBSA);
			stringBSA += String.format("%1.4f, ",BSA);
			strindBSA += String.format("%1.4f, ",dBSA);
			stringNpos += String.format("%1.0f, ",Np);
			stringNneg += String.format("%1.0f, ",Nm);
		}
		stringBSA += "};";
		strindBSA += "};";
		stringNpos += "};";
		stringNneg += "};";
		System.out.println(stringBSA);
		System.out.println(strindBSA);
		System.out.println(stringNpos);
		System.out.println(stringNneg);
		System.out.println("<xB>="+H_dvcs_Q2xB.projectionX().getMean()+" , <Q2>="+H_dvcs_Q2xB.projectionY().getMean());
		System.out.println("<-t>="+H_dvcs_tphi.projectionY().getMean());
                EmbeddedCanvas can_dvcs_asym  = new EmbeddedCanvas();
                can_dvcs_asym.setSize(1500,1500);
                can_dvcs_asym.divide(2,2);
                can_dvcs_asym.setAxisTitleSize(18);
                can_dvcs_asym.setAxisFontSize(18);
                can_dvcs_asym.setTitleSize(18);
		H_dvcs_phi_plus.setLineColor(2);H_dvcs_phi_minus.setLineColor(4);
		can_dvcs_asym.cd(0);can_dvcs_asym.draw(H_dvcs_phi);
		can_dvcs_asym.cd(1);can_dvcs_asym.draw(H_dvcs_phi_plus);can_dvcs_asym.draw(H_dvcs_phi_minus,"same");
		can_dvcs_asym.cd(2);can_dvcs_asym.draw(g_dvcs_bsa);
		can_dvcs_asym.save(String.format("plots"+runNum+"/dvcsasym.png"));
		System.out.println(String.format("saved plots"+runNum+"/dvcsasym.png"));

                EmbeddedCanvas can_pi0_asym  = new EmbeddedCanvas();
                can_pi0_asym.setSize(1500,1500);
                can_pi0_asym.divide(2,2);
                can_pi0_asym.setAxisTitleSize(18);
                can_pi0_asym.setAxisFontSize(18);
                can_pi0_asym.setTitleSize(18);
		H_pi0_phi_plus.setLineColor(2);H_pi0_phi_minus.setLineColor(4);
		can_pi0_asym.cd(0);can_pi0_asym.draw(H_pi0_phi);
		can_pi0_asym.cd(1);can_pi0_asym.draw(H_pi0_phi_plus);can_pi0_asym.draw(H_pi0_phi_minus,"same");
		can_pi0_asym.save(String.format("plots"+runNum+"/pi0asym.png"));
		System.out.println(String.format("saved plots"+runNum+"/pi0asym.png"));

		EmbeddedCanvas can_elast = new EmbeddedCanvas();
		//can_elast.setSize(3500,5500);
		can_elast.setSize(3500,6000);
		can_elast.divide(7,12);
		can_elast.setAxisTitleSize(18);
		can_elast.setAxisFontSize(18);
		can_elast.setTitleSize(18);
		for(int s=0;s<7;s++){
			can_elast.cd(s);can_elast.draw(H_elast_e_th_mom[s]);
			can_elast.cd(7+s);can_elast.draw(H_elast_e_th_phi[s]);
			can_elast.cd(14+s);can_elast.draw(H_elast_e_phi_mom[s]);
			can_elast.cd(21+s);can_elast.draw(H_elast_p_th_mom[s]);
			can_elast.cd(28+s);can_elast.draw(H_elast_p_th_phi[s]);
			can_elast.cd(35+s);can_elast.draw(H_elast_p_phi_mom[s]);
			can_elast.cd(42+s);can_elast.draw(H_elast_th_th[s]);
			can_elast.cd(49+s);can_elast.draw(H_elast_ph_ph[s]);
			can_elast.cd(56+s);can_elast.draw(H_elast_EB_bb[s]);
			can_elast.cd(63+s);can_elast.draw(H_elast_W_th[s]);
			can_elast.cd(70+s);can_elast.draw(H_elast_W_phi[s]);
			can_elast.cd(77+s);can_elast.draw(H_elast_vz_vz[s]);
		}
		if(runNum>0){
			can_elast.save(String.format("plots"+runNum+"/elast.png"));
		}
		//else{
		//	can_elast.save(String.format("plots/elast.png"));
		//}
	}
        public void write() {
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/parts/");
                dirout.cd("/parts/");
		dirout.addDataSet(H_elec_theta_mom,H_prot_theta_mom,H_g_theta_E,H_g1_theta_E,H_g2_theta_E);
		dirout.addDataSet(H_elec_theta_phi,H_prot_theta_phi,H_g_theta_phi,H_g1_theta_phi,H_g2_theta_phi);
		dirout.addDataSet(H_elec_phi_mom,H_prot_phi_mom,H_g_phi_E,H_g1_phi_E,H_g2_phi_E);
		dirout.addDataSet(H_prot_miss_phi,H_prot_pi0_miss_phi,H_phot_elec_angle,H_pi0_elec_angle,H_phot_virt_angle,H_pi0_virt_angle);
                dirout.addDataSet(H_prot_beta_p,H_prot_e_vz,H_phot_b_p);
		dirout.mkdir("/dvcs/");
                dirout.cd("/dvcs/");
		dirout.addDataSet(H_dvcs_Q2xB,H_dvcs_tphi,H_dvcs_copl,H_dvcs_MM_eg,H_dvcs_pT,H_dvcs_cone);
		dirout.addDataSet(H_dvcs_ME,H_dvcs_MM_ep,H_dvcs_MM_epg,H_dvcs_Phi);
		dirout.addDataSet(g_dvcs_bsa);
                dirout.mkdir("/pi0/");
                dirout.cd("/pi0/");
		dirout.addDataSet(H_IM_pi0,H_pi0_Q2xB,H_pi0_tphi,H_pi0_copl,H_pi0_MM_epi0,H_pi0_pT,H_pi0_cone);
		dirout.addDataSet(H_pi0_ME,H_pi0_MM_ep,H_pi0_MM_eppi0,H_pi0_Phi);
		if(runNum>0)dirout.writeFile("plots"+runNum+"/dvcspi0_"+runNum+".hipo");
                else dirout.writeFile("plots/dvcspi0.hipo");
	}

////////////////////////////////////////////////
        public static void main(String[] args) {
                System.setProperty("java.awt.headless", "true");
		GStyle.setPalette("kRainBow");
                int count = 0;
		int runNum = 0;
		float EBreq = 7;
		int reqNphi=12;
		int maxevents = 500000000;
		int CaloChoice = 0;// 0=both 1=FTCAL 2=PCAL
		String filelist = "list_of_files.txt";
		if(args.length>0)runNum = Integer.parseInt(args[0]);
		if(args.length>1)filelist = args[1];
		if(args.length>2)EBreq = Float.parseFloat(args[2]);
		if(args.length>3)reqNphi = Integer.parseInt(args[3]);
		if(args.length>4)maxevents = Integer.parseInt(args[4]);
		if(args.length>5)CaloChoice = Integer.parseInt(args[5]);
		if( CaloChoice!=0 && CaloChoice!=1 && CaloChoice!=2){
			System.out.println("Wrong calorimeter choice! "+CaloChoice+" ... aborting");
			return;
		}
		dvcspi0_hipo4 ana = new dvcspi0_hipo4(runNum,EBreq,reqNphi,CaloChoice);
                List<String> toProcessFileNames = new ArrayList<String>();
                File file = new File(filelist);
                Scanner read;
                try {
                        read = new Scanner(file);
                        do { 
                                String filename = read.next();
                                toProcessFileNames.add(filename);

                        }while (read.hasNext());
                        read.close();
                }catch(IOException e){
                        e.printStackTrace();
                }
		int progresscount=0;int filetot = toProcessFileNames.size();
		for (String runstrg : toProcessFileNames)if(count<maxevents) {
			progresscount++;
			System.out.println(String.format(">>>>>>>>>>>>>>>> %s",runstrg));
			if(runstrg.indexOf("phys2_") != -1) {
				int runnum = Integer.parseInt(runstrg.substring(runstrg.indexOf("phys2_")+6,runstrg.indexOf("phys2_")+10));
				System.out.println("Extracted run num "+runnum);
				if(runnum>3850 && runnum<4248)ana.setHWP(-1);
				if(runnum>4247)ana.setHWP(1);
			}
			File varTmpDir = new File(runstrg);
			if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
			System.out.println("READING NOW "+runstrg);
			HipoDataSource reader = new HipoDataSource();
			reader.open(runstrg);
			int filecount = 0;
			while( reader.hasEvent()&& count<maxevents && count<maxevents){
				DataEvent event = reader.getNextEvent();
				ana.processEvent(event);
				filecount++;count++;
				if(count%10000 == 0) System.out.println(count/1000 + "k events (this is ep(g) analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot);
			}
			reader.close();
		}
		System.out.println("Total : " + count + " events ; N dvcs=" + ana.Ndvcs + " ; N pi0 = " + ana.Npi0 );
		ana.plot();
		ana.write();
        }
}
