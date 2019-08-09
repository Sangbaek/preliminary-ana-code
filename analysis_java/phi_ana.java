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

public class phi_ana {
        public int CHAN;
	public double EB;
	public int Nexclu, runNum;
	public int e_sect, p_sect, haz_elec, haz_prot, haz_kp, haz_km;
        public double elec_mom, elec_theta, elec_phi, elec_vz, elec_phi_sect;
        public double prot_mom, prot_theta, prot_phi, prot_vz, prot_beta, prot_phi_sect;
	public double kp_mom, kp_theta, kp_phi, kp_vz, kp_beta;
	public double km_mom, km_theta, km_phi, km_vz, km_beta;
	public double coplanarity;

        public LorentzVector VB, VT, VE, VGS, VPROT, Vmand, VKP, VKM, VPHI, VMISS, VmissP, VmissPHI, VmissKP, VmissKM, VLambda, VmissLambda;
	public Vector3 Vlept, Vhadr, Vhad2;
        public int oldInd, NOLD;
	public int doMix;
	public LorentzVector[] VOLDGS, VOLDKP, VOLDKM;

        public H2F H_elec_theta_mom, H_elec_theta_phi, H_elec_phi_mom, H_elec_vz_theta;
        public H2F H_prot_theta_mom, H_prot_theta_phi, H_prot_phi_mom, H_prot_vz_theta;
        public H2F H_kp_theta_mom, H_kp_theta_phi, H_kp_phi_mom, H_kp_vz_theta;
        public H2F H_km_theta_mom, H_km_theta_phi, H_km_phi_mom, H_km_vz_theta;
	public H2F H_prot_beta_mom, H_kp_beta_mom, H_km_beta_mom;

	public H1F H_phi_nocut_IM;
	
	public H2F H_CM_angles;
	public H1F H_phi_IM, H_phi_copl, H_phi_IM_back;
	public H2F H_Q2_xB, H_t_phi;
	public H2F H_prot_miss_azimuth, H_misspT_missE, H_IMphi_MMep, H_MMepkp_MMepkm;
       	public H2F H_IMphi_MMepphi, H_angleKp_angle_Km, H_angleKp_angleP, H_KKopen_phiE;

	public H2F H_Q2_xB_cut, H_t_phi_cut;

	public H1F H_IM_KmP;
	public H2F H_IM_KmP_MMeKp, H_IM_KmP_IMphi;

	public phi_ana(int reqrunNum, float EBreq) {
                CHAN = 1;
		runNum = reqrunNum;
                EB = 10.6f;//6.42f;
                if(EBreq<8f)EB = 6.42f;
                Nexclu = 0;
                VB = new LorentzVector(0,0,EB,EB);
                VT = new LorentzVector(0,0,0,0.938);
		
		NOLD = 2;
		VOLDGS = new LorentzVector[NOLD];
		VOLDKP = new LorentzVector[NOLD];
		VOLDKM = new LorentzVector[NOLD];
		for(oldInd=0;oldInd<NOLD;oldInd++){
			VOLDGS[oldInd] = new LorentzVector(0,0,0,0);
			VOLDKP[oldInd] = new LorentzVector(0,0,0,0);
			VOLDKM[oldInd] = new LorentzVector(0,0,0,0);
		}
		oldInd = 0;
		doMix = 0;
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
                H_elec_vz_theta = new H2F("H_elec_vz_theta","H_elec_vz_theta",100,0,40,100,-25,25);
                H_elec_vz_theta.setTitle("elec vz vs #theta");
                H_elec_vz_theta.setTitleX("#theta (^o)");
                H_elec_vz_theta.setTitleY("vz (cm)");
                H_prot_theta_mom = new H2F("H_prot_theta_mom","H_prot_theta_mom",100,0,4,100,0,90);
                H_prot_theta_mom.setTitle("prot #theta vs p");
                H_prot_theta_mom.setTitleX("p (GeV)");
                H_prot_theta_mom.setTitleY("#theta (^o)");
                H_prot_theta_phi = new H2F("H_prot_theta_phi","H_prot_theta_phi",100,-180,180,100,0,90);
                H_prot_theta_phi.setTitle("prot #theta vs #phi");
                H_prot_theta_phi.setTitleX("#phi (^o)");
                H_prot_theta_phi.setTitleY("#theta (^o)");
                H_prot_phi_mom = new H2F("H_prot_phi_mom","H_prot_phi_mom",100,0,4,100,-180,180);
                H_prot_phi_mom.setTitle("prot #phi vs p");
                H_prot_phi_mom.setTitleX("p (GeV)");
                H_prot_phi_mom.setTitleY("#phi (^o)");
                H_prot_vz_theta = new H2F("H_prot_vz_theta","H_prot_vz_theta",100,0,90,100,-25,25);
                H_prot_vz_theta.setTitle("prot vz vs #theta");
                H_prot_vz_theta.setTitleX("#theta (^o)");
                H_prot_vz_theta.setTitleY("vz (cm)");
                H_kp_theta_mom = new H2F("H_kp_theta_mom","H_kp_theta_mom",100,0,EB/2,100,0,50);
                if(CHAN==1)H_kp_theta_mom.setTitle("K^+ #theta vs p");
                if(CHAN==2)H_kp_theta_mom.setTitle("#pi^+ #theta vs p");
                H_kp_theta_mom.setTitleX("p (GeV)");
                H_kp_theta_mom.setTitleY("#theta (^o)");
                H_kp_theta_phi = new H2F("H_kp_theta_phi","H_kp_theta_phi",100,-180,180,100,0,50);
                if(CHAN==1)H_kp_theta_phi.setTitle("K^+ #theta vs #phi");
                if(CHAN==2)H_kp_theta_phi.setTitle("#pi^+ #theta vs #phi");
                H_kp_theta_phi.setTitleX("#phi (^o)");
                H_kp_theta_phi.setTitleY("#theta (^o)");
                H_kp_phi_mom = new H2F("H_kp_phi_mom","H_kp_phi_mom",100,0,EB/2,100,-180,180);
                if(CHAN==1)H_kp_phi_mom.setTitle("K^+ #phi vs p");
                if(CHAN==2)H_kp_phi_mom.setTitle("#pi^+ #phi vs p");
                H_kp_phi_mom.setTitleX("p (GeV)");
                H_kp_phi_mom.setTitleY("#phi (^o)");
                H_kp_vz_theta = new H2F("H_kp_vz_theta","H_kp_vz_theta",100,0,50,100,-25,25);
                if(CHAN==1)H_kp_vz_theta.setTitle("K^+ vz vs #theta");
                if(CHAN==2)H_kp_vz_theta.setTitle("#pi^+ vz vs #theta");
                H_kp_vz_theta.setTitleX("#theta (^o)");
                H_kp_vz_theta.setTitleY("vz (cm)");
                H_km_theta_mom = new H2F("H_km_theta_mom","H_km_theta_mom",100,0,EB/2,100,0,50);
                if(CHAN==1)H_km_theta_mom.setTitle("K^- #theta vs p");
                if(CHAN==2)H_km_theta_mom.setTitle("#pi^- #theta vs p");
                H_km_theta_mom.setTitleX("p (GeV)");
                H_km_theta_mom.setTitleY("#theta (^o)");
                H_km_theta_phi = new H2F("H_km_theta_phi","H_km_theta_phi",100,-180,180,100,0,50);
                if(CHAN==1)H_km_theta_phi.setTitle("K^- #theta vs #phi");
                if(CHAN==2)H_km_theta_phi.setTitle("#pi^- #theta vs #phi");
                H_km_theta_phi.setTitleX("#phi (^o)");
                H_km_theta_phi.setTitleY("#theta (^o)");
                H_km_phi_mom = new H2F("H_km_phi_mom","H_km_phi_mom",100,0,EB/2,100,-180,180);
                if(CHAN==1)H_km_phi_mom.setTitle("K^- #phi vs p");
                if(CHAN==2)H_km_phi_mom.setTitle("#pi^- #phi vs p");
                H_km_phi_mom.setTitleX("p (GeV)");
                H_km_phi_mom.setTitleY("#phi (^o)");
                H_km_vz_theta = new H2F("H_km_vz_theta","H_km_vz_theta",100,0,50,100,-25,25);
                if(CHAN==1)H_km_vz_theta.setTitle("K^- vz vs #theta");
                if(CHAN==2)H_km_vz_theta.setTitle("#pi^- vz vs #theta");
                H_km_vz_theta.setTitleX("#theta (^o)");
                H_km_vz_theta.setTitleY("vz (cm)");
		H_prot_beta_mom = new H2F("H_prot_beta_mom","H_prot_beta_mom",100,0,4,100,0,1.2);
		H_prot_beta_mom.setTitle("prot #beta vs p");
		H_prot_beta_mom.setTitleX("p (GeV)");
		H_prot_beta_mom.setTitleY("#beta");
		H_kp_beta_mom = new H2F("H_kp_beta_mom","H_kp_beta_mom",100,0,EB/2,100,0,1.2);
		if(CHAN==1)H_kp_beta_mom.setTitle("K^+ #beta vs p");
		if(CHAN==2)H_kp_beta_mom.setTitle("#pi^+ #beta vs p");
		H_kp_beta_mom.setTitleX("p (GeV)");
		H_kp_beta_mom.setTitleY("#beta");
		H_km_beta_mom = new H2F("H_km_beta_mom","H_km_beta_mom",100,0,EB/2,100,0,1.2);
		if(CHAN==1)H_km_beta_mom.setTitle("K^- #beta vs p");
		if(CHAN==2)H_km_beta_mom.setTitle("#pi^- #beta vs p");
		H_km_beta_mom.setTitleX("p (GeV)");
		H_km_beta_mom.setTitleY("#beta");
		H_phi_IM = new H1F("H_phi_IM","H_phi_IM",100,0.8,2.8);
		if(CHAN==1)H_phi_IM.setTitle("K^+K^- IM");
		if(CHAN==1)H_phi_IM.setTitleX("K^+K^- IM (GeV)");
		if(CHAN==2)H_phi_IM.setTitle("#pi^+#pi^- IM");
		if(CHAN==2)H_phi_IM.setTitleX("#pi^+#pi^- IM (GeV)");
		H_phi_nocut_IM = new H1F("H_phi_nocut_IM","H_phi_nocut_IM",100,0.8,2.8);
		if(CHAN==1)H_phi_nocut_IM.setTitle("K^+K^- IM before cuts");
		if(CHAN==1)H_phi_nocut_IM.setTitleX("K^+K^- IM (GeV)");
		if(CHAN==2)H_phi_nocut_IM.setTitle("#pi^+#pi^- IM before cuts");
		if(CHAN==2)H_phi_nocut_IM.setTitleX("#pi^+#pi^- IM (GeV)");
		H_phi_IM_back = new H1F("H_phi_IM_back","H_phi_IM_back",100,0.8,2.8);
		if(CHAN==1)H_phi_IM_back.setTitle("Mixing K^+K^- IM");
		if(CHAN==1)H_phi_IM_back.setTitleX("K^+K^- IM (GeV)");
		if(CHAN==2)H_phi_IM_back.setTitle("#pi^+#pi^- IM");
		if(CHAN==2)H_phi_IM_back.setTitleX("#pi^+#pi^- IM (GeV)");
		H_phi_copl = new H1F("H_phi_copl","H_phi_copl",100,0,180);
		H_phi_copl.setTitle("coplanarity angle");
		H_phi_copl.setTitleX("#Delta#phi (^o)");
		
		H_Q2_xB = new H2F("H_Q2_xB","H_Q2_xB",100,0,1,100,0,EB);
		H_Q2_xB.setTitle("Q^2 vs x_B");
		H_Q2_xB.setTitleX("x_B");
		H_Q2_xB.setTitleY("Q^2");
		H_t_phi = new H2F("H_t_phi","H_t_phi",100,-180,180,100,0,5);
		H_t_phi.setTitle("-t vs #phi");
		H_t_phi.setTitleX("#phi (^o)");
		H_t_phi.setTitleY("-t (GeV^2)");

		H_Q2_xB_cut = new H2F("H_Q2_xB_cut","H_Q2_xB_cut",100,0,1,100,0,EB);
		H_Q2_xB_cut.setTitle("Q^2 vs x_B");
		H_Q2_xB_cut.setTitleX("x_B");
		H_Q2_xB_cut.setTitleY("Q^2");
		H_t_phi_cut = new H2F("H_t_phi_cut","H_t_phi_cut",100,-180,180,100,0,5);
		H_t_phi_cut.setTitle("-t vs #phi");
		H_t_phi_cut.setTitleX("#phi (^o)");
		H_t_phi_cut.setTitleY("-t (GeV^2)");

		H_misspT_missE = new H2F("H_misspT_missE","H_misspT_missE",100,-3,5,100,0,2.5);
		if(CHAN==1)H_misspT_missE.setTitle("DV#phi miss pT vs miss E");
		if(CHAN==2)H_misspT_missE.setTitle("DV#rho miss pT vs miss E");
		H_misspT_missE.setTitleX("miss E (GeV)");
		H_misspT_missE.setTitleY("miss pT (GeV)");
		H_prot_miss_azimuth = new H2F("H_prot_miss_azimuth","H_prot_miss_azimuth",100,-180,180,100,-180,180);
		if(CHAN==1)H_prot_miss_azimuth.setTitle("Proton vs #phi");
		if(CHAN==2)H_prot_miss_azimuth.setTitle("Proton vs #rho");
		H_prot_miss_azimuth.setTitleX("proton #phi (^o)");
		if(CHAN==1)H_prot_miss_azimuth.setTitleY("#phi #phi (^o)");
		if(CHAN==2)H_prot_miss_azimuth.setTitleY("#rho #phi (^o)");
		H_IMphi_MMep = new H2F("H_IMphi_MMep","H_IMphi_MMep",100,0,3,100,0,3);
		if(CHAN==1)H_IMphi_MMep.setTitle("IM K^+K^- vs MM ep");
		if(CHAN==2)H_IMphi_MMep.setTitle("IM #pi^+#pi^- vs MM ep");
		H_IMphi_MMep.setTitleX("MM ep (GeV)");
		if(CHAN==1)H_IMphi_MMep.setTitleY("IM K^+K^- (GeV)");
		if(CHAN==2)H_IMphi_MMep.setTitleY("IM #pi^+#pi^- (GeV)");
		H_MMepkp_MMepkm = new H2F("H_MMepkp_MMepkm","H_MMepkp_MMepkm",100,-2,5,100,-2,5);
		if(CHAN==1)H_MMepkp_MMepkm.setTitle("MM epK^+ vs MM epK^-");
		if(CHAN==1)H_MMepkp_MMepkm.setTitleX("MM epK^- (GeV)");
		if(CHAN==1)H_MMepkp_MMepkm.setTitleY("MM epK^+ (GeV)");
		if(CHAN==2)H_MMepkp_MMepkm.setTitle("MM ep#pi^+ vs MM ep#pi^-");
		if(CHAN==2)H_MMepkp_MMepkm.setTitleX("MM ep#pi^- (GeV)");
		if(CHAN==2)H_MMepkp_MMepkm.setTitleY("MM ep#pi^+ (GeV)");
		H_IMphi_MMepphi = new H2F("H_IMphi_MMepphi","H_IMphi_MMepphi",100,-2,4,100,0,3);
		if(CHAN==1)H_IMphi_MMepphi.setTitle("IM K^+K^- vs MM epK^+K^-");
		if(CHAN==2)H_IMphi_MMepphi.setTitle("IM #pi^+#pi^- vs MM ep#pi^+#pi^-");
		if(CHAN==1)H_IMphi_MMepphi.setTitleX("MM epK^+K^- (GeV)");
		if(CHAN==2)H_IMphi_MMepphi.setTitleX("MM ep#pi^+#pi^- (GeV)");
		if(CHAN==1)H_IMphi_MMepphi.setTitleY("IM K^+K^- (GeV)");
		if(CHAN==2)H_IMphi_MMepphi.setTitleY("IM #pi^+#pi^- (GeV)");
		H_angleKp_angle_Km = new H2F("H_angleKp_angle_Km","H_angleKp_angle_Km",100,0,180,100,0,180);
		if(CHAN==1)H_angleKp_angle_Km.setTitle("Colinearity K^+ vs K^-");
		if(CHAN==2)H_angleKp_angle_Km.setTitle("Colinearity #pi^+ vs #pi^-");
		if(CHAN==1)H_angleKp_angle_Km.setTitleX("colinear K^- (^o)");
		if(CHAN==2)H_angleKp_angle_Km.setTitleX("colinear #pi^- (^o)");
		if(CHAN==1)H_angleKp_angle_Km.setTitleY("colinear K^+ (^o)");
		if(CHAN==2)H_angleKp_angle_Km.setTitleY("colinear #pi^+ (^o)");
		H_angleKp_angleP = new H2F("H_angleKp_angleP","H_angleKp_angleP",100,0,180,100,0,180);
		if(CHAN==1)H_angleKp_angleP.setTitle("Colinearity K^+ vs prot");
		if(CHAN==2)H_angleKp_angleP.setTitle("Colinearity #pi^+ vs prot");
		H_angleKp_angleP.setTitleX("colinear prot (^o) (^o)");
		if(CHAN==1)H_angleKp_angleP.setTitleX("colinear #pi^+ (^o)");
		if(CHAN==2)H_angleKp_angleP.setTitleX("colinear #pi^+ (^o)");
		H_KKopen_phiE = new H2F("H_KKopen_phiE","H_KKopen_phiE",100,0,EB,100,0,180);
		if(CHAN==1)H_KKopen_phiE.setTitle("K^+K^- open angle vs #phi E");
		if(CHAN==2)H_KKopen_phiE.setTitle("#pi^+#pi^- open angle vs #rho E");
		if(CHAN==1)H_KKopen_phiE.setTitleX("#phi E (GeV)");
		if(CHAN==2)H_KKopen_phiE.setTitleX("#rho E (GeV)");
		if(CHAN==1)H_KKopen_phiE.setTitleY("K^+K^- open angle (^o)");
		if(CHAN==2)H_KKopen_phiE.setTitleY("#pi^+#pi^- open angle (^o)");
		
		H_IM_KmP = new H1F("H_IM_KmP","H_IM_KmP",100,1.25,2.75);
		H_IM_KmP.setTitle("invariant mass PK^-");
		H_IM_KmP.setTitleX("IM PK^- (GeV)");
		H_IM_KmP_MMeKp = new H2F("H_IM_KmP_MMeKp","H_IM_KmP_MMeKp",100,1.25,2.75,100,0,5);
		H_IM_KmP_MMeKp.setTitle("MM eK^+ vs IM PK^-");
		H_IM_KmP_MMeKp.setTitleX("IM PK^- (GeV)");
		H_IM_KmP_MMeKp.setTitleY("MM eK^+ (GeV)");
		H_IM_KmP_IMphi = new H2F("H_IM_KmP_IMphi","H_IM_KmP_IMphi",100,1.25,2.75,100,1.8,2.8);
		H_IM_KmP_IMphi.setTitle("IM K^+K^- vs vs IM PK^-");
		H_IM_KmP_IMphi.setTitleX("IM PK^- (GeV)");
		H_IM_KmP_IMphi.setTitleY("IM K^+K^- (GeV");
	}
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
                        if(recPart.getInt("pid",p)==11 ){
                                float px = recPart.getFloat("px", p); 
                                float py = recPart.getFloat("py", p); 
                                float pz = recPart.getFloat("pz", p); 
                                float vz = recPart.getFloat("vz", p); 
                                float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
                                float theta = (float)Math.toDegrees(Math.acos(pz/mom));
                                float phi = (float)Math.toDegrees(Math.atan2(py,px));
                                phi+=20;
                                if(phi<0)phi+=180;
                                if(mom>1.75 && theta>7 && Math.abs(vz)<20 && theta>17*(1-mom/7) ){
                                        e_sect = (int)phi/60;
                                        haz_elec=p;
                                }   
                        }   
                }   
        }   
        public void MakeProton(DataBank recPart){
		for(int p=0;p<recPart.rows();p++){
			if(recPart.getInt("pid",p)==2212 ){
				float px = recPart.getFloat("px", p); 
				float py = recPart.getFloat("py", p); 
				float pz = recPart.getFloat("pz", p); 
				float vz = recPart.getFloat("vz", p); 
				prot_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
				prot_theta = (float)Math.toDegrees(Math.acos(pz/prot_mom));
				//if(prot_mom>0.4 && prot_mom<4 && prot_theta>15 && prot_theta<35 && Math.abs(vz)<20){}
				if(prot_mom>0.4 && prot_mom<4 && prot_theta>15 && prot_theta<75 && Math.abs(vz)<20){
					haz_prot=p;
				}   
                        }   
                }
        }
	public void MakeKP(DataBank recPart){
		for(int p=0;p<recPart.rows();p++){
			//if( recPart.getInt("pid",p)==321 || recPart.getInt("pid",p)==211 ){}
			if( ( CHAN==1 && recPart.getInt("pid",p)==321) ||( CHAN==2 && recPart.getInt("pid",p)==211) ){
				float px = recPart.getFloat("px", p); 
				float py = recPart.getFloat("py", p); 
				float pz = recPart.getFloat("pz", p); 
				float vz = recPart.getFloat("vz", p); 
				kp_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
				kp_theta = (float)Math.toDegrees(Math.acos(pz/kp_mom));
				//if(kp_mom>0.4 && kp_theta<75 && Math.abs(vz)<20){}
				if(kp_mom>1.1 && kp_theta<35 && Math.abs(vz)<20){
					haz_kp=p;
				}
			}
		}
	}
	public void MakeKM(DataBank recPart){
		for(int p=0;p<recPart.rows();p++){
			//if( recPart.getInt("pid",p)==-321 || recPart.getInt("pid",p)==-211 ){}
			if( ( CHAN==1 && recPart.getInt("pid",p)==-321) ||( CHAN==2 && recPart.getInt("pid",p)==-211) ){
				float px = recPart.getFloat("px", p); 
				float py = recPart.getFloat("py", p); 
				float pz = recPart.getFloat("pz", p); 
				float vz = recPart.getFloat("vz", p); 
				km_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
				km_theta = (float)Math.toDegrees(Math.acos(pz/km_mom));
				//if(km_mom>0.4 && km_theta<75 && Math.abs(vz)<20){}
				if(km_mom>1.1 && km_theta<35 && Math.abs(vz)<20){
					haz_km=p;
				}
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
                VmissPHI = new LorentzVector(0,0,0,0);
                VmissPHI.add(VB);
                VmissPHI.add(VT);
                VmissPHI.sub(VE);
                VmissPHI.sub(VPROT);
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

		p = haz_kp;
                px = recPart.getFloat("px", p); 
                py = recPart.getFloat("py", p); 
                pz = recPart.getFloat("pz", p); 
                if(CHAN==1)E = (float)Math.sqrt(px*px+py*py+pz*pz + 0.4937*0.4937);
                if(CHAN==2)E = (float)Math.sqrt(px*px+py*py+pz*pz + 0.1396*0.1396);
		VKP = new LorentzVector(px,py,pz,E);
		kp_beta = recPart.getFloat("beta", p); 
                kp_vz = recPart.getFloat("vz", p); 

		p = haz_km;
		px = recPart.getFloat("px", p); 
		py = recPart.getFloat("py", p); 
		pz = recPart.getFloat("pz", p); 
		if(CHAN==1)E = (float)Math.sqrt(px*px+py*py+pz*pz + 0.4937*0.4937);
                if(CHAN==2)E = (float)Math.sqrt(px*px+py*py+pz*pz + 0.1396*0.1396);
		VKM = new LorentzVector(px,py,pz,E);
		km_beta = recPart.getFloat("beta", p); 
                km_vz = recPart.getFloat("vz", p); 

		VmissP = new LorentzVector(0,0,0,0);
		VmissP.add(VB);
		VmissP.add(VT);
		VmissP.sub(VE);
		VmissP.sub(VKP);
		VmissP.sub(VKM);

		VMISS = new LorentzVector(0,0,0,0);
		VMISS.add(VB);
		VMISS.add(VT);
		VMISS.sub(VE);
		VMISS.sub(VPROT);
		VMISS.sub(VKP);
		VMISS.sub(VKM);
		VPHI = new LorentzVector(0,0,0,0);
		VPHI.add(VKP);
		VPHI.add(VKM);
		Vhad2 = (VGS.vect()).cross(VPHI.vect());
		VmissKP = new LorentzVector(0,0,0,0);
	       	VmissKP.add(VB);
	       	VmissKP.add(VT);
	       	VmissKP.sub(VE);
	       	VmissKP.sub(VPROT);
	       	VmissKP.sub(VKM);
		VmissKM = new LorentzVector(0,0,0,0);
		VmissKM.add(VB);
		VmissKM.add(VT);
		VmissKM.sub(VE);
		VmissKM.sub(VPROT);
		VmissKM.sub(VKP);

		//CM vectors
                //Vector3 CMphi = new Vector3(VPHI.boostVector());
                //CMphi.negative();
                //LorentzVector CMKp = new LorentzVector(pKp);
                //LorentzVector CMKm = new LorentzVector(pKm);
                //CMKp.boost(CMphi);
                //CMKm.boost(CMphi);
                // Math.cos( CMKp.theta() ) ;
		
		VLambda = new LorentzVector(VPROT);
	       	VLambda.add(VKM);
		VmissLambda = new LorentzVector(VT);
		VmissLambda.add(VGS);
		VmissLambda.sub(VKP);
	}
	public boolean KineCut(){
		boolean res = false;
		float MissMesCut = 0.75f;
		if(CHAN==1)MissMesCut += 1.5;
		float MissMassCut = 0.3f;
		if(CHAN==1)MissMassCut += 0.3;
		float backback = (float)Math.toDegrees(VPROT.phi() - VmissP.phi());
		while(backback<-180)backback += 360;
		while(backback>180)backback -= 360;
		coplanarity = (float)Vangle(Vhadr,Vhad2);
		while(coplanarity<-180)coplanarity += 360;
		while(coplanarity>180)coplanarity -= 360;
		if(true
		  && (VKM.p()<2.5 && VKP.p()<2.5)
		  //&& VmissKP.mass()+VmissKM.mass() > -0.8
		  //&& VmissKP.mass()+VmissKM.mass() <  1.6
		  //&& VPHI.mass() < 2.0
		  && VPHI.mass() < 3.0
		  && VmissPHI.mass() < 3.0
		  //&& coplanarity< 90 
		  //&& backback < 90
		  && VMISS.e() < 2.5 && Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py())<0.5
		  && Math.abs(VmissKP.mass2()) < MissMesCut
		  && Math.abs(VmissKM.mass2()) < MissMesCut
		  && Math.abs(VMISS.mass2()) < MissMassCut
		  && Vangle(VPROT.vect(),VmissP.vect()) < 40
		  && Vangle(VKM.vect(),VmissKM.vect()) < 30
		  && Vangle(VKP.vect(),VmissKP.vect()) < 30
		)res = true;
		return res;
	}
	public void FillNoCutsHists(){
		H_phi_nocut_IM.fill(VPHI.mass());
	}
        public void FillHists(){
                H_elec_theta_mom.fill(  VE.p(),Math.toDegrees(VE.theta()));
                H_prot_theta_mom.fill(VPROT.p(),Math.toDegrees(VPROT.theta()));
                H_kp_theta_mom.fill( VKP.p(),Math.toDegrees(VKP.theta()));
                H_km_theta_mom.fill( VKM.p(),Math.toDegrees(VKM.theta()));

                H_elec_theta_phi.fill(Math.toDegrees(VE.phi())  ,Math.toDegrees(VE.theta()));
                H_prot_theta_phi.fill(Math.toDegrees(VPROT.phi()),Math.toDegrees(VPROT.theta()));
                H_kp_theta_phi.fill(Math.toDegrees(VKP.phi()) ,Math.toDegrees(VKP.theta()));
                H_km_theta_phi.fill(Math.toDegrees(VKM.phi()) ,Math.toDegrees(VKM.theta()));

                H_elec_phi_mom.fill(  VE.p(),Math.toDegrees(VE.phi()));
                H_prot_phi_mom.fill(VPROT.p(),Math.toDegrees(VPROT.phi()));
                H_kp_phi_mom.fill( VKP.p(),Math.toDegrees(VKP.phi()));
                H_km_phi_mom.fill( VKM.p(),Math.toDegrees(VKM.phi()));

		H_elec_vz_theta.fill( Math.toDegrees(VE.theta()) , elec_vz);
		H_prot_vz_theta.fill( Math.toDegrees(VPROT.theta()) , prot_vz);
		H_kp_vz_theta.fill( Math.toDegrees(VKP.theta()) , kp_vz);
		H_km_vz_theta.fill( Math.toDegrees(VKM.theta()) , km_vz);

                H_prot_beta_mom.fill(VPROT.p(),prot_beta);
                H_kp_beta_mom.fill( VKP.p(),kp_beta);
                H_km_beta_mom.fill( VKM.p(),km_beta);

                float TrentoAng = (float)Vangle(Vlept,Vhadr);
                if((VPROT.vect()).dot(Vlept)<0)TrentoAng=-TrentoAng;

		//should only mix events that pass exclusivity cuts
		if(doMix>NOLD-1)for(int io=0;io<NOLD;io++){
			LorentzVector dumKP = new LorentzVector(VOLDKP[io]);
			LorentzVector dumKM = new LorentzVector(VOLDKM[io]);
			LorentzVector dumGS = new LorentzVector(VOLDGS[io]);
			double polarComp = VGS.phi() - dumGS.phi();
			dumKP.rotateZ(polarComp);
			dumKM.rotateZ(polarComp);
			dumGS.rotateZ(polarComp);
			LorentzVector boost4Vect = new LorentzVector(VOLDGS[io]);
			boost4Vect.add(VT);
			Vector3 boostVect = boost4Vect.boostVector();
			VOLDKM[io].boost(boostVect);
			VOLDKP[io].boost(boostVect);
			boost4Vect = new LorentzVector(VGS);
			boost4Vect.add(VT);
			boostVect = boost4Vect.boostVector();
			boostVect.negative();
			VOLDKM[io].boost(boostVect);
			VOLDKP[io].boost(boostVect);
			//done mixing
			LorentzVector VMIX1PHI = new LorentzVector(0,0,0,0);
			VMIX1PHI.add(VKP);
			VMIX1PHI.add(dumKM);
			LorentzVector VMIX2PHI = new LorentzVector(0,0,0,0);
			VMIX2PHI.add(dumKP);
			VMIX2PHI.add(VKM);
			H_phi_IM_back.fill(VMIX1PHI.mass());
			H_phi_IM_back.fill(VMIX2PHI.mass());
		}
		VOLDGS[oldInd] = VGS;
		VOLDKP[oldInd] = VKP;
		VOLDKM[oldInd] = VKM;
		doMix++;
		oldInd = doMix%NOLD;
		//
		H_phi_IM.fill(VPHI.mass());
		H_phi_copl.fill(coplanarity);
		
		H_Q2_xB.fill(-VGS.mass2()/(2*0.938*VGS.e()),-VGS.mass2());
		H_t_phi.fill(TrentoAng,-Vmand.mass2());
		H_prot_miss_azimuth.fill(Math.toDegrees(VPROT.phi()),Math.toDegrees(VmissP.phi()));
		H_IMphi_MMep.fill(VmissPHI.mass(),VPHI.mass());
		H_MMepkp_MMepkm.fill(VmissKP.mass2(),VmissKM.mass2()); 
		H_misspT_missE.fill(VMISS.e() , Math.sqrt(VMISS.px()*VMISS.px()+VMISS.py()*VMISS.py()) );

		H_IMphi_MMepphi.fill(VMISS.mass2() , VPHI.mass() );
		H_angleKp_angle_Km.fill( Vangle(VKM.vect(),VmissKM.vect()) , Vangle(VKP.vect(),VmissKP.vect()) );
		H_angleKp_angleP.fill( Vangle(VPROT.vect(),VmissP.vect()) , Vangle(VKP.vect(),VmissKP.vect()) );
		H_KKopen_phiE.fill( VPHI.e() , Vangle(VKP.vect(),VKM.vect()) );


		      H_IM_KmP.fill(VLambda.mass() );
		H_IM_KmP_MMeKp.fill(VLambda.mass() , VmissLambda.mass());
		H_IM_KmP_IMphi.fill(VLambda.mass() , VPHI.mass());

		if(VPHI.mass()<1.15){
			H_Q2_xB_cut.fill(-VGS.mass2()/(2*0.938*VGS.e()),-VGS.mass2());
			H_t_phi_cut.fill(TrentoAng,-Vmand.mass2());
		}
	}
	public void processEvent(DataEvent event) {
		haz_elec=-1;haz_prot=-1;haz_kp=-1;haz_km=-1;
                if(event.hasBank("REC::Particle"))HasElectron(event.getBank("REC::Particle"));
                if(haz_elec>-1)MakeProton(event.getBank("REC::Particle"));
                if(haz_elec>-1 && haz_prot>-1 )MakeKP(event.getBank("REC::Particle"));
		if(haz_elec>-1 && haz_prot>-1 && haz_kp>-1 )MakeKM(event.getBank("REC::Particle"));
		if(haz_elec>-1 && haz_prot>-1 && haz_kp>-1 && haz_km>-1){
			boolean IsPion = false;
			boolean IsKaon = false;
			CHAN=2;
			MakeParticles(event.getBank("REC::Particle"));
			IsPion = KineCut();
			CHAN=1;
			MakeParticles(event.getBank("REC::Particle"));
			IsKaon = KineCut();
			FillNoCutsHists();
			if( !IsPion && IsKaon ){
				FillHists();
				Nexclu++;
			}
		}
	}
	public void plot(){
                EmbeddedCanvas can_phi_parts  = new EmbeddedCanvas();
                can_phi_parts.setSize(2000,2500);
                can_phi_parts.divide(4,5);
                can_phi_parts.setAxisTitleSize(18);
                can_phi_parts.setAxisFontSize(18);
                can_phi_parts.setTitleSize(18);
                can_phi_parts.cd(0);can_phi_parts.draw(H_elec_theta_mom);
                can_phi_parts.cd(1);can_phi_parts.draw(H_prot_theta_mom);
                can_phi_parts.cd(2);can_phi_parts.draw(H_kp_theta_mom);
                can_phi_parts.cd(3);can_phi_parts.draw(H_km_theta_mom);
    
                can_phi_parts.cd(4);can_phi_parts.draw(H_elec_theta_phi);
                can_phi_parts.cd(5);can_phi_parts.draw(H_prot_theta_phi);
                can_phi_parts.cd(6);can_phi_parts.draw(H_kp_theta_phi);
                can_phi_parts.cd(7);can_phi_parts.draw(H_km_theta_phi);
    
                can_phi_parts.cd(8);can_phi_parts.draw(H_elec_phi_mom);
                can_phi_parts.cd(9);can_phi_parts.draw(H_prot_phi_mom);
                can_phi_parts.cd(10);can_phi_parts.draw(H_kp_phi_mom);
                can_phi_parts.cd(11);can_phi_parts.draw(H_km_phi_mom);
    
                can_phi_parts.cd(12);can_phi_parts.draw(H_elec_vz_theta);
                can_phi_parts.cd(13);can_phi_parts.draw(H_prot_vz_theta);
                can_phi_parts.cd(14);can_phi_parts.draw(H_kp_vz_theta);
                can_phi_parts.cd(15);can_phi_parts.draw(H_km_vz_theta);

		can_phi_parts.cd(16);can_phi_parts.draw(H_prot_beta_mom);
		can_phi_parts.cd(17);can_phi_parts.draw(H_kp_beta_mom);
		can_phi_parts.cd(18);can_phi_parts.draw(H_km_beta_mom);
    
                can_phi_parts.save(String.format("plots"+runNum+"/phi_parts.png"));
                System.out.println(String.format("saved plots"+runNum+"/phi_parts.png"));

		EmbeddedCanvas can_phi_exclu = new EmbeddedCanvas();
		can_phi_exclu.setSize(2000,1500);
		can_phi_exclu.divide(4,3);
		can_phi_exclu.setAxisTitleSize(18);
		can_phi_exclu.setAxisFontSize(18);
		can_phi_exclu.setTitleSize(18);
		int ipad=0;
		//L1
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_Q2_xB);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_t_phi);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_phi_IM);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_phi_copl);ipad++;
		//L2
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_prot_miss_azimuth);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_IMphi_MMep);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_MMepkp_MMepkm);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_misspT_missE);ipad++;
		//L3
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_IMphi_MMepphi);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_angleKp_angle_Km);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_angleKp_angleP);ipad++;
		can_phi_exclu.cd(ipad);can_phi_exclu.draw(H_KKopen_phiE);ipad++;
		can_phi_exclu.save(String.format("plots"+runNum+"/phi_exclu.png"));
		System.out.println(String.format("saved plots"+runNum+"/phi_exclu.png"));
		
		EmbeddedCanvas can_phi_other = new EmbeddedCanvas();
		can_phi_other.setSize(2000,1000);
		can_phi_other.divide(4,2);
		can_phi_other.setAxisTitleSize(18);
		can_phi_other.setAxisFontSize(18);
		can_phi_other.setTitleSize(18);
		ipad=0;
		can_phi_other.cd(ipad);can_phi_other.draw(H_Q2_xB_cut);ipad++;
		can_phi_other.cd(ipad);can_phi_other.draw(H_t_phi_cut);ipad++;
		can_phi_other.cd(ipad);can_phi_other.draw(H_phi_nocut_IM);ipad++;
		can_phi_other.cd(ipad);can_phi_other.draw(H_phi_IM_back);ipad++;
		can_phi_other.cd(ipad);can_phi_other.draw(H_IM_KmP);ipad++;
		can_phi_other.cd(ipad);can_phi_other.draw(H_IM_KmP_MMeKp);ipad++;
		can_phi_other.cd(ipad);can_phi_other.draw(H_IM_KmP_IMphi);
		can_phi_other.save(String.format("plots"+runNum+"/phi_other.png"));
		System.out.println(String.format("saved plots"+runNum+"/phi_other.png"));

	}
	public void write(){
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/particles/");
                dirout.cd("/particles/");
                dirout.addDataSet(H_elec_theta_mom,H_elec_theta_phi,H_elec_phi_mom,H_elec_vz_theta);
		dirout.addDataSet(H_prot_theta_mom,H_prot_theta_phi,H_prot_phi_mom,H_prot_vz_theta,H_prot_beta_mom);
		dirout.addDataSet(H_kp_theta_mom,H_kp_theta_phi,H_kp_phi_mom,H_kp_vz_theta,H_kp_beta_mom);
		dirout.addDataSet(H_km_theta_mom,H_km_theta_phi,H_km_phi_mom,H_km_vz_theta,H_km_beta_mom);
                dirout.mkdir("/exclu/");
                dirout.cd("/exclu/");
		dirout.addDataSet(H_phi_nocut_IM);
                dirout.addDataSet(H_Q2_xB,H_t_phi,H_phi_IM,H_phi_copl,H_phi_IM_back);
		dirout.addDataSet(H_prot_miss_azimuth,H_IMphi_MMep,H_MMepkp_MMepkm,H_misspT_missE,H_IMphi_MMepphi,H_angleKp_angle_Km,H_angleKp_angleP,H_KKopen_phiE);
		dirout.addDataSet(H_Q2_xB_cut,H_t_phi_cut);
                if(runNum>0)dirout.writeFile("plots"+runNum+"/phi_"+runNum+".hipo");
                else dirout.writeFile("plots/phi.hipo");

	}
	void testBoost(){
		H2F H_pi0_open_E = new H2F("H_pi0_open_E","H_pi0_open_E",200,0.5,6.5,200,0,25);
		for(int e=0;e<100000;e++){
			//labe pi0
			double theta = Math.random() * 3.1415927f;
			double phi = 2 * Math.random() * 3.1415927f;
			double E = 1 + Math.random() * 5;
			double p = Math.sqrt(E*E-0.135f*0.135f);
			double px = p * Math.cos(Math.toDegrees(phi)) * Math.sin(Math.toDegrees(theta));
			double py = p * Math.sin(Math.toDegrees(phi)) * Math.sin(Math.toDegrees(theta));
			double pz = p * Math.cos(Math.toDegrees(theta));
			LorentzVector VPI0 = new LorentzVector(px,py,pz,E);
			Vector3 boostVect = VPI0.boostVector();
			boostVect.negative();
			//cm pi0
			theta = Math.random() * 3.1415927f;
			phi = 2 * Math.random() * 3.1415927f;
			p = 0.135f/2f;
			px = p * Math.cos(Math.toDegrees(phi)) * Math.sin(Math.toDegrees(theta));
			py = p * Math.sin(Math.toDegrees(phi)) * Math.sin(Math.toDegrees(theta));
			pz = p * Math.cos(Math.toDegrees(theta));
			LorentzVector VG1 = new LorentzVector(px,py,pz,p);
			LorentzVector VG2 = new LorentzVector(-px,-py,-pz,p);
			VG1.boost(boostVect);
			VG2.boost(boostVect);
			H_pi0_open_E.fill(VPI0.e() , Vangle(VG1.vect(),VG2.vect()) );
		}
                EmbeddedCanvas can_pi0  = new EmbeddedCanvas();
                can_pi0.setSize(800,800);
		can_pi0.draw(H_pi0_open_E);
		F1D funMin = new F1D("funMin","2*asin(0.135/x)*180/3.1415927",1.5,6.5);
		funMin.setLineWidth(3);
		can_pi0.draw(funMin,"same");
		can_pi0.save("pi0.png");

	}
////////////////////////////////////////////////
        public static void main(String[] args) {
                System.setProperty("java.awt.headless", "true");
                GStyle.setPalette("kRainBow");
                int count = 0;
                int runNum = 0;
                float EBreq = 7;
                int maxevents = 500000000;
                String filelist = "list_of_files.txt";
                if(args.length>0)runNum = Integer.parseInt(args[0]);
                if(args.length>1)filelist = args[1];
                if(args.length>2)EBreq = Float.parseFloat(args[2]);
                if(args.length>3)maxevents = Integer.parseInt(args[3]);
                phi_ana ana = new phi_ana(runNum,EBreq);
		System.out.println("testing boost");
		ana.testBoost();
		System.out.println("Done");
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
				if(count%100000 == 0){
					String message = "> " + count/1000 + "k events (this is phi analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot + " ; Nphi=" + ana.Nexclu;
					System.out.println(message);
				}
                                //if(count%10000 == 0) System.out.println(count/1000 + "k events (this is phi analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot);
                        }   
                        reader.close();
                }   
                System.out.println("Total : " + count + " events ; N phi = " + ana.Nexclu );
                ana.plot();
		ana.write();
        }   
}

