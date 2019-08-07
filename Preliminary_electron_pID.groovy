import java.io.*;
import java.util.*;
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.GraphErrors
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.graphics.EmbeddedCanvas;

def run = args[0].toInteger()

TDirectory out = new TDirectory()
out.mkdir('/'+run)
out.cd('/'+run)

float EB = 10.6f
if(run>6607) EB=10.2f
// H1F H_elec_mom[6];
// H2F H_elec_vz_mom[6], H_elec_theta_mom[6], H_elec_theta_phi[6];
H_elec_mom =(0..5).collect{
	def h1 = new H1F("H_elec_mom_S"+(it+1), "H_elec_mom_S"+(it+1),100, 0, EB);
	return h1
}

H_elec_vz_mom =(0..<6).collect{
	def h1 = new H2F("H_elec_vz_mom_S"+(it+1), "H_elec_vz_mom_S"+(it+1),100,0,EB,100,-25,25);
	return h1
}

H_elec_theta_mom =(0..<6).collect{
	def h1 = new H2F("H_elec_theta_mom_S"+(it+1), "H_elec_theta_mom_S"+(it+1),100,0,EB,100,0,40);
	return h1
}

H_elec_phi_mom = (0..<6).collect{
	def h1 = new H2F("H_elec_phi_mom_S"+(it+1), "H_elec_phi_mom_S"+(it+1),100,0,EB,100,-180,180);
	return h1
}

H_elec_theta_phi =(0..<6).collect{
	def h1 = new H2F("H_elec_theta_phi_S"+(it+1), "H_elec_theta_phi_S"+(it+1),100,-180,180,100,0,40);
	return h1
}

H_elec_HTCC_nphe =(0..<6).collect{
	def h1 = new H1F("H_elec_HTCC_nphe_S"+(it+1), "H_elec_HTCC_nphe_S"+(it+1),100,0,100);
	return h1
}

H_elec_EC_Sampl =(0..<6).collect{
	def h1 = new H1F("H_elec_EC_Sampl_S"+(it+1), "H_elec_EC_Sampl_S"+(it+1),100,0,1);
	return h1
}

//H_elec_mom=new H1F("H_elec_mom", "H_elec_mom_S",100, 0, EB);


int e_index, e_sect, e_nphe
float e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_ecal_E, e_Sampl_frac
LorentzVector Ve = new LorentzVector()

filenum=-1
for (arg in args){
	filenum=filenum+1
	if (filenum==0) continue
	HipoDataSource reader = new HipoDataSource();
	reader.open(arg);

	while( reader.hasEvent()){
		DataEvent event = reader.getNextEvent();
		processEvent(event)
	}
}

(0..<6).each{
	out.addDataSet(H_elec_mom[it])
	out.addDataSet(H_elec_vz_mom[it])
	out.addDataSet(H_elec_theta_mom[it])
	out.addDataSet(H_elec_phi_mom[it])
	out.addDataSet(H_elec_theta_phi[it])
	out.addDataSet(H_elec_HTCC_nphe[it])
	out.addDataSet(H_elec_EC_Sampl[it])
}

//out.addDataSet(H_elec_mom)
out.writeFile('electron_pID_'+run+'.hipo')

public void processEvent(DataEvent event) {
	if(!event.hasBank("REC::Particle")) return 
	if(!event.hasBank("REC::Calorimeter")) return 
	if(!event.hasBank("REC::Cherenkov")) return
	DataBank partBank = event.getBank("REC::Particle");
	e_index=-1
	if (!hasElectron(partBank, event)) return
	if (e_index>-1){
		makeElectron(partBank)
		fillHists()
		//H_elec_mom.fill(e_mom)
	}
	else return;
}


public void fillHists(){
	H_elec_mom[e_sect-1].fill(e_mom);
	H_elec_vz_mom[e_sect-1].fill(e_mom,e_vz);
	H_elec_theta_mom[e_sect-1].fill(e_mom,e_theta)
	H_elec_phi_mom[e_sect-1].fill(e_mom,e_phi)
	H_elec_theta_phi[e_sect-1].fill(e_phi,e_theta);
	H_elec_HTCC_nphi[e_sect-1].fill(e_nphe)
	H_elec_EC_Sampl[e_sect-1].fill(e_Sampl_frac)
}

public void makeElectron(DataBank recPart, DataEvent cur_event){
		int ei=e_index
		float px = recPart.getFloat("px",ei)
		float py = recPart.getFloat("py",ei)
		float pz = recPart.getFloat("pz",ei)
		e_mom = (float)Math.sqrt(px*px+py*py+pz*pz)
		e_vz = recPart.getFloat("vz",ei)
		e_vx = recPart.getFloat("vx",ei)
		e_vy = recPart.getFloat("vy",ei)
		e_phi = (float)Math.toDegrees(Math.atan2(py,px))
		if(e_phi<0) e_phi+=360;
		e_phi=360-e_phi;
		e_phi=e_phi-140;
		if (e_phi<0) e_phi+=360;
		e_sect = (int) Math.ceil(e_phi/60);
		Ve = new LorentzVector(px,py,pz,e_mom)
		e_phi = (float) Math.toDegrees(Ve.phi())
		e_theta = (float) Math.toDegrees(Ve.theta())
		DataBank HTCCbank = cur_event.getBank("REC::Cherenkov");
		for(int l = 0; l < HTCCbank.rows(); l++){
			if(HTCCbank.getShort("pindex",l)==ei && HTCCbank.getInt("detector",l)==15){
				//HTCCnphe = HTCCbank.getInt("nphe",l);
				e_nphe = (int) HTCCbank.getFloat("nphe",l);
			}
		}
		DataBank ECALbank = cur_event.getBank("REC::Calorimeter")
		for(int l = 0; l < ECALbank.rows(); l++){
			if(ECALbank.getShort("pindex",l)==p){
				e_ecal_E += ECALbank.getFloat("energy",l);
			}
		}
		float e_Sampl_frac = e_ecal_E/mom

}

public boolean hasElectron(DataBank recPart, DataEvent cur_event){
	boolean found = false
	for(int p=0;p<recPart.rows();p++){
		if (isElectron(recPart, cur_event, p)){
			if (found) System.out.println ("Error, two or more electrons found!")
			found=true
		}
	}
	return found
}

public boolean isElectron(DataBank recPart, DataEvent cur_event, int p){
	if (CC_nphe_cut(recPart,cur_event,p) && EC_sampling_fraction_cut(recPart,cur_event,p)&& ele_default_PID_cut(recPart,p)&& ele_charge_cut(recPart,p) && ele_kine_cut(recPart,p)&& inDC(recPart,p)){
		//System.out.println("Electron Found!")
		e_index=p
		return true
	}
	else return false
}

// pID function
// Cherenkov, Calo cut needed
// edep proportional to e_mom from DC
// ~30% ratio cut, Sampling fraction

public boolean EC_sampling_fraction_cut(DataBank recPart, DataEvent cur_event, int p){
	float px = recPart.getFloat("px", p);
	float py = recPart.getFloat("py", p);
	float pz = recPart.getFloat("pz", p);
	float vz = recPart.getFloat("vz", p);
	float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
	DataBank ECALbank = cur_event.getBank("REC::Calorimeter")
	float ecal_E=0
	for(int l = 0; l < ECALbank.rows(); l++){
		if(ECALbank.getShort("pindex",l)==p){
			ecal_E += ECALbank.getFloat("energy",l);
		}
	}
	float Sampl_frac = ecal_E/mom
	if (Sampl_frac>0.18) return true
	else return false
}

public boolean CC_nphe_cut(DataBank recPart, DataEvent cur_event, int p){
	DataBank HTCCbank = cur_event.getBank("REC::Cherenkov");
	HTCCnphe=-1
	for(int l = 0; l < HTCCbank.rows(); l++){
		if(HTCCbank.getShort("pindex",l)==p && HTCCbank.getInt("detector",l)==15){
			//HTCCnphe = HTCCbank.getInt("nphe",l);
			HTCCnphe = HTCCbank.getFloat("nphe",l);
		}
	}
	if (HTCCnphe>1) return true
	else return false
}
public boolean inDC(DataBank recPart, int p){
	int status = recPart.getShort("status", p);
	if (status<0) status = -status; // if negative, trigger electron
	boolean inDC = (status>=2000 && status<4000); // checks FD or CD
}


public boolean ele_default_PID_cut(DataBank recPart, int p){
  if(recPart.getInt("pid",p)==11) return true;
  else return false;
}

public boolean ele_charge_cut(DataBank recPart, int p){
  if(recPart.getInt("charge",p)==-1) return true;
  else return false;
}

public boolean ele_kine_cut(DataBank recPart, int p){
	float px = recPart.getFloat("px", p);
	float py = recPart.getFloat("py", p);
	float pz = recPart.getFloat("pz", p);
	float vz = recPart.getFloat("vz", p);
	float mom = (float)Math.sqrt(px*px+py*py+pz*pz);
	float theta = (float)Math.toDegrees(Math.acos(pz/mom));
	float phi = (float)Math.toDegrees(Math.atan2(py,px));
	if(mom>1.75 && theta>7 && Math.abs(vz)<15 && theta>17*(1-mom/7) ){
		return true;
	}
	else return false;
}
