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

HistoDef()

int e_part_ind, e_sect, e_nphe
float e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_ecal_E, e_Sampl_frac
LorentzVector Ve = new LorentzVector()

filenum=-1

for (arg in args){
	filenum++
	if (filenum==0) continue
	HipoDataSource reader = new HipoDataSource();
	reader.open(arg);
	while( reader.hasEvent()){
		DataEvent event = reader.getNextEvent();
		processEvent(event)
	}
}

(0..<6).each{
	out.addDataSet(H_elec_vz[it])
	out.addDataSet(H_elec_HTCC_nphe[it])
	out.addDataSet(H_elec_EC_Sampl[it])
	out.addDataSet(H_neg_vz[it])
	out.addDataSet(H_neg_HTCC_nphe[it])
	out.addDataSet(H_neg_EC_Sampl[it])
}

out.writeFile('electron_pID_'+run+'.hipo')

public void processEvent(DataEvent event) {
	if(!event.hasBank("REC::Particle")) return 
	if(!event.hasBank("REC::Calorimeter")) return 
    def evc = event.getBank("REC::Calorimeter")
    def secs = [evc.getShort('pindex')*.toInteger(), evc.getByte('sector')].transpose().collectEntries()
    def evp = event.getBank("REC::Particle")
    evp.getInt("pid").eachWithIndex{pid, ind ->
    	if(secs[ind]==null) return
	    int charge = evp.getInt("charge",ind)
	    vz = evp.getFloat("vz",ind)
	    px = evp.getFloat("px",ind)
	    py = evp.getFloat("py",ind)
	    pz = evp.getFloat("pz",ind)
	    mom = (float) Math.sqrt(px*px+py*py+pz*pz)
	    energy = 0
	    evc.getInt("pindex").eachWithIndex{ pindex, ind_c ->
	    	if (pindex==ind) energy+=evc.getFloat("energy",ind_c)
	    }
   	    sampl_frac = energy/mom
	    if (charge<0){
	    	H_neg_vz[secs[ind]-1].fill(vz) //vz
	    	if (pid==11) H_elec_vz[secs[ind]-1].fill(vz) //electron
	    	H_neg_EC_Sampl[secs[ind]-1].fill(sampl_frac) // sampling Fraction
	    	if (pid==11) H_elec_EC_Sampl[secs[ind]-1].fill(sampl_frac) //electron
	    	if(!event.hasBank("REC::Cherenkov")) return
	    	def evh = event.getBank("REC::Cherenkov")
	    	evh.getInt("pindex").eachWithIndex{pindex, ind_h ->
	    		if(evh.getInt("detector",ind_h)!=15) return
    			if(pindex==ind) H_neg_HTCC_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h))
    			if(pid==11 && pindex==ind) H_elec_HTCC_nphe[secs[ind]-1].fill(evh.getFloat("nphe",ind_h))
	    	}
    	}
    }
}

public void HistoDef(){
	H_elec_vz =(0..<6).collect{
		def h1 = new H1F("H_elec_vz_S"+(it+1), "H_elec_vz_S"+(it+1),100,-25,25);
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

	H_neg_vz =(0..<6).collect{
		def h1 = new H1F("H_neg_vz_S"+(it+1), "H_neg_vz_S"+(it+1),100,-25,25);
		return h1
	}

	H_neg_HTCC_nphe =(0..<6).collect{
		def h1 = new H1F("H_neg_HTCC_nphe_S"+(it+1), "H_neg_HTCC_nphe_S"+(it+1),100,0,100);
		return h1
	}

	H_neg_EC_Sampl =(0..<6).collect{
		def h1 = new H1F("H_neg_EC_Sampl_S"+(it+1), "H_neg_EC_Sampl_S"+(it+1),100,0,1);
		return h1
	}

}