import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.GraphErrors
// import ROOTFitter

def grtl = (1..6).collect{
  def gr = new GraphErrors('elec sec'+it)
  gr.setTitle("vz")
  gr.setTitleY("vz sigma (cm)")
  gr.setTitleX("run number")
  return gr
}

def grtl2 = (1..6).collect{
  def gr = new GraphErrors('neg sec'+it)
  gr.setTitle("vz")
  gr.setTitleY("vz sigma (cm)")
  gr.setTitleX("run number")
  return gr
}


TDirectory out = new TDirectory()

for(arg in args) {
  TDirectory dir = new TDirectory()
  dir.readFile(arg)

  def name = arg.split('/')[-1]
  def m = name =~ /\d\d\d\d/
  def run = m[0].toInteger()

  (0..<6).each{
    def h1 = dir.getObject('/elec/H_elec_vz_S'+(it+1))
    def h2 = dir.getObject('/elec/H_neg_vz_S'+(it+1))

    grtl[it].addPoint(run, h1.getRMS(), 0, 0)
    grtl2[it].addPoint(run, h2.getRMS(), 0, 0)

    out.mkdir('/'+run)
    out.cd('/'+run)
    out.addDataSet(h1)
    out.addDataSet(h2)
  }
}

out.mkdir('/timelines')
out.cd('/timelines')
grtl.each{ out.addDataSet(it) }
grtl2.each{ out.addDataSet(it) }

out.writeFile('vz_sigma.hipo')