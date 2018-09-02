from shelxfile.shelx import ShelXFile
import numpy as np

def grow(shx: 'ShelXFile'):
  """
  Completes covalent bound molecule fragments.
  """
  showatoms = []
  for at in range(shx.atoms):
    showatoms.append(at)
  showbonds = []
  brauchSymm = []
  bs = ""
  prime = np.array(3.1)
  dp = np.array(3.1)
  D = np.array(3.1)
  floorD = np.array(3.1)
  dk, dddd = 0.0, 0.0
  for (int k =0; k<sdm.size();k++){
    if ((sdm.at(k).covalent)||(growQPeak)){
      if (!(growQPeak)&&(asymm[sdm.at(k).a1].molindex<1)) continue;
      for (int n=0;n<cell.symmops.size();  n++){
        if (((asymm[sdm.at(k).a1].part!=0)&&(asymm[sdm.at(k).a2].part!=0)&&(asymm[sdm.at(k).a1].part!=asymm[sdm.at(k).a2].part)))continue;
        if ((asymm[sdm.at(k).a1].an==asymm[sdm.at(k).a2].an)&&(asymm[sdm.at(k).a1].an==0)) continue;
        prime=cell.symmops.at(n) * asymm[sdm.at(k).a1].frac + cell.trans.at(n);
        D=prime - asymm[sdm.at(k).a2].frac+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        dp=D - floorD - V3(0.5,0.5,0.5);
        if ((n==0)&&(V3(0,0,0)==floorD)) continue;
        dk=fl(dp.x,dp.y,dp.z);
        //printf ("%f n%d\n",dk,n);
        dddd=(sdm.at(k).covalent)?(sdm.at(k).d+0.2):0;
        // an<0 means hydrogen atom (this is a hydrogen contact):
        if ((asymm[sdm.at(k).a1].an<0)&&(asymm[sdm.at(k).a2].an<0)) dddd=1.8;//||
        if ( (dk>0.001)&&(dddd>=dk)) {
          bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-(int)floorD.x).arg(5-(int)floorD.y).arg(5-(int)floorD.z).arg(asymm[sdm.at(k).a1].molindex);
          //	     printf("%s %s %g %g\n",asymm[sdm.at(k).a1].Label.toStdString().c_str(),asymm[sdm.at(k).a2].Label.toStdString().c_str(),sdm.at(k).d,dddd);
          //             printf("(grow)%s\n",bs.toStdString().c_str());
          if  ((!brauchSymm.contains(bs))) {
            brauchSymm.append(bs);
          }
        }
      }
    }
  }