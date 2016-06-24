#include <iostream>
#include <fstream>
#include <vector>
#include "Pythia8/Pythia.h"
using namespace std;
using namespace Pythia8;

double invariant(double,double,double,double,double,double,double,double,double,double);
void drcalc(double phi1, double phi2, double y1, double y2, double *pdphi, double *pdeta, double *pdR);

int main() {
    const char *path = "/Users/Yuichi/root/macros/mytext25.txt";
    const char *pathb = "/Users/Yuichi/root/macros/mytext25b.txt";
    const char *pathc = "/Users/Yuichi/root/macros/mytext25c.txt";
    const char *pathd = "/Users/Yuichi/root/macros/mytext25d.txt";
    
    ofstream myfile(path);
    ofstream myfileb(pathb);
    ofstream myfilec(pathc);
    ofstream myfiled(pathd);
    
    
  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("HiggsSM:ffbar2HZ = on");
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfMatch = 5 -5");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfMatch = 13 -13");
  pythia.readString("PhaseSpace:pTHatMin = 0.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();
    
    // Common parameters for the two jet finders.
    double etaMax   = 2.4;
    double radius   = 0.4;
    double pTjetMin = 20.;
    // Exclude neutrinos (and other invisible) from study.
    int    nSel     = 2;
    
    // Set up SlowJet jet finder, with anti-kT clustering
    // and pion mass assumed for non-photons..
    
    SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);
    
    //Histogram

    
    //Variables
    
    int id;
    int mom1,mom2;
    int mom1id,mom2id;
    int jetsize,jetsizeaf;
    double hm;
    double zm;
    double invm;
    double dpT,dphi,deta,dR;
    double px1,py1,pz1,e1,m1;
    double px2,py2,pz2,e2,m2;
    double px3,py3,pz3,e3,m3;
    double px4,py4,pz4,e4,m4;
    double pT1,pT2,pT3,pT4;
    double y1,y2,y3,y4,phi1,phi2,phi3,phi4;
    double detam,dphim,dRm;
    double detambar,dphimbar,dRmbar;
    double detab,dphib,dRb;
    double detabbar,dphibbar,dRbbar;
    Vec4 jet1,jet2;
    
    // Begin event loop. Generate event. Skip if error. List first one.
    for (int iEvent = 0; iEvent < 1000; ++iEvent)
    {
        //Defining foo for for loop in jet removal
        int foo=0;
        
        if(!pythia.next()) continue;
        
        cout << "This is " << iEvent << " th event" << string(114,'=') << endl;
        
        //if (iEvent < 5) {pythia.info.list(); pythia.event.list();}    //print particle lists
        
        // Analyze Slowet jet properties. Every jet variables should be defined after this process.
        slowJet.analyze( pythia.event );      //<== IMPORTANT!!!
        
        //defining jet size before removal
        jetsize = slowJet.sizeJet();

        //Particle loop =========================================================
        for(int i = 0; i < pythia.event.size();++i)
        {
            //declear id, mom1, mom2, mom1id, mom2id
            id = pythia.event[i].id();
            mom1 = pythia.event[i].mother1();
            mom2 = pythia.event[i].mother2();
            mom1id = pythia.event[mom1].id();
            mom2id = pythia.event[mom2].id();
            
            //defining px1,px2,... for b and bbar coming from Higgs
            if(id == 5 && (mom1id == 25 || mom2id == 25))
            {
                pT1 = pythia.event[i].pT();
                y1  = pythia.event[i].y();
                phi1= pythia.event[i].phi();
                px1 = pythia.event[i].px();
                py1 = pythia.event[i].py();
                pz1 = pythia.event[i].pz();
                e1  = pythia.event[i].e();
                m1  = pythia.event[i].m();
            }
            
            if(id == -5 && (mom1id == 25 || mom2id == 25))
            {
                pT2 = pythia.event[i].pT();
                y2  = pythia.event[i].y();
                phi2= pythia.event[i].phi();
                px2 = pythia.event[i].px();
                py2 = pythia.event[i].py();
                pz2 = pythia.event[i].pz();
                e2  = pythia.event[i].e();
                m2  = pythia.event[i].m();
            }
            
            //defining px3,px4,... for muon
            if(id == 13 && (mom1id == 23 || mom2id == 23))
            {
                pT3 = pythia.event[i].pT();
                y3  = pythia.event[i].y();
                phi3= pythia.event[i].phi();
                px3 = pythia.event[i].px();
                py3 = pythia.event[i].py();
                pz3 = pythia.event[i].pz();
                e3 = pythia.event[i].e();
                m3 = pythia.event[i].m();
            }
            
            if(id == -13 && (mom1id == 23 || mom2id == 23))
            {
                pT4 = pythia.event[i].pT();
                y4  = pythia.event[i].y();
                phi4= pythia.event[i].phi();
                px4 = pythia.event[i].px();
                py4 = pythia.event[i].py();
                pz4 = pythia.event[i].pz();
                e4 = pythia.event[i].e();
                m4 = pythia.event[i].m();
            }
            
        }//end of particle loop
        //=======================================================================
        
        //slowJet.list();     //Initial jet list
        
        //Remove jet by discarding jet with dR value
        for(int j=0;j < jetsize;++j)
        {
            //calculation of dR of b + jet
            {
                //dR of b + jet calcualation
                drcalc(slowJet.phi(foo),phi1,slowJet.y(foo),y1,&dphib,&detab,&dRb);
                
                //dR of bbar + jet calcualation
                drcalc(slowJet.phi(foo),phi2,slowJet.y(foo),y2,&dphibbar,&detabbar,&dRbbar);
            
            }
            //calculation of dR of mu + jet
            {
                //dR of mu- + jet calcualation
                drcalc(slowJet.phi(foo),phi3,slowJet.y(foo),y3,&dphim,&detam,&dRm);
                
                //dR of mu+ + jet calcualation
                drcalc(slowJet.phi(foo),phi4,slowJet.y(foo),y4,&dphimbar,&detambar,&dRmbar);
            }
            
            //removal of jet
            {
                /*
                if(dRm < 0.1 || dRmbar < 0.1)
                {
                    slowJet.removeJet(foo);--foo;
                }
                if(!(dRb < 0.01 || dRbbar < 0.01) && !(dRm < 0.1 || dRmbar < 0.1))
                {
                    slowJet.removeJet(foo);--foo;
                }*/
            }
            
            //cout << "dRb = " << dRb << ", dRbbar = " << dRbbar << ", dRm = " << dRm << ", dRmbar = " << dRmbar << endl;
            
            //plots dR distribution of muons and jets
            myfilec << setw(10) << dRm;
            myfilec << setw(10) << dRmbar;
            myfilec << setw(10) << dRb;
            myfilec << setw(10) << dRbbar;
            myfilec << endl;
            
            //counting number of times program went through
            ++foo;
       }
        
        slowJet.list();     //Final jet list
        
        //counting multiplicities of jets
        for(int p=0;p < slowJet.sizeJet();++p)
        {
            cout << "jet mult = " << slowJet.multiplicity(p) << endl;
            myfiled << setw(10) << slowJet.multiplicity(p);
            myfiled << endl;
        }
        
        //Calculation involving dijets
        if(slowJet.sizeJet() > 1)
        {
            //declearing number of jets used to calculate
            jetsizeaf = slowJet.sizeJet();
            
            //selecting two different jets
            for(int n = 0; n < slowJet.sizeJet()-1; ++n)
            {
                for(int l = n+1; l < slowJet.sizeJet(); ++l)
                {
                    //dpT,dphi,deta,dR calculation
                    {
                        dpT = abs( slowJet.pT(n) - slowJet.pT(l) );
                        
                        drcalc(slowJet.phi(n),slowJet.phi(l),slowJet.y(n),slowJet.y(l),&dphi,&deta,&dR);
                    }
                    
                    //Higgs invariant mass
                    {
                        jet1 = slowJet.p(n);
                        jet2 = slowJet.p(l);
                        
                        invm = invariant(slowJet.m(n),slowJet.m(l),jet1.e(),jet2.e(),jet1.px(),jet2.px(),jet1.py(),jet2.py(),jet1.pz(),jet2.pz());
                    }
                    
                    //prints out variable for any JET COMBINATIONS
                    {
                        myfileb << setw(10) << dpT;
                        myfileb << setw(10) << dphi;
                        myfileb << setw(10) << deta;
                        myfileb << setw(10) << dR;
                        myfileb << setw(10) << invm;
                        myfileb << endl;
                    }
                }
            }
        }
        
        //Calling function for Higgs/Z0 mass
        {
            hm = invariant(m1,m2,e1,e2,px1,px2,py1,py2,pz1,pz2);
            zm = invariant(m3,m4,e3,e4,px3,px4,py3,py4,pz3,pz4);
            
            cout << "px of -mu is " << px3 << " and px of +mu is " << px4 << endl;
            cout << "higgs mass from b bbar: " << hm << endl;
            cout << "z mass from mu mubar: " << zm << endl;
        }
        
        //prints out variable for EACH EVENT
        {
            myfile << setw(10) << hm;
            myfile << setw(10) << zm;
            myfile << setw(10) << pT1;
            myfile << setw(10) << pT2;
            myfile << setw(10) << pT3;
            myfile << setw(10) << pT4;
            myfile << setw(10) << y1;
            myfile << setw(10) << y2;
            myfile << setw(10) << y3;
            myfile << setw(10) << y4;
            myfile << setw(10) << jetsizeaf;
            myfile << endl;
        }
        
        cout << endl;
    }//end of iEvent
    
    pythia.stat();
    return 0;
}

//function for calculation of invariant mass
double invariant(double m1,double m2,double e1,double e2,double px1,double px2,double py1,double py2,double pz1,double pz2){
    
    double a;
    
    a = sqrt(m1*m1 + m2*m2 + 2*(e1*e2 - px1*px2 - py1*py2 - pz1*pz2));
    
    return a;
}

void drcalc(double phi1, double phi2, double y1, double y2, double *pdphi, double *pdeta, double *pdR){
    
    *pdphi = abs( phi1 - phi2 );
    *pdeta = y1 - y2;
    if(*pdphi > M_PI) *pdphi = 2. * M_PI - *pdphi;
    *pdR = sqrt( (*pdeta) * (*pdeta) + (*pdphi) * (*pdphi) );
    
}






























