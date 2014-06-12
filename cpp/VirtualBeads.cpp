extern config CFG;

template <typename T> void Swap_elem(std::vector<T> &X, unsigned int E1, unsigned int E2) {
  T buf=X[E1];
  X[E1]=X[E2];
  X[E2]=buf;
}

void BubbleSort_IVEC_using_associated_dVals(vector<int> &X, vector<double> &val) {

  if (X.size()>=2) {

    bool Swapped;
    unsigned int N=X.size();
    unsigned int I;

    do {
      Swapped=false;
      for (I=0;I<=(N-2);I++) {
        if (val[I]>val[I+1]) {
          Swap_elem(X,I,I+1);
          Swap_elem(val,I,I+1);
          Swapped=true;
        }//if
      }//I
      N--;
    } while((Swapped==true)&&(N>1));
  }
}

class VirtualBeads{

    public:

    std::vector<double> S_0;

    std::vector<vec3D> U_found;
    std::vector<vec3D> R_found;
    std::vector<vec3D> U_guess;

    std::vector< std::map<size_t,double> > I;
    std::vector< std::map<size_t,double> > Itrans;
    std::vector< std::map<size_t,double> > ItransI;
    std::vector< vec3D > ItransUfound;

    std::vector<bool> vbead;
    std::vector<bool> outofstack;

    std::vector< std::vector< mat3D > > A_sp;
    std::vector< std::vector< size_t > > A_j;

    std::vector<double> u_x;
    std::vector<double> u_y;
    std::vector<double> u_z;

    std::vector<double> localweight;

    std::vector< std::list<size_t> > conconnections;
    std::vector< std::list<size_t> > oldconconnections;

    std::vector< vec3D > b;

    vec3D Drift;

    int sX,sY,sZ;
    double dX,dY,dZ;

    void correctHorizontalDrift(FiniteBodyForces& M);

    void loadVbeads(std::string VBname);
    void findBeads(const stack3D& stackr);
    void loadUfound(std::string Uname, std::string Sname);
    void loadGuess(FiniteBodyForces& M);
    void loadScatteredRfound(std::string Rname);
    void computeInterpolationMatrix(const FiniteBodyForces& M);
    void allBeads(const FiniteBodyForces&);
    void computeOutOfStack(const FiniteBodyForces&);
    void computeConconnections(FiniteBodyForces& M);
    void computeConconnections_Laplace(FiniteBodyForces& M);

    void mulA(std::vector<double>& u, std::vector<double>& f, const FiniteBodyForces& M);
    double solve_CG(FiniteBodyForces& M, bool local);
    vec3D relax(FiniteBodyForces& M, bool local);

    double testDrift(const stack3D& stack1, const stack3D& stack2, vec3D D);

    void findDisplacements(const stack3D& stack_r, const stack3D& stack_a, FiniteBodyForces& M,double lambda);
    void findDisplacementsScattered(const stack3D& stack_r, const stack3D& stack_a, double lambda);
    void refineDisplacements(const stack3D& stack_r, const stack3D& stack_a, FiniteBodyForces& M, double lambda);
    vec3D findDrift(const stack3D& stack1, const stack3D& stack2);
    vec3D findDriftCoarse(const stack3D& stack1, const stack3D& stack2, double range, double step);

    void computeAAndb(FiniteBodyForces& M, double lambda, double lambda2);
    void recordRelaxationStatus(FiniteBodyForces& M, DRec3D& relrec);

    void updateLocalWeigth(const FiniteBodyForces& M);

    void storeUfound(std::string Uname, std::string Sname);
    void storeRfound(std::string Rname);
    void storeRRandDD(std::string Rname, std::string Uname);
    void storeLocalweights(std::string wname);

    VirtualBeads(int ssX, int ssY, int ssZ, double ddX, double ddY, double ddZ);
    VirtualBeads();

};


VirtualBeads::VirtualBeads(){};

VirtualBeads::VirtualBeads(int ssX, int ssY, int ssZ, double ddX, double ddY, double ddZ){

    sX=ssX;
    sY=ssY;
    sZ=ssZ;
    dX=ddX;
    dY=ddY;
    dZ=ddZ;

}

void VirtualBeads::correctHorizontalDrift(FiniteBodyForces& M){


    std::vector<int> indicesx;
    std::vector<int> indicesy;
    std::vector<int> indicesz;
    std::vector<double> valuesx;
    std::vector<double> valuesy;
    std::vector<double> valuesz;

    vec3D Umean=vec3D();

    int c,vbeadcount=0;

    for(c=0; c<M.N_c; c++) if(vbead[c]){

        if(norm(U_found[c])!=0){

            indicesx.push_back(c);
            indicesy.push_back(c);
            indicesz.push_back(c);
            valuesx.push_back(U_found[c].x);
            valuesy.push_back(U_found[c].y);
            valuesz.push_back(U_found[c].z);

            vbeadcount++;

        }

    }


    if(vbeadcount>0){

        BubbleSort_IVEC_using_associated_dVals(indicesx,valuesx);
        BubbleSort_IVEC_using_associated_dVals(indicesy,valuesy);
        BubbleSort_IVEC_using_associated_dVals(indicesz,valuesz);

        double Umx=valuesx[floor(vbeadcount/2)];
        double Umy=valuesy[floor(vbeadcount/2)];
        double Umz=valuesz[floor(vbeadcount/2)];

        for(c=0; c<M.N_c; c++) if(vbead[c]){

            U_found[c]-=vec3D(Umx,Umy,Umz);

        }

    }

    for(c=0; c<M.N_c; c++) if(vbead[c]){

        U_found[c]-=Umean;

    }

}

void VirtualBeads::allBeads(const FiniteBodyForces& M){

    double thresh=int(CFG["VB_SX"])*dX+2.0*double(CFG["DRIFT_STEP"]);

    vbead.assign(M.N_c,false);

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D scaledShift=vec3D(sX/2,sY/2,sZ/2);

    mat3D Trans=Scale;

    vec3D R;

    for(int t=0; t<M.N_c; t++){

        R=(Trans*((M.R[t])))+scaledShift;

        if( R.x>thresh && R.x<(sX-thresh) && R.y>thresh && R.y<(sY-thresh) && R.z>thresh && R.z<(sZ-thresh) ) vbead[t]=true;

    }

}

void VirtualBeads::computeOutOfStack(const FiniteBodyForces& M){

    double thresh=int(CFG["VB_SX"])*dX+2.0*double(CFG["DRIFT_STEP"]);

    outofstack.assign(M.N_c,true);

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D scaledShift=vec3D(sX/2,sY/2,sZ/2);

    mat3D Trans=Scale;

    vec3D R;

    for(int t=0; t<M.N_c; t++){

        R=(Trans*((M.R[t])))+scaledShift;

        if( R.x>thresh && R.x<(sX-thresh) && R.y>thresh && R.y<(sY-thresh) && R.z>thresh && R.z<(sZ-thresh) ) outofstack[t]=false;

    }

}

void VirtualBeads::loadVbeads(std::string VBname){

    std::vector< std::vector<double> > vbead_temp;

    readFromDatFile(VBname.c_str(),vbead_temp);

    int i=0;
    int iend=vbead_temp.size();

    vbead.resize(iend);

    for(i=0; i<iend; i++){

        vbead[i]=(vbead_temp[i][0]>0.5);

    }

}

void VirtualBeads::findBeads(const stack3D& stackr){

    std::multiset<std::tuple<unsigned char,int,int,int>> vlist=findLocalMinima(stackr,floor(int(CFG["VB_SX"])*0.5),floor(int(CFG["VB_SY"])*0.5),floor(int(CFG["VB_SZ"])*0.5),int(CFG["VB_COUNT"]));

    DRec3D mrec=DRec3D();

    int i,j,k;

    mat3D T=mat3D();

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D scaledShift=vec3D(sX/2,sY/2,sZ/2);

    mat3D Trans=Scale;

    mat3D Transinv=Trans.inv();

    int sgX=int(CFG["VB_SX"]);
    int sgY=int(CFG["VB_SY"]);
    int sgZ=int(CFG["VB_SZ"]);

    std::vector< std::vector< std::vector<double> > > substack;

    substack.assign(sgX,
        std::vector< std::vector< double > >(sgY,
            std::vector< double >(sgZ, 0.0 )
        )
    );

    for(auto n : vlist){

        i=std::get<1>(n);
        j=std::get<2>(n);
        k=std::get<3>(n);

        mrec.record(i,j,k);

        std::vector< std::vector< std::vector<double> > > substacktemp=getSubstack(stackr,vec3D(i,j,k));

        for(auto ii=0; ii<sgX; ii++) for(auto jj=0; jj<sgY; jj++) for(auto kk=0; kk<sgZ; kk++){

            substack[ii][jj][kk]+=substacktemp[ii][jj][kk];


        }

    }

    std::string outdir=std::string(CFG["DATAOUT"]);
    std::string mrecname=outdir+std::string("/")+std::string(CFG["VB_MREC"]);
    mrec.store(mrecname.c_str());


    DRec ssrec=DRec();

    double var=0.0;

    for(auto ii=0; ii<sgX; ii++) for(auto jj=0; jj<sgY; jj++) for(auto kk=0; kk<sgZ; kk++){

        var+=substack[ii][jj][kk]*substack[ii][jj][kk];


    }

    var=sqrt(var);

    for(auto ii=0; ii<sgX; ii++) for(auto jj=0; jj<sgY; jj++) for(auto kk=0; kk<sgZ; kk++){

        substack[ii][jj][kk]=substack[ii][jj][kk]/var;

        ssrec.record(substack[ii][jj][kk]);

    }

    std::string ssrecname=outdir+std::string("/ssrec.dat");
    ssrec.store(ssrecname.c_str());

    DRec Srec=DRec();

    R_found.clear();

    int nend=vlist.size();

    int count=0;

    for(auto n : vlist){

        i=std::get<1>(n);
        j=std::get<2>(n);
        k=std::get<3>(n);

        vec3D Rf=vec3D(i,j,k)+findLocalDisplacement(substack,stackr,vec3D(i,j,k),vec3D(0.0,0.0,0.0),Srec,double(CFG["VB_REGPARA"]));

        if(Srec.data.back()>double(CFG["VB_MINMATCH"])) R_found.push_back( Transinv*(Rf-scaledShift) );

        std::cout<<"examining candidates "<<(floor((count/(nend+0.0))*1000)+0.0)/10.0<<"% S="<<Srec.data.back()<<"      \r";

        count++;

    }


}

void VirtualBeads::loadGuess(FiniteBodyForces& M){

    std::string indir=std::string(CFG["DATAIN"]);

    M.loadConfiguration((indir+std::string("/Uguess.dat")).c_str());

    M.computeConnections();

    computeInterpolationMatrix(M);

    U_guess.clear();
    U_guess.assign(R_found.size(),vec3D(0.0,0.0,0.0));

    for(auto i=0; i<R_found.size(); i++) for(auto itIT=I[i].begin(); itIT!=I[i].end(); itIT++) {

        //auto j=itIT->first;

        U_guess[i] += M.U[itIT->first] * itIT->second;

    }

    M.U.clear();
    M.U.assign(M.N_c,vec3D(0.0,0.0,0.0));

}

void VirtualBeads::updateLocalWeigth(const FiniteBodyForces& M){


    std::vector<int> indices;
    std::vector<double> Fvalues;

    indices.clear();
    Fvalues.clear();

    double k=1.345;

    if(std::string(CFG["ROBUSTMETHOD"])==std::string("bisquare")) k=4.685;
    if(std::string(CFG["ROBUSTMETHOD"])==std::string("cauchy")) k=2.385;

    int c;
    int scount=0;

    for(c=0; c<M.N_c; c++) if(M.var[c]){

        indices.push_back(c);
        Fvalues.push_back(sqrt(M.f_glo[3*c]*M.f_glo[3*c]+M.f_glo[3*c+1]*M.f_glo[3*c+1]+M.f_glo[3*c+2]*M.f_glo[3*c+2]));
        scount++;

    }

    BubbleSort_IVEC_using_associated_dVals(indices,Fvalues);

    int Nrenew=floor(scount/2);
    //int Nrenew=50;

    localweight.assign(M.N_c,1.0);
    //for(c=0; c<M.N_c; c++) if(outofstack[c]) localweight[c]=double(CFG["LAMBDA2"]);

    double counter=0.0,counterall=0.0;

    c=indices[scount-Nrenew];

    //size_t cmax=indices[scount-1];

    double Flocal=0.0;

    double Fmedian=sqrt(M.f_glo[3*c]*M.f_glo[3*c]+M.f_glo[3*c+1]*M.f_glo[3*c+1]+M.f_glo[3*c+2]*M.f_glo[3*c+2]);

    if(std::string(CFG["ROBUSTMETHOD"])==std::string("singlepoint")) localweight[int(CFG["REG_FORCEPOINT"])]=1.0e-10;

    for(int i=0; i<scount; i++){

        c=indices[i];

        //localweight[c]=exp(-1.0*(M.f_glo[3*c]*M.f_glo[3*c]+M.f_glo[3*c+1]*M.f_glo[3*c+1]+M.f_glo[3*c+2]*M.f_glo[3*c+2])/1000.0)+0.001;

        Flocal=sqrt(M.f_glo[3*c]*M.f_glo[3*c]+M.f_glo[3*c+1]*M.f_glo[3*c+1]+M.f_glo[3*c+2]*M.f_glo[3*c+2]);

        if(std::string(CFG["ROBUSTMETHOD"])==std::string("bisquare")){

            if(Flocal<(k*Fmedian)) localweight[c]*=( 1 - (Flocal/k/Fmedian)*(Flocal/k/Fmedian) )*( 1 - (Flocal/k/Fmedian)*(Flocal/k/Fmedian) );
            else localweight[c]*=1e-10;

        }

        if(std::string(CFG["ROBUSTMETHOD"])==std::string("cauchy")){

            if(Fmedian>0) localweight[c]*=1.0/(1.0+pow((Flocal/k/Fmedian),2.0));
            else localweight[c]*=1.0;

        }

        if(std::string(CFG["ROBUSTMETHOD"])==std::string("huber")){

            if(Flocal>(k*Fmedian)) localweight[c]*=k*Fmedian/Flocal;
            else localweight[c]*=1.0;

        }

        if(localweight[c]<1e-10) localweight[c]=1e-10;

        counter+=(1.0-localweight[c]);
        counterall++;

    }

    std::cout<<"total weight: "<<counter<<"/"<<counterall<<" \n";

}

void VirtualBeads::loadUfound(std::string Uname, std::string Sname){

    DRec3D Urec=DRec3D();
    std::vector< std::vector<double> > Srec;
    //std::cout<<"check inside loadUfound    \n";

    Urec.read(Uname.c_str());
    //Srec.read(Sname.c_str());
    readFromDatFile(Sname.c_str(),Srec);

    std::vector< std::vector<double> > Ufound_temp;

    int i=0;
    int iend=Urec.data.size();

    int js=Srec[0].size()-1;

    //int outcount=0;

    U_found.resize(iend);
    S_0.resize(iend);
    vbead.resize(iend);

    //std::cout<<"check inside loadUfound    \n";
    for(i=0; i<iend; i++){

        //std::cout<<"check inside loadUfound in loop at i="<<i<<" S_0=="<<Srec[i][js]<<"    \r";
        U_found[i]=Urec.data[i];
        S_0[i]=Srec[i][js];

        if(S_0[i]>0.7) vbead[i]=true;
        else vbead[i]=false;

    }
    //std::cout<<"check inside loadUfound    \n";

    Urec.data.clear();
    Srec.clear();

}

void VirtualBeads::loadScatteredRfound(std::string Rname){

    DRec3D Rrec=DRec3D();

    Rrec.read(Rname.c_str());

    int i=0;
    int iend=Rrec.data.size();

    //int outcount=0;

    R_found.resize(iend);

    for(i=0; i<iend; i++){

        R_found[i]=Rrec.data[i];

    }

    Rrec.data.clear();

}

void VirtualBeads::computeInterpolationMatrix(const FiniteBodyForces& M){

    int i=0;
    int iend=R_found.size();

    mat3D indentity=mat3D(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);

    I.resize(iend);
    for(int i=0; i<(iend); i++) I[i].clear();

    Itrans.resize(M.N_c);
    for(int i=0; i<(M.N_c); i++) Itrans[i].clear();


    int foundcount=0;


    for(i=0; i<iend; i++){

        if( (i%10)==0 ) std::cout<<"localizing scattered input in mesh "<<(floor((i/(iend+0.0))*1000)+0.0)/10.0<<"%                \r";

        //std::cout<<"Rf=("<<R_found[i].x<<","<<R_found[i].y<<","<<R_found[i].z<<") \n";

        for(int t=0; t<M.N_T; t++){

            //std::cout<<"check  check               \n";

//            vec3D a=M.R[M.T[t][0]];
//            vec3D b=M.R[M.T[t][1]];
//            vec3D c=M.R[M.T[t][2]];
//            vec3D d=M.R[M.T[t][3]];


            mat3D B=mat3D(

                M.R[M.T[t][1]].x-M.R[M.T[t][0]].x, M.R[M.T[t][2]].x-M.R[M.T[t][0]].x, M.R[M.T[t][3]].x-M.R[M.T[t][0]].x,
                M.R[M.T[t][1]].y-M.R[M.T[t][0]].y, M.R[M.T[t][2]].y-M.R[M.T[t][0]].y, M.R[M.T[t][3]].y-M.R[M.T[t][0]].y,
                M.R[M.T[t][1]].z-M.R[M.T[t][0]].z, M.R[M.T[t][2]].z-M.R[M.T[t][0]].z, M.R[M.T[t][3]].z-M.R[M.T[t][0]].z

            );


            //std::cout<<B.det()<<"       \n";

//            B=mat3D(
//                R[T[t][1]].x-R[T[t][0]].x, R[T[t][1]].y-R[T[t][0]].y, R[T[t][1]].z-R[T[t][0]].z,
//                R[T[t][2]].x-R[T[t][0]].x, R[T[t][2]].y-R[T[t][0]].y, R[T[t][2]].z-R[T[t][0]].z,
//                R[T[t][3]].x-R[T[t][0]].x, R[T[t][3]].y-R[T[t][0]].y, R[T[t][3]].z-R[T[t][0]].z
//            );

            if(B.det()>0){

                vec3D q=B.inv()*(R_found[i]- M.R[M.T[t][0]]);

                /*

                if(abs((R_found[i]- M.R[M.T[t][0]]))<double(CFG["BM_GRAIN"])){



                    std::cout<<"q=("<<q.x<<","<<q.y<<","<<q.z<<") \n"
                              <<"Ra=("<<M.R[M.T[t][0]].x<<","<<M.R[M.T[t][0]].y<<","<<M.R[M.T[t][0]].z<<") \n"
                             <<"Rb=("<<M.R[M.T[t][1]].x<<","<<M.R[M.T[t][1]].y<<","<<M.R[M.T[t][1]].z<<") \n"
                            <<"Rc=("<<M.R[M.T[t][2]].x<<","<<M.R[M.T[t][2]].y<<","<<M.R[M.T[t][2]].z<<") \n"
                            <<"Rd=("<<M.R[M.T[t][3]].x<<","<<M.R[M.T[t][3]].y<<","<<M.R[M.T[t][3]].z<<") \n\n";


                }*/

                if( q.x>0.0 && q.y>0.0 && q.z>0.0 && (q.x+q.y+q.z)<1.0 ){

                    //auto Xi=M.Phi[t]*(R_found[i] - M.R[TT[0]]);
                    //Xi[0]+=1;

                    //int debug=M.T[t][0];
                    //double debug2=1.0-(q.x+q.y+q.z);

                    //if(i==230) std::cout<<"i="<<i<<" \n"
                    //           <<"q=("<<q.x<<","<<q.y<<","<<q.z<<") \n"
                    //           <<"T=("<<M.T[t][0]<<","<<M.T[t][1]<<","<<M.T[t][2]<<","<<M.T[t][3]<<","<<")\n\n";

                    I[i][M.T[t][0]]+=1.0-(q.x+q.y+q.z);
                    I[i][M.T[t][1]]+=q.x;
                    I[i][M.T[t][2]]+=q.y;
                    I[i][M.T[t][3]]+=q.z;

                    Itrans[M.T[t][0]][i]+=1.0-(q.x+q.y+q.z);
                    Itrans[M.T[t][1]][i]+=q.x;
                    Itrans[M.T[t][2]][i]+=q.y;
                    Itrans[M.T[t][3]][i]+=q.z;

                    foundcount++;

                    //std::cout<<"found a bead!   \n";

                    break;

                }

            }


        }



    }

    if(foundcount<iend) std::cout<<"\nWARNING: "<<iend-foundcount<<" of "<<iend<<" beads could not be localized in mesh!!!        \n";

    if(false){

        DRec3D Irec=DRec3D();

        for(int i=0; i<iend; i++) for(auto itj=I[i].begin(); itj!=I[i].end(); itj++) Irec.record(i+1,itj->first+1,itj->second);

        std::string outdir=std::string(CFG["DATAOUT"]);
        std::string relrecname=outdir+std::string("/I.dat");
        Irec.store(relrecname.c_str());

    }

    ItransI.resize(M.N_c);
    for(int i=0; i<(M.N_c); i++) ItransI[i].clear();

    for(int i1=0; i1<M.N_c; i1++){

        if( (i%100)==0 ) std::cout<<"computing interpolation matrix "<<(floor((i/(iend+0.0))*1000)+0.0)/10.0<<"%                \r";

        for(auto itc=M.connections[i1].begin(); itc!=M.connections[i1].end(); itc++){

            int i2=*itc;

            auto it1 = Itrans[i1].begin();
            auto it2 = Itrans[i2].begin();
            while ( it1 != Itrans[i1].end() && it2 != Itrans[i2].end() )
            {
                //auto debug1=it1->first;
                //auto debug2=it2->first;

                if (it1->first < it2->first){

                    ++it1;

                }
                else if (it2->first < it1->first) {

                    ++it2;

                }
                else {

                    ItransI[i1][i2] += it1->second * it2->second;

                    ++it1;
                    ++it2;

                }

            }

        }

    }

    if(false){

        DRec3D ItransIrec=DRec3D();

        for(int i=0; i<M.N_c; i++) for(auto itj=ItransI[i].begin(); itj!=ItransI[i].end(); itj++) ItransIrec.record(i+1,itj->first+1,itj->second);

        std::string outdir=std::string(CFG["DATAOUT"]);
        std::string relrecname=outdir+std::string("/ItransI.dat");
        ItransIrec.store(relrecname.c_str());

    }

    ItransUfound.clear();
    ItransUfound.assign(M.N_c,vec3D(0.0,0.0,0.0));

    if(U_found.size()>0){

        std::cout<<"check ItransUfound size of U_found is "<<U_found.size()<<"   \n";

        for(i=0; i<M.N_c; i++) for(auto itIT=Itrans[i].begin(); itIT!=Itrans[i].end(); itIT++) {

            //auto j=itIT->first;

            ItransUfound[i] += U_found[itIT->first] * itIT->second;
        }

    }


    if(false){

        DRec3D ItransUfoundrec=DRec3D();

        for(int i=0; i<M.N_c; i++) ItransUfoundrec.record(ItransUfound[i]);

        std::string outdir=std::string(CFG["DATAOUT"]);
        std::string relrecname=outdir+std::string("/ItransUfound.dat");
        ItransUfoundrec.store(relrecname.c_str());

    }

    //I.clear();
    Itrans.clear();

}

void VirtualBeads::computeConconnections(FiniteBodyForces& M){

    int c1,c2,c3; //tt,t1,t2,

    conconnections.resize(M.N_c);

    for(c1=0;c1<M.N_c;c1++){

        conconnections[c1].clear();

    }

    std::list<size_t>::iterator it,it1,it2;

    bool found=false;

    for(c1=0; c1<M.N_c; c1++){

        if( (c1%100)==0 ) std::cout<<"computing conconnections "<<(floor((c1/(M.N_c+0.0))*1000)+0.0)/10.0<<"%          \r";

        for(it1=M.connections[c1].begin(); (it1!=M.connections[c1].end()); it1++){

            c2=*it1;

            for(it2=M.connections[c2].begin(); (it2!=M.connections[c2].end()); it2++){

                c3=*it2;

                found=false;

                for(it=conconnections[c1].begin(); (it!=conconnections[c1].end() && !found); it++){

                    if(*it==c3) found=true;

                }

                if(!found) conconnections[c1].push_back(c3);

            }

        }

    }


    for(c1=0;c1<M.N_c;c1++) conconnections[c1].sort();

}

void VirtualBeads::computeConconnections_Laplace(FiniteBodyForces& M){

    int c1,c2,c3; //tt,t1,t2,

    oldconconnections=std::vector< std::list<size_t> >(conconnections);

    conconnections.resize(M.N_c);

    for(c1=0;c1<M.N_c;c1++){

        conconnections[c1].clear();

    }

    std::list<size_t>::iterator it,it1,it2;

    bool found=false;

    for(c1=0; c1<M.N_c; c1++){

        if( (c1%100)==0 ) std::cout<<"computing conconnections "<<(floor((c1/(M.N_c+0.0))*1000)+0.0)/10.0<<"%          \r";

        for(it1=oldconconnections[c1].begin(); (it1!=oldconconnections[c1].end()); it1++){

            c2=*it1;

            for(it2=oldconconnections[c2].begin(); (it2!=oldconconnections[c2].end()); it2++){

                c3=*it2;

                found=false;

                for(it=conconnections[c1].begin(); (it!=conconnections[c1].end() && !found); it++){

                    if(*it==c3) found=true;

                }

                if(!found) conconnections[c1].push_back(c3);

            }

        }

    }


    for(c1=0;c1<M.N_c;c1++) conconnections[c1].sort();

}

vec3D VirtualBeads::findDriftCoarse(const stack3D& stackr, const stack3D& stacka, double range, double step){

    double Stemp;
    double Smax=-1.0;
    double dxmax=0.0,dymax=0.0,dzmax=0.0;

    //DRec3D recxyz=DRec3D();
    //DRec recS=DRec();

    for(double dx=-range; dx<range; dx+=step) for(double dy=-range; dy<range; dy+=step) for(double dz=-range; dz<range; dz+=step){

        Stemp=testDrift(stackr,stacka,vec3D(dx,dy,dz));

        if(Stemp>Smax){
            dxmax=dx;
            dymax=dy;
            dzmax=dz;
            Smax=Stemp;

            //std::cout<<dx<<" "<<dy<<" "<<dz<<" : "<<Smax<<"   \n";
        }

        //recxyz.record(vec3D(dx,dy,dz));
        //recS.record(Stemp);

    }

    //recxyz.store("coarsdriftxyz.dat");
    //recS.store("coarsdriftS.dat");

    return vec3D(dxmax,dymax,dymax);

}

vec3D VirtualBeads::findDrift(const stack3D& stackr, const stack3D& stacka){


    double lambda=0.0;

    std::vector<vec3D> P;

    std::vector<double> S;

    double hinit=double(CFG["DRIFT_STEP"]);
    double subpixel=double(CFG["SUBPIXEL"]);

    double vsx=double(CFG["VOXELSIZEX"]);
    double vsy=double(CFG["VOXELSIZEY"]);
    double vsz=double(CFG["VOXELSIZEZ"]);

    P.resize(4);
    S.resize(4);

    P[0]=vec3D(1.0,1.0,1.0)*hinit;
    P[1]=vec3D(-1.0,-1.0,1.0)*hinit;
    P[2]=vec3D(-1.0,1.0,-1.0)*hinit;
    P[3]=vec3D(1.0,-1.0,-1.0)*hinit;

    S[0]=testDrift(stackr,stacka,P[0]+Drift)-lambda*norm(P[0]);
    S[1]=testDrift(stackr,stacka,P[1]+Drift)-lambda*norm(P[1]);
    S[2]=testDrift(stackr,stacka,P[2]+Drift)-lambda*norm(P[2]);
    S[3]=testDrift(stackr,stacka,P[3]+Drift)-lambda*norm(P[3]);

    //std::cout<<S[0]<<" "<<S[1]<<" "<<S[2]<<" "<<S[3]<<" \n";

    size_t mini=minimal(S);
    size_t maxi=maximal(S);

    bool done=false;

    vec3D P_ref;
    vec3D P_exp;
    vec3D P_con;
    vec3D P_plane;

    double S_ref;
    double S_exp;
    double S_con;

    double mx,my,mz,stdx,stdy,stdz;

    int i;

    int ii;

    for(ii=0; (ii<100000 && !done); ii++){

        mini=minimal(S);
        maxi=maximal(S);

        //std::cout<<"| new cycle mini = "<<mini<<std::endl;

        //reflect

        P_plane=vec3D(0.0,0.0,0.0);
        for(i=0; i<4; i++) if(i!=mini) P_plane+=P[i];
        P_plane=P_plane*0.333333333;

        P_ref=P[mini]+(P_plane-P[mini])*2.0;

        S_ref=testDrift(stackr,stacka,P_ref+Drift)-lambda*norm(P_ref);

        if(S_ref>S[maxi]){

            //expand

            //std::cout<<"| expanding "<<std::endl;

            P_exp=P[mini]+(P_plane-P[mini])*3.0;

            S_exp=testDrift(stackr,stacka,P_exp+Drift)-lambda*norm(P_exp);

            if(S_exp>S_ref){

                //std::cout<<"| took expanded "<<std::endl;

                P[mini]=P_exp;
                S[mini]=S_exp;

            }else{

                //std::cout<<"| took reflected (expanded was worse)"<<std::endl;

                P[mini]=P_ref;
                S[mini]=S_ref;

            }

        }else{

            bool bsw=false;
            for(i=0;i<4;i++) if(i!=mini) if(S_ref>S[i]) bsw=true;

            if(bsw){

                //std::cout<<"| took reflected (better than second worst)"<<std::endl;

                P[mini]=P_ref;
                S[mini]=S_ref;

            }else{

                if(S_ref>S[maxi]){

                    //std::cout<<"| took reflected (not better than second worst)"<<std::endl;

                    P[mini]=P_ref;
                    S[mini]=S_ref;

                }else{

                    P_con=P[mini]+(P_plane-P[mini])*0.5;
                    S_con=testDrift(stackr,stacka,P_con+Drift)-lambda*norm(P_con);

                    if(S_con>S[mini]){

                        //std::cout<<"| took contracted"<<std::endl;

                        P[mini]=P_con;
                        S[mini]=S_con;

                    }else{

                        //std::cout<<"| contracting myself"<<std::endl;

                        for(i=0; i<4; i++) if(i!=maxi) {

                            P[i]=(P[maxi]+P[i])*0.5;
                            S[i]=testDrift(stackr,stacka,P[i]+Drift)-lambda*norm(P[i]);

                        }

                    }

                }

            }

        }

        //std::cout<<" S_ref = "<<S_ref<<std::endl;

        mx=(P[0].x+P[1].x+P[2].x+P[3].x)*0.25;
        my=(P[0].y+P[1].y+P[2].y+P[3].y)*0.25;
        mz=(P[0].z+P[1].z+P[2].z+P[3].z)*0.25;

        stdx=sqrt( (P[0].x-mx)*(P[0].x-mx) + (P[1].x-mx)*(P[1].x-mx) + (P[2].x-mx)*(P[2].x-mx) + (P[3].x-mx)*(P[3].x-mx) )*0.25;
        stdy=sqrt( (P[0].y-my)*(P[0].y-my) + (P[1].y-my)*(P[1].y-my) + (P[2].y-my)*(P[2].y-my) + (P[3].y-my)*(P[3].y-my) )*0.25;
        stdz=sqrt( (P[0].z-mz)*(P[0].z-mz) + (P[1].z-mz)*(P[1].z-mz) + (P[2].z-mz)*(P[2].z-mz) + (P[3].z-mz)*(P[3].z-mz) )*0.25;

        //std::cout<<"stdx = "<<stdx<<" ; stdy = "<<stdy<<" ; stdz = ";

        if(stdx<subpixel*vsx && stdy<subpixel*vsy && stdz<subpixel*vsz) done=true;


    }

    mini=minimal(S);
    maxi=maximal(S);


    return P[maxi]+Drift;


}

double VirtualBeads::testDrift(const stack3D& stack1, const stack3D& stack2, vec3D D){

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D U=(Scale*(D));

    return crosscorrelateStacks(stack1,stack2,U);

}

void VirtualBeads::findDisplacements(const stack3D& stack_r, const stack3D& stack_a, FiniteBodyForces& M, double lambda){

    DRec Srec=DRec();

    vec3D R,U;

    mat3D T=mat3D();

    U_found.assign(M.N_c,vec3D(0.0,0.0,0.0));
    S_0.assign(M.N_c,0.0);

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D scaledShift=vec3D(sX/2,sY/2,sZ/2);

    mat3D Trans=Scale;

    mat3D Transinv=Trans.inv();

    int t;
    int tend=M.R.size();

    vec3D U_new,Umean;

    int count,tt;

    double Stemp=0.0;

    int n=0;

    std::list<double> Sd;
    std::list<vec3D> Ud;

    double dx,dy,dz;

    for(t=0; t<tend; t++) if(vbead[t]){

        std::cout<<"finding Displacements "<<(floor((t/(tend+0.0))*1000)+0.0)/10.0<<"% S="<<Stemp<<" n="<<n<<"      \r";

        if(int(CFG["VB_N"])==1){

            R=(Trans*((M.R[t]))+scaledShift);

            Umean=Drift;

            Umean=Trans*Umean;

            substack substackr=getSubstack(stack_r,R);
            U_new=findLocalDisplacement(substackr,stack_a,R,Umean,Srec,lambda);
            S_0[t]=Srec.data.back();
            Stemp=Srec.data.back();

            U_new=Transinv*U_new;

            U_found[t]=U_new;

        }else{

            Ud.clear();
            Sd.clear();

            int imax=int(CFG["VB_N"])-1;

            for(int i=0; i<(imax+1); i++) for(int j=0; j<(imax+1); j++) for(int k=0; k<(imax+1); k++){

                dx=(i-(imax*0.5))/(imax+1.0)*double(CFG["VB_SX"]);
                dy=(j-(imax*0.5))/(imax+1.0)*double(CFG["VB_SY"]);
                dz=(k-(imax*0.5))/(imax+1.0)*double(CFG["VB_SZ"]);

                R=(Trans*((M.R[t]))+scaledShift+vec3D(dx,dy,dz));

                Umean=Drift;

                Umean=Trans*Umean;

                substack substackr=getSubstack(stack_r,R);
                U_new=findLocalDisplacement(substackr,stack_a,R,Umean,Srec,lambda);
                Stemp=Srec.data.back();

                U_new=Transinv*U_new;

                if(Stemp>double(CFG["VB_MINMATCH"])){

                    Ud.push_back(U_new);
                    Sd.push_back(Stemp);

                }

            }

            if(Sd.size()>0){

                U_new=vec3D(0.0,0.0,0.0);
                Stemp=0.0;

                for(std::list<double>::iterator it=Sd.begin(); it!=Sd.end(); it++) Stemp+=*it;
                for(std::list<vec3D>::iterator it=Ud.begin(); it!=Ud.end(); it++) U_new+=*it;

                U_new=U_new*(1.0/Sd.size());
                Stemp=Stemp*(1.0/Sd.size());

                n=Sd.size();

            }else{

                U_new=vec3D(0.0,0.0,0.0);
                Stemp=-1.0;
                n=0;

            }

            S_0[t]=Srec.data.back();
            U_found[t]=U_new;

        }

    }

}

void VirtualBeads::findDisplacementsScattered(const stack3D& stack_r, const stack3D& stack_a, double lambda){

    DRec Srec=DRec();

    vec3D R,U;

    mat3D T=mat3D();

    int t;
    int tend=R_found.size();

    U_found.assign(tend,vec3D(0.0,0.0,0.0));
    S_0.assign(tend,0.0);

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D scaledShift=vec3D(sX/2,sY/2,sZ/2);

    mat3D Trans=Scale;

    mat3D Transinv=Trans.inv();

    vec3D U_new,Umean;

    int count,tt;

    double Stemp=0.0;

    int n=0;

    std::list<double> Sd;
    std::list<vec3D> Ud;

    double dx,dy,dz;

    bool guesspresent=(U_guess.size()>0);

    for(t=0; t<tend; t++){

        std::cout<<"finding Displacements "<<(floor((t/(tend+0.0))*1000)+0.0)/10.0<<"% S="<<Stemp<<" n="<<n<<"      \r";

        R=(Trans*((R_found[t]))+scaledShift);

        Umean=Drift;

        if(guesspresent) Umean+=U_guess[t];

        Umean=Trans*Umean;

        substack substackr=getSubstack(stack_r,R);
        U_new=findLocalDisplacement(substackr,stack_a,R,Umean,Srec,lambda);
        S_0[t]=Srec.data.back();
        Stemp=Srec.data.back();

        U_new=Transinv*U_new;

        U_found[t]=U_new;

    }

    Umean=vec3D(0.0,0.0,0.0);

    int vbcount=0;

    for(t=0; t<U_found.size(); t++){

        if(S_0[t]>float(CFG["VB_MINMATCH"])) {

        Umean+=U_found[t];
        vbcount++;

        }else{

            S_0.erase(S_0.begin()+t);
            U_found.erase(U_found.begin()+t);
            R_found.erase(R_found.begin()+t);
            t--;

        }

    }

    Umean=Umean/(vbcount+0.0);

    for(t=0; t<U_found.size(); t++) U_found[t]-=Umean;

}

void VirtualBeads::refineDisplacements(const stack3D& stack_r, const stack3D& stack_a, FiniteBodyForces& M, double lambda){

    DRec Srec=DRec();

    mat3D Scale=mat3D(1.0/dX,0.0,0.0,0.0,1.0/dY,0.0,0.0,0.0,1.0/dZ);

    vec3D scaledShift=vec3D(sX/2,sY/2,sZ/2);

    mat3D Trans=Scale;

    mat3D Transinv=Trans.inv();

    std::vector<int> indices;
    std::vector<double> Svalues;

    indices.clear();
    Svalues.clear();

    int c;
    int vbeadcount=0;

    for(c=0; c<M.N_c; c++) if(vbead[c]){

        indices.push_back(c);
        Svalues.push_back(S_0[c]);
        vbeadcount++;

    }

    BubbleSort_IVEC_using_associated_dVals(indices,Svalues);

    int Nrenew=floor(vbeadcount/4);

    std::list<size_t>::iterator it;

    int cc,cccount;

    vec3D Umean,U_new,R;

    double Sold;

    for(int i=0; i<Nrenew; i++){

        std::cout<<"refining displacements: "<<(floor((i/(Nrenew+0.0))*1000)+0.0)/10.0<<"%     S="<<Svalues[i]<<"            \r";

        c=indices[i];
        cccount=0;

        Umean=vec3D(0.0,0.0,0.0);

        for(it=M.connections[c].begin(); it!=M.connections[c].end(); it++){

            cc=*it;

            if(vbead[cc]){

                Umean+=U_found[cc];
                cccount++;

            }

        }

        if(cccount>0){

            Umean=Umean/(cccount+0.0);

            R=(Trans*((M.R[c]))+scaledShift);
            Umean=(Trans*((Umean)));

            substack substackr=getSubstack(stack_r,R);
            U_new=findLocalDisplacement(substackr,stack_a,R,Umean,Srec,lambda);
            S_0[c]=Srec.data.back();
            substack substacka=getSubstack(stack_a,R);
            U_new-=findLocalDisplacement(substacka,stack_r,R+Umean,Umean*(-1.0),Srec,lambda);
            S_0[c]+=Srec.data.back();

            S_0[c]*=0.5;

            U_new=Transinv*U_new*0.5;

            U_found[c]=U_new;

        }


    }


}

void VirtualBeads::computeAAndb(FiniteBodyForces& M, double lambda, double lambda2){

    bool uselaplace=(std::string(CFG["REGMETHOD"])=="laplace");

    bool scatteredRfound=bool(CFG["SCATTEREDRFOUND"]);

    double lagrain=double(CFG["REG_LAPLACEGRAIN"]);
    lagrain=lagrain*lagrain;

    double llambda_z=double(CFG["REG_SIGMAZ"]);

    A_sp.resize(M.N_c);
    A_j.resize(M.N_c);
    b.resize(M.N_c);

    double mod;

    int i,j,k,saddlecount=0;

    mat3D A_ijtemp,KA_ik,AK_kj,KL_ik,LK_kj,KApL2_ij;

    bool hasentry,ifound,jfound;

    std::list<size_t>::iterator it,itj,it2;

    for(i=0; i<M.N_c; i++) if(M.var[i]){

        if( (i%100)==0 ) std::cout<<"computing A "<<(floor((i/(M.N_c+0.0))*1000)+0.0)/10.0<<"%               \r";

        A_j[i].clear();
        A_sp[i].clear();

        for(itj=conconnections[i].begin(); itj!=conconnections[i].end(); itj++){

            j=*itj;

            A_ijtemp=mat3D(0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0 );

            hasentry=false;

            if(M.var[j]){

                std::list<size_t> intersection=std::list<size_t>();
                intersection.clear();


                if(uselaplace){

                    std::set_intersection (oldconconnections[i].begin(), oldconconnections[i].end(),
                                              oldconconnections[j].begin(), oldconconnections[j].end(),
                                               std::back_inserter(intersection));

                }else{

                    std::set_intersection (M.connections[i].begin(), M.connections[i].end(),
                                              M.connections[j].begin(), M.connections[j].end(),
                                               std::back_inserter(intersection));

                }


                /*
                std::cout<<"oldcon[i]: ";
                for(std::list<size_t>::iterator it2=oldconconnections[i].begin(); it2!=oldconconnections[i].end(); it2++) std::cout<<" "<<*it2<<",";
                std::cout<<"\noldcon[j]: ";
                for(std::list<size_t>::iterator it2=oldconconnections[j].begin(); it2!=oldconconnections[j].end(); it2++) std::cout<<" "<<*it2<<",";
                std::cout<<"\nintersection: ";
                for(std::list<size_t>::iterator it2=intersection.begin(); it2!=intersection.end(); it2++) std::cout<<" "<<*it2<<",";
                std::cout<<"\n";
                */

                for(std::list<size_t>::iterator it2=intersection.begin(); it2!=intersection.end(); it2++) if(M.var[*it2]){

                    k=*it2;

                    //std::cout<<"in cAAndB at i="<<i<<" j="<<j<<" k="<<k<<"      \n";

                    hasentry=true;

                    if(uselaplace){

                        KL_ik=dot(M.K_glo[i],M.Laplace[k]);
                        LK_kj=dot(M.Laplace[k],M.K_glo[j]);

                        if(M.K_glo[i].find(k)!=M.K_glo[i].end()) KA_ik=M.K_glo[i][k];
                        else KA_ik=mat3D(0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0);
                        if(M.K_glo[k].find(j)!=M.K_glo[k].end()) AK_kj=M.K_glo[k][j];
                        else AK_kj=mat3D(0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0);


                        A_ijtemp+=KA_ik*AK_kj*lambda2*localweight[k];
                        //A_ijtemp+=AK_ik*LK_kj*lambda2*sqrt(localweight[k]*lagrain);
                        //A_ijtemp+=LK_ik*AK_kj*lambda2*sqrt(localweight[j]*lagrain);
                        A_ijtemp+=KL_ik*LK_kj*lambda2*lagrain;


                    }else{

                        KA_ik=M.K_glo[i][k];
                        AK_kj=M.K_glo[k][j];

                        A_ijtemp+=KA_ik*AK_kj*lambda2*localweight[k];

                    }


                }

            }

            if( scatteredRfound && ItransI[i].find(j)!=ItransI[i].end() ) {

                //std::cout<<"ItransI has entry i="<<i<<" j="<<j<<" : "<<ItransI[i][j]<<"   \n";

                hasentry=true;

                A_ijtemp.xx+=ItransI[i][j];
                A_ijtemp.yy+=ItransI[i][j];
                A_ijtemp.zz+=ItransI[i][j]*llambda_z;

            }else if(i==j && vbead[i]){

                hasentry=true;

                A_ijtemp.xx+=1.0;
                A_ijtemp.yy+=1.0;
                A_ijtemp.zz+=llambda_z;

            }

            if(hasentry){

                A_j[i].push_back(j);
                A_sp[i].push_back(A_ijtemp);

            }

        }

    }

    vec3D btemp,utemp;

    b.assign(M.N_c,vec3D(0.0,0.0,0.0));

    std::vector< vec3D > f=std::vector< vec3D >();
    f.assign(M.N_c,vec3D(0.0,0.0,0.0));

    for(i=0; i<M.N_c; i++) f[i]=vec3D( M.f_glo[3*i], M.f_glo[3*i+1], M.f_glo[3*i+2] );

    std::vector< vec3D > lf=std::vector< vec3D >();
    std::vector< vec3D > llf=std::vector< vec3D >();

    if(uselaplace){

        lf.assign(M.N_c,vec3D(0.0,0.0,0.0));
        mul(M.Laplace, f, lf);
        llf.assign(M.N_c,vec3D(0.0,0.0,0.0));
        mul(M.Laplace, lf, llf);

    }

    for(i=0; i<M.N_c; i++) if(M.var[i]){

        if( (i%100)==0 ) std::cout<<"computing b "<<(floor((i/(M.N_c+0.0))*1000)+0.0)/10.0<<"%               \r";

        //if(H_S[i].xx<0 && H_S[i].yy<0 && H_S[i].zz<0 && vbead[i]) b[i]-=H_S[i].transpose()*grad_S[i];

        if( scatteredRfound ) {

            vec3D utemp=vec3D(0.0,0.0,0.0);

            for(auto itIT=ItransI[i].begin(); itIT!=ItransI[i].end(); itIT++) utemp += M.U[itIT->first] * itIT->second;

            btemp=ItransUfound[i]-utemp;
            btemp.z*=llambda_z;
            b[i]+=btemp;


        }else if(vbead[i]){

            //utemp=U_found[i];
            //utemp.z=llambda_z*utemp.z;
            btemp=(U_found[i]-M.U[i]);
            btemp.z*=llambda_z;
            b[i]+=btemp;

        }

        for(it=M.connections[i].begin(); it!=M.connections[i].end(); it++){

            j=*it;

            if(M.var[j]){

                b[i]+=M.K_glo[i][j]*f[j]*lambda2*localweight[j];
                if(uselaplace) b[i]+=M.K_glo[i][j]*llf[j]*lambda2*lagrain;

            }

        }

    }

}

void VirtualBeads::mulA(std::vector<double>& u, std::vector<double>& f, const FiniteBodyForces& M){

    int c1,c2,jend,j;

    vec3D ff,uu;
    mat3D A;

    for(c1=0; c1<M.N_c; c1++){

        ff=vec3D(0.0,0.0,0.0);

        jend=A_j[c1].size();

        for(j=0; j<jend; j++){

            c2=A_j[c1][j];

            uu=vec3D( u[3*c2],u[3*c2+1],u[3*c2+2] );

            A=A_sp[c1][j];

            ff+=A*uu;

        }

        f[3*c1]=ff.x;
        f[3*c1+1]=ff.y;
        f[3*c1+2]=ff.z;

    }

}



void VirtualBeads::recordRelaxationStatus(FiniteBodyForces& M, DRec3D& relrec){

    double Uf=0.0,ff=0.0,L=0.0,u2=0.0,uuf2=0.0,suuf=0.0;

    bool scatteredInput=bool(CFG["SCATTEREDRFOUND"]);

    double lagrain=double(CFG["REG_LAPLACEGRAIN"]);
    lagrain=lagrain*lagrain;

    int b,ii;
    int bcount=0;

    vec3D btemp;

    double lambda=double(CFG["LAMBDA"]);

    double lambda2=lambda;
    double llambda_z=double(CFG["REG_SIGMAZ"]);


    // Displ. part of penalty function

    if(!scatteredInput){

        for(b=0; b<M.N_c; b++) if(M.var[b] && vbead[b]){

            //std::cout<<b<<"\n";

            //Uf+=abs(U_found[b]);
            btemp=(U_found[b]-M.U[b]);
            btemp.z*=llambda_z;
            uuf2+=norm(btemp);
            suuf+=abs(btemp);
            bcount++;

        }

    }else{

        //std::vector< vec3D > u=std::vector< vec3D >();
        //u.assign(M.N_c,vec3D(0.0,0.0,0.0));

        for(auto i=0; i<U_found.size(); i++){

            bcount++;

            btemp=U_found[i]*-1.0;

            for(auto it=I[i].begin(); it!=I[i].end(); it++){

                btemp+=M.U[it->first]*it->second;

            }

            btemp.z*=llambda_z;

            uuf2+=norm(btemp);
            suuf+=abs(btemp);

        }

    }

    // Total displacements

    for(b=0; b<M.N_c; b++) if(M.var[b]) u2+=norm(M.U[b]);

    //Force part of penalty function

    std::vector< vec3D > f=std::vector< vec3D >();
    f.assign(M.N_c,vec3D(0.0,0.0,0.0));

    for(ii=0; ii<M.N_c; ii++) if(M.var[ii]) f[ii]=vec3D( M.f_glo[3*ii], M.f_glo[3*ii+1], M.f_glo[3*ii+2] );


    if(std::string(CFG["REGMETHOD"])=="laplace") {

        std::vector< vec3D > lf=std::vector< vec3D >();
        lf.assign(M.N_c,vec3D(0.0,0.0,0.0));

        mul(M.Laplace, f, lf);

        for(b=0; b<M.N_c; b++) if(M.var[b]) f[b]+=lf[b]*lagrain;

    }

    ff=0.0;

    for(b=0; b<M.N_c; b++) if(M.var[b])  ff+=norm(f[b])*localweight[b];

    L=lambda2*ff+uuf2;

    std::cout<<"|u-uf|^2 = "<<uuf2<<"\t\tperbead="<<suuf/(bcount+0.0)<<"   \n";
    std::cout<<"|w*f|^2  = "<<ff<<"\t\t|u|^2 = "<<u2<<"        \n";
    std::cout<<"L = |u-uf|^2 + lambda*|w*f|^2 = "<<L<<"          \n";

    relrec.record(L,uuf2,ff);

    std::string outdir=std::string(CFG["DATAOUT"]);
    std::string relrecname=outdir+std::string("/")+std::string(CFG["REG_RELREC"]);
    relrec.store(relrecname.c_str());

}


vec3D VirtualBeads::relax(FiniteBodyForces& M, bool local=false){

    double lagrain=double(CFG["REG_LAPLACEGRAIN"]);
    lagrain=lagrain*lagrain;

    DRec3D relrec=DRec3D();

    double lambda=double(CFG["LAMBDA"]);
    double lambda2=lambda;

    localweight.assign(M.N_c,1.0);
    //for(int c=0; c<M.N_c; c++) if(outofstack[c]) localweight[c]=double(CFG["LAMBDA2"]);

    std::cout<<"going to update glo f and K\n";

    M.updateGloFAndK();

    recordRelaxationStatus(M,relrec);

    int i_max=int(CFG["REG_ITERATIONS"]);

    std::cout<<"check before relax !       \n";

    for(int i=0; i<i_max; i++){

        if(std::string(CFG["REGMETHOD"])!=std::string("normal")) updateLocalWeigth(M);
        computeAAndb(M,1.0,lambda2);
        double uu=solve_CG(M,local);

        M.updateGloFAndK();

        std::cout<<"Round "<<i+1<<" |du|="<<uu<<"                 \n";
        recordRelaxationStatus(M,relrec);

        if(i>6){

            int s=relrec.data.size();

            double ufmean=0.0,Lmean=0.0,cc=0.0;

            for(int ii=(s-1); ii>s-6; ii--){

                ufmean+=relrec.data[ii][0];
                Lmean+=relrec.data[ii][1];
                cc=cc+1.0;

            }

            ufmean=ufmean/cc;
            Lmean=Lmean/cc;

            double ufstd=0.0,Lstd=0.0;
            cc=0.0;

            for(int ii=(s-1); ii>s-6; ii--){

                ufstd+=(ufmean-relrec.data[ii][0])*(ufmean-relrec.data[ii][0]);
                Lstd+=(Lmean-relrec.data[ii][1])*(Lmean-relrec.data[ii][1]);
                cc=cc+1.0;

            }

            ufstd=sqrt(ufstd)/cc;
            Lstd=sqrt(Lstd)/cc;

            if( (Lstd/Lmean)<double(CFG["REG_CONV_CRIT"]) ) i=i_max;

        }

    }

    return relrec.data.back();

}


double VirtualBeads::solve_CG(FiniteBodyForces& M, bool local=false){

    //b=f_glo
    //A=K_glo

    double tol=(M.N_c+0.0)*double(CFG["REG_SOLVER_PRECISION"]);

    int maxiter=25*floor(pow(M.N_c,0.33333)+0.5);

    double resid,alpha_0,beta_0,rho_0,rho1_0,alpha,rsnew;

    std::vector<double> pp, zz, qq, rr, kk, uu, ff,fff, Ap;

    pp.assign(3*M.N_c,0.0);
    Ap.assign(3*M.N_c,0.0);
    zz.assign(3*M.N_c,0.0);
    qq.assign(3*M.N_c,0.0);
    rr.assign(3*M.N_c,0.0);
    kk.assign(3*M.N_c,0.0);
    uu.assign(3*M.N_c,0.0);
    ff.assign(3*M.N_c,0.0);
    fff.assign(3*M.N_c,0.0);

    //plusminus(f_glo,1.0,f_ext,ff);

    int i;

    for(i=0; i<M.N_c; i++){

        ff[3*i]=b[i].x;
        ff[3*i+1]=b[i].y;
        ff[3*i+2]=b[i].z;

    }

    double normb=imul(ff,ff);

    double norm0=normb;
    //std::cout<<ff[0]<<" "<<ff[1]<<" "<<ff[2]<<" "<<ff[3]<<" "<<ff[4]<<" "<<ff[5]<<" "<<ff[6]<<std::endl;

    if(normb>tol){

        //std::cout<<normb<<" IS BIGGER THAN "<<tol<<"       \n";

        mulA(uu,kk,M);

        plusminus(ff,-1.0,kk,rr);

        if(normb==0.0) normb=1.0;

        //std::cout<<"computin resid fist stime"<<std::endl;

        setequal(rr,pp);

        resid=imul(pp,pp);

        int i=0;
        //std::cout<<i<<": "<<resid<<std::endl;
        for( i = 1; i <= maxiter; i++){

            //precond_J3(rr,zz,M);

            //rho_0=imul(rr,zz);

            mulA(pp,Ap,M);

            alpha=resid/imul(pp,Ap);

            plusminus(uu,alpha,pp,uu);
            plusminus(rr,-alpha,Ap,rr);

            rsnew=imul(rr,rr);

            if( rsnew<tol ) break;

            plusminus(rr,rsnew/resid,pp,pp);

            resid=rsnew;

            if( (i%100)==0 ) std::cout<<i<<": "<<resid<<" alpha="<<alpha<<"\r";


            /*

            //std::cout<<rho_0<<std::endl;

            if(i==1) setequal(zz,pp);
            else{
                beta_0=rho_0/rho1_0;
                plusminus(zz,beta_0,pp,pp);
            }

            mulA(pp,qq,M);

            alpha_0=rho_0/imul(pp,qq);

            plusminus(uu,alpha_0,pp,uu);
            plusminus(rr,-1.0*alpha_0,qq,rr);

            resid=imul(pp,pp)*alpha_0*alpha_0;

            if(resid <= tol) i=maxiter+1;
            else if( (i%100)==0 ){

                std::cout<<i<<": "<<resid<<" alpha="<<alpha_0<<std::endl;

            }

            rho1_0=rho_0;

            */

        }

        int c;

        double du=0.0;

        double stepper=double(CFG["REG_SOLVER_STEP"]);

        for(c=0; c<M.N_c; c++){

            M.U[c].x+=stepper*uu[3*c];
            M.U[c].y+=stepper*uu[3*c+1];
            M.U[c].z+=stepper*uu[3*c+2];

            du+=stepper*stepper*(uu[3*c]*uu[3*c] + uu[3*c+1]*uu[3*c+1] + uu[3*c+2]*uu[3*c+2]);

        }

        return sqrt(du/(M.N_c+0.0));

    }else{

        //std::cout<<normb<<"IS SMALLER THAN "<<tol<<"       \n";

        return 0.0;

    }



}

void VirtualBeads::storeUfound(std::string Uname, std::string Sname){

    DRec3D Urec=DRec3D();

    DRec Srec=DRec();

    int c;

    int cend=U_found.size();

    for(c=0; c<cend; c++){
        Urec.record(U_found[c]);
        Srec.record(S_0[c]);
    }

    Urec.store(Uname.c_str());
    Srec.store(Sname.c_str());


}

void VirtualBeads::storeRfound(std::string Rname){

    DRec3D Rrec=DRec3D();

    int c;

    int cend=R_found.size();

    for(c=0; c<cend; c++){
        Rrec.record(R_found[c]);
    }

    Rrec.store(Rname.c_str());

}

void VirtualBeads::storeLocalweights(std::string wname){

    DRec Wrec=DRec();

    int c;

    int cend=localweight.size();

    for(c=0; c<cend; c++){

        Wrec.record(localweight[c]);
    }

    Wrec.store(wname.c_str());

}

