extern config CFG;

void saveEpsilon(std::vector<double> epsilon, const char* fname){

    DRecXY epsrec=DRecXY();

    int imax=ceil( (double(CFG["EPSMAX"])+1.0)/double(CFG["EPSSTEP"]) );

    double lambda;

    for(int i=0; i<imax; i++){

        lambda=(i*double(CFG["EPSSTEP"]))-1.0;

        epsrec.record(lambda,epsilon[i]);

    }

    epsrec.store(fname);

}

std::vector<double> buildEpsilon(std::vector<double>& epsilon, double k1, double k2, double ds0, double s1, double ds1){

    //std::vector<double> epsilon;
    std::vector<double> epsbar;
    std::vector<double> epsbarbar;

    int imax=ceil( (double(CFG["EPSMAX"])+1.0)/double(CFG["EPSSTEP"]) );

    epsilon.assign(imax,0.0);
    epsbar.assign(imax,0.0);
    epsbarbar.assign(imax,0.0);

    double lambda;

    int i;

    double k00=(k2-k1)*(erf(((s1)/ds1))*0.5+0.5)+k1;

    double ds00=k00/ds0;

    std::cout<<"k00: "<<k00<<" ds00: "<<ds00<<"\n";

    for(i=0; i<imax; i++){

        lambda=(i*double(CFG["EPSSTEP"]))-1.0;

        epsbarbar[i]=0;

        if(lambda>0) {

            //epsbarbar[i]+=k1-k0;

            epsbarbar[i]+=(k1);

            epsbarbar[i]+=(k2-k1)*(erf(((lambda-s1)/ds1))*0.5+0.5);

        }else{

            //epsbarbar[i]+=(k1)*exp( ((lambda)/ds0) );

            epsbarbar[i]+=(k00)*exp( ((lambda)*ds00) );

        }

        //if(lambda<1.0) std::cout<<i<<": lambda="<<lambda<<"   epsbarbar="<<epsbarbar[i]<<std::endl;

    }

    double sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbarbar[i]*double(CFG["EPSSTEP"]);
        epsbar[i]=sum;
    }

    int imid=ceil(1.0/double(CFG["EPSSTEP"]));

    double off=epsbar[imid];

    for(i=0; i<imax; i++) epsbar[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsbar="<<epsbar[i]<<std::endl;

    sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbar[i]*double(CFG["EPSSTEP"]);
        epsilon[i]=sum;
    }

    off=epsilon[imid];
    for(i=0; i<imax; i++) epsilon[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsilon="<<epsilon[i]<<std::endl;

    saveEpsilon(epsbarbar,"epsbarbar.dat");
    saveEpsilon(epsbar,"epsbar.dat");

    return epsilon;

}

std::vector<double> buildEpsilon2(std::vector<double>& epsilon, double u_b, double A, double kappa){

    //std::vector<double> epsilon;
    std::vector<double> epsbar;
    //std::vector<double> epsbarbar;

    int imax=ceil( (double(CFG["EPSMAX"])+1.0)/double(CFG["EPSSTEP"]) );

    epsilon.assign(imax,0.0);
    epsbar.assign(imax,0.0);
    //epsbarbar.assign(imax,0.0);

    double lambda;

    int imid=ceil(1.0/double(CFG["EPSSTEP"]));
    int i;

    for(i=imid; i<imax; i++){

        lambda=(i*double(CFG["EPSSTEP"]))-1.0;

        epsbar[i]= -(A - lambda + kappa*u_b - pow(A*A + 2*A*kappa*u_b - 2*A*lambda + kappa*kappa*u_b*u_b + 2*kappa*lambda*u_b + lambda*lambda,0.5))/(2*kappa);

    }


    double kmid=(epsbar[imid+1]-epsbar[imid])/double(CFG["EPSSTEP"]);

    std::cout<<"kmid = "<<kmid<<"\n";

    for(i=0; i<imid; i++){

        lambda=(i*double(CFG["EPSSTEP"]))-1.0;

        epsbar[i]= 38.0*exp(38.0*lambda/kmid)-38.0;
    }

    double sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbar[i]*double(CFG["EPSSTEP"]);
        epsilon[i]=sum;
    }

    double off=epsilon[imid];
    for(i=0; i<imax; i++) epsilon[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsilon="<<epsilon[i]<<std::endl;

    saveEpsilon(epsbar,"epsbar.dat");

    return epsilon;

}

void buildEpsilon3(std::vector<double>& epsilon, std::vector<double>& epsbar, std::vector<double>& epsbarbar, double k1, double ds0, double s1, double ds1){

    std::cout<<k1<<" "<<ds0<<" "<<s1<<" "<<ds1<<std::endl;

    int imax=ceil( (double(CFG["EPSMAX"])+1.0)/double(CFG["EPSSTEP"]) );

    epsilon.assign(imax,0.0);
    epsbar.assign(imax,0.0);
    epsbarbar.assign(imax,0.0);

    double lambda;

    int i;

    for(i=0; i<imax; i++){

        lambda=(i*double(CFG["EPSSTEP"]))-1.0;

        epsbarbar[i]=0;

        if(lambda>0) {

            //epsbarbar[i]+=k1-k0;

            epsbarbar[i]+=(k1);

            if(lambda>s1 && ds1>0.0){

                epsbarbar[i]+=(k1)*(exp((lambda-s1)/ds1)-1.0);

            }

        }else{

            //epsbarbar[i]+=(k1)*exp( ((lambda)/ds0) );

            if(ds0!=0.0) epsbarbar[i]+=(k1)*exp( ((lambda)/ds0) );
            else epsbarbar[i]+=(k1);

        }

        if(epsbarbar[i]>10e10) epsbarbar[i]=10e10;

        //if(lambda<1.0) std::cout<<i<<": lambda="<<lambda<<"   epsbarbar="<<epsbarbar[i]<<std::endl;

    }

    double sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbarbar[i]*double(CFG["EPSSTEP"]);
        epsbar[i]=sum;
    }

    int imid=ceil(1.0/double(CFG["EPSSTEP"]));

    double off=epsbar[imid];

    for(i=0; i<imax; i++) epsbar[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsbar="<<epsbar[i]<<std::endl;

    sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbar[i]*double(CFG["EPSSTEP"]);
        epsilon[i]=sum;
    }

    off=epsilon[imid];
    for(i=0; i<imax; i++) epsilon[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsilon="<<epsilon[i]<<std::endl;

    //saveEpsilon(epsbarbar,"epsbarbar.dat");
    //saveEpsilon(epsbar,"epsbar.dat");
    //saveEpsilon(epsilon,"epsilon.dat");

}

void buildEpsilon3new(std::vector<double>& epsilon, std::vector<double>& epsbar, std::vector<double>& epsbarbar, double k1, double ds0, double s1, double ds1){

    //std::vector<double> epsilon;
    //std::vector<double> epsbar;
    //std::vector<double> epsbarbar;

    int imax=ceil( (double(CFG["EPSMAX"])+1.0)/double(CFG["EPSSTEP"]) );

    epsilon.assign(imax,0.0);
    epsbar.assign(imax,0.0);
    epsbarbar.assign(imax,0.0);

    double lambda;

    int i;

    for(i=0; i<imax; i++){

        lambda=(i*double(CFG["EPSSTEP"]))-1.0;

        epsbarbar[i]=0;

        if(lambda>0) {

            //epsbarbar[i]+=k1-k0;

            epsbarbar[i]+=(k1);

            if(lambda>s1){

                epsbarbar[i]+=(k1)*(exp((lambda-s1)*ds1)-1.0);
                //epsbarbar[i]+=(k1)*(lambda-s1)*ds1;
                //epsbarbar[i]+=(k1)*(lambda-s1)*ds1*(lambda-s1)*ds1*0.5;
                //epsbarbar[i]+=(k1)*(lambda-s1)*ds1*(lambda-s1)*ds1*(lambda-s1)*ds1/6.0;
                //epsbarbar[i]+=(k1)*(lambda-s1)*ds1*(lambda-s1)*ds1*(lambda-s1)*ds1*(lambda-s1)*ds1/24.0;

            }

        }else{

            //epsbarbar[i]+=(k1)*exp( ((lambda)/ds0) );

            epsbarbar[i]+=(k1)*exp( ((lambda)/ds0) );

        }

        if(epsbarbar[i]>10e7) epsbarbar[i]=10e7;

        //if(lambda<1.0) std::cout<<i<<": lambda="<<lambda<<"   epsbarbar="<<epsbarbar[i]<<std::endl;

    }

    double sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbarbar[i]*double(CFG["EPSSTEP"]);
        epsbar[i]=sum;
    }

    int imid=ceil(1.0/double(CFG["EPSSTEP"]));

    double off=epsbar[imid];

    for(i=0; i<imax; i++) epsbar[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsbar="<<epsbar[i]<<std::endl;

    sum=0.0;

    for(i=0; i<imax; i++){
        sum+=epsbar[i]*double(CFG["EPSSTEP"]);
        epsilon[i]=sum;
    }

    off=epsilon[imid];
    for(i=0; i<imax; i++) epsilon[i]-=off;

    //for(i=0; i<(imax/5*2); i=i+100) std::cout<<i<<": lambda="<<(i*double(CFG["EPSSTEP"]))-1.0<<"   epsilon="<<epsilon[i]<<std::endl;

    //saveEpsilon(epsbarbar,"epsbarbar.dat");
    //saveEpsilon(epsbar,"epsbar.dat");
    saveEpsilon(epsilon,"epsilon.dat");

}

