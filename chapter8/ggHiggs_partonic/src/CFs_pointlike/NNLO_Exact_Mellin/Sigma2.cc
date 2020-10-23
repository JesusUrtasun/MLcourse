#include "Sigma2.hh"

dcomplex sigma_NNLO_gg_N(dcomplex NN){
  dcomplex N=NN-1.;
  dcomplex N1=N+1.;
  dcomplex N2=N+2.;
  dcomplex N3=N+3.;
  dcomplex N4=N+4.;
  MellinFunc MF;
  long double zeta2=gsl_sf_zeta(2.);
  long double zeta3=gsl_sf_zeta(3.);
  long double zeta2q=zeta2*zeta2;
  long double W1=11./2.;
  
  
  dcomplex sigmaA_gg_N=(837./16.+67./2.*zeta2-9./20.*zeta2q-165./4.*zeta3)+(-101./3.+33.*zeta2+351.*zeta3/2.)*MF.D0(N)+(67.-90.*zeta2)*MF.D1(N)
  -33.*MF.D2(N)+72.*MF.D3(N)+9.*(5.*MF.Li3zplus(N)+14.*MF.Li3zplus(N1)+52.*MF.Li3zplus(N2)+32.*MF.Li3zplus(N3)-7.*MF.Li3zplus(N4))
  -18.*(MF.S12z2plus(N)+2.*MF.S12z2plus(N1)+3.*MF.S12z2plus(N2)+2.*MF.S12z2plus(N3)+MF.S12z2plus(N4))
  +9.*(4.*MF.S12mzplus(N4)+8.*MF.S12mzplus(N3)+21.*MF.S12mzplus(N2)+14.*MF.S12mzplus(N1)+7.*MF.S12mzplus(N))
  +0.5*9.*(10.*MF.S12zplus(N)-50.*MF.S12zplus(N1)-103.*MF.S12zplus(N2)-46.*MF.S12zplus(N3)+5.*MF.S12zplus(N4))
  -9./2.*(8.*MF.Li3mzplus(N4)+8.*MF.Li3mzplus(N3)-3.*MF.Li3mzplus(N2)-2.*MF.Li3mzplus(N1)-MF.Li3mzplus(N))
  -9./2.*(16.*MF.Li2zLogzplusminus(N)+13.*MF.Li2zLogzplusminus(N+5.)-40.*MF.Li2zLogzplusminus(N3)-67.*MF.Li2zLogzplusminus(N4)+64.*MF.Li2zLogzplusminus(N2)
  +36.*MF.Li2zLogzplusminus(N1))+9./2.*(2.*MF.Li2mzLogzplus(N4)-15.*MF.Li2mzLogzplus(N2)-10.*MF.Li2mzLogzplus(N1)-5.*MF.Li2mzLogzplus(N))
  -9./4.*(59.*MF.LogzLogminus2minus(N)+177.*MF.LogzLogminus2minus(N2)-116.*MF.LogzLogminus2minus(N3)+59.*MF.LogzLogminus2minus(N4)-118.*MF.LogzLogminus2minus(N1))
  +27.*(3.*MF.Li2mzLogplusplus(N2)+2.*MF.Li2mzLogplusplus(N1)+MF.Li2mzLogplusplus(N))
  +9.*(6.*MF.Logz2Logminusminus(N)-11.*MF.Logz2Logminusminus(N3)+18.*MF.Logz2Logminusminus(N2)-12.*MF.Logz2Logminusminus(N1)+6.*MF.Logz2Logminusminus(N4))
  +9./2.*(2.*MF.Li2zLogminus(N)-4.*MF.Li2zLogminus(N1)+5.*MF.Li2zLogminus(N2)-3.*MF.Li2zLogminus(N3))
  -3./2.*(7.*MF.Logz3minusplus(N1)-7.*MF.Logz3minusplus(N3)+4.*MF.Logz3minusplus(N)+18.*MF.Logz3minusplus(N2)-17.*MF.Logz3minusplus(N4)+9.*MF.Logz3minusplus(N+5.))
  +9./2.*(8.*MF.Logplusplus(N4)+16.*MF.Logplusplus(N3)+33.*MF.Logplusplus(N2)+22.*MF.Logplusplus(N1)+11.*MF.Logplusplus(N))*zeta2
  -36.*(MF.Li2zLogplusplus(N)+2.*MF.Li2zLogplusplus(N1)+3.*MF.Li2zLogplusplus(N2)+2.*MF.Li2zLogplusplus(N3)+MF.Li2zLogplusplus(N4))
  -9./4.*(4.*MF.Logz2Logplusplus(N4)+8.*MF.Logz2Logplusplus(N3)+27.*MF.Logz2Logplusplus(N2)+18.*MF.Logz2Logplusplus(N1)+9.*MF.Logz2Logplusplus(N))
  +(-21.*MF.LogzLogplus(N)+63./2.*MF.LogzLogplus(N2)-18.*MF.LogzLogplus(N1)+33./2.*MF.LogzLogplus(N3))
  +27./2.*(3.*MF.LogzLogplus2plus(N2)+2.*MF.LogzLogplus2plus(N1)+MF.LogzLogplus2plus(N))+3./4.*(278.*MF.Li2z(N)-116.*MF.Li2z(N1)-137.*MF.Li2z(N2)+143.*MF.Li2z(N3))
  +(-21.*MF.Li2mz(N)+63./2.*MF.Li2mz(N2)-18.*MF.Li2mz(N1)+33./2.*MF.Li2mz(N3))
  +(-2559./4.*MF.Logminus(N3)+1079./2.*MF.Logminus(N2)-2687./4.*MF.Logminus(N1)+2027./4.*MF.Logminus(N))
  -3./8.*(374.*MF.Logz2minus(N4)-389.*MF.Logz2minus(N1)+154.*MF.Logz2minus(N)+699.*MF.Logz2minus(N2)-827.*MF.Logz2minus(N3))
  +(330.*MF.Logminus2(N3)-348.*MF.Logminus2(N2)+381.*MF.Logminus2(N1)-297.*MF.Logminus2(N))
  +3./4.*(-1180.*MF.LogzLogminusminus(N3)+641.*MF.LogzLogminusminus(N)-1238.*MF.LogzLogminusminus(N1)+1227.*MF.LogzLogminusminus(N2)+605.*MF.LogzLogminusminus(N4))
  -72.*(2.*MF.Logminus3(N1)-MF.Logminus3(N2)+MF.Logminus3(N3))
  -1./8.*(4318.*MF.Logzminus(N4)-6955.*MF.Logzminus(N3)+6447.*MF.Logzminus(N2)-5611.*MF.Logzminus(N1)+2333.*MF.Logzminus(N))
  -3./4.*(-27./N+173./(N1)-391./(N2)+495./N3)*zeta2
  +9.*(6.*MF.Logzminusplus(N1)+18.*MF.Logzminusplus(N2)+2.*MF.Logzminusplus(N)+10.*MF.Logzminusplus(N+5.)-6.*MF.Logzminusplus(N3)-19.*MF.Logzminusplus(N4))*zeta2
  +9./2.*(-2.*MF.Logminus(N)+44.*MF.Logminus(N1)-25.*MF.Logminus(N2)+23.*MF.Logminus(N3))*zeta2
  -0.5*9.*(17.*MF.plus(N)+50.*MF.plus(N1)+31.*MF.plus(N2)+18.*MF.plus(N3)+33.*MF.plus(N4))*zeta3+7539./16./(N3)-24107./48./N2+22879./48./(N1)-18157./48./N
  +162./2.*MF.Li3zregminus(N)+9./2.*MF.S12zregminus(N)+33./4.*MF.Li2zregminus(N)+9./2.*MF.Li2zregLogminusminus(N)
  -W1*(-12.*(-MF.Logminus(N2)+MF.Logminus(N3)+2.*MF.Logminus(N1))-6.*(MF.Logzminus(N)-2.*MF.Logzminus(N1)+3.*MF.Logzminus(N2)-2.*MF.Logzminus(N3)+MF.Logzminus(N4))
  -33./(N*(6.+11.*N+6.*N*N+N*N*N)));
  
  dcomplex sigmaF_gg_N=1./36.*(-247.-60.*zeta2+30.*zeta3)+(14./9.-2.*zeta2)*MF.D0(N)-10./3.*MF.D1(N)+2.*MF.D2(N)
  +(31./6.*MF.S12z(N1)+1./6.*MF.S12z(N)+65./12.*MF.S12z(N2))+(-31./12.*MF.Li3z(N2)+1./6.*MF.Li3z(N)-17./6.*MF.Li3z(N1))
  +(47./12.*MF.Li2zLogz(N2)+25./6.*MF.Li2zLogz(N1)-1./6.*MF.Li2zLogz(N))+(-1./12.*MF.Logminus(N2)+1./6.*MF.Logminus(N1)-1./6.*MF.Logminus(N))*zeta2
  -4.*zeta2*(MF.Logz(N1)+MF.Logz(N2))+(-1./6.*MF.Li2zLogminus(N1)+1./6.*MF.Li2zLogminus(N)+1./12.*MF.Li2zLogminus(N2))
  +(1./12.*MF.LogzLogminus2(N)-1./12.*MF.LogzLogminus2(N1)+1./24.*MF.LogzLogminus2(N2))-(1./12.*MF.Logz2Logminus(N)-1./12.*MF.Logz2Logminus(N1)+1./24.*MF.Logz2Logminus(N2))
  +5./9.*(MF.Logz3(N1)+MF.Logz3(N2))+(-17./6./(N2)-7./3./N1-1./3./N)*zeta3+(-34./9.*MF.Logminus2(N3)+2./3.*MF.Logminus2(N2)-8./3.*MF.Logminus2(N1)+16./9.*MF.Logminus2(N))
  -(-34./9./N3+2./3./N2-8./3./N1+16./9./N)*zeta2-2./9.*(21.*MF.LogzLogminusminus(N2)+7.*MF.LogzLogminusminus(N1)+25.*MF.LogzLogminusminus(N4)+17.*MF.LogzLogminusminus(N)
  -61.*MF.LogzLogminusminus(N3))+(785./54.*MF.Logminus(N3)-83./36.*MF.Logminus(N2)+49./18.*MF.Logminus(N1)-461./54.*MF.Logminus(N))
  +1./72.*(-351.*MF.Logz2minus(N3)+117.*MF.Logz2minus(N2)+68.*MF.Logz2minus(N)+132.*MF.Logz2minus(N4)+52.*MF.Logz2minus(N1))
  +1./36.*(227.*MF.Li2minusminus(N3)+68.*MF.Li2minusminus(N)+4.*MF.Li2minusminus(N4)-302.*MF.Li2minusminus(N1)+21.*MF.Li2minusminus(N2))
  +1./216.*(333.*MF.Logzminus(N2)+2384.*MF.Logzminus(N4)-598.*MF.Logzminus(N1)-3041.*MF.Logzminus(N3)+1282.*MF.Logzminus(N))
  -8887./648./N3+1267./432./N2-497./216./N1+12923./1296./N;
  
  long double Nf=5.;
  return(sigmaA_gg_N+Nf*sigmaF_gg_N);
}

dcomplex sigma_NNLO_qg_N(dcomplex NN){
  dcomplex N=NN-1.;
  dcomplex N1=N+1.;
  dcomplex N2=N+2.;
  dcomplex N3=N+3.;
  //dcomplex N4=N+4.;
  MellinFunc MF;
  long double zeta2=gsl_sf_zeta(2.);
  long double zeta3=gsl_sf_zeta(3.);
  //long double zeta2q=zeta2*zeta2;
  long double W1=11./2.;
  
  
  dcomplex sigmaA_qg_N=(170./3.*MF.Li3z(N1)+338./9.*MF.Li3z(N)+119./3.*MF.Li3z(N2))+(4.*MF.Li3mz(N1)+4.*MF.Li3mz(N)+2.*MF.Li3mz(N2))
  +(16.*MF.S12mz(N)+8.*MF.S12mz(N2)+16.*MF.S12mz(N1))+(-614./9.*MF.S12z(N1)-269./9.*MF.S12z(N2)-74./9.*MF.S12z(N))
  +(-2.*MF.S12z2(N2)-4.*MF.S12z2(N1)-4.*MF.S12z2(N))+(367./27.*MF.Logminus3(N)+367./54.*MF.Logminus3(N2)-367./27.*MF.Logminus3(N1))
  +(2.*MF.Li2zLogminus(N)+MF.Li2zLogminus(N2)-2.*MF.Li2zLogminus(N1))-(446./9.*MF.Li2zLogz(N1)+214./9.*MF.Li2zLogz(N)+281./9.*MF.Li2zLogz(N2))
  -(8.*MF.Li2zLogplus(N)+4.*MF.Li2zLogplus(N2)+8.*MF.Li2zLogplus(N1))+(8.*MF.Li2mzLogplus(N)+8.*MF.Li2mzLogplus(N1)+4.*MF.Li2mzLogplus(N2))
  -(8.*MF.Li2mzLogz(N)+8.*MF.Li2mzLogz(N1)+4.*MF.Li2mzLogz(N2))+(-115./9.*MF.LogzLogminus2(N2)+230./9.*MF.LogzLogminus2(N1)-230./9.*MF.LogzLogminus2(N))
  +(107./9.*MF.Logz2Logminus(N)+107./18.*MF.Logz2Logminus(N2)-107./9.*MF.Logz2Logminus(N1))+(-145./54.*MF.Logz3(N2)-71./27.*MF.Logz3(N1)-2.*MF.Logz3(N))
  +(-3.*MF.Logz2Logplus(N2)-6.*MF.Logz2Logplus(N)-6.*MF.Logz2Logplus(N1))+(4.*MF.LogzLogplus2(N1)+4.*MF.LogzLogplus2(N)+2.*MF.LogzLogplus2(N2))
  +(-4./27.*MF.Li2mz(N3)-74./9.*MF.Li2mz(N1)-11./9.*MF.Li2mz(N2)-166./27.*MF.Li2mz(N))+(2605./54.*MF.Li2z(N)-146./9.*MF.Li2z(N1)+74./27.*MF.Li2z(N3)-79./6.*MF.Li2z(N2))
  +(1139./18.*MF.Logminus2(N1)+37./12.*MF.Logminus2(N2)+8.*MF.Logminus2(N3)-72.*MF.Logminus2(N))
  +(-121./18.*MF.LogzLogminus(N2)-326./27.*MF.LogzLogminus(N3)-826./9.*MF.LogzLogminus(N1)+5935./54.*MF.LogzLogminus(N))
  +(113./27.*MF.Logz2(N3)+244./9.*MF.Logz2(N1)-13./3*MF.Logz2(N2)-31./2.*MF.Logz2(N))+(-4./27.*MF.LogzLogplus(N3)-74./9.*MF.LogzLogplus(N1)-11./9.*MF.LogzLogplus(N2)-166./27.*MF.LogzLogplus(N))
  +zeta2*(-59./9.*MF.Logminus(N2)+118./9.*MF.Logminus(N1)-118./9.*MF.Logminus(N))+zeta2*(140./9.*MF.Logz(N1)+128./9.*MF.Logz(N2)+52./9.*MF.Logz(N))
  +zeta2*(12.*MF.Logplus(N)+12.*MF.Logplus(N1)+6.*MF.Logplus(N2))+(-392./81.*MF.Logminus(N3)-49./3.*MF.Logminus(N2)+23671./162.*MF.Logminus(N)-106.*MF.Logminus(N1))
  +(1985./108.*MF.Logz(N2)+800./9.*MF.Logz(N1)-12209./162.*MF.Logz(N)+616./81.*MF.Logz(N3))
  +(-292./27./(N3)-82./3./N1+16./3./N2+221./27./N)*zeta2+(-18./N1+10./N2+92./9./N)*zeta3-210115./1944./N+1537./486./N3+16465./162./N1+2393./648./N2
  -W1*(-2./3.*(2.*MF.Logz(N)-2.*MF.Logz(N1)+MF.Logz(N2))+4./3.*(2.*MF.Logminus(N)-2.*MF.Logminus(N1)+MF.Logminus(N2))-1./N+2./N1-1./3./N2);
  
  dcomplex sigmaF_qg_N=(1./18.*MF.Logminus2(N2)-1./9.*MF.Logminus2(N1)+1./9.*MF.Logminus2(N))+(-38./27.*MF.Logz(N1)+19./27.*MF.Logz(N2)+29./27.*MF.Logz(N))
  -209./81./N1+265./162./N+(-4./9.*MF.LogzLogminus(N)+4./9.*MF.LogzLogminus(N1)-2./9.*MF.LogzLogminus(N2))
  +(-MF.Logminus(N2)+16./9.*MF.Logminus(N1)-13./9.*MF.Logminus(N))+179./162./N2+(1./9.*MF.Logz2(N2)-2./9.*MF.Logz2(N1)+2./9.*MF.Logz2(N));
  
  long double Nf=5.;
  return(sigmaA_qg_N+Nf*sigmaF_qg_N);
}

dcomplex sigma_NNLO_qqbar_N(dcomplex NN){
  dcomplex N=NN-1.;
  dcomplex N1=N+1.;
  dcomplex N2=N+2.;
  dcomplex N3=N+3.;
  //dcomplex N4=N+4.;
  MellinFunc MF;
  long double zeta2=gsl_sf_zeta(2.);
  long double zeta3=gsl_sf_zeta(3.);
  //long double zeta2q=zeta2*zeta2;
  long double W1=11./2.;
  
  dcomplex sigmaA_qqbar_N=(-16./9.*MF.Li3mz(N)-16./9.*MF.Li3mz(N1)-8./9.*MF.Li3mz(N2))+(-16./27.*MF.S12mz(N2)-32./27.*MF.S12mz(N)-32./27.*MF.S12mz(N1))
  +32./9.*(MF.Li3z(N2)+4.*MF.Li3z(N1)+4.*MF.Li3z(N))-32./9.*(MF.S12z(N2)+4.*MF.S12z(N1)+4.*MF.S12z(N))-4./27.*(MF.Logz3(N2)+4.*MF.Logz3(N1)+4.*MF.Logz3(N))
  +4./9.*(2.*MF.Logz2Logplus(N)+2.*MF.Logz2Logplus(N1)+MF.Logz2Logplus(N2))-8./27.*(2.*MF.LogzLogplus2(N)+2.*MF.LogzLogplus2(N1)+MF.LogzLogplus2(N2))
  -8./3.*(MF.Li2zLogz(N2)+4.*MF.Li2zLogz(N1)+4.*MF.Li2zLogz(N))+8./9.*(2.*MF.Li2mzLogz(N)+2.*MF.Li2mzLogz(N1)+MF.Li2mzLogz(N2))
  -16./27.*(2.*MF.Li2mzLogplus(N)+2.*MF.Li2mzLogplus(N1)+MF.Li2mzLogplus(N2))+32./81.*(-14.*MF.Logminus2(N)-21.*MF.Logminus2(N1)+48.*MF.Logminus2(N2)-13.*MF.Logminus2(N3))
  -16./81.*(-44.*MF.LogzLogminus(N)-57.*MF.LogzLogminus(N1)+138.*MF.LogzLogminus(N2)-37.*MF.LogzLogminus(N3))
  -8./81.*(44.*MF.Logz2(N3)+39.*MF.Logz2(N1)-81.*MF.Logz2(N2)+27.*MF.Logz2(N))+16./27.*(MF.LogzLogplus(N2)+6.*MF.LogzLogplus(N3)+2.*MF.LogzLogplus(N1))
  +8./81.*(42.*MF.Li2z(N1)-87.*MF.Li2z(N2)+12.*MF.Li2z(N)+10.*MF.Li2z(N3))+16./27.*(MF.Li2mz(N2)+6.*MF.Li2mz(N3)+2.*MF.Li2mz(N1))
  -4./81.*(-75.*MF.Logminus(N)-892.*MF.Logminus(N1)+1351.*MF.Logminus(N2)-384.*MF.Logminus(N3))+(-16./27./N2-32./27./N-32./27./N1)*zeta3
  +8./9*(MF.Logz(N2)+4.*MF.Logz(N1)+4.*MF.Logz(N))*zeta2+4222./81.*MF.Logz(N2)-2896./81.*MF.Logz(N1)-512./27.*MF.Logz(N3)-10./3.*MF.Logz(N)
  -8./27.*(2.*MF.Logplus(N)+2.*MF.Logplus(N1)+MF.Logplus(N2))*zeta2+(752./81./N3-544./27./N2+80./81./N+400./27./N1)*zeta2
  +4./81.*(373./N-2298./N1+2708./N2-783/N3)-W1*64./(9.*N*(6.+11.*N+6.*N*N+N*N*N));
  
  dcomplex sigmaF_qqbar_N=32./81.*(MF.Logminus(N)-3.*MF.Logminus(N1)+3.*MF.Logminus(N2)-MF.Logminus(N3))
  +(-64./27.*MF.Logz(N2)+64./81.*MF.Logz(N3)-16./27.*MF.Logz(N)+80./27.*MF.Logz(N1))-8./243.*(23./N-111./(N1)+129./N2-41./N3);
  
  long double Nf=5.;
  return(sigmaA_qqbar_N+Nf*sigmaF_qqbar_N);
}

dcomplex sigma_NNLO_qq_N(dcomplex NN){
  dcomplex N=NN-1.;
  dcomplex N1=N+1.;
  dcomplex N2=N+2.;
  //dcomplex N3=N+3.;
  //dcomplex N4=N+4.;
  MellinFunc MF;
  long double zeta2=gsl_sf_zeta(2.);
  long double zeta3=gsl_sf_zeta(3.);
  //long double zeta2q=zeta2*zeta2;
  //long double W1=11./2.;
  
  dcomplex sigmaA_qq_N=(368./27.*MF.Li3z(N1)+104./27.*MF.Li3z(N2)+400./27.*MF.Li3z(N))-32./9.*(MF.S12z(N2)+4.*MF.S12z(N1)+4.*MF.S12z(N))
  -4./27.*(2.*MF.Logz2Logminus(N)+MF.Logz2Logminus(N2)-2.*MF.Logz2Logminus(N1))-4./27.*(MF.Logz3(N2)+4.*MF.Logz3(N1)+4.*MF.Logz3(N))
  -16./27.*(19.*MF.Li2zLogz(N)+5.*MF.Li2zLogz(N2)+17.*MF.Li2zLogz(N1))-32./9.*(-MF.Logminus2(N2)-2.*MF.Logminus2(N1)+3.*MF.Logminus2(N))
  +16./3.*(-MF.LogzLogminus(N2)-2.*MF.LogzLogminus(N1)+3.*MF.LogzLogminus(N))+4./27.*(26.*MF.Logz2(N1)-18.*MF.Logz2(N)+9.*MF.Logz2(N2))
  -8./9.*(-6.*MF.Li2z(N)+MF.Li2z(N2)+4.*MF.Li2z(N1))+4./3.*(-5.*MF.Logminus(N2)-12.*MF.Logminus(N1)+17.*MF.Logminus(N))
  +(8./9.*(MF.Logz(N2)+4.*MF.Logz(N1)+4.*MF.Logz(N))*zeta2-118./9.*MF.Logz(N)+248./27.*MF.Logz(N1)+46./9.*MF.Logz(N2))
  +(-8./27./N2-16./27./N+16./27./N1)*zeta3+(16./3./N-32./9./N1-8./3./N2)*zeta2-4./27.*(-27./N2-133./N1+160./N);
  
  dcomplex sigmaF_qq_N=0.;
  
  long double Nf=5.;
  return(sigmaA_qq_N+Nf*sigmaF_qq_N);
}

dcomplex sigma_NNLO_qqprime_N(dcomplex NN){
  dcomplex N=NN-1.;
  dcomplex N1=N+1.;
  dcomplex N2=N+2.;
  //dcomplex N3=N+3.;
  //dcomplex N4=N+4.;
  MellinFunc MF;
  long double zeta2=gsl_sf_zeta(2.);
  //long double zeta3=gsl_sf_zeta(3.);
  //long double zeta2q=zeta2*zeta2;
  //long double W1=11./2.;

  dcomplex sigmaA_qqprime_N=32./9.*(MF.Li3z(N2)+4.*MF.Li3z(N1)+4.*MF.Li3z(N))-32./9.*(MF.S12z(N2)+4.*MF.S12z(N1)+4.*MF.S12z(N))
  -8./3.*(MF.Li2zLogz(N2)+4.*MF.Li2zLogz(N1)+4.*MF.Li2zLogz(N))-4./27.*(MF.Logz3(N2)+4.*MF.Logz3(N1)+4.*MF.Logz3(N))-8./9.*(4.*MF.Li2z(N1)-6.*MF.Li2z(N)+MF.Li2z(N2))
  -32./9.*(-MF.Logminus2(N2)-2.*MF.Logminus2(N1)+3.*MF.Logminus2(N))+16./3.*(-MF.LogzLogminus(N2)-2.*MF.LogzLogminus(N1)+3.*MF.LogzLogminus(N))
  +8./9.*(MF.Logz2(N2)+4.*MF.Logz2(N1)-3.*MF.Logz2(N))+8./9.*zeta2*(MF.Logz(N2)+4.*MF.Logz(N1)+4.*MF.Logz(N))
  +4./3.*(-5.*MF.Logminus(N2)-12.*MF.Logminus(N1)+17.*MF.Logminus(N))+2./9.*(29.*MF.Logz(N2)+44.*MF.Logz(N1)-59.*MF.Logz(N))
  +(16./3./N-32./9./N1-8./3./N2)*zeta2-2./9.*(-11./N2-94./N1+105./N);
  
  dcomplex sigmaF_qqprime_N=0.;
  
  long double Nf=5.;
  return(sigmaA_qqprime_N+Nf*sigmaF_qqprime_N);
}





//Table of Important Mellin of logarithms and harmonic polylogarithms
MellinFunc::MellinFunc():H(false,false,false){
  zeta2=gsl_sf_zeta(2.);
  zeta3=gsl_sf_zeta(3.);
  Li4=0.5174790616738993863;
  log2=std::log(2);
  log2q=log2*log2;
  log2c=log2*log2*log2;
  zeta2q=zeta2*zeta2;
  EulerGamma=0.577215664901532860606512090082; 
}

dcomplex MellinFunc::Li3zplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-pow(-1.,n)*(H.HS(-3,1,N-1.)-zeta2*H.HS(-2,N-1.)+zeta3*H.HS(-1,N-1.)+3./5.*zeta2q-2.*Li4-3./4.*zeta3*log2+0.5*zeta2*log2q-1./12.*log2q*log2q));
}
dcomplex MellinFunc::S12z2plus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(37./20.*zeta2q+2.*(H.HS(-2,-1,-1,N-1.)+H.HS(-2,1,1,N-1.)+H.HS(2,-1,1,N-1.)+H.HS(2,1,-1,N-1.))
  +2.*(H.HS(-2,-1,N-1.)-H.HS(-2,1,N-1.)-H.HS(2,-1,N-1.)+H.HS(2,1,N-1))*log2-4.*Li4-H.HS(-1,N-1.)*zeta3-H.HS(-2,N-1.)*(zeta2-2.*log2q)
  +H.HS(2,N-1.)*(zeta2-2.*log2q)+zeta2*log2q-1./6.*log2q*log2q-9./2.*zeta3*log2));
}
dcomplex MellinFunc::S12mzplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return (pow(-1.,n)*(H.HS(2,1,-1,N-1.)+(H.HS(2,1,N-1.)-H.HS(2,-1,N-1.))*log2-0.5*(H.HS(2,N-1.)-H.HS(-2,N-1.))*log2q-1./8.*zeta3*H.HS(-1,N-1.)-3.*Li4+6./5.*zeta2q
  -11./4.*zeta3*log2+3./4.*zeta2*log2q-1./8.*log2q*log2q));
}

dcomplex MellinFunc::S12zplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(H.HS(-2,1,1,N-1.)-zeta3*H.HS(-1,N-1.)+Li4-1./8.*zeta2q-1./8.*zeta3*log2-1./4.*zeta2*log2q+1./24.*log2q*log2q));
}
dcomplex MellinFunc::Li3mzplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-pow(-1.,n)*(H.HS(3,-1,N-1.)+(H.HS(3,N-1.)-H.HS(-3,N-1.))*log2+0.5*zeta2*H.HS(-2,N-1.)-3./4.*zeta3*H.HS(-1,N-1.)+1./8.*zeta2q-3./4.*zeta3*log2));
}
dcomplex MellinFunc::Li2zLogminus(dcomplex N){
  return(1./N*(-2.*zeta3-zeta2*H.HS(1,N)+1./N*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N))+H.HS(2,1,N)));
}

dcomplex MellinFunc::Li2zLogplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./(2.*N*N)*(-zeta2+H.HS(-1,N)*H.HS(-1,N)+H.HS(2,N)+2.*H.HS(-1,N)*log2-2.*H.HS(1,N)*log2+log2q+2.*N*zeta2*log2
  +pow(-1.,n)*(2.*N*H.HS(-2,1,N)-2.*N*zeta2*(H.HS(-1,N)+log2)+(zeta2+2.*H.HS(-1,1,N)-log2q)+5./4.*N*zeta3)));
}
dcomplex MellinFunc::Li2mzLogplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./(N*N)*(log2*(-0.5*zeta2*N+log2)+pow(-1.,n)*(2.*H.HS(1,-1,N)+N*H.HS(2,-1,N)+H.HS(-1,N)*(0.5*N*zeta2-2.*log2)
  +0.5*N*zeta2*log2-log2*(N*H.HS(-2,N)-2.*H.HS(1,N)-N*H.HS(2,N)+log2)-1./4.*N*zeta3)));
}
dcomplex MellinFunc::Li2mzLogz(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./(N*N)*0.5*zeta2-2./(N*N*N)*log2+pow(-1.,n)*(1./(N*N)*(0.5*zeta2+H.HS(-2,N))+1./(N*N*N)*(2.*H.HS(-1,N)+2.*log2)));
}
dcomplex MellinFunc::Li2zLogzplusminus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./960.*(160.*6.*zeta2*(pow(-1.,n)*H.HS(-2,N-1.)+H.HS(2,N-1))+pow(-1.,n)*(36.*zeta2q-960.*H.HS(-3,1,N-1.)-480.*H.HS(-2,2,N-1.))
  -4.*(36.*zeta2q+60.*(H.HS(2,N-1.)*H.HS(2,N-1.)+H.HS(4,N-1.))+240.*H.HS(3,1,N-1.))));
}
dcomplex MellinFunc::Li2mzLogzplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-pow(-1.,n)*(2.*H.HS(3,-1,N-1.)+H.HS(2,-2,N-1.)-2.*H.HS(-3,N-1.)*log2+2.*H.HS(3,N-1.)*log2+0.5*zeta2*H.HS(2,N-1.)
  +0.5*zeta2*H.HS(-2,N-1.)-4.*Li4+13./8.*zeta2q-7./2.*zeta3*log2+zeta2*log2q-1./6.*log2q*log2q)
  );
}
dcomplex MellinFunc::LogzLogminus2minus(dcomplex N){
  return(H.HS(2,N-1.)*H.HS(2,N-1.)-zeta2*H.HS(2,N-1.)+2.*H.HS(4,N-1.)+H.HS(1,N-1.)*H.HS(1,N-1.)*(H.HS(2,N-1.)-zeta2)+2.*H.HS(1,N-1.)*(H.HS(3,N-1.)-zeta3)-4./5.*zeta2q);
}
dcomplex MellinFunc::Li2mzLogplusplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-pow(-1.,n)*(H.HS(1,2,-1,N-1.)+2.*H.HS(2,1,-1,N-1.)+(H.HS(2,1,N-1.)-H.HS(1,-2,N-1.)-2.*H.HS(2,-1,N-1.)+H.HS(1,N-1.)*H.HS(2,N-1.)+H.HS(3,N-1.)-0.5*zeta2*H.HS(-1,N-1.))*log2
    +0.5*zeta2*H.HS(1,-1,N-1.)-(H.HS(2,N-1.)-H.HS(-2,N-1.))*log2q-(1./4.*zeta3-0.5*zeta2*log2)*H.HS(1,N-1.)-3.*Li4+6./5.*zeta2q-21./8.*zeta3*log2+0.5*zeta2*log2q-1./8.*log2q*log2q)
  );
}
dcomplex MellinFunc::Logz2Logminusminus(dcomplex N){
  return(-2.*zeta3*H.HS(1,N-1.)+2.*H.HS(1,N-1.)*H.HS(3,N-1.)-2.*zeta2*H.HS(2,N-1.)+H.HS(2,N-1.)*H.HS(2,N-1.)+3.*H.HS(4,N-1.)-zeta2q/5.);
}
dcomplex MellinFunc::Liz2Logminus(dcomplex N){
  return(1./N*(-2.*zeta3-zeta2*H.HS(1,N)+1./N*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N))+H.HS(2,1,N)));
}
dcomplex MellinFunc::Logz3minusplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(21./20.*zeta2q+3.*H.HS(-4,N-1.))-6./5.*zeta2q+3.*H.HS(4,N-1.));
}
dcomplex MellinFunc::Logplusplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(H.HS(1,-1,N-1.)-0.5*log2q+(-H.HS(-1,N-1.)+H.HS(1,N-1.))*log2));
}
dcomplex MellinFunc::Li2zLogplusplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(3./40.*zeta2q-H.HS(-2,-1,-1,N-1.)-H.HS(1,-2,1,N-1.)-H.HS(2,-1,1,N-1.)+0.5*H.HS(-2,N-1.)*(zeta2-log2q)
    -0.5*H.HS(2,N-1.)*(zeta2-log2q)-(H.HS(-2,-1,N-1.)-H.HS(-2,1,N-1.))*log2-1./4.*zeta2*(-4.*H.HS(1,-1,N-1.)+2.*log2q
    +4.*(H.HS(-1,N-1.)-H.HS(1,N-1.))*log2)-5./8.*H.HS(1,N-1.)*zeta3)
  );
}
dcomplex MellinFunc::Logz2Logplusplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(2.*H.HS(3,-1,N-1.)+2.*H.HS(1,-3,N-1.)+2.*H.HS(2,-2,N-1.)+zeta2*H.HS(2,N-1.)+3./2.*zeta3*H.HS(1,N-1.)+2.*H.HS(3,N-1)*log2
	-2.*H.HS(-3,N-1.)*log2-4.*Li4+3./2.*zeta2q-7./2.*zeta3*log2+zeta2*log2q-1./6.*log2q*log2q));
}
dcomplex MellinFunc::LogzLogplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./12./(N*N)*(-12.*log2+pow(-1.,n)*(6.*N*zeta2+12.*N*H.HS(-2,N)+12.*H.HS(-1,N)+12.*log2)));
}
dcomplex MellinFunc::Logz2Logplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./(N*N*N)*2.*log2-pow(-1.,n)*(1./(N*N)*(zeta2+2.*H.HS(-2,N))+1./(N*N*N)*(2.*H.HS(-1,N)+2.*log2)+1./N*(2.*H.HS(-3,N)+3./2.*zeta3)));
}
dcomplex MellinFunc::LogzLogplus2(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-1./(12.*N*N)*(12.*log2q+pow(-1.,n)*(2.*H.HS(1,N)*(6.*N*zeta2+12.*log2)+3.*(8.*N*H.HS(1,-2,N)+8.*H.HS(1,-1,N)+8.*N*H.HS(2,-1,N)
  -4.*log2*(2.*N*H.HS(-2,N)+2.*H.HS(-1,N)-2.*N*H.HS(2,N)+log2)-N*zeta3))));
}
dcomplex MellinFunc::LogzLogplus2plus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(2.*pow(-1.,n)*(H.HS(1,1,-2,N-1.)+H.HS(1,2,-1,N-1.)+H.HS(2,1,-1,N-1.)-H.HS(1,-2,N-1.)*log2
    -H.HS(2,-1,N-1)*log2+H.HS(1,N-1.)*H.HS(2,N-1.)*log2+0.25*zeta2*H.HS(1,N-1.)*H.HS(1,N-1.)+H.HS(3,N-1.)*log2
    +0.25*zeta2*H.HS(2,N-1.)-0.5*H.HS(2,N-1.)*log2q+0.5*H.HS(-2,N-1.)*log2q-1./8.*zeta3*H.HS(1,N-1.)-Li4+
    2./5.*zeta2q-7./8.*zeta3*log2+0.25*zeta2*log2q-1./24.*log2q*log2q)
  );
}
dcomplex MellinFunc::Li2z(dcomplex N){
  return(1./N*zeta2-1./(N*N)*H.HS(1,N));
}
dcomplex MellinFunc::Li2mz(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-0.5*zeta2/N-1./(N*N)*(pow(-1.,n)*(H.HS(-1,N)+log2)-log2));
}
dcomplex MellinFunc::S12z(dcomplex N){
  return(-1./(2.*N*N)*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N)-2.*N*zeta3));
}
dcomplex MellinFunc::S12mz(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./(8.*N*N)*(-8.*pow(-1.,n)*H.HS(1,-1,N)+8.*pow(-1.,n)*H.HS(-1,N)*log2
  -8.*pow(-1.,n)*H.HS(1,N)*log2-4.*log2q+4.*pow(-1.,n)*log2q+N*zeta3));
}
dcomplex MellinFunc::S12z2(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-1./(6.*N*N)*((-1.+pow(-1.,n))*6.*zeta2+12.*pow(-1.,n)*H.HS(-2,N)+12.*H.HS(2,N)+6.*(H.HS(-1,N)*H.HS(-1,N)+H.HS(1,N)*H.HS(1,N)+2.*log2q
  +2.*pow(-1.,n)*(H.HS(1,N)-log2)*(H.HS(-1,N)+log2)+(H.HS(-1,N)-H.HS(1,N))*2.*log2)-6.*N*zeta3));
}
dcomplex MellinFunc::Li3z(dcomplex N){
  return(1./(N*N*N)*H.HS(1,N)-1./(N*N)*zeta2+1./N*zeta3);
}
dcomplex MellinFunc::Li3mz(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*1./(N*N*N)*(log2+H.HS(-1,N))-log2/(N*N*N)+0.5/(N*N)*zeta2-3./4./N*zeta3);
}
dcomplex MellinFunc::Li2zLogz(dcomplex N){
  return(2./(N*N*N)*H.HS(1,N)+1./(N*N)*H.HS(2,N)-2./(N*N)*zeta2);
}
dcomplex MellinFunc::Logminus(dcomplex N){
  return(-H.HS(1,N)/N);
}
dcomplex MellinFunc::Logplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(1./N*(log2-pow(-1.,n)*(H.HS(-1,N)+log2)));
}
dcomplex MellinFunc::Logz(dcomplex N){
  return(-1./(N*N));
}
dcomplex MellinFunc::Logz2(dcomplex N){
  return(2./(N*N*N));
}
dcomplex MellinFunc::Logz3(dcomplex N){
  return(-6./(N*N*N*N));
}
dcomplex MellinFunc::LogzLogminus(dcomplex N){
  return(-1./N*zeta2+1./N*H.HS(2,N)+1./(N*N)*H.HS(1,N));
}
dcomplex MellinFunc::LogzLogminus2(dcomplex N){
  return(1./(N*N)*(-H.HS(1,N)*H.HS(1,N)+2.*N*H.HS(1,N)*(zeta2-H.HS(2,N))-H.HS(2,N)-2.*N*H.HS(3,N)+2.*N*zeta3));
}
dcomplex MellinFunc::Logz2Logminus(dcomplex N){
  return(-2./(N*N*N)*H.HS(1,N)+1./(N*N)*(2.*zeta2-2.*H.HS(2,N))-2./N*(H.HS(3,N)-zeta3));
}
dcomplex MellinFunc::Logz2minus(dcomplex N){
  return(2.*(-H.HS(3,N-1.)+zeta3));
}
dcomplex MellinFunc::Logminus2(dcomplex N){
  return(1/N*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N)));
}
dcomplex MellinFunc::Logminus3(dcomplex N){
  return(-1./N*(H.HS(1,N)*H.HS(1,N)*H.HS(1,N)+3.*H.HS(1,N)*H.HS(2,N)+2.*H.HS(3,N)));
}
dcomplex MellinFunc::LogzLogminusminus(dcomplex N){
  return(zeta3+zeta2*H.HS(1,N-1.)-H.HS(1,N-1.)*H.HS(2,N-1.)-H.HS(3,N-1.));
}
dcomplex MellinFunc::Logzminus(dcomplex N){
  return(-zeta2+H.HS(2,N-1.));
}
dcomplex MellinFunc::Logzminusplus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(pow(-1.,n)*(0.25*zeta2+0.5*H.HS(-2,N-1.))+0.5*H.HS(2,N-1.)-0.5*zeta2);
}
dcomplex MellinFunc::plus(dcomplex N){
  long double frac,n=0.;
  frac=modf(real(N),&n);
  return(-pow(-1.,n)*(H.HS(-1,N-1.)+log2));
}
dcomplex MellinFunc::Li2minusminus(dcomplex N){
  return(-zeta2*H.HS(1,N-1.)+H.HS(1,2,N-1.)+zeta3);
}
dcomplex MellinFunc::S12zregminus(dcomplex N){
  return(-6./5.*zeta2q+H.HS(2,1,1,N-1.));
}
dcomplex MellinFunc::Li2zregminus(dcomplex N){
  return(H.HS(2,1,N-1.)-2.*zeta3);
}
dcomplex MellinFunc::Li3zregminus(dcomplex N){
  return(-0.5*zeta2q+zeta2*H.HS(2,N-1.)-H.HS(3,1,N-1.));
}
dcomplex MellinFunc::Li2zregLogminusminus(dcomplex N){
  return(6./5.*zeta2q-H.HS(1,2,1,N-1.)-2.*H.HS(2,1,1,N-1.)+2.*H.HS(1,N-1.)*zeta3);
}
dcomplex MellinFunc::D0(dcomplex N){
  return(-H.HS(1,N-1.));
}
dcomplex MellinFunc::D1(dcomplex N){
  return(H.HS(1,1,N-1.));
}
dcomplex MellinFunc::D2(dcomplex N){
  return(-2.*H.HS(1,1,1,N-1.));
}
dcomplex MellinFunc::D3(dcomplex N){
  return(6.*H.HS(1,1,1,1,N-1.));
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
