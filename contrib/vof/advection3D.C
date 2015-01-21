

// ======================================
// ======================================
//  Advection
// ======================================
// ======================================
// Split advection along the x-direction:
void  MGSolCC::lagrangeX(const unsigned int Level,const double dt){

  const unsigned int nxl=_nx[Level]+1;  
  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1;  int fd=(nxl-2)/nxc;
  double st=dt/hx;
  // First line
  for (unsigned int ix=0; ix< nxl*9; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }
  /*c1=0.*/
  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for (int k=0; k<nz ; k++) {
    for (int j=0; j<ny ; j++) {
      for (int i=0; i<nxl-2 ; i++) {


        int ixyz[3];double utemp[6];
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;
        int indl=i+((k%3)*3+(j%3))*nxl;
        double vof1 = 0.;double  vof2 = 0.;double vof3 = 0.;
        const double  c_0=_tempc1[indl+1];
        // c=1.
        if ( c_0 > 1.-MIN_VAL) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[0]*st;double a2=utemp[1]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          vof1 =c_0* MAX(-a1,0.);
          vof2 =c_0* ( 1.- MAX(a1,0.)-MAX(-a2,0.));
          vof3 =c_0* MAX(a2,0.);
        }
        // 0.< c <1.
        else if (c_0 > MIN_VAL ) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[0]*st;double a2=utemp[1]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          // unit normal
          double mxt=_tempmx1[indl+1];double myt=_tempmy1[indl+1];double mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          mxt /= mm; myt /= mm; mzt /= mm;
          // alpha
          double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          double mmx=mxt;double mmy=myt;double mmz=mzt;double tnot=alphat;
          // splitting correction factor
          mmx=mmx/(1.-a1+a2);tnot=tnot+mmx*a1;
          if (a1<0.) vof1=get_vol3D(mmx,mmy,mmz,tnot,a1,-a1);
          if (a2>0.) vof3=get_vol3D(mmx,mmy,mmz,tnot,1.,a2);
          vof2=get_vol3D(mmx,mmy,mmz,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
        }
        // temporary storage
        _tempc2[indl+1]+=vof2; _tempc2[indl]+=vof1; _tempc2[indl+2]+=vof3;
      }
#ifdef PERIODIC
      // periodic boundary along x-axis
      _tempc2[((k%3)*3+(j%3))*nxl+1] +=_tempc2[((k%3)*3+(j%3))*nxl+nxl-1];
      // _tempc2[((k%3)*3+(j%3))*nxl+nxl-2] +=_tempc2[((k%3)*3+(j%3))*nxl];
#endif
      // Contracting _tempc2 -> c1
      CtrRow(_tempc2,&c1[Level],Level,j+k*ny,((j%3)+(k%3)*3)*nxl);
      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for (int k2=-1;k2<2;k2++) {
        int indj=j+2+(k+k2)*ny; int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
      }

    }
    // new block for new k line
    // storage fo the last line in j
    // CtrRow(_tempc2,c1[Level],Level,ny-1+(k)*ny,((ny-1)%3+((k)%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl *3*3; ix++) {
      _tempc1[ix]=0.;_tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=k; k0<k+3; k0++) for (unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }
  return;
}

// ---------------------------------------------------------------
// Split advection along the y-direction
void  MGSolCC::lagrangeY(const unsigned int Level,const double dt) {

  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1; int fd=(nxl-2)/nxc;
  double st=dt/(hy);
  // First line
  for (unsigned int ix=0; ix< nxl *3*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }

  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for (int k=0 ; k<nz ; k++) {
    for (int j=0 ; j<ny ; j++) {
      for (int i=0 ; i<nxl-2 ; i++) {
        int ixyz[3];double utemp[6];
        //int indv = i+(j+k*(ny+1))*(nxl-1);
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;
        int indl=i+( (k%3)*3+(j%3) )*nxl;

        const double  c_0=_tempc1[indl+1];
        double vof1 = 0.;double  vof2 = 0.;double vof3 = 0.;
        // c=1.
        if ( c_0 > 1.-MIN_VAL) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[2]*st; double a2=utemp[3]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          vof1 =c_0* MAX(-a1,0.);
          vof2 =c_0* ( 1. - MAX(a1,0.) - MAX(-a2,0.));
          vof3 =c_0* MAX(a2,0.);
        }
        // 0.< c <1.
        else if (c_0 > MIN_VAL) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[2]*st;double a2=utemp[3]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}

          double mxt=_tempmx1[indl+1];double myt=_tempmy1[indl+1];double mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          mxt /= mm; myt /= mm; mzt /= mm;
          double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          double mmx=mxt;double mmy=myt; double mmz=mzt; double tnot=alphat;

          mmy=mmy/(1.-a1+a2);tnot=tnot+mmy*a1;
          if (a1<0.) vof1=get_vol3D(mmy,mmx,mmz,tnot,a1,-a1);
          if (a2>0.) vof3=get_vol3D(mmy,mmx,mmz,tnot,1.,a2);
          vof2=get_vol3D(mmy,mmx,mmz,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
        }
        // new values of c1
        _tempc2[indl+1]+=vof2;
        _tempc2[i+((j-1+3)%3+(k%3)*3)*nxl+1]+=vof1;
        _tempc2[i+((j+1)%3+(k%3)*3)*nxl+1] +=vof3;
      }

      // contracting  tempc2 -> c2
      if (j>0) CtrRow(_tempc2,&c1[Level],Level,j-1+k*ny,((j-1+3)%3+(k%3)*3)*nxl);

      //Expanding matrix c1_old ->line  _tempc1
      for (int k2=-1;k2<2;k2++) {
        int indjt=j+2+(k+k2)*ny;   int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indjt,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indjt,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indjt,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indjt,indjl);
      }
    }
    // new block for new k line
    // storage fo the last line in j
    CtrRow(_tempc2,&c1[Level],Level,ny-1+k*ny,( (ny-1)%3+(k%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.; _tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=k; k0<k+3; k0++) for (unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }
  return;
}
// ---------------------------------------------------------------
// Split advection along the z-direction
// ---------------------------------------------------------------
void  MGSolCC::lagrangeZ(const unsigned int Level,
                         const double dt) {
  // Set up
  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1;  int fd=(nxl-2)/nxc;
  double vof1,vof2,vof3,mmx,mmy,mmz,tnot;
  double mxt,myt,mzt,alphat;
  double a1,a2;int ixyz[3];double utemp[6];
  double st=dt/(hz);
  // First line
  for (unsigned int ix=0; ix< nxl*3*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int ind=k0*ny+j0; int indl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,ind,indl);
      ExpRow(&_mx1[Level],_tempmx1,Level,ind,indl);
      ExpRow(&_my1[Level],_tempmy1,Level,ind,indl);
      ExpRow(&_mz1[Level],_tempmz1,Level,ind,indl);
    }

  // Domain cartesian loop
  for (int j=0 ; j<ny ; j++) {
    for (int k=0 ; k<nz ; k++) {
      for (int i=0 ; i<nxl-2 ; i++) {
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;

        // c=0. and fluxes=0.
        int indl=i+((k%3)*3+(j%3))*nxl;
        vof1 = 0.; vof2 = 0.; vof3 = 0.;
        const double  c_0=_tempc1[indl+1];
        // c=1.
        if ( c_0 >= 1.) {
          get_vel(ixyz,fd,utemp);
          a1=utemp[4]*st;a2=utemp[5]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          vof1 =c_0*MAX(-a1,0.);
          vof2 =c_0*(1.-MAX(a1,0.)-MAX(-a2,0.));
          vof3 =c_0*MAX(a2,0.);
        }
        //
        else if (c_0 > 0.) {
          get_vel(ixyz,fd,utemp);
          a1=utemp[4]*st;a2=utemp[5]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          // unit normal
          mxt=_tempmx1[indl+1];myt=_tempmy1[indl+1];mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          mxt /= mm; myt /= mm; mzt /= mm;

          alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          mmx=mxt;mmy=myt;;mmz=mzt;tnot=alphat;

          mmz=mmz/(1.-a1+a2); tnot=tnot+mmz*a1;
          if (a1<0.)   vof1=get_vol3D(mmz,mmx,mmy,tnot,a1,-a1);
          if (a2>0.)   vof3=get_vol3D(mmz,mmx,mmy,tnot,1.,a2);
          vof2=get_vol3D(mmz,mmx,mmy,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
        }
        // new values of c1 (in temp2)
        // ---------------------------------------
        _tempc2[indl+1]+=vof2;
        _tempc2[i+((j%3)+((k-1+3)%3)*3)*nxl+1]+=vof1; _tempc2[i+(j%3+((k+1)%3)*3)*nxl+1]+=vof3;
      }

      // contraction temp2 -> c1
      if (k>0) CtrRow(_tempc2,&c1[Level],Level,j+(k-1)*ny,((j%3)+((k-1)%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for (int j2=-1;j2<2;j2++) {
        int indk=j+j2+(k+2)*ny; int indkl=((j+j2+3)%3+((k+2)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indkl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indkl);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indkl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indk,indkl);
      }
    }
    // new block for new j line
    // storage fo the last line in k
    CtrRow(_tempc2,&c1[Level],Level,j+(nz-1)*ny,(j%3+((nz-1)%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.;_tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=j; j0<j+3; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0%3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }

  return;
}

// ---------------------------------------------------------------
// ---------------------------------------------------------------
// Split advection along the x-direction:
void  MGSolCC::eulerianoX(const unsigned int Level,const double dt) {

  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1;  int fd=(nxl-2)/nxc;
  double st=dt/hx;
  // First line
  for (unsigned int ix=0; ix< nxl*9; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }
  /*c1=0.*/
  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for (int k=1; k<nz ; k++) {
    for (int j=1; j<ny ; j++) {

      for (int i=0; i<nxl-1 ; i++) {

        int ixyz[3]; ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;

        // c=0. and fluxes=0.
        int indl=i+((k%3)*3+(j%3))*nxl;
        const double  c_0=_tempc1[indl+1];

        if (c_0+_tempc1[indl]+_tempc1[indl+2] > MIN_VAL){
          double vof1 = 0.; double vof2 = 0.; double  vof3 = 0.;
          double utemp[6];
          //double implicito;
          get_vel(ixyz,fd,utemp);
          double a1=utemp[0]*st; double a2=utemp[1]*st;
          // compression factor
          double implicito=(a2-a1); if (fabs(implicito) < MIN_VAL) implicito=0.;
          _tempc2[i+(((k+1)%3)*3+((j+1)%3))*nxl+1] =implicito;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          // c=1.
          if ( c_0 > 1.-MIN_VAL) {vof1 = MAX(-a1,0.); vof3 =MAX(a2,0.);}
          // 0.< c <1.
          else if (c_0 > MIN_VAL) {
            // unit normal
            double mxt=_tempmx1[indl+1]; double myt=_tempmy1[indl+1]; double mzt=_tempmz1[indl+1];
            double mm=fabs(mxt)+fabs(myt)+fabs(mzt)+MIN_VAL;
            mxt /= mm; myt /= mm; mzt /= mm;
            // alpha
            double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
            alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
            double mmx=mxt; double mmy=myt; double mmz=mzt; double tnot=alphat;
            // splitting correction factor
            //  mmx=mmx/(1.-a1+a2);tnot=tnot+mmx*a1;
            if (a1<0.) vof1=get_vol3D(mmx,mmy,mmz,tnot,0.,-a1);
            //    vof1=get_vol3D(mmx,mmy,mmz,tnot,a1,-a1);
            if (a2>0.) vof3=get_vol3D(mmx,mmy,mmz,tnot,1.-a2,a2);
            //  vof3=get_vol3D(mmx,mmy,mmz,tnot,1.,a2);
            // vof2=get_vol3D(mmx,mmy,mmz,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));
          }
          _tempc2[indl+1] +=(c_0-vof1-vof3 );
          _tempc2[indl]+=vof1;
          _tempc2[indl+2]+=vof3;
        }

      }
      // compression factor      i+((k%3)*3+(j%3))*nxl;
      for (int i=0; i<nxl-1 ; i++)  if (_tempc2[i+((k%3)*3+(j%3))*nxl+1] >MIN_VAL){
          _tempc2[i+((k%3)*3+(j%3))*nxl+1] /=(1.-_tempc2[i+(((k+1)%3)*3+((j+1)%3))*nxl+1]);
        }
      // Contracting _tempc2 -> c1
      CtrRow(_tempc2,&c1[Level],Level,j+k*ny,((j%3)+(k%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for (int k2=-1;k2<2;k2++) {
        int indj=j+2+(k+k2)*ny; int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
      }
    }
    // new block for new k line
    // storage fo the last line in j
    // CtrRow(_tempc2,c1[Level],Level,ny-1+(k)*ny,((ny-1)%3+((k)%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl *3*3; ix++) {
      _tempc1[ix]=0.;_tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=k; k0<k+3; k0++) for (unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }
  return;
}


// ---------------------------------------------------------------
// Split advection along the y-direction
void  MGSolCC::eulerianoY(const unsigned int Level,const double dt) {



  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1; int fd=(nxl-2)/nxc;
  double st=dt/(hy);
  // First line
  for (unsigned int ix=0; ix< nxl *3*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }
  int ixyz[3];double utemp[6];
  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for (int k=0 ; k<nz ; k++) {
    for (int j=0 ; j<ny ; j++) {
      for (int i=0 ; i<nxl-2 ; i++) {

        //int indv = i+(j+k*(ny+1))*(nxl-1);
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;
        int indl=i+( (k%3)*3+(j%3) )*nxl;

        const double  c_0=_tempc1[indl+1];
        double vof1 = 0.;double  vof2 = 0.;double vof3 = 0.;
        // c=1.
        if ( c_0 >= 1.) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[2]*st; double a2=utemp[3]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          vof1 =c_0* MAX(-a1,0.);
          vof2 =c_0* ( 1. - MAX(a1,0.) - MAX(-a2,0.));
          vof3 =c_0* MAX(a2,0.);
        }
        // 0.< c <1.
        else if (c_0 > 0.) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[2]*st;double a2=utemp[3]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          double mxt=_tempmx1[indl+1];double myt=_tempmy1[indl+1];double mzt=_tempmz1[indl+1];
          // mxt=0.;myt=0.;mzt=1.;
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          if (mm >0.) { mxt /= mm; myt /= mm; mzt /= mm;  }
          double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          double mmx=mxt;double mmy=myt; double mmz=mzt; double tnot=alphat;

          // mmy=mmy/(1.-a1+a2);tnot=tnot+mmy*a1;
          if (a1<0.)// vof1=get_vol3D(mmy,mmz,mmx,tnot,a1,-a1);
            vof1=get_vol3D(mmy,mmz,mmx,tnot,0.,-a1);
          if (a2>0.) // vof3=get_vol3D(mmy,mmz,mmx,tnot,1.,a2);
            vof3=get_vol3D(mmy,mmz,mmx,tnot,1.-a2,a2);
          //  vof2=get_vol3D(mmy,mmz,mmx,tnot,MAX(0.,a1),1-MAX(0.,-a2)-MAX(0.,a1));


        }
        // new values of c1
        _tempc2[indl+1]+=c_0-vof1-vof3;//
//_tempc2[indl+1]+=vof2;
        _tempc2[i+((j-1+3)%3+(k%3)*3)*nxl+1]+=vof1;
        _tempc2[i+((j+1)%3+(k%3)*3)*nxl+1] +=vof3;
      }
      //eulerian correction
      for (unsigned int ix=0; ix< nxl; ix++) if (_tempc2[((j-1+3)%3+(k%3)*3)*nxl+ix+1] >  MIN_VAL){
          ixyz[0]=ix; ixyz[1]=j-1;ixyz[2]=k;
          get_vel(ixyz,fd,utemp);
          double implicito=(utemp[3]-utemp[2])*st;
          if (fabs(implicito) > MIN_VAL)
            _tempc2[((j-1+3)%3+(k%3)*3)*nxl+ix+1] /=(1.-implicito);
        }
      // contracting  tempc2 -> c2
      if (j>0) CtrRow(_tempc2,&c1[Level],Level,j-1+k*ny,((j-1)%3+(k%3)*3)*nxl);

      //Expanding matrix c1_old ->line  _tempc1
      for (int k2=-1;k2<2;k2++) {
        int indjt=j+2+(k+k2)*ny;   int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indjt,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indjt,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indjt,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indjt,indjl);
      }
    }
    // new block for new k line
    // storage fo the last line in j
    CtrRow(_tempc2,&c1[Level],Level,ny-1+k*ny,( (ny-1)%3+(k%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.; _tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=k; k0<k+3; k0++) for (unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }
  return;
}

// ---------------------------------------------------------------
// Split advection along the y-direction
void  MGSolCC::eulerianoZ(const unsigned int Level,
                          const double dt) {
  // Set up
  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1;  int fd=(nxl-2)/nxc;
  double st=dt/(hz);
  // First line
  for (unsigned int ix=0; ix< nxl*3*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int ind=k0*ny+j0; int indl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,ind,indl);
      ExpRow(&_mx1[Level],_tempmx1,Level,ind,indl);
      ExpRow(&_my1[Level],_tempmy1,Level,ind,indl);
      ExpRow(&_mz1[Level],_tempmz1,Level,ind,indl);
    }

  int ixyz[3]; double utemp[6];
  // Domain cartesian loop
  for (int j=0 ; j<ny ; j++) {
    for (int k=0 ; k<nz ; k++) {
      for (int i=0 ; i<nxl-2 ; i++) {
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;

        int indl=i+((k%3)*3+(j%3))*nxl;
        double  vof1 = 0.; double  vof3 = 0.;
        const double  c_0=_tempc1[indl+1];
        // c=1.
        if ( c_0 >= 1.) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[4]*st;double a2=utemp[5]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          vof1 =c_0*MAX(-a1,0.);   vof3 =c_0*MAX(a2,0.);
        }
        //
        else if (c_0 > 0.) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[4]*st;double a2=utemp[5]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          // unit normal
          double mxt=_tempmx1[indl+1];double myt=_tempmy1[indl+1];double mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          mxt /= mm; myt /= mm; mzt /= mm;

          double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          double mmx=mxt;double mmy=myt; double mmz=mzt;double tnot=alphat;
          if (a1<0.)  vof1=get_vol3D(mmz,mmx,mmy,tnot,0,-a1);
          if (a2>0.)  vof3=get_vol3D(mmz,mmx,mmy,tnot,1.-a2,a2);
        }
        // new values of c1 (in temp2)
        // ---------------------------------------
        _tempc2[indl+1]+=c_0-vof1-vof3;//vof2;
        _tempc2[i+((j%3)+((k-1+3)%3)*3)*nxl+1]+=vof1; _tempc2[i+(j%3+((k+1)%3)*3)*nxl+1]+=vof3;
      }
      // eulerian correction
      for (unsigned int ix=0; ix< nxl; ix++) if (_tempc2[((j%3)+((k-1+3)%3)*3)*nxl+ix+1]>MIN_VAL){
          ixyz[0]=ix; ixyz[1]=j;ixyz[2]=k-1;
          get_vel(ixyz,fd,utemp);
          double implicito=(utemp[5]-utemp[4])*st;
          if (fabs(implicito) > MIN_VAL)  _tempc2[((j%3)+((k-1+3)%3)*3)*nxl+ix+1] /=(1.-implicito);
        }

      // contraction temp2 -> c1
      if (k>0) CtrRow(_tempc2,&c1[Level],Level,j+(k-1)*ny,((j%3)+((k-1)%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for (int j2=-1;j2<2;j2++) {
        int indk=j+j2+(k+2)*ny; int indkl=((j+j2+3)%3+((k+2)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indkl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indkl);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indkl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indk,indkl);
      }
    }
    // new block for new j line
    // storage fo the last line in k
    CtrRow(_tempc2,&c1[Level],Level,j+(nz-1)*ny,(j%3+((nz-1)%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.;_tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=j; j0<j+3; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0%3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }

  return;
}

// ---------------------------------------------------------------
// ---------------------------------------------------------------
// Split advection along the x-direction:
void  MGSolCC::eulerbalotX(const unsigned int Level,const double dt,const int idir2) {

  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1;  int fd=(nxl-2)/nxc;
  double st=dt/hx;
  // First line
  for (unsigned int ix=0; ix< nxl*9; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }
  /*c1=0.*/
  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for (int k=1; k<nz ; k++) {
    for (int j=1; j<ny ; j++) {

      for (int i=0; i<nxl-1 ; i++) {

        int ixyz[3]; ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;

        // c=0. and fluxes=0.
        int indl=i+((k%3)*3+(j%3))*nxl;
        const double  c_0=_tempc1[indl+1];

        if (c_0+_tempc1[indl]+_tempc1[indl+2] > MIN_VAL){
          double vof1 = 0.; double vof2 = 0.; double  vof3 = 0.;
          double utemp[6];
          //double implicito;
          get_vel(ixyz,fd,utemp);
          double a1=utemp[0]*st; double a2=utemp[1]*st;
          double b1=utemp[idir2*2]*st;double b2=utemp[idir2*2+1]*st;
          // compression factor
          double implicito=(a2-a1)+(b2-b1); if (fabs(implicito) < MIN_VAL) implicito=0.;
          _tempc2[i+(((k+1)%3)*3+((j+1)%3))*nxl+1] =implicito;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          // c=1.
          if ( c_0 > 1.-MIN_VAL) {vof1 = MAX(-a1,0.); vof3 =MAX(a2,0.);}
          // 0.< c <1.
          else if (c_0 > MIN_VAL) {
            // unit normal
            double mxt=_tempmx1[indl+1]; double myt=_tempmy1[indl+1]; double mzt=_tempmz1[indl+1];
            double mm=fabs(mxt)+fabs(myt)+fabs(mzt)+MIN_VAL;
            mxt /= mm; myt /= mm; mzt /= mm;
            // alpha
            double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
            alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
            double mmx=mxt; double mmy=myt; double mmz=mzt; double tnot=alphat;
            // splitting correction factor
            if (a1<0.) vof1=get_vol3D(mmx,mmy,mmz,tnot,0.,-a1);
            if (a2>0.) vof3=get_vol3D(mmx,mmy,mmz,tnot,1.-a2,a2);
          }
          _tempc2[indl+1] +=(c_0*(1-(b2-b1))-vof1-vof3 );
          _tempc2[indl]+=vof1;
          _tempc2[indl+2]+=vof3;
        }

      }
      // compression factor      i+((k%3)*3+(j%3))*nxl;
      for (int i=0; i<nxl-1 ; i++)  if (_tempc2[i+((k%3)*3+(j%3))*nxl+1] >MIN_VAL){
          _tempc2[i+((k%3)*3+(j%3))*nxl+1] /=(1.-_tempc2[i+(((k+1)%3)*3+((j+1)%3))*nxl+1]);
        }
      // Contracting _tempc2 -> c1
      CtrRow(_tempc2,&c1[Level],Level,j+k*ny,((j%3)+(k%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for (int k2=-1;k2<2;k2++) {
        int indj=j+2+(k+k2)*ny; int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
      }
    }
    // new block for new k line
    // storage fo the last line in j
    // CtrRow(_tempc2,c1[Level],Level,ny-1+(k)*ny,((ny-1)%3+((k)%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl *3*3; ix++) {
      _tempc1[ix]=0.;_tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=k; k0<k+3; k0++) for (unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }
  return;
}


// ---------------------------------------------------------------
// Split advection along the y-direction
void  MGSolCC::eulerbalotY(const unsigned int Level,const double dt,const int idir2) {



  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1; int fd=(nxl-2)/nxc;
  double st=dt/(hy);
  // First line
  for (unsigned int ix=0; ix< nxl *3*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int indj=k0*ny+j0; int indjl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,indj,indjl);
      ExpRow(&_mx1[Level],_tempmx1,Level,indj,indjl);
      ExpRow(&_my1[Level],_tempmy1,Level,indj,indjl);
      ExpRow(&_mz1[Level],_tempmz1,Level,indj,indjl);
    }
  int ixyz[3];double utemp[6];
  //for (i=0 ; i<n_nodes ; i++)  c1[i]=0.;
  for (int k=0 ; k<nz ; k++) {
    for (int j=0 ; j<ny ; j++) {
      for (int i=0 ; i<nxl-2 ; i++) {

        //int indv = i+(j+k*(ny+1))*(nxl-1);
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;
        int indl=i+( (k%3)*3+(j%3) )*nxl;

        const double  c_0=_tempc1[indl+1];
        double vof1 = 0.;double vof3 = 0.;
        double b1,b2;
        // c=1.
        if ( c_0 > 1.-MIN_VAL) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[2]*st; double a2=utemp[3]*st;
          b1=utemp[2*idir2]*st;b2=utemp[2*idir2+1]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          vof1 =c_0* MAX(-a1,0.);
          vof3 =c_0* MAX(a2,0.);
        }
        // 0.< c <1.
        else if (c_0 > MIN_VAL) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[2]*st;double a2=utemp[3]*st;
          b1=utemp[2*idir2]*st;b2=utemp[2*idir2+1]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 \n ");}
          double mxt=_tempmx1[indl+1];double myt=_tempmy1[indl+1];double mzt=_tempmz1[indl+1];
          // mxt=0.;myt=0.;mzt=1.;
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          mxt /= mm; myt /= mm; mzt /= mm;
          double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          double mmx=mxt;double mmy=myt; double mmz=mzt; double tnot=alphat;
          if (a1<0.)  vof1=get_vol3D(mmy,mmz,mmx,tnot,0.,-a1);
          if (a2>0.)  vof3=get_vol3D(mmy,mmz,mmx,tnot,1.-a2,a2);


        }
        // new values of c1
        _tempc2[indl+1]+=c_0*(1.-b2+b1)-vof1-vof3;
        _tempc2[i+((j-1+3)%3+(k%3)*3)*nxl+1]+=vof1;
        _tempc2[i+((j+1)%3+(k%3)*3)*nxl+1] +=vof3;
      }
      //eulerian correction
      for (unsigned int ix=0; ix< nxl; ix++) if (_tempc2[((j-1+3)%3+(k%3)*3)*nxl+ix+1] >  MIN_VAL){
          ixyz[0]=ix; ixyz[1]=j-1;ixyz[2]=k;
          get_vel(ixyz,fd,utemp);
          double implicito=(utemp[3]-utemp[2]+utemp[2*idir2+1]-utemp[2*idir2])*st;
          if (fabs(implicito) > MIN_VAL)
            _tempc2[((j-1+3)%3+(k%3)*3)*nxl+ix+1] /=(1.-implicito);
        }
      // contracting  tempc2 -> c2
      if (j>0) CtrRow(_tempc2,&c1[Level],Level,j-1+k*ny,((j-1)%3+(k%3)*3)*nxl);

      //Expanding matrix c1_old ->line  _tempc1
      for (int k2=-1;k2<2;k2++) {
        int indjt=j+2+(k+k2)*ny;   int indjl=((j+2)%3+((k+k2+3)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indjt,indjl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indjt,indjl);
        ExpRow(&_my1[Level],_tempmy1,Level,indjt,indjl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indjt,indjl);
      }
    }
    // new block for new k line
    // storage fo the last line in j
    CtrRow(_tempc2,&c1[Level],Level,ny-1+k*ny,( (ny-1)%3+(k%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.; _tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=k; k0<k+3; k0++) for (unsigned int j0=0; j0<2; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }
  return;
}

// ---------------------------------------------------------------
// Split advection along the z-direction
void  MGSolCC::eulerbalotZ(const unsigned int Level,
                           const double dt,const int idir2) {
  // Set up
  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const double hx=LX/(nxl-2);const double hy=LY/ny;const double hz=LZ/nz;
  unsigned int nxc=_nx[0]-1;  int fd=(nxl-2)/nxc;
  double st=dt/(hz);
  // First line
  for (unsigned int ix=0; ix< nxl*3*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;
    _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
  }
  for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=0; j0<2; j0++) {
      int ind=k0*ny+j0; int indl=(k0*3+j0)*nxl;
      ExpRow(&c1_old[Level],_tempc1,Level,ind,indl);
      ExpRow(&_mx1[Level],_tempmx1,Level,ind,indl);
      ExpRow(&_my1[Level],_tempmy1,Level,ind,indl);
      ExpRow(&_mz1[Level],_tempmz1,Level,ind,indl);
    }

  int ixyz[3]; double utemp[6];
  // Domain cartesian loop
  for (int j=0 ; j<ny ; j++) {
    for (int k=0 ; k<nz ; k++) {
      for (int i=0 ; i<nxl-2 ; i++) {
        ixyz[0]=i;ixyz[1]=j;ixyz[2]=k;

        int indl=i+((k%3)*3+(j%3))*nxl;
        double  vof1 = 0.; double  vof3 = 0.;
        double b1;double b2;
        const double  c_0=_tempc1[indl+1];
        // c=1.
        if ( c_0 >= 1.) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[4]*st;double a2=utemp[5]*st;
          b1=utemp[idir2*2]*st;b2=utemp[idir2*2+1]*st;
          // compression factor
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          vof1 =c_0*MAX(-a1,0.);   vof3 =c_0*MAX(a2,0.);
        }
        //
        else if (c_0 > 0.) {
          get_vel(ixyz,fd,utemp);
          double a1=utemp[4]*st;double a2=utemp[5]*st;
          b1=utemp[idir2*2]*st; b2=utemp[idir2*2+1]*st;
          if (fabs(a1)>1. || fabs(a2)>1.){ printf("\n cfl>1 %e %e \n ",a1,a2);}
          // unit normal
          double mxt=_tempmx1[indl+1];double myt=_tempmy1[indl+1];double mzt=_tempmz1[indl+1];
          double mm=fabs(mxt)+fabs(myt)+fabs(mzt);
          mxt /= mm; myt /= mm; mzt /= mm;

          double alphat=get_alpha(fabs(mxt),fabs(myt),fabs(mzt),c_0) ;
          alphat += MIN(mxt,0.)+MIN(myt,0.)+MIN(mzt,0.);
          double mmx=mxt;double mmy=myt; double mmz=mzt;double tnot=alphat;
          if (a1<0.)  vof1=get_vol3D(mmz,mmx,mmy,tnot,0,-a1);
          if (a2>0.)  vof3=get_vol3D(mmz,mmx,mmy,tnot,1.-a2,a2);
        }
        // new values of c1 (in temp2)
        // ---------------------------------------
        _tempc2[indl+1]+=c_0*(1.-(b2-b1))-vof1-vof3;//vof2;
        _tempc2[i+((j%3)+((k-1+3)%3)*3)*nxl+1]+=vof1; _tempc2[i+(j%3+((k+1)%3)*3)*nxl+1]+=vof3;
      }
      // eulerian correction
      for (unsigned int ix=0; ix< nxl; ix++) if (_tempc2[((j%3)+((k-1+3)%3)*3)*nxl+ix+1]>MIN_VAL){
          ixyz[0]=ix; ixyz[1]=j;ixyz[2]=k-1;
          get_vel(ixyz,fd,utemp);
          double implicito=(utemp[5]-utemp[4]+utemp[idir2*2+1]-utemp[idir2*2])*st;
          if (fabs(implicito) > MIN_VAL)  _tempc2[((j%3)+((k-1+3)%3)*3)*nxl+ix+1] /=(1.-implicito);
        }

      // contraction temp2 -> c1
      if (k>0) CtrRow(_tempc2,&c1[Level],Level,j+(k-1)*ny,((j%3)+((k-1)%3)*3)*nxl);

      //Expanding matrix c1_old (jy+memh)->line (jy+memh)%_mem _tempc1
      for (int j2=-1;j2<2;j2++) {
        int indk=j+j2+(k+2)*ny; int indkl=((j+j2+3)%3+((k+2)%3)*3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indkl);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indkl);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indkl);
        ExpRow(&_mz1[Level],_tempmz1,Level,indk,indkl);
      }
    }
    // new block for new j line
    // storage fo the last line in k
    CtrRow(_tempc2,&c1[Level],Level,j+(nz-1)*ny,(j%3+((nz-1)%3)*3)*nxl);
    // zero the temp storage
    for (unsigned int ix=0; ix< nxl*3*3; ix++) {
      _tempc1[ix]=0.;_tempc2[ix]=0.;
      _tempmx1[ix]=0.;_tempmy1[ix]=0.;_tempmz1[ix]=0.;
    }
    // set up the temp stprage for new k line
    for (unsigned int k0=0; k0<2; k0++) for (unsigned int j0=j; j0<j+3; j0++) {
        int indk= k0*ny+j0;  int indlk=((k0%3)*3+j0%3)*nxl;
        ExpRow(&c1_old[Level],_tempc1,Level,indk,indlk);
        ExpRow(&_mx1[Level],_tempmx1,Level,indk,indlk);
        ExpRow(&_my1[Level],_tempmy1,Level,indk,indlk);         ExpRow(&_mz1[Level],_tempmz1,Level,indk,indlk);
      }
  }

  return;
}

