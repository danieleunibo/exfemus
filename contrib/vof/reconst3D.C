
// ---------------------
// RECONSTRUCTION
// ---------------------

// generation of the cell recontruction domain
void MGSolCC::rec_3x3x3(const unsigned int nxl,
                        const unsigned int ix,const unsigned int jy,const unsigned int kz,
                        double ccl[]) {

  for (unsigned int k0=0;k0<3;k0++){
    for (unsigned int j0=0;j0<3;j0++){
      for (unsigned int i0=0;i0<3;i0++) {
        int ind=((ix+1)+i0-1+nxl)%nxl+((jy+j0-1+3)%3+((kz+k0-1+3)%3)*3)*nxl;
        ccl[i0+3*j0+9*k0]=_tempc1[ind];
      }
    }
  }
}

// --------------------------------------------------
// generation of the modified cell recontruction domain
// --------------------------------------------------
void MGSolCC::rec_3x3x3mod(const unsigned int nxl,
                           const unsigned int ix,const unsigned int jy,const unsigned int kz,
                           double ccl[]) {

  double mmx,mmy,mmz,alp;
//   unsigned int nxl=_nx[Level]+1;unsigned int ny=_ny[Level]-1;
//   unsigned int nz=_nz[Level]-1;
  int ind[3][3][3]; //double ccl[3][3][3];

  // initial  ccl
  for (unsigned int k0=0;k0<3;k0++)
    for (unsigned int j0=0;j0<3;j0++)
      for (unsigned int i0=0;i0<3;i0++) {
        ind[i0][j0][k0]=((ix+1)+i0-1+nxl)%nxl+((jy+j0-1+3)%3+((kz+k0-1+3)%3)*3)*nxl;
        ccl[i0+3*j0+9*k0]=_tempc1[ind[i0][j0][k0]];
      }

  // key values
  double ccc=_tempc1[ind[1][1][1]];
  double amx=_tempmx1[ind[1][1][1]];
  double amy=_tempmy1[ind[1][1][1]];
  double amz=_tempmz1[ind[1][1][1]];

  if (fabs(amx)>fabs(amy)) {
    // amx > (amy,amz)
    if (fabs(amx)>fabs(amz)) {
      int i12 = (amx > 0.)? 1:-1;
      for (int i1=-1;i1<2;i1++)
        for (int j1=-1;j1<2;j1++)
          for (int k1=-1;k1<2;k1++) {
            if (i1 !=i12 || j1 !=1 || k1 !=1) {
              //int indx=ind+i1+(j1+k1*ny)*nx;
              double mxl=_tempmx1[ind[i1+1][j1+1][k1+1]];
              if ( mxl*amx < 0.) {
                ccl[i1+1+3*(j1+1)+9*(k1+1)]=1.;
                if (fabs(i1-i12) < 2) {
                  //int indx1=ind+i1-i12+(j1+k1*ny)*nx;
                  double mxl2=_tempmx1[ind[i1-i12+1][j1+1][k1+1]];
                  if (mxl2*amx < 0. ) ccl[i1-i12+1+3*(j1+1)+9*(k1+1)]=1.;
                }
              }
            }
          }
    } // end amx > amy &&    amx > amz)
    else {
      int k12 = (amz > 0.)? 1:-1;
      for (int i1=-1;i1<2;i1++)
        for (int j1=-1;j1<2;j1++)
          for (int k1=-1;k1<2;k1++) {
            if (i1 !=1  || j1 !=1 || k1 !=k12) {
              //int indx=ind+i1+(j1+k1*ny)*nx;
              double mzl=_tempmz1[ind[i1+1][j1+1][k1+1]];
              if (mzl*amz < 0.) {
                ccl[i1+1+3*(j1+1)+9*(k1+1)]=1.;
                if (fabs(k1-k12) < 2) {
                  //int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
                  double mzl2=_tempmx1[ind[i1+1][j1+1][k1-k12+1]];
                  if (mzl2*amz < 0.) ccl[i1+1+3*(j1+1)+9*(k1-k12+1)]=1.;
                }
              }
            }
          }
    }
  } else {
    if (fabs(amy)>fabs(amz)) {
      int j12 = (amy > 0.)? 1:-1;
      for (int i1=-1;i1<2;i1++)
        for (int j1=-1;j1<2;j1++)
          for (int k1=-1;k1<2;k1++) {
            if (i1 !=1 || j1 !=j12 || k1 !=1) {
              //  int indx=ind+i1+(j1+k1*ny)*nx;
              double myl=_tempmy1[ind[i1+1][j1+1][k1+1]];
              if (myl*mmy < 0.) {
                ccl[i1+1+3*(j1+1)+9*(k1+1)]=1.;
                if (fabs(j1-j12) < 2) {
                  // int indx1=ind+i1+(j1-j12+k1*ny)*nx;
                  double myl2=_tempmy1[ind[i1+1][j1-j12+1][k1+1]];
                  if (myl2*mmy <0) ccl[i1+1+3*(j1-j12+1)+9*(k1+1)]=1.;
                }
              }
            }
          }
    } else {
      int k12 = (amz > 0.)? 1:-1;
      for (int i1=-1;i1<2;i1++)
        for (int j1=-1;j1<2;j1++)
          for (int k1=-1;k1<2;k1++) {
            if (i1 !=1 || j1 !=1 || k1 !=k12) {
              // int indx=ind+i1+(j1+k1*ny)*nx;
              double mzl=_tempmz1[ind[i1+1][j1+1][k1+1]];
              if (mzl*amz <0) {
                ccl[i1+1+3*(j1+1)+9*(k1+1)]=1.;
                if (fabs(k1-k12) < 2) {
                  // int indx1=ind+i1+(j1+(k1-k12)*ny)*nx;
                  double mzl2=_tempmz1[ind[i1+1][j1+1][k1-k12+1]];
                  if (mzl2*amz <0) ccl[i1+1+3*(j1+1)+9*(k1-k12+1)]=1.;
                }
              }
            }
          }
    }
  }
  return;
}


// -----------------------------------------------------------------
void MGSolCC::rec_Cent(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  const unsigned int nxl=_nx[Level]+1; const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const int ix=ixyz[0]; const int jy=ixyz[1];const int kz=ixyz[2];

  //local ccl
  double ccl[27];  REC3x3x3(nxl,ix,jy,kz,ccl);

  // sum of columns in colccx colccy colccz
  double colccx[3][3];double colccy[3][3];double colccz[3][3];
  for (int l1=0;l1<3;l1++){
    for (int l2=0;l2<3;l2++){
      colccx[l1][l2]=0.;colccy[l1][l2]=0.;colccz[l1][l2]=0.;
      for (int l=0;l<3;l++)  {
        colccx[l1][l2] += ccl[l+3*l1+9*l2];
        colccy[l1][l2] += ccl[l1+3*l+9*l2];
        colccz[l1][l2] += ccl[l1+3*l2+9*l];
      }
    }
  }
  double ccc=ccl[13];
  double g0[9];
  double mmyy,mmzz,mmxx,mmx,mmy,mmz,alphat; double m1,m2;
  double s0[9];

  // direction derivatives  top-bottom
  double  dirxtb=0; double dirytb  =0; double  dirztb =0;
  for (int l3=0;l3<3;l3++) {
    dirxtb += (colccy[2][l3]-colccy[0][l3]);
    dirytb += (colccz[l3][2]-colccz[l3][0]);
    dirztb += (colccx[l3][2]-colccx[l3][0]);
  }
  // test case
  int icase=0;
  if (fabs(dirxtb) > fabs(dirytb)){
    if (fabs(dirxtb) > fabs(dirztb))   icase=0;
    else  icase=2;
  }
  else {
    if (fabs(dirytb) > fabs(dirztb)) icase=1;
    else icase=2;
  }
  // weights
  g0[0]=1.;g0[1]=1.;g0[2]=1.;
  g0[3]=1.;g0[4]=1.;g0[5]=1.;
  g0[6]=1.;g0[7]=1.;g0[8]=1.;

  switch (icase){

    case 0: //  case normal x
      //  color sums
      s0[0]=colccx[0][0]; s0[1]=colccx[1][0]; s0[2]=colccx[2][0];
      s0[3]=colccx[0][1]; s0[4]=colccx[1][1]; s0[5]=colccx[2][1];
      s0[6]=colccx[0][2]; s0[7]=colccx[1][2]; s0[8]=colccx[2][2];


      rec_Cent0(dirxtb,g0,s0,mxyz);
      //  order y z x
      mmx=mxyz[0]; mmy=mxyz[1]; mmz=mxyz[2];
      break;

    case 1: // case y normal ---------------
      //  color sums
      s0[0]=colccy[0][0]; s0[1]=colccy[1][0]; s0[2]=colccy[2][0];
      s0[3]=colccy[0][1]; s0[4]=colccy[1][1]; s0[5]=colccy[2][1];
      s0[6]=colccy[0][2]; s0[7]=colccy[1][2]; s0[8]=colccy[2][2];

      rec_Cent0(dirytb,g0,s0,mxyz); // order yzx
      // order x z y
      mmx=mxyz[1]; mmy=mxyz[0]; mmz=mxyz[2];
      break;

    case 2://z direction --
      //  color sums
      s0[0]=colccz[0][0]; s0[1]=colccz[1][0]; s0[2]=colccz[2][0];
      s0[3]=colccz[0][1]; s0[4]=colccz[1][1]; s0[5]=colccz[2][1];
      s0[6]=colccz[0][2]; s0[7]=colccz[1][2]; s0[8]=colccz[2][2];

      rec_Cent0(dirztb,g0,s0,mxyz); // order yzx
      // order x  y z
      mmx=mxyz[1]; mmy=mxyz[2]; mmz=mxyz[0];
      break;

    default:
      printf("Error switch central reconstr");
      break;
  }
  // get alpha
  alphat=get_alpha(fabs(mmx),fabs(mmy),fabs(mmz),ccc) ;
  alphat += MIN(0.,mmx) + MIN(0.,mmy)+ MIN(0.,mmz);

  //   now update the minimum value and final line
  mxyz[3] =alphat; mxyz[0] = mmx; mxyz[1] =mmy ;mxyz[2] = mmz;

  return;
}


// -------------------------------------
/// reconstruction by using minimizations
void MGSolCC::rec_Cent0(double dirxtb,double g0[],double s0[],double mxyz[]) {

  double mmxx, mmyy,mmzz; double mmx, mmy,mmz;
  // solution sistem  mmx mmy (mmz) alph  3x3
  double a00=0.125*(g0[0]+g0[1]+9*g0[2]+g0[3]+g0[4]+9*g0[5]+g0[6]+g0[7]+9*g0[8]);
  double a01=0.125*(g0[0]-g0[1]-3*g0[2]-g0[3]+g0[4]+3*(g0[5]-g0[6] +g0[7]+3*g0[8]));
  double a02=0.25*(g0[0]-g0[1]-3*g0[2]+g0[3]-g0[4]-3*g0[5]+g0[6]-g0[7]-3*g0[8]);

  double a10=a01;
  double a11=0.125*(g0[0]+g0[1]+g0[2]+g0[3]+g0[4]+g0[5]+9*(g0[6]+g0[7]+g0[8]));
  double a12=0.25*(g0[0]+g0[1]+g0[2]-g0[3]-g0[4]-g0[5]-3*g0[6]-3*g0[7]-3*g0[8]);

  double a20=0.25*(-g0[0]+g0[1]+3*g0[2]-g0[3]+g0[4]+3*g0[5]-g0[6]+g0[7]+3*g0[8]);
  double a21=0.25*(-g0[0]-g0[1]-g0[2]+g0[3]+g0[4]+g0[5]+3*(g0[6] + g0[7] + g0[8]));
  double a22=0.5*(-g0[0]-g0[1]-g0[2]-g0[3]-g0[4]-g0[5]-g0[6]-g0[7]-g0[8]);

  double det=1./(-a02*a11*a20 + a01*a12*a20 + a02*a10*a21 -
                 a00*a12*a21 - a01*a10*a22 + a00*a11*a22);
  double c00=0.25*(-g0[0]*s0[0]+g0[1]*s0[1]+ 3*g0[2]*s0[2]-g0[3]*s0[3]+g0[4]*s0[4]+
                   3*g0[5]*s0[5]-g0[6]*s0[6]+g0[7]*s0[7]+3*g0[8]*s0[8]);
  double c10=0.25*(-g0[0]*s0[0]-g0[1]*s0[1]-g0[2]*s0[2]+g0[3]*s0[3]+g0[4]*s0[4]+
                   g0[5]*s0[5]+3*g0[6]*s0[6]+3*g0[7]*s0[7]+3*g0[8]*s0[8]);
  double c20=0.5*(g0[0]*s0[0]+g0[1]*s0[1]+g0[2]*s0[2]+g0[3]*s0[3]+
                  g0[4]*s0[4]+g0[5]*s0[5]+g0[6]*s0[6]+g0[7]*s0[7] +g0[8]*s0[8]);

  mmxx = (dirxtb >0.)?-1.:1.;
  mmyy = -det*(-a12*a21*c00+a11*a22*c00 + a02*a21*c10 - a01*a22*c10 - a02*a11*c20 +
               a01*a12*c20);
  mmzz = det*(-a12*a20*c00 +  a10*a22*c00 + a02*a20*c10
              - a00*a22*c10 - a02*a10*c20 +  a00*a12*c20);


// normalization  mmx+mmy+mmz=1
  mmz =mmzz;mmx =mmxx;mmy =mmyy;
  double mm=1.+fabs(mmy)+fabs(mmz);
  mmx /= mm; mmy /= mm; mmz /= mm;
  // storage case
  mxyz[0]=mmx; mxyz[1]=mmy; mxyz[2]=mmz;
  return;
}

// =================================================================
void MGSolCC::rec_Elv1(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  const unsigned int nxl=_nx[Level]+1; const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;double min_err = 16.;
  const int ix=ixyz[0]; const int jy=ixyz[1];const int kz=ixyz[2];
  // initial  local ccl

  double ccl[27]; REC3x3x3(nxl,ix,jy,kz,ccl);

  // sum of columns in colccx colccy colccz
  double colccx[3][3];double colccy[3][3];double colccz[3][3];
  for (int l1=0;l1<3;l1++)
    for (int l2=0;l2<3;l2++){
      colccx[l1][l2]=0.;colccy[l1][l2]=0.;colccz[l1][l2]=0.;
      for (int l=0;l<3;l++)  {
        colccx[l1][l2] += ccl[l+3*l1+9*l2];
        colccy[l1][l2] += ccl[l1+3*l+9*l2];
        colccz[l1][l2] += ccl[l1+3*l2+9*l];
      }
    }

  double ccc=ccl[13];
  double g0[9];
  double mmyy,mmzz,mmxx,mmx,mmy,mmz,alphat; double m1,m2;
  double s0[9];

  // direction derivatives  top-bottom
  double  dirxtb=0; double dirytb  =0; double  dirztb =0;
  for (int l3=0;l3<3;l3++) {
    dirxtb += (colccy[2][l3]-colccy[0][l3]);
    dirytb += (colccz[l3][2]-colccz[l3][0]);
    dirztb += (colccx[l3][2]-colccx[l3][0]);
  }

  // ----------------------------------------------------
  // minimization cases
  // ----------------------------------------------------
  // candidates in
  double mxc[10],myc[10],mzc[10],alc[10],sum_area_diff[10];
  //  case 0  height meth in x
  //  case 1  height meth in y
  //  case 2  height  meth in z
  //  case 3 Young
  // ------------------------------------
  //  case with normal x ========================
  int icase=0;
  //  color sums
  s0[0]=colccx[0][0]; s0[1]=colccx[1][0]; s0[2]=colccx[2][0];
  s0[3]=colccx[0][1]; s0[4]=colccx[1][1]; s0[5]=colccx[2][1];
  s0[6]=colccx[0][2]; s0[7]=colccx[1][2]; s0[8]=colccx[2][2];
  // weights
  g0[0]=1.;g0[1]=2.;g0[2]=1.;
  g0[3]=2.;g0[4]=4.;g0[5]=2.;
  g0[6]=1.;g0[7]=2.;g0[8]=1.;
  for (int ki=0;ki<9;ki++){if (s0[ki]<0.5)   g0[ki]=0.;}
  // recontruction
  rec_Cent0(dirxtb,g0,s0,mxyz);
  //  storage (order y z x)
  mxc[icase]=mxyz[0]; myc[icase]=mxyz[1]; mzc[icase]=mxyz[2];

  // case y with normal =========================
  icase=1;
  //  color sums
  s0[0]=colccy[0][0]; s0[1]=colccy[1][0]; s0[2]=colccy[2][0];
  s0[3]=colccy[0][1]; s0[4]=colccy[1][1]; s0[5]=colccy[2][1];
  s0[6]=colccy[0][2]; s0[7]=colccy[1][2]; s0[8]=colccy[2][2];
// weights
  g0[0]=1.;g0[1]=2.;g0[2]=1.;
  g0[3]=2.;g0[4]=4.;g0[5]=2.;
  g0[6]=1.;g0[7]=2.;g0[8]=1.;
  for (int ki=0;ki<9;ki++){if (s0[ki]<0.5)   g0[ki]=0.;}
  // recontruction
  rec_Cent0(dirytb,g0,s0,mxyz); // order yzx
  // storage (but order x z y)
  mxc[icase]=mxyz[1]; myc[icase]=mxyz[0]; mzc[icase]=mxyz[2];

  // case with z normal =========================
  icase=2;
  //  color sums
  s0[0]=colccz[0][0]; s0[1]=colccz[1][0]; s0[2]=colccz[2][0];
  s0[3]=colccz[0][1]; s0[4]=colccz[1][1]; s0[5]=colccz[2][1];
  s0[6]=colccz[0][2]; s0[7]=colccz[1][2]; s0[8]=colccz[2][2];
  // weights
  g0[0]=1.;g0[1]=2.;g0[2]=1.;
  g0[3]=2.;g0[4]=4.;g0[5]=2.;
  g0[6]=1.;g0[7]=2.;g0[8]=1.;

  for (int ki=0;ki<9;ki++){if (s0[ki]<0.5)   g0[ki]=0.;}
  // recontruction
  rec_Cent0(dirztb,g0,s0,mxyz); // order yzx
  // storage (but order x  y z)
  mxc[icase]=mxyz[1]; myc[icase]=mxyz[2]; mzc[icase]=mxyz[0];

  //  Young case ==========================
  icase=3;
  // recontruction
  rec_Young0(mxyz,ccl);
  // storage
  mxc[icase]=mxyz[0]; myc[icase]=mxyz[1]; mzc[icase]=mxyz[2];

  // ============================================
  //       test (minimization of the directions
  // ============================================
  int gcase;  // best case
  // get area diff. with the surrounding cells
  for (int lcase=0; lcase<=3; lcase++) {
    sum_area_diff[lcase] = 0.; double mmx=mxc[lcase];
    double mmy =myc[lcase];    double mmz=mzc[lcase];
    alphat=get_alpha(fabs(mmx),fabs(mmy),fabs(mmz),ccc) ;
    alphat += MIN(0.,mmx) + MIN(0.,mmy)+ MIN(0.,mmz);

    for (int kk=-1; kk<=1; kk++) {
      for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
          double diffarea =get_vol3D(mmx,mmy,mmz,alphat-ii*mmx-jj*mmy-kk*mmz,0.,1.)
                           -ccl[ii+1+3*(jj+1)+9*(kk+1)];
          sum_area_diff[lcase] += diffarea* diffarea;
        }
      }
    }
    //   the best test
    if (sum_area_diff[lcase] < min_err-MIN_VAL) {
      min_err = sum_area_diff[lcase];gcase=lcase;
    }
  }
//   now update the minimum value and final line
  mxyz[3] = alc[gcase]; mxyz[0] = mxc[gcase]; mxyz[1] = myc[gcase];mxyz[2] = mzc[gcase];

  return;
}


// =================================================================
void MGSolCC::rec_Height(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  const unsigned int nxl=_nx[Level]+1; const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;double min_err = 16.;
  const int ix=ixyz[0]; const int jy=ixyz[1];const int kz=ixyz[2];
  // initial  local ccl

  double ccl[27]; REC3x3x3(nxl,ix,jy,kz,ccl);

  // sum of columns in colccx colccy colccz
  double colccx[3][3];double colccy[3][3];double colccz[3][3];
  for (int l1=0;l1<3;l1++)
    for (int l2=0;l2<3;l2++){
      colccx[l1][l2]=0.;colccy[l1][l2]=0.;colccz[l1][l2]=0.;
      for (int l=0;l<3;l++)  {
        colccx[l1][l2] += ccl[l+3*l1+9*l2];
        colccy[l1][l2] += ccl[l1+3*l+9*l2];
        colccz[l1][l2] += ccl[l1+3*l2+9*l];
      }
    }

  double ccc=ccl[13];
  double g0[9];
  double mmyy,mmzz,mmxx,mmx,mmy,mmz,alphat; double m1,m2;
  double s0[9];

  // direction derivatives  top-bottom
  double  dirxtb=0; double dirytb  =0; double  dirztb =0;
  for (int l3=0;l3<3;l3++) {
    dirxtb += (colccy[2][l3]-colccy[0][l3]);
    dirytb += (colccz[l3][2]-colccz[l3][0]);
    dirztb += (colccx[l3][2]-colccx[l3][0]);
  }

  // candidates in
  double mxyzc[10];
  // ------------------------------------
  //  case with normal x ========================
 if ( fabs(dirxtb) > fabs(dirytb) && fabs(dirxtb) > fabs(dirztb)){
  //  color sums
  s0[0]=colccx[0][0]; s0[1]=colccx[1][0]; s0[2]=colccx[2][0];
  s0[3]=colccx[0][1]; s0[4]=colccx[1][1]; s0[5]=colccx[2][1];
  s0[6]=colccx[0][2]; s0[7]=colccx[1][2]; s0[8]=colccx[2][2];
  // weights
  g0[0]=1.;g0[1]=2.;g0[2]=1.;
  g0[3]=2.;g0[4]=4.;g0[5]=2.;
  g0[6]=1.;g0[7]=2.;g0[8]=1.;

  for (int ki=0;ki<9;ki++){if (s0[ki]<0.5)   g0[ki]=0.;}
  // recontruction
  rec_Cent0(dirxtb,g0,s0,mxyzc);
  //  storage (order y z x)
  mxyz[0]=mxyzc[0]; mxyz[1]=mxyzc[1]; mxyz[2]=mxyzc[2];
}
  // case y with normal =========================
else if ( fabs(dirytb)> fabs(dirxtb) && fabs(dirytb)> fabs(dirztb)){

  //  color sums
  s0[0]=colccy[0][0]; s0[1]=colccy[1][0]; s0[2]=colccy[2][0];
  s0[3]=colccy[0][1]; s0[4]=colccy[1][1]; s0[5]=colccy[2][1];
  s0[6]=colccy[0][2]; s0[7]=colccy[1][2]; s0[8]=colccy[2][2];
// weights
  g0[0]=1.;g0[1]=2.;g0[2]=1.;
  g0[3]=2.;g0[4]=4.;g0[5]=2.;
  g0[6]=1.;g0[7]=2.;g0[8]=1.;
  for (int ki=0;ki<9;ki++){if (s0[ki]<0.5)   g0[ki]=0.;}
  // recontruction
  rec_Cent0(dirytb,g0,s0,mxyzc); // order yzx
  // storage (but order x z y)
mxyz[0]=mxyzc[1]; mxyz[1]=mxyzc[0]; mxyz[2]=mxyzc[2];
}
else{
  // case with z normal =========================
  //  color sums
  s0[0]=colccz[0][0]; s0[1]=colccz[1][0]; s0[2]=colccz[2][0];
  s0[3]=colccz[0][1]; s0[4]=colccz[1][1]; s0[5]=colccz[2][1];
  s0[6]=colccz[0][2]; s0[7]=colccz[1][2]; s0[8]=colccz[2][2];
  // weights
  g0[0]=1.;g0[1]=2.;g0[2]=1.;
  g0[3]=2.;g0[4]=4.;g0[5]=2.;
  g0[6]=1.;g0[7]=2.;g0[8]=1.;

  for (int ki=0;ki<9;ki++){if (s0[ki]<0.5)   g0[ki]=0.;}
  // recontruction
  rec_Cent0(dirztb,g0,s0,mxyzc); // order yzx
  // storage (but order x  y z)
mxyz[0]=mxyzc[1]; mxyz[1]=mxyzc[2]; mxyz[2]=mxyzc[0];
} 

  return;
}



// ============================================
void MGSolCC::rec_ElvR(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  const unsigned int nxl=_nx[Level]+1; const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;double min_err = 16.;
  const int ix=ixyz[0]; const int jy=ixyz[1];const int kz=ixyz[2];
  // initial  local ccl

  double ccl[27]; REC3x3x3(nxl,ix,jy,kz,ccl);

  // sum of columns in colccx colccy colccz
  double colccx[3][3];double colccy[3][3];double colccz[3][3];
  for (int l1=0;l1<3;l1++)
    for (int l2=0;l2<3;l2++){
      colccx[l1][l2]=0.;colccy[l1][l2]=0.;colccz[l1][l2]=0.;
      for (int l=0;l<3;l++)  {
        colccx[l1][l2] += ccl[l+3*l1+9*l2];
        colccy[l1][l2] += ccl[l1+3*l+9*l2];
        colccz[l1][l2] += ccl[l1+3*l2+9*l];
      }
    }

  double ccc=ccl[13];
 // double g0[9];
  double mmyy,mmzz,mmxx,mmx,mmy,mmz,alphat; double m1,m2;
  double s0[9];

  // direction derivatives  top-bottom
  double  dirxtb=0; double dirytb  =0; double  dirztb =0;
  for (int l3=0;l3<3;l3++) {
    dirxtb += (colccy[2][l3]-colccy[0][l3]);
    dirytb += (colccz[l3][2]-colccz[l3][0]);
    dirztb += (colccx[l3][2]-colccx[l3][0]);
  }

  // ----------------------------------------------------
  // minimization cases
  // ----------------------------------------------------
  // candidates in
  double mxc[40],myc[40],mzc[40],alc[40],sum_area_diff[30];
  double mmm[9][2];
  int ncase=-1;

  //  case with normal x ========================
  if ( fabs(dirxtb) > fabs(dirytb) && fabs(dirxtb) > fabs(dirztb)){

    s0[0]=colccx[0][0]; s0[1]=colccx[1][0]; s0[2]=colccx[2][0];
    s0[3]=colccx[0][1]; s0[4]=colccx[1][1]; s0[5]=colccx[2][1];
    s0[6]=colccx[0][2]; s0[7]=colccx[1][2]; s0[8]=colccx[2][2];
    //  horizontal height function	(c_t,c_1,c_b)
    double c_t =s0[6]+s0[7]+s0[8];
    double c_1 =s0[3]+s0[4]+s0[5];
    double c_b =s0[0]+s0[1]+s0[2];
    //  vertical  height function	(c_r,c_2,c_l)
    double c_r =s0[2]+s0[5]+s0[8];
    double c_2 =s0[1]+s0[4]+s0[7];
    double c_l =s0[0]+s0[3]+s0[6];

    mmxx = (dirxtb>0.)? -3.:3.;
    // from backward centered and forward finite differences
    mmm[0][0]=c_l-c_2;      mmm[0][1]=c_b-c_1;
    mmm[1][0]=0.5*(c_l-c_r);mmm[1][1]=c_b-c_1;
    mmm[2][0]=(c_2-c_r);    mmm[2][1]=c_b-c_1;
    mmm[3][0]=c_l-c_2;      mmm[3][1]=0.5*(c_b-c_t);
    mmm[4][0]=0.5*(c_l-c_r);mmm[4][1]=0.5*(c_b-c_t);
    mmm[5][0]=(c_2-c_r);    mmm[5][1]=0.5*(c_b-c_t);
    mmm[6][0]=c_l-c_2;      mmm[6][1]=c_1-c_t;
    mmm[7][0]=0.5*(c_l-c_r);mmm[7][1]=c_1-c_t;
    mmm[8][0]=(c_2-c_r);    mmm[8][1]=c_1-c_t;

    for (int km=0;km<9;km++){ ncase++;
      double mm = 1./(fabs(mmm[km][0])+fabs(mmm[km][1])+fabs(mmxx));
      double mmx = mmxx*mm;double mmy = mmm[km][0]*mm;
      double mmz = mmm[km][1]*mm;
      mxc[ncase]=mmx; myc[ncase]=mmy; mzc[ncase]=mmz;
     
    }
  }
  // case y with normal =========================
  else if ( fabs(dirytb)> fabs(dirxtb) && fabs(dirytb)> fabs(dirztb)){

    //  color sums
    s0[0]=colccy[0][0]; s0[1]=colccy[1][0]; s0[2]=colccy[2][0];
    s0[3]=colccy[0][1]; s0[4]=colccy[1][1]; s0[5]=colccy[2][1];
    s0[6]=colccy[0][2]; s0[7]=colccy[1][2]; s0[8]=colccy[2][2];

    //  horizontal height function	(c_t,c_1,c_b)
    double c_t =s0[6]+s0[7]+s0[8];
    double c_1 =s0[3]+s0[4]+s0[5];
    double c_b =s0[0]+s0[1]+s0[2];
    //  vertical  height function	(c_r,c_2,c_l)
    double c_r =s0[2]+s0[5]+s0[8];
    double c_2 =s0[1]+s0[4]+s0[7];
    double c_l =s0[0]+s0[3]+s0[6];

    mmxx = (dirytb>0.)? -3.:3.;

    // from backward centered and forward finite differences
    mmm[0][0]=c_l-c_2;      mmm[0][1]=c_b-c_1;
    mmm[1][0]=0.5*(c_l-c_r);mmm[1][1]=c_b-c_1;
    mmm[2][0]=(c_2-c_r);    mmm[2][1]=c_b-c_1;
    mmm[3][0]=c_l-c_2;      mmm[3][1]=0.5*(c_b-c_t);
    mmm[4][0]=0.5*(c_l-c_r);mmm[4][1]=0.5*(c_b-c_t);
    mmm[5][0]=(c_2-c_r);    mmm[5][1]=0.5*(c_b-c_t);
    mmm[6][0]=c_l-c_2;      mmm[6][1]=c_1-c_t;
    mmm[7][0]=0.5*(c_l-c_r);mmm[7][1]=c_1-c_t;
    mmm[8][0]=(c_2-c_r);    mmm[8][1]=c_1-c_t;

    for (int km=0;km<9;km++){ ncase++;
      double mm = 1./(fabs(mmm[km][0])+fabs(mmm[km][1])+fabs(mmxx));
      double mmx = mmxx*mm;double mmy = mmm[km][0]*mm;
      double mmz = mmm[km][1]*mm;
      mxc[ncase]=mmy; myc[ncase]=mmx; mzc[ncase]=mmz;
    
    }
  }
  // case with z normal =========================
  //if ( fabs(dirztb)> fabs(dirxtb) && fabs(dirztb)> fabs(dirytb)){
else{
//  color sums
    s0[0]=colccz[0][0]; s0[1]=colccz[1][0]; s0[2]=colccz[2][0];
    s0[3]=colccz[0][1]; s0[4]=colccz[1][1]; s0[5]=colccz[2][1];
    s0[6]=colccz[0][2]; s0[7]=colccz[1][2]; s0[8]=colccz[2][2];

//  horizontal height function	(c_t,c_1,c_b)
    double c_t =s0[6]+s0[7]+s0[8];
    double c_1 =s0[3]+s0[4]+s0[5];
    double c_b =s0[0]+s0[1]+s0[2];
    //  vertical  height function	(c_r,c_2,c_l)
    double c_r =s0[2]+s0[5]+s0[8];
    double c_2 =s0[1]+s0[4]+s0[7];
    double c_l =s0[0]+s0[3]+s0[6];

//  double ccc=s0[4];
    double mm;double mmx;double mmy;double mmz;
    mmxx = (dirztb>0.)? -3.:3.;
    // from backward centered and forward finite differences
    mmm[0][0]=c_l-c_2;      mmm[0][1]=c_b-c_1;
    mmm[1][0]=0.5*(c_l-c_r);mmm[1][1]=c_b-c_1;
    mmm[2][0]=(c_2-c_r);    mmm[2][1]=c_b-c_1;
    mmm[3][0]=c_l-c_2;      mmm[3][1]=0.5*(c_b-c_t);
    mmm[4][0]=0.5*(c_l-c_r);mmm[4][1]=0.5*(c_b-c_t);
    mmm[5][0]=(c_2-c_r);    mmm[5][1]=0.5*(c_b-c_t);
    mmm[6][0]=c_l-c_2;      mmm[6][1]=c_1-c_t;
    mmm[7][0]=0.5*(c_l-c_r);mmm[7][1]=c_1-c_t;
    mmm[8][0]=(c_2-c_r);    mmm[8][1]=c_1-c_t;

    for (int km=0;km<9;km++){ ncase++;
      double mm = 1./(fabs(mmm[km][0])+fabs(mmm[km][1])+fabs(mmxx));
      double mmx = mmxx*mm;double mmy = mmm[km][0]*mm;
      double mmz = mmm[km][1]*mm;
      mxc[ncase]=mmy; myc[ncase]=mmz; mzc[ncase]=mmx;
     
    }

  }

  //  Young case ==========================
 ncase++;
//recontruction
 rec_Young0(mxyz,ccl);
//  storage
  mxc[ncase]=mxyz[0]; myc[ncase]=mxyz[1]; mzc[ncase]=mxyz[2];

  // ============================================
  //       test (minimization of the directions
  // ============================================
  int gcase;  // best case
  // get area diff. with the surrounding cells
  for (int lcase=0; lcase<=ncase; lcase++) {
    sum_area_diff[lcase] = 0.; double mmx=mxc[lcase];
    double mmy =myc[lcase];    double mmz=mzc[lcase];
    alphat=get_alpha(fabs(mmx),fabs(mmy),fabs(mmz),ccc) ;
    alphat += MIN(0.,mmx) + MIN(0.,mmy)+ MIN(0.,mmz);

    for (int kk=-1; kk<=1; kk++) {
      for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
          double diffarea =get_vol3D(mmx,mmy,mmz,alphat-ii*mmx-jj*mmy-kk*mmz,0.,1.)
                           -ccl[ii+1+3*(jj+1)+9*(kk+1)];
          sum_area_diff[lcase] += diffarea*diffarea;
        }
      }
    }
    //   the best test
    if (sum_area_diff[lcase] < min_err-MIN_VAL) {
      min_err = sum_area_diff[lcase];gcase=lcase;
    }
  }
//   now update the minimum value and final line
  mxyz[3] = alc[gcase]; mxyz[0] = mxc[gcase]; mxyz[1] = myc[gcase];mxyz[2] = mzc[gcase];

  return;
}


// -------------------------------------
/// Wrapper for Young reconsruction (see Young0)
void MGSolCC::rec_Young(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  const unsigned int nxl=_nx[Level]+1;
  const unsigned int ny=_ny[Level]-1;
  const unsigned int nz=_nz[Level]-1;
  const int ix=ixyz[0]; const int jy=ixyz[1];const int kz=ixyz[2];

// initial  local ccl
  double ccl[27]; REC3x3x3(nxl,ix,jy,kz,ccl);
//   unit normal by Young ->mxyz
  rec_Young0(mxyz,ccl);

  return;
}

// -------------------------------------
/// reconstruction by using Young stensil,
///.brutal finite difference across the cell with weight
/// {[1,2,1][2,4,2][1,2,1]} over each plane
void MGSolCC::rec_Young0(double mxyz[],double ccl[]) {

  // unit normal  mmx,mmy,mmz;
  double m1,m2,mm;
// z direction
  m1=(ccl[0]+ccl[6]+ccl[2]+ccl[8])+
     (ccl[3]+ccl[5]+ccl[1]+ccl[7])*2.+
     ccl[4]*4.;
  m2=(ccl[18]+ccl[24]+ccl[20]+ccl[26])+
     (ccl[21]+ccl[23]+ccl[19]+ccl[25])*2.+
     ccl[22]*4.;
  double mmz=m1-m2;
  // y-direction
  m1=(ccl[0]+ccl[18]+ccl[2]+ccl[20])+
     (ccl[9]+ccl[11]+ccl[1]+ccl[19])*2.+
     ccl[10]*4.;
  m2=(ccl[6]+ccl[24]+ccl[8]+ccl[26])+
     (ccl[7]+ccl[15]+ccl[17]+ccl[25])*2.+
     ccl[16]*4.;
  double mmy=m1-m2;
  // x-direction
  m1=(ccl[0]+ccl[6]+ccl[18]+ccl[24])+
     (ccl[3]+ccl[9]+ccl[15]+ccl[21])*2.+
     ccl[12]*4.;
  m2=(ccl[2]+ccl[8]+ccl[20]+ccl[26])+
     (ccl[5]+ccl[11]+ccl[17]+ccl[23])*2.+
     ccl[14]*4.;
  double mmx=m1-m2;
  // normalization  mmx+mmy+mmz=1
  mm=fabs(mmx)+fabs(mmy)+fabs(mmz)+MIN_VAL;
  mmx /= mm; mmy /= mm; mmz /= mm;
  // storage
  mxyz[0]=mmx; mxyz[1]=mmy; mxyz[2]=mmz;

  return;
}





