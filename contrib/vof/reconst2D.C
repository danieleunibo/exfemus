
// -----------------------------------------
// RECONSTRUCTION
// ------------------------------------------

// ================================
/// Metodo Centrale Modificato
void MGSolCC::rec_Cent(const unsigned int Level,int ixy[],
                       double mxy[]) {

  // set up
  int ix=ixy[0];int jy=ixy[1];
  const unsigned int nxl=_nx[Level]+1;
  const unsigned int ny=_ny[Level];
  double min_err = 16.;

  // index def
  const int ind= ix+((jy)   %3)*nxl+1;
  const int indt=ix+((jy+1) %3)*nxl+1;
  const int indb=ix+((jy+2) %3)*nxl+1;

  double ccl[9];  
  ccl[0]=_tempc1[indb-1];ccl[1]=_tempc1[indb];ccl[2]=_tempc1[indb+1];
  ccl[3]=_tempc1[ind-1]; ccl[4]=_tempc1[ind]; ccl[5]=_tempc1[ind+1];
  ccl[6]=_tempc1[indt-1];ccl[7]=_tempc1[indt];ccl[8]=_tempc1[indt+1];

  //  horizontal height function (c_t,c_1,c_b)
  double c_t = ccl[6]+ccl[7]+ ccl[8];
  double c_b = ccl[0]+ccl[1]+ ccl[2];
  //  vertical  height function	(c_r,c_2,c_l)
  double c_r = ccl[2]+ccl[5]+ ccl[8];
  double c_l = ccl[0]+ccl[3]+ ccl[6]; 

  //  case with normal x ---------------
  if (fabs(c_l-c_r) > fabs(c_t-c_b)){
    // from backward centered and forward finite differences
    double myy=0.5*(c_b-c_t);
    double mm = 1./(fabs(myy)+1.);
    double mmx =(c_l>c_r)?mm:-mm;double mmy = myy*mm;
    mxy[0]=mmx; mxy[1]=mmy;
  }
  else{//  case with normal y -------------
    // from backward centered and forward finite differences
    double mxx=0.5*(c_l-c_r);
    double mm = 1./(fabs(mxx)+1.);
    double mmx = (c_t-c_b>0.)? -mm:mm;double mmy = mxx*mm;
    mxy[0]=mmy; mxy[1]=mmx;
  }

  return;
}
// ============================
/// Youngs local reconstruction
void MGSolCC::rec_Young(const unsigned int Level ,int ixy[],
                        double mxy[]) {

  int ix=ixy[0];int jy=ixy[1];
  const unsigned int nxl=_nx[Level]+1;
  const unsigned int ny=_ny[Level];

  // index def
  const int ind= ix+((jy)   %3)*nxl+1;
  const int indt=ix+((jy+1) %3)*nxl+1;
  const int indb=ix+((jy+2) %3)*nxl+1;

  double ccl[9];
  ccl[0]=_tempc1[indb-1];ccl[1]=_tempc1[indb];ccl[2]=_tempc1[indb+1];
  ccl[3]=_tempc1[ind-1]; ccl[4]=_tempc1[ind]; ccl[5]=_tempc1[ind+1];
  ccl[6]=_tempc1[indt-1];ccl[7]=_tempc1[indt];ccl[8]=_tempc1[indt+1];
  double ccc= ccl[4];

  // Young reconstruction
  rec_Young0(mxy,ccl);
  mxy[2] = get_alpha2D(mxy[0],mxy[1],ccc);
}

// ===================================================
void MGSolCC::rec_Young0(double mxy[],double ccl[]) {
  // unit normal  mmx,mmy,mmz;
  double mmx,mmy,mm;
  mmy=ccl[0]+ccl[2]+2.*ccl[1]-(ccl[6]+ccl[8]+ccl[7]*2.);
  mmx=ccl[0]+ccl[6]+ccl[3]*2.-(ccl[2]+ccl[8]+ccl[5]*2.);
  // normalization
  mm=fabs(mmx)+fabs(mmy);
  mmx /= mm; mmy /= mm;
  // storage
  mxy[0]=mmx; mxy[1]=mmy;
  return;
}


// ===============================
// 2D Elvira reconstruction
/// elvira reconstruction by using contracing matrices
/// c_old and mx1,my1,mz1
// --------------------------------------------------
// --------------------------------------------------
void MGSolCC::rec_Elv2D(unsigned int Level,int ixyz[],double mxy[]) {
  // set up
  int ix=ixyz[0];int jy=ixyz[1];
  const unsigned int nxl=_nx[Level]+1;
  const unsigned int ny=_ny[Level];
  double m[2][3],mmx,mmy,alp; double min_err = 16.;

  // index def
  const int ind= ix+((jy)   %3)*nxl+1;
  const int indt=ix+((jy+1) %3)*nxl+1;
  const int indb=ix+((jy+2) %3)*nxl+1;

  //  horizontal height function	(c_t,c_1,c_b)
  double c_t = _tempc1[indt-1]  + _tempc1[indt] + _tempc1[indt+1];
  double c_1 = _tempc1[ind-1]   + _tempc1[ind]  + _tempc1[ind+1];
  double c_b = _tempc1[indb-1]  + _tempc1[indb] + _tempc1[indb+1];
  //  vertical  height function	(c_r,c_2,c_l)
  double c_r = _tempc1[indb+1]  + _tempc1[ind+1]  + _tempc1[indt+1];
  double c_2 = _tempc1[indb]    + _tempc1[ind]    + _tempc1[indt];
  double c_l = _tempc1[indb-1]  + _tempc1[ind-1]  + _tempc1[indt-1];

  double ccc=_tempc1[ind];

  // angular coefficients from backward centered and forward differences
  m[0][0]=c_l-c_2; m[0][1]=0.5*(c_l-c_r);m[0][2]=c_2-c_r;
  m[1][0]=c_b-c_1; m[1][1]=0.5*(c_b-c_t);m[1][2]=c_1-c_t;
  int iyx = 1; if (c_t>c_b)iyx = -1;
  int ixy = 1; if (c_r>c_l)ixy = -1;

  // k=0: y=y(x), k=1: x = x(y)
  for (int k=0; k<=1; k++) {

    int invy = (1-k) *iyx + k*ixy;
    for (int kk=0; kk<=2; kk++) {

      int invx =1;mmy = 1.;mmx = fabs(m[k][kk]);
      if (m[k][kk] < 0.) invx = -1;
      // mx and my are now positive, set  mx + my = 1
      double mm = mmx+ mmy; double mmx1 = mmx/mm; double mmy1 = mmy/mm;
      mmx = (k) *invy*mmy1 + (1-k) *invx*mmx1;
      mmy = (1-k) *invy*mmy1 + (k) *invx*mmx1;
      alp = get_alpha2D(mmx,mmy,ccc);

      // get area diff. with the surrounding cells
      double sum_area_diff = 0.;
      for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
          int indl=(ix+ii+1)+((jy+jj+3)%3)*nxl;
          double area = get_vol2D(mmx,mmy,alp-mmx*ii-mmy*jj,0.,0.,1.,1.);
          sum_area_diff +=(area-_tempc1[indl])*(area-_tempc1[indl]);
        }
      }
      // now update the minimum value and final line
      if (sum_area_diff < min_err-MIN_VAL) {
        min_err = sum_area_diff;mxy[2] = alp; mxy[0] = mmx; mxy[1] = mmy;
      } else if (fabs(sum_area_diff-min_err) < MIN_VAL){ mxy[0] = 0.5*(mxy[0]+mmx); mxy[1] = 0.5*(mxy[1]+mmy);      }
    }
  }
// if (ix == 0 || ix== nxl-2) {mxy[0] = 0.; mxy[1] = 1.;}
  return;
}

// ===============================
/// elvira reconstruction with a restricted direction
void MGSolCC::rec_ElvR(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  int ix=ixyz[0];int jy=ixyz[1];
  const unsigned int nxl=_nx[Level]+1;
  const unsigned int ny=_ny[Level];
  double min_err = 16.;

  // index def
  const int ind= ix+((jy)   %3)*nxl+1;
  const int indt=ix+((jy+1) %3)*nxl+1;
  const int indb=ix+((jy+2) %3)*nxl+1;

  double ccl[9];
  ccl[0]=_tempc1[indb-1];ccl[1]=_tempc1[indb];ccl[2]=_tempc1[indb+1];
  ccl[3]=_tempc1[ind-1]; ccl[4]=_tempc1[ind]; ccl[5]=_tempc1[ind+1];
  ccl[6]=_tempc1[indt-1];ccl[7]=_tempc1[indt];ccl[8]=_tempc1[indt+1];
  double ccc= ccl[4];

  //  horizontal height function (c_t,c_1,c_b)
  double c_t = ccl[6]+ccl[7]+ ccl[8];
  double c_1 = ccl[3]+ccl[4]+ ccl[5];
  double c_b = ccl[0]+ccl[1]+ ccl[2];
  //  vertical  height function	(c_r,c_2,c_l)
  double c_r = ccl[2]+ccl[5]+ ccl[8];
  double c_2 = ccl[1]+ccl[4]+ ccl[7];
  double c_l = ccl[0]+ccl[3]+ ccl[6];

  // candidates in
  double mxc[7],myc[7],sum_area_diff;double mmm[3];
  int ncase=-1;

  //  case with normal x ---------------
  if (fabs(c_l-c_r) > fabs(c_t-c_b)){
    // from backward centered and forward finite differences
    mmm[0]=c_1-c_t; mmm[1]=0.5*(c_b-c_t); mmm[2]=(c_b-c_1);
    for (int km=0;km<3;km++){
      double mm = 1./(fabs(mmm[km])+1.);
      double mmx =(c_l>c_r)?mm:-mm;double mmy = mmm[km]*mm;
      ncase++;mxc[ncase]=mmx; myc[ncase]=mmy;
    }
  }

  else{//  case with normal y -------------
    // from backward centered and forward finite differences
    mmm[0]=c_l-c_2; mmm[1]=0.5*(c_l-c_r); mmm[2]=(c_2-c_r);
    for (int km=0;km<3;km++){
      double mm = 1./(fabs(mmm[km])+1.);
      double mmx = (c_t-c_b>0.)? -mm:mm;double mmy = mmm[km]*mm;
      ncase++; mxc[ncase]=mmy; myc[ncase]=mmx;
    }
  }
  //  Young case ---------------------
  rec_Young0(mxyz,ccl);
  //   storage
  ncase++; mxc[ncase]=mxyz[0]; myc[ncase]=mxyz[1];

  // ============================================
  //              minimization
  int gcase;  // best case
  // get area diff. with the surrounding cells
  for (int lcase=0; lcase<=ncase; lcase++) {
    sum_area_diff = 0.;
    double mmx=mxc[lcase]; double mmy =myc[lcase];
    double alp=get_alpha2D(mmx,mmy,ccc) ;

    for (int jj=-1; jj<=1; jj++) {
      for (int ii=-1; ii<=1; ii++) {
        int indl=(ix+ii+1)+((jy+jj+3)%3)*nxl;
        double area = get_vol2D(mmx,mmy,alp-mmx*ii-mmy*jj,0.,0.,1.,1.);
        sum_area_diff +=(area-_tempc1[indl])*(area-_tempc1[indl]);
      }
    }
    //   the best test
    if (sum_area_diff < min_err) {min_err = sum_area_diff;gcase=lcase;}
  }
  //   now update the minimum value and final line
  mxyz[0] = mxc[gcase]; mxyz[1] = myc[gcase];

  return;
}


// ===============================
/// reconstruction with an height method
void MGSolCC::rec_Height(const unsigned int Level,int ixyz[],double mxyz[]) {

  // set up
  int ix=ixyz[0];int jy=ixyz[1];
  const unsigned int nxl=_nx[Level]+1;
  const unsigned int ny=_ny[Level];
  double min_err = 16.;

  // index def
  const int ind= ix+((jy)   %3)*nxl+1;
  const int indt=ix+((jy+1) %3)*nxl+1;
  const int indb=ix+((jy+2) %3)*nxl+1;

  double ccl[9];
  ccl[0]=_tempc1[indb-1];ccl[1]=_tempc1[indb];ccl[2]=_tempc1[indb+1];
  ccl[3]=_tempc1[ind-1]; ccl[4]=_tempc1[ind]; ccl[5]=_tempc1[ind+1];
  ccl[6]=_tempc1[indt-1];ccl[7]=_tempc1[indt];ccl[8]=_tempc1[indt+1];
  double ccc= ccl[4];

  //  horizontal height function	(c_t,c_1,c_b)
  double c_t = ccl[6]+ccl[7]+ ccl[8];
  double c_1 = ccl[3]+ccl[4]+ ccl[5];
  double c_b = ccl[0]+ccl[1]+ ccl[2];
  //  vertical  height function	(c_r,c_2,c_l)
  double c_r = ccl[2]+ccl[5]+ ccl[8];
  double c_2 = ccl[1]+ccl[4]+ ccl[7];
  double c_l = ccl[0]+ccl[3]+ ccl[6];

  // weight
  double g0,g1,g2;
  //  case with normal x ---------------
  if (fabs(c_l-c_r) < fabs(c_t-c_b)){

    g0=(c_l>.5 )? 1.:0.; g1=1.;g2=(c_r>.5 )? 1.:0.;
    if (g0+g2 ==0.) {g0=1.;g2=1.;}
    // minimization
    double mmxx=(-c_r*(2*g0 +g1)*g2+ c_2*g1*(-g0+g2)+c_l*g0*(g1+2*g2))/(g1*g2+g0*(g1+4*g2));
    double mm = 1./(fabs(mmxx)+ 1.); double mmx = mmxx*mm; double mmy = (c_t > c_b)? -mm:mm;
    mxyz[0] = mmx; mxyz[1] = mmy;
  }

  else{//  case with normal y -------------
    // weight factor
    g0=(c_b>.5 )? 1.:0.; g1=1.;g2=(c_t>.5 )? 1.:0.;
    if (g0+g2 ==0.) {g0=1.;g2=1.;}
    // minimization
    double mmyy=(-c_t*(2*g0 +g1)*g2+ c_1*g1*(-g0+g2)+c_b*g0*(g1+2*g2))/(g1*g2+g0*(g1+4*g2));
    double mm = 1./(fabs(mmyy)+ 1.);double mmx = (c_r > c_l)? -mm:mm;double mmy = mmyy*mm;
    mxyz[0] = mmx; mxyz[1] = mmy;
  }
  return;
}



// --------------------------------------------------
// 2D alpha reconstruction
///  alpha reconstruction (elvira modified for low resolutions)
///  by using contracing matrices c_old and mx1,my1,mz1
// --------------------------------------------------
void MGSolCC::rec_3x3x3mod(const unsigned int nxl,
                           const unsigned int ix,const unsigned int jy,const unsigned int kz,
                           double ccl[]) {
// // void MGSolCC::rec_melv1(unsigned int Level,int ixyz[],
// //                         double mxy[]) {
//
// //   int ix=ixyz[0];int jy=ixyz[1];
// //   double m[2][3],mmx,mmy,alp; double min_err = 16.;
// //   unsigned int nxl=_nx[Level]+1;unsigned int ny=_ny[Level]-1;
// //   int ind[3][3];double crw[3][3];  double ccl[3][3];
//
// //   int tab1[9];int invtab1[9]; tab1[0]=0; tab1[1]=1; tab1[2]=2;
// //   tab1[3]=7; tab1[4]=8; tab1[5]=3;tab1[6]=6; tab1[7]=5; tab1[8]=4;
// //   invtab1[0]=0; invtab1[1]=1; invtab1[2]=2;
// //   invtab1[7]=3; invtab1[8]=4; invtab1[3]=5; invtab1[6]=6; invtab1[5]=7; invtab1[4]=8;
//
//
//   // initial crw and ccl
//   for (unsigned int j0=0;j0<3;j0++)
//     for (unsigned int i0=0;i0<3;i0++) {
//       int ind=ix+i0-1+((jy+j0+2)%3) *nxl;
//       crw[i0][j0]= ccl[i0][j0]=_tempc1[ind];
//     }
//
//   double ccc=_tempc1[ind[1][1]];
//   double amx=_tempmx1[ind[1][1]];
//   double amy=_tempmy1[ind[1][1]];
//
// // correction to color function *********************************
//   int kl[3];kl[1]=1;
//   double valfk[3]={1.,1.,0.};
//   if (3.*fabs(amx) < fabs(amy)) {
//     kl[0]= (amy >=0) ? 0:2;  kl[2]=2-kl[0];
//     // corrected ccl matrix (h column function)
//     for (int j1=0;j1<3;j1++) {
//       for (int i1=0;i1<3;i1++) {
//         double myl=_tempmy1[ind[i1][kl[j1]]]; double mxl=_tempmx1[ind[i1][kl[j1]]];
//         if (myl*amy+ mxl*amx<0.) {
//           ccl[i1][kl[j1]]=valfk[j1];
//           if (j1 == 1 && amy*_tempmy1[ind[i1][kl[0]]]+_tempmx1[ind[i1][kl[0]]]*amx < MIN_VAL)
//             ccl[i1][kl[0]]=1.;
//         }
//       }
//     }
//   } else if (fabs(amx) >fabs(amy) *3.) {
//     kl[0]= (amx >=0) ? 0:2; kl[1]=1; kl[2]=2-kl[0];
//     //corrected crw matrix (h row function)
//     for (int j1=0;j1<3;j1++) {
//       for (int i1=0;i1<3;i1++) {
//         double myl=_tempmy1[ind[kl[i1]][j1]];
//         double mxl=_tempmx1[ind[kl[i1]][j1]];
//         if (mxl*amx+myl*amy < 0.) {
//           crw[kl[i1]][j1]=valfk[i1];
//           if (i1 == 1 && amy*_tempmy1[ind[kl[0]][j1]]+_tempmx1[ind[kl[0]][j1]]*amx < MIN_VAL)
//             crw[kl[0]][j1]=1.;
//         }
//       }
//     }
//   }
//
//
//
//
//
//
// //   else {
// //     int j1= (amx >=0) ? 0:2; int i1= (amy >=0) ? 0:2;
// //     int k1=j1+3*i1; int kp=invtab1[(tab1[k1]+1) %8];
// //     int  km=invtab1[(tab1[k1]+7) %8];
// //
// //     // bottom  line
// //     int ind0=ind-1+j1+ (i1-1) * (nx+1);
// //     double mx0=_mx.Cmp[ind0+1];double my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy < 0.) {
// //       crw[i1][j1]=1.; ccl[i1][j1]=1.;
// //     }
// //
// //     ind0=ind-1+ (kp%3) + (kp/3-1) * (nx+1);
// //     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy <0.) {
// //       crw[kp/3][kp%3]=1.; ccl[kp/3][kp%3]=1.;
// //     }
// //
// //     ind0=ind-1+ (km%3) + (km/3-1) * (nx+1);
// //     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy <0.) {
// //       crw[km/3][km%3]=1.; ccl[km/3][km%3]=1.;
// //     }
// //
// //     // center line
// //     kp=invtab1[(tab1[k1]+2) %8];km=invtab1[(tab1[k1]+6) %8];
// //     ind0=ind-1+ (kp%3) + (kp/3-1) * (nx+1);
// //     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy <0.) {
// //       crw[kp/3][kp%3]=1.; ccl[kp/3][kp%3]=1.;
// //     }
// //
// //     ind0=ind-1+ (km%3) + (km/3-1) * (nx+1);
// //     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy <0.) {
// //       crw[km/3][km%3]=1.; ccl[km/3][km%3]=1.;
// //     }
// //
// //     // top line
// //     int ind2=ind-1+2-j1+ (2-i1-1) * (nx+1);
// //     double mx2=_mx.Cmp[ind2+1];double my2=_my.Cmp[ind2+1];
// //     if (mx2*amx +my2*amy < 0.) {
// //       crw[2-i1][2-j1]=0.;ccl[2-i1][2-j1]=0.;
// //     }
// //
// //     kp=invtab1[(tab1[k1]+5) %8];km=invtab1[(tab1[k1]+3) %8];
// //     ind0=ind-1+ (kp%3) + (kp/3-1) * (nx+1);
// //     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy <0.) {
// //       crw[kp/3][kp%3]=0.; ccl[kp/3][kp%3]=0.;
// //     }
// //
// //     ind0=ind-1+ (km%3) + (km/3-1) * (nx+1);
// //     mx0=_mx.Cmp[ind0+1];my0=_my.Cmp[ind0+1];
// //     if (mx0*amx+my0*amy <0.) {
// //       crw[km/3][km%3]=0.;ccl[km/3][km%3]=0.;
// //     }
// //   }
//   // *****************************
//
//   // height function
//   double c_r = ccl[2][0]+ccl[2][1]+ccl[2][2];
//   double c_2 = ccl[1][0]+ccl[1][1]+ccl[1][2];
//   double c_l = ccl[0][0]+ccl[0][1]+ccl[0][2];
//   double c_t = crw[0][2]+crw[1][2]+crw[2][2];
//   double c_1 = crw[0][1]+crw[1][1]+crw[2][1];
//   double c_b = crw[0][0]+crw[1][0]+crw[2][0];
//
//   /* angular coefficients from backward centered and forward finite  differences */
//   m[0][0] = c_l-c_2; m[0][1] = 0.5* (c_l-c_r); m[0][2] = c_2 - c_r;
//   m[1][0] = c_b - c_1; m[1][1] = 0.5* (c_b-c_t);  m[1][2] = c_1 - c_t;
// //   int iyx = 1; if (c_t > c_b)    iyx = -1;
// //   int ixy = 1; if (c_r > c_l)    ixy = -1;
//   int iyx = 1; if (ccl[0][2]+ccl[1][2]+ccl[2][2] > ccl[0][0]+ccl[1][0]+ccl[2][0])    iyx = -1;
//   int ixy=  1; if (crw[2][0]+crw[2][1]+crw[2][2] > crw[0][0]+crw[0][1]+crw[0][2])    ixy = -1;
//
// // for (int k=0; k<=1; k++) {
//   int k=1; if (fabs(amx) <fabs(amy)) k=0;   {
//
//     // unit normal
//     int invy = (1-k) *iyx + k*ixy;
//     for (int kk=0; kk<=2; kk++) {
//       int invx = 1;   mmy = 1.;     mmx = fabs(m[k][kk]);
//       if (m[k][kk] < 0.) invx = -1;
//       /* mx and my are now positive, set  mx + my = 1 */
//       double mm = mmx+ mmy; double mmx1 = mmx/mm; double mmy1 = mmy/mm;
//       /* get alpha for the equation of the interface */
//       mm = MIN(mmx1,mmy1); alp = get_alpha(mm,1.-mm,0.,ccc);
//       /* now back to the original line */
//       mmx = (k) *invy*mmy1 + (1-k) *invx*mmx1; mmy = (1-k) *invy*mmy1 + (k) *invx*mmx1;
//       alp += MIN(0.,mmx) + MIN(0.,mmy);
//
//       /* get area diff. with the surrounding cells */
//       double sum_area_diff = 0.;     double sum1=0.;
//       if (k==0) {
//         for (int ii=-1; ii<=1; ii++) for (int jj=-1; jj<=1; jj++) {
//             sum1 += ccl[ii+1][jj+1]*ccl[ii+1][jj+1];
//             double area = get_area(mmx,mmy,alp, (double) ii, (double) jj);
//             sum_area_diff += (area-ccl[ii+1][jj+1]) * (area - ccl[ii+1][jj+1]);
//           }
//       } else { // k !=0
//         for (int ii=-1;ii<=1;ii++) for (int jj=-1;jj<=1;jj++) {
//             sum1 +=crw[ii+1][jj+1]*crw[ii+1][jj+1];
//             double area = get_area(mmx,mmy,alp, (double) ii, (double) jj);
//             sum_area_diff += (area-crw[ii+1][jj+1]) * (area-crw[ii+1][jj+1]);
//           }
//       }
//       /* now update the minimum value and final line */
//       sum_area_diff =  sum_area_diff;
//       if (sum_area_diff < min_err) {
//         min_err = sum_area_diff;
//         mxy[2] = alp; mxy[0] = mmx; mxy[1] = mmy;
//       }
//     }
//   }

  return;
}
