

// ------------------------------------------------------------------
/// \f$  QVector \right  QVector : \quad \vec{u}(t+dt) =Adv( \vec{u}(t) \f$
/// Advection by fluxing the volume with  VOF bubble unsplit method:
/// xy direction
void MGSolCC::AdvVelLXY(unsigned int Level,Matrix &savemat,double dt) {
//void MGSolCC::AdvVelLXY(unsigned int Level, double dt) {
  // geometry level
  unsigned int nxl=_nx[Level]+1;unsigned int ny=_ny[Level]-1;
  unsigned int nxc=_nx[0]-1;unsigned int nyc=_ny[0]-1;
  double hx=LX/(nxl-2);  double hy=LY/ny;
  int fd=(nxl-2)/nxc;
  double stx=dt/hx; double sty=dt/hy;
  const int *dofn=_invnode_dof;
  double uu[3];double utemp[4];int ixy[3];

  // First line
  for (unsigned int ix=0; ix< nxl *3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.;_tempmx1[ix]=0.;_tempmy1[ix]=0.;
  }
  for (unsigned int k=0; k<2; k++) {
    ExpRow(&c1_old[Level],_tempc1,Level,k,k*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,k,k*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,k,k*nxl);
  }

  // Cartesian domain loop
  for (uint jy=0; jy<ny; jy++) {

    // Update  temporary line index (jy-1)
    int indup=((jy+2)%3)*nxl;
    // test for interface
    if (M__GetLen(&c1_old[Level],jy) +
        M__GetLen(&c1_old[Level],(jy+1)%(ny+1)) +
        M__GetLen(&c1_old[Level],(jy+2)%(ny+1)) > 0) {
      // column ix

      for (uint ix=0; ix<nxl-2; ix++)  {
        ixy[0]=ix;ixy[1]=jy;
        get_vel(ixy,fd,utemp);
        uu[0] =utemp[0]*stx; uu[2] =-utemp[1]*stx; // ul and -ur
        double vb=utemp[2]*sty; double vt =utemp[3]*sty;
        // xy  Advection and volume computation *********************
        // geometric parameters
        double du = -uu[2]-uu[0]; double aa = 1.-du;double ai=1./aa;
        double dyt = MAX(0.,vt); double dyb = -MIN(0.,vb);
        double y0 = MAX(0.,vb);  double  y1 = MIN(1.,1.+vt);
        double dyc = y1 - y0;
        double x0 = MAX(0.,-uu[0]);   double x1 = MIN(1.,1.+uu[2]);
        uu[1] = x1 - x0;
        double x6 = MAX(0.0,uu[0]*ai); //double x7 = MIN(1.0,1.0 + (-uu[2]*ai));

        // loop left(0) center(1) right(2) cell
        for (int ii=0;ii<3;ii++) {
          // color function
          int indl= ix+ii-1+(jy%3)*nxl; double cc1=_tempc1[indl+1];
          if (cc1 > 0. && uu[ii]>0.) {
            // set up volume computations
            double vcold = uu[ii]; double dx = vcold*ai;
            double vct = dx*dyt; double vcb = dx*dyb;
            double vcc = dx*dyc;
            if (cc1 < 1.) {
              // unit normal and alpha
              double mx=_tempmx1[indl+1];  double my=_tempmy1[indl+1];
              double alpha=get_alpha2D(mx,my,cc1);
              // computation of the volume fraction
              double w1=0.5*(1.-uu[0])*(2-ii)*(1-ii)+x0*(2-ii)*(ii);
              vcold = get_vol2D(mx,my,alpha,w1,0.,uu[ii],1.0);
              vct = vcb = vcc = 0.;
              if (vcold > 0.) {
                double w2=0.5*(1.-dx)*(ii-1)*(ii)+x6*(2-ii)*(ii);
                my *= ai;  alpha += mx*(uu[0]+ii-1.) + my*vb;  mx *= aa;
                if (vt > 0.) vct += get_vol2D(mx,my,alpha,w2,1.,dx,dyt);
                if (vb < 0.) vcb += get_vol2D(mx,my,alpha,w2,vb,dx,dyb);
                vcc +=  get_vol2D(mx,my,alpha,w2,y0,dx,dyc);
              }
            }
            // storage of the volume fraction
            _tempc2[ix+((jy+1)%3)*nxl+1] += vct;_tempc2[ix+((jy)%3)*nxl+1] += vcc;
            _tempc2[ix+((jy+2)%3)*nxl+1] += vcb;
          }
        }

      } // global grid ix
    } // if len>0
    // Contracting tempc2 -> c1
    if (jy>0) {
      CtrRow(_tempc2,&savemat,Level,jy-1,indup);
    }
    // Expanding matrix c1_old jy+2 ->line indup  _tempc1
    ExpRow(&c1_old[Level],_tempc1,Level,jy+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,jy+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,jy+2,indup);
  }// global grid jy
  // Last line contraction
  CtrRow(_tempc2,&savemat,Level,ny-1,((ny-1)%3)*nxl);

  return;
}

// ------------------------------------------------------------------
/// \f$  QVector \right  QVector : \quad \vec{u}(t+dt) =Adv( \vec{u}(t) \f$
/// Advection by fluxing the volume with  VOF bubble unsplit method:
/// xy direction
void MGSolCC::AdvVelLXY1(unsigned int Level,Matrix &savemat,double dt) {
  // geometry level
  unsigned int nx=_nx[Level]-1;unsigned int nxl=nx+2;
  unsigned int ny=_ny[Level]-1;
  unsigned int nxc=_nx[0]-1;unsigned int nyc=_ny[0]-1;
  double hx=LX/(nxl-2);  double hy=LY/ny;
  int n_nodes= (nxc+1) * (nyc+1);
  int fd=(nxl-2)/nxc;
  double stx=dt/hx; double sty=dt/hy;
  const int *dofn=_invnode_dof;double utemp[4];int ixy[3];

  // First line
  for (unsigned int ix=0; ix< nxl *3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.;_tempmx1[ix]=0.;_tempmy1[ix]=0.;
  }
  for (unsigned int k=0; k<2; k++) {
    ExpRow(&c1_old[Level],_tempc1,Level,k,k*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,k,k*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,k,k*nxl);
  }
  // Cartesian domain loop
  for (uint jy=0; jy<ny; jy++) {

    // Update  temporary line index (jy-1)
    int indup=((jy+2)%3)*nxl;
    // test for interface
    if (M__GetLen(&c1_old[Level],jy) +
        M__GetLen(&c1_old[Level],(jy+1)%(ny+1)) +
        M__GetLen(&c1_old[Level],(jy+2)%(ny+1)) > 0) {
      // column ix
      for (uint ix=0; ix<nx; ix++)  {
        // Color function stensil  left top center bottom right
        int indc=ix+(jy%3)*nxl;    double ccc=_tempc1[indc+1];
        int indt=ix+((jy+1)%3)*nxl;int indb=ix+indup; // not use
        int indl=ix-1+(jy%3)*nxl;  double ccl=_tempc1[indl+1];
        int indr=ix+1+(jy%3)*nxl;  double ccr=_tempc1[indr+1];

        // velocity field  ul(left),ur(right),vt(top),vb(bot)
        ixy[0]=ix;ixy[1]=jy;
        get_vel(ixy,fd,utemp);
        double ul = utemp[0]*stx;
        double ur = utemp[1]*stx; // ul and ur
        double vb = utemp[2]*sty;
        double vt = utemp[3]*sty;

        // xy  Advection and volume computation *********************
        double vcold, dx,vct,vcb,vcc, x0,x1,dxold ;
        double du = ur - ul;  double    aa = 1. - du;   double   ai = 1./aa;
        double dyt = MAX(0.,vt);  double    dyb = - MIN(0.,vb);
        double  y0 = MAX(0.,vb);   double    y1 = MIN(1.,1.+vt);     double  dyc = y1 - y0;

        // left side ----------------------------
        if (ccl > MIN_VAL && ul > MIN_VAL) {
          vcold = ul;	 dx = ul*ai;
          vct = dx*dyt; vcb = dx*dyb;vcc = dx*dyc;
          if (ccl < 1.-MIN_VAL) {
            // unit normal
            double mx=_tempmx1[indl+1]; double my=_tempmy1[indl+1];
            double alpha=get_alpha2D(mx,my,ccl);
            // Volume
            vcold = get_vol2D(mx,my,alpha,1.-ul,0.,ul,1.);
            vct = vcb = vcc = 0.;
            if (vcold > MIN_VAL) {
              my *= ai;  alpha += mx*(ul-1.) + my*vb;  mx *= aa;
              if (vt > MIN_VAL)  vct = get_vol2D(mx,my,alpha,0.,1.,dx,dyt);
              if (vb < -MIN_VAL)  vcb = get_vol2D(mx,my,alpha,0.,vb,dx,dyb);
              vcc =  get_vol2D(mx,my,alpha,0.,y0,dx,dyc);
            }
          }
          _tempc2[indt+1] += vct; _tempc2[indc+1] += vcc;_tempc2[indb+1] += vcb;
        }
        // right side ---------------------------------
        if (ccr > MIN_VAL && ur < -MIN_VAL) {
          vcold = -ur;	dx = -ur*ai;
          vct = dx*dyt;	vcb = dx*dyb;	vcc = dx*dyc;
          if (ccr < 1.-MIN_VAL) {
            double mx=_tempmx1[indr+1]; double my=_tempmy1[indr+1];
            double alpha=get_alpha2D(mx,my,ccr);
            // Volume
            vcold = get_vol2D(mx,my,alpha,0.,0.,-ur,1.);
            vct = vcb = vcc = 0.;
            if (vcold > MIN_VAL) {
              my *= ai; alpha += mx*(ul+1.) + my*vb;   mx *= aa;
              if (vt > MIN_VAL) vct = get_vol2D(mx,my,alpha,1.-dx,1.,dx,dyt);
              if (vb < -MIN_VAL) vcb = get_vol2D(mx,my,alpha,1.-dx,vb,dx,dyb);
              vcc =  get_vol2D(mx,my,alpha,1.-dx,y0,dx,dyc);
            }
          }
          _tempc2[indt+1] += vct; _tempc2[indc+1] += vcc;_tempc2[indb+1] += vcb;
        }
        // center ----------------------------------------
        if (ccc > MIN_VAL) { // c > 0.
          x0 = MAX(0.,-ul);x1 = MIN(1.,1.-ur);
          dxold = x1 - x0; vcold = dxold;	dx = dxold*ai;
          vct = dyt*dx;	vcb = dyb*dx;	vcc = dyc*dx;
          if (ccc < 1.-MIN_VAL) { // c < 1.
            double mx=_tempmx1[indc+1]; double my=_tempmy1[indc+1];
            double alpha=get_alpha2D(mx,my,ccc);
            // Volume
            vcold = get_vol2D(mx,my,alpha,x0,0.0,dxold,1.0);
            vct = vcb = vcc = 0.0;
            if (vcold > MIN_VAL) {
              my *= ai; alpha += mx*ul + my*vb; mx *= aa;
              x0 = MAX(0.0,ul*ai); x1 = MIN(1.0,1.0 + (ur*ai));
              if (vt > MIN_VAL) vct = get_vol2D(mx,my,alpha,x0,1.0,dx,dyt);
              if (vb < -MIN_VAL) vcb = get_vol2D(mx,my,alpha,x0,vb,dx,dyb);
              vcc =  get_vol2D(mx,my,alpha,x0,y0,dx,dyc);
            }
          }
          _tempc2[indt+1] += vct;_tempc2[indc+1] += vcc;_tempc2[indb+1] += vcb;
        }
      } // global grid ix
    } // if
    // Contraction and expansion  in  indup line
    if (jy>0)  CtrRow(_tempc2,&savemat,Level,jy-1,indup);
    ExpRow(&c1_old[Level],_tempc1,Level,jy+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,jy+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,jy+2,indup);
  }// global grid jy
  // final line
  CtrRow(_tempc2,&savemat,Level,ny-1,((ny-1)%3)*nxl);

  return;
}

// ------------------------------------------------------------------
/// \f$  QVector \right  QVector : \quad \vec{u}(t+dt) =Adv( \vec{u}(t) \f$
/// Advection by fluxing the volume with  VOF bubble unsplit method:
/// yx direction
void MGSolCC::AdvVelLYX1(unsigned int Level,Matrix &savemat, double dt) {
  // geometry level
  unsigned int nxl=_nx[Level]+1;unsigned int ny=_ny[Level]-1;
  unsigned int nxc=_nx[0]-1;unsigned int nyc=_ny[0]-1;
  double hx=LX/(nxl-2);  double hy=LY/ny;
  int n_nodes= (nxc+1) * (nyc+1);
  int fd=(nxl-2)/nxc;
  double stx=dt/hx; double sty=dt/hy;
  const int *dofn=_invnode_dof;double utemp[4];int ixy[3];

  // First line
  for (unsigned int ix=0; ix< nxl *3; ix++) {
    _tempc1[ix]=0.; _tempc2[ix]=0.;_tempmx1[ix]=0.;_tempmy1[ix]=0.;
  }
  for (unsigned int k=0; k<2; k++) {
    ExpRow(&c1_old[Level],_tempc1,Level,k,k*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,k,k*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,k,k*nxl);
  }

  // Cartesian domain
  for (uint jy=0; jy<ny; jy++) {
    // Update  temporary line index (jy-1)
    int indup=((jy+2)%3)*nxl;
    // test for interface
    if (M__GetLen(&c1_old[Level],jy) +
        M__GetLen(&c1_old[Level],(jy+1)%(ny+1)) +
        M__GetLen(&c1_old[Level],(jy+2)%(ny+1)) > 0) {
      // column ix
      for (uint ix=0; ix<nxl-2; ix++)  {
        // Color function stensil  left top center bottom right
        int indc=ix+(jy%3)*nxl;    double ccc=_tempc1[indc+1];
        int indt=ix+((jy+1)%3)*nxl;double cct=_tempc1[indt+1];
        int indb=ix+indup;         double ccb=_tempc1[indb+1];
        int indl=ix-1+(jy%3)*nxl;  int indr=ix+1+(jy%3)*nxl;

        // velocity field  ul(left),ur(right),vt(top),vb(bot)
        ixy[0]=ix;ixy[1]=jy;
        get_vel(ixy,fd,utemp);
        double ul =utemp[0]*stx; double ur =utemp[1]*stx; // ul and -ur
        double vb=utemp[2]*sty; double vt =utemp[3]*sty;

        // Volume computations
        double du,dv,aa,ai,y0,y1,x0,x1,dx,dy;
        double dyc,dxc,dyt,dyb,dxr,dxl,vcold,vct,vcb,vcr,vcl,vcc;
        double dxold,dyold;

        // Advection yx
        dv = vt - vb;      aa = 1. - dv;     ai = 1./aa;
        dxr = MAX(0.0,ur);     dxl = - MIN(0.0,ul);
        x0 = MAX(0.0,ul);      x1 = MIN(1.0,1.0+ur);    dxc = x1 - x0;
        // bottom side -------------------------------
        if (ccb > 0. && vb > 0.) {
          vcold = vb;	dy = vb*ai;
          vcr = dxr*dy;	vcl = dxl*dy;	vcc = dxc*dy;

          if (ccb < 1.) {
            // unit normal
            double mx=_tempmx1[indb+1]; double my=_tempmy1[indb+1];
            double alpha=get_alpha2D(mx,my,ccb);
            // Volume
            vcold = get_vol2D(mx,my,alpha,0.,1.-vb,1.,vb);
            vcr = vcl = vcc = 0.;
            if (vcold > 0.) {
              mx *= ai;  alpha += my*(vb-1.) + mx*ul;  my *= aa;
              if (ur > 0.) vcr = get_vol2D(mx,my,alpha,1.,0.,dxr,dy);
              if (ul < 0.) vcl = get_vol2D(mx,my,alpha,ul,0.,dxl,dy);
              vcc =  get_vol2D(mx,my,alpha,x0,0.,dxc,dy);
            }
          }
          _tempc2[indl+1] += vcl; _tempc2[indc+1] += vcc;
          _tempc2[indr+1] += vcr;
        }
        // top -----------------------------
        if (cct > 0. && vt < 0.) {
          vcold = -vt;	dy = -vt*ai;
          vcr = dxr*dy;	vcl = dxl*dy;	vcc = dxc*dy;
          if (cct < 1.) {
            // unit normal
            double mx=_tempmx1[indt+1]; double my=_tempmy1[indt+1];
            double alpha=get_alpha2D(mx,my,cct);
            // Volume
            vcold = get_vol2D(mx,my,alpha,0.,0.,1.,-vt);
            vcr = vcl = vcc = 0.;
            if (vcold > 0.) {
              mx *= ai;  alpha += my*(vb+1.0) + mx*ul;   my *= aa;
              if (ur > 0.) vcr = get_vol2D(mx,my,alpha,1.,1.-dy,dxr,dy);
              if (ul < 0.) vcl = get_vol2D(mx,my,alpha,ul,1.-dy,dxl,dy);
              vcc =  get_vol2D(mx,my,alpha,x0,1.-dy,dxc,dy);
            }
          }
          _tempc2[indl+1] += vcl; _tempc2[indc+1] += vcc;_tempc2[indr+1] += vcr;
        }
        // center  --------------------------------------
        if (ccc > 0.) {
          y0 = MAX(0.,-vb);	y1 = MIN(1.,1.-vt);
          dyold = y1 - y0;	vcold = dyold;	dy = dyold*ai;
          vcr = dxr*dy;	vcl = dxl*dy;	vcc = dxc*dy;
          if (ccc < 1.) {
            // unit normal
            double mx=_tempmx1[indc+1]; double my=_tempmy1[indc+1];
            double alpha=get_alpha2D(mx,my,ccc);
            // volume
            vcold = get_vol2D(mx,my,alpha,0.,y0,1.,dyold);
            vcr = vcl = vcc = 0.;
            if (vcold > 0.) {
              mx *= ai;	    alpha += my*vb + mx*ul;    my *= aa;
              y0 = MAX(0.,vb*ai);	    y1 = MIN(1.,1. + (vt*ai));
              if (ur > 0.)     vcr = get_vol2D(mx,my,alpha,1.,y0,dxr,dy);
              if (ul < 0.)    vcl = get_vol2D(mx,my,alpha,ul,y0,dxl,dy);
              vcc =  get_vol2D(mx,my,alpha,x0,y0,dxc,dy);
            }
          }
          _tempc2[indl+1] += vcl; _tempc2[indc+1] += vcc; _tempc2[indr+1] += vcr;
        }
      }
    }
    if (jy>0) {
      // storage in indup line
      CtrRow(_tempc2,&savemat,Level,jy-1,indup);
    }
    ExpRow(&c1_old[Level],_tempc1,Level,jy+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,jy+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,jy+2,indup);
  }// global grid jy
  // last line
  CtrRow(_tempc2,&savemat,Level,ny-1,((ny-1)%3)*nxl);

  return;
}

// ---------------------------------------------------------------
// Split advection along the x-direction:
/// Standard advection in one direction by using the normal from
/// (mx1 my1) previously computed by Young function
void  MGSolCC::lagrangeX(const unsigned int Level,
                         const double dt) {

  // Set up
  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nxc=_nx[0]-1;  const unsigned int nyc=_ny[0]-1;
  const unsigned int n_nodes=(nxl-1)*(ny+1);
  const double hx=LX/(nxl-2);const double hy=LY/ny;
  const int fd=(nxl-2)/nxc;
  double stx=dt/hx;double utemp[4];int ixy[2];

  // First line
  for (unsigned int ix=0; ix< nxl *3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.; _tempmx1[ix]=0.;_tempmy1[ix]=0.;
  }
  for (unsigned int j0=0; j0<2; j0++) {
    int indl=j0*nxl;
    ExpRow(&c1_old[Level],_tempc1,Level,j0,indl);
    ExpRow(&_mx1[Level],_tempmx1,Level,j0,indl);
    ExpRow(&_my1[Level],_tempmy1,Level,j0,indl);
  }

  // Domain loop (i,j)=[0,nz]x[0,ny]
  for (uint j=0; j<ny ; j++) {

    // Update temporary line -> j-1 && solved line -> j
    int indup=((j+2)%3)*nxl; int indout=(j%3)*nxl;
    // index i loop
    for (uint i=0; i<nxl-2 ; i++) {
      // color function
      int indc=i+(j%3)*nxl; const double  ccc=_tempc1[indc+1];

      // velocity field  ul(left),ur(right),vt(top),vb(bot)
      ixy[0]=i;ixy[1]=j;
      get_vel(ixy,fd,utemp);
      double a1 =utemp[0]*stx; double a2 =utemp[1]*stx; // ul and -ur

      // Volume computation ******************************
      double vof1 = 0.;double  vof2 = 0.;double vof3 = 0.;// c=0.
      if ( ccc > 1.-MIN_VAL) { // c=1.
        vof1 =ccc* MAX(-a1,0.);    vof3 =ccc* MAX(a2,0.);
        vof2 =ccc* ( 1.- MAX(a1,0.)-MAX(-a2,0.));
      } else if (ccc > 0. ) { // 0.< c <1.
        // unit normal
        double mx=_tempmx1[indc+1];double my=_tempmy1[indc+1];
        // intersection alpha
        double alpha=get_alpha2D(mx,my,ccc);
        // Advected Volume
        mx=mx/(1.-a1+a2);alpha +=mx*a1;
        if (a1<0.) vof1=get_vol2D(mx,my,alpha,a1,0.,-a1,1.);
        if (a2>0.) vof3=get_vol2D(mx,my,alpha,1.,0.,a2,1.);
        vof2=get_vol2D(mx,my,alpha,MAX(0.,a1),0.,1-MAX(0.,-a2)-MAX(0.,a1),1.);
      }
      // temporary storage
      _tempc2[indc+1]+=vof2;_tempc2[indc]+=vof1; _tempc2[indc+2]+=vof3;
    }
    // Contracting _tempc2 -> c1
    CtrRow(_tempc2,&c1[Level],Level,j,indout);
    //Expanding matrix c1_old j+2->line indup _tempc1
    ExpRow(&c1_old[Level],_tempc1,Level,j+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,j+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,j+2,indup);
  }
  return;
}

// ---------------------------------------------------------------
// Split advection along the y-direction
/// Standard advection in one direction by using the normal from
/// (mx1,my1) previously computed by Young function.
/// The matrices mx1,my1 are in contracted form
void  MGSolCC::lagrangeY(const unsigned int Level,
                         const double dt) {

  const unsigned int nxl=_nx[Level]+1;  const unsigned int ny=_ny[Level]-1;
  const unsigned int nxc=_nx[0]-1;  const unsigned int nyc=_ny[0]-1;
  const unsigned int n_nodes=(nxc+1)*(nyc+1);
  const double hx=LX/(nxl-2);const double hy=LY/ny;
  const int fd=(nxl-2)/nxc;
  double sty=dt/(hy);double utemp[4];int ixy[2];

  // First line
  for (unsigned int ix=0; ix< nxl*3; ix++) {
    _tempc1[ix]=0.;_tempc2[ix]=0.;_tempmx1[ix]=0.;_tempmy1[ix]=0.;
  }
  for (unsigned int j0=0; j0<2; j0++) {
    ExpRow(&c1_old[Level],_tempc1,Level,j0,j0*nxl);
    ExpRow(&_mx1[Level],_tempmx1,Level,j0,j0*nxl);
    ExpRow(&_my1[Level],_tempmy1,Level,j0,j0*nxl);
  }

  // Cartesian domain loop  (i,j)=[0,nx]x[0,ny]
  for (uint j=0 ; j<ny ; j++) {

    // Update temporary line -> j-1 && solved line -> j
    int indup=((j+2)%3)*nxl;  int indout=((j-1)%3)*nxl;
    // index i loop
    for (uint i=0; i<nxl-2 ; i++) {
      // color function
      const int indc=i+(j%3)*nxl; const double  ccc=_tempc1[indc+1];
      const int indb=i+((j+2)%3)*nxl;const int indt=i+((j+1)%3)*nxl;
      ixy[0]=i; ixy[1]=j;
      // velocity field  ul(left),ur(right),vt(top),vb(bot) -> utemp
      get_vel(ixy,fd,utemp);
      double a1=utemp[2]*sty; double a2 =utemp[3]*sty;  //vb vt

      // advected volume ---------------------------------
      double vof1=0.;double vof2=0.;double vof3=0.; // c=0.
      if ( ccc >= 1.) {     // c=1.
        vof1 =ccc* MAX(-a1,0.); vof3 =ccc* MAX(a2,0.);
        vof2 =ccc* ( 1. - MAX(a1,0.) - MAX(-a2,0.));
      } else if (ccc > 0.) { // 0.< c <1.
        // unit normal
        double mx=_tempmx1[indc+1];double my=_tempmy1[indc+1];
        // intersect alpha
        double alpha=get_alpha2D(mx,my,ccc) ;
        // advected volume
        my=my/(1.-a1+a2);alpha +=my*a1;
        if (a1<0.)  vof1=get_vol2D(my,mx,alpha,a1,0.,-a1,1.);
        if (a2>0.)  vof3=get_vol2D(my,mx,alpha,1.,0.,a2,1.);
        vof2=get_vol2D(my,mx,alpha,MAX(0.,a1),0.,1-MAX(0.,-a2)-MAX(0.,a1),1.);
      }
      // new values of c1
      _tempc2[indb+1]+=vof1; _tempc2[indc+1]+=vof2; _tempc2[indt+1] +=vof3;
    }

    // Contracting  tempc2[((j-1)%3)*nxl] -> c2[j-1] --------------
    if (j>0) CtrRow(_tempc2,&c1[Level],Level,j-1,indout);
    // Expanding matrix c1_old[j+2] ->line  _tempc1[((j+2)%3)*nxl]
    ExpRow(&c1_old[Level],_tempc1,Level,j+2,indup);
    ExpRow(&_mx1[Level],_tempmx1,Level,j+2,indup);
    ExpRow(&_my1[Level],_tempmy1,Level,j+2,indup);
  }
  // storage of the last line
  CtrRow(_tempc2,&c1[Level],Level,ny-1,((ny-1)%3)*nxl);
  return;
}

