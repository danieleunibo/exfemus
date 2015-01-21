/****************************************************************************/
/*                                itersolv.c                                */
/****************************************************************************/
/*                                                                          */
/* ITERative SOLVers for systems of linear equations                        */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#include <stddef.h>

#ifdef __cplusplus
/* # include <cmath> */
# include <math.h>
#else
# include <math.h>
#endif

#include "itersolv.h"
#include "elcmp.h"
#include "errhandl.h"
#include "operats.h"
#include "rtc.h"
#include "copyrght.h"

/* number of GMRES steps bevore restart */
static int GMRESSteps = 10;


QVector *RK4Iter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                 PrecondProcType Dummy, _LPDouble dt) {
    int Iter;
    _LPReal bNorm;
    size_t Dim,i;
    QVector krk4;
    QVector krk3;
    QVector krk2;
    QVector krk1;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);

    V_Constr(&krk1,"krk1",Dim,Normal, _LPTrue);
    V_Constr(&krk2,"krk2",Dim,Normal, _LPTrue);
    V_Constr(&krk3,"krk3",Dim,Normal, _LPTrue);
    V_Constr(&krk4,"krk4",Dim,Normal, _LPTrue);

    Asgn_VV(&krk1,Add_VV(Mul_QV(A,x),b));
    /*  Asgn_VV(&krk2,Add_VV(Mul_QV(A,Add_VV(x,Mul_SV(0.5*dt,&krk1))),b)); */
    /*   Asgn_VV(&krk3,Add_VV(Mul_QV(A,Add_VV(x,Mul_SV(0.5*dt,&krk2))),b)); */
    /*   Asgn_VV(&krk4,Add_VV(Mul_QV(A,Add_VV(x,Mul_SV(dt,&krk3))),b)); */




    Asgn_VV(&krk2,Add_VV(&krk1,Mul_SV(0.5*dt,Mul_QV(A,&krk1))));
    Asgn_VV(&krk3,Add_VV(&krk1,Mul_SV(0.5*dt,Mul_QV(A,&krk2))));
    Asgn_VV(&krk4,Add_VV(&krk1,Mul_SV(dt,Mul_QV(A,&krk3))));

    for (i=1; i<=Dim; i++) {
        double val=dt* (V__GetCmp(&krk1,i) +2.*V__GetCmp(&krk2,i) +
                        2.*V__GetCmp(&krk3,i) +V__GetCmp(&krk4,i)) /6.;
        V__SetCmp(b,i,val);
    }

    V_Destr(&krk1);
    V_Destr(&krk2);
    V_Destr(&krk3);
    V_Destr(&krk4);
    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);


    // printf("Lamp \n");
    // exit(0);


}

/* ************************************************************* */
QVector *LumpIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                  PrecondProcType Dummy, _LPDouble Omega) {
    int Iter,i12;
    _LPReal bNorm;
    size_t Dim,i;
    QVector r;
    QVector diagl;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&diagl, "diagl", Dim, Normal, _LPTrue);
    V_SetAllCmp(&r,1.);
    Asgn_VV(&diagl, Mul_QV(A,&r));


//   printf(" diag ");
//    for(i12=1;i12<=Dim;i12++) printf(" %g ",V_GetCmp(&diagl,i12));
    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        while (!RTCResult(Iter, l2Norm_V(&r), bNorm, LumpIterId)
                && Iter < MaxIter) {
            Iter++;
            for (i=1;i<=Dim;i++) V_SetCmp(&r,i,V_GetCmp(&r,i) /V_GetCmp(&diagl,i));
            /* x(i+1) = x(i) + Omega * D^(-1) * r */
            AddAsgn_VV(x, Mul_SV(Omega,&r));
            if (Q_KerDefined(A))     OrthoRightKer_VQ(x, A);
            /* r = b - A * x(i) */
            if (Iter < MaxIter)
                Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        }
    }

    V_Destr(&r);
    V_Destr(&diagl);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);


    // printf("Lamp \n");
    // exit(0);


}

/* ============================================================= */
QVector *VankaNSIterD(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                     PrecondProcType Dummy, _LPDouble Omega) {
    int Iter,i,Len,Lenj,j,k,lk,lj, dof_u,lsj,kk;
    _LPReal bNorm;
    double res,aa_j,val,bb_i,bb_j,udd,bb,ppp,utemp,ftot,btot,p_i_old;
    size_t Dim;
//   double ax[10000];double bx[10000];int iu[10000];int *iju;
    double *ax;
    double *bx;/*int iu[10000];*/
    int *iju;
    QVector dd;
    QVector r;
    Omega=0.;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);    // dimension
    iju = (int *) malloc((Dim+2) *sizeof(int));
    ax = (double *) malloc((Dim+2) *sizeof(double));
    bx = (double *) malloc((Dim+2) *sizeof(double));
    // temporary x
    V_Constr(&dd, "dd", Dim, Normal,_LPTrue);
    V_SetAllCmp(&dd,1.);
    // ax=Diag_Q(A)
    Asgn_VV(&dd, Mul_QV(Diag_Q(A),&dd));          //dd= diag(A)
    dof_u=Dim;
    while (fabs(V__GetCmp(&dd,dof_u)) == 0.) dof_u--;
    printf(" pdof %d \n",dof_u);
    //   V_SetAllCmp(&dd,0.);

//   Iter = 0;
    // r = b - A * x(i)
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_SetAllCmp(&r,0.);
    /* r = b - A * x(i) */
    bNorm=1.;
    Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
    Iter = 0;
    int MaxIter1=8;

    while (!RTCResult(Iter, l2Norm_V(&r), bNorm, VankaIterId) && Iter < MaxIter1) {
        Iter++;
        /*


        while (Iter < 8) {
          Iter++;*/

        // pressure Vertex loop
        for (i=dof_u+1;i<=Dim;i++) {
            //if(fabs(Q__GetVal(A,i,i))< 1.e-20){
            //vertex system ++++++++++++++++++++++++++++++++++
            Len=Q__GetLen(A,i);    // dim local system  //row length
            //  local system for pressure p(i)
            //  x(i+1)=D^(-1) [(D-A)x(i)+b]
            p_i_old= V__GetCmp(x,i);
            kk=0;
            ftot=0.;
            btot=0.;
            // loop over the term in eq p(i)
            for (lj=0;lj<Len;lj++) {
                j=Q__GetPos(A,i,lj);

                if (j<=dof_u) {   // no pressure terms
                    bb_i= Q__GetVal(A,i,lj);    // divergence term B^T(i,j)

                    // iu[kk]=lj;
                    kk++;
                    Lenj=Q__GetLen(A,j);

                    // extraction line j
                    for (lsj=0;lsj<Lenj;lsj++)  iju[Q__GetPos(A,j,lsj)]=lsj;
                    aa_j= Q__GetVal(A,j,iju[j]);   // extraction A(j,j)
                    bb_j= Q__GetVal(A,j,iju[i]);   // extraction B(j,i)

                    // residual
                    // res(j) = b(j)+D(j,k)*u(k)+B(j,k)p(k=i)
                    //          -A(j,lk)*u(lk)-B(j,lk)p(lk!=i)
                    res= V__GetCmp(b,j) + aa_j*V__GetCmp(x,j) +bb_j*p_i_old;
                    for (lk=0;lk<Lenj;lk++) res -= Q__GetVal(A,j,lk) * V__GetCmp(x,Q__GetPos(A,j,lk));

                    // storage
                    ax[lj]=res/aa_j; // ax = D^(-1)*res
                    bx[lj]=bb_j/aa_j; //  bx = D^-1*b1
                    // divergence
//                     ftot +=ax[lj]*bb_i;
//                     btot += bx[lj]*bb_i; // b1^T*px
                    
                    ftot +=ax[lj]*bb_i;
                    btot += bx[lj]*bb_i; // b1^T*px
                }

            }
            // for (lj=0;lj<kk;lj++)
            // pressure   ppp=Omega*p_i_old+(1-Omega)*ftot/btot;
            ppp=ftot/btot;
            V__SetCmp(x,i,ppp);
	       printf("  ppp= %g ft=%g btot=%g \n",ppp,ftot,btot);
            // velocity
            for (lj=0;lj<kk;lj++)
                V_SetCmp(x,Q__GetPos(A,i,lj),ax[lj]-bx[lj]*ppp);
        }
        if (Iter < MaxIter)
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
//     }
    }
    free(iju);
    free(ax);
    free(bx);
    V_Destr(&dd);
    V_Destr(&r);
    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);
    // printf("I am VankaNS \n");
    // exit(0);
    return (x);
}
/* ============================================================= */
QVector *VankaNSIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                      PrecondProcType Dummy, _LPDouble Omega) {
    int Iter,i,Len,Lenj,j,k,il,lk,lj,jl,Lenl,k1, dof_u,lsj,kk,i1,i2,j1,j2;
    _LPReal bNorm;
    int pos1[1000],loc_idx[1000], count;int *indx;
    double res,aa_j,val,val2,bb_i,bb_j,udd,ppp,sum,utemp,ftot,btot,p_i_old;
    size_t Dim;
//   double ax[10000];double bx[10000];int iu[10000];int *iju;

    QVector xb;
    QVector xf;
    QVector bb;
    QVector bf;
    QMatrix Mat;
    QVector r;
    QVector dd;
    Omega=0.;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);    // dimension
//   iju = (int *) malloc((Dim+2) *sizeof(int));
//   ax = (double *) malloc((Dim+2) *sizeof(double));
//   bx = (double *) malloc((Dim+2) *sizeof(double));
    // temporary x
    V_Constr(&dd, "dd", Dim, Normal,_LPTrue); V_SetAllCmp(&dd,1.);
   Asgn_VV(&dd, Mul_QV(Diag_Q(A),&dd)); 
    // ax=Diag_Q(A)
             //dd= diag(A)
    dof_u=Dim;
    while (fabs(V__GetCmp(&dd,dof_u)) < 1.e-17) dof_u--;
    printf(" pdof %d \n",dof_u);
    //   V_SetAllCmp(&dd,0.);

//   Iter = 0;
    // r = b - A * x(i)
    V_Constr(&r, "r", Dim, Normal, _LPTrue);  Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
    
    bNorm=1.;  Iter = 0;   int MaxIter1=1;
    indx = (int *) malloc(Dim * sizeof(int)); 
    
//     int indd[300];
   
    while (!RTCResult(Iter, l2Norm_V(&r), bNorm, VankaIterId) && Iter < MaxIter1) {
        Iter++;
        /*


        while (Iter < 8) {
          Iter++;*/

        // pressure Vertex loop
        for (i=dof_u+1;i<=Dim;i++) {
            //if(fabs(Q__GetVal(A,i,i))< 1.e-20){
            //vertex system ++++++++++++++++++++++++++++++++++
            count=Q__GetLen(A,i);    // dim local system  //row length
	   for (j1=0;j1<Dim;j1++) indx[j1]=-1;
//           pos1 = (int *) malloc(count * sizeof(int)); 
// 	  loc_idx= (int *) malloc(count * sizeof(int));
            //  local system for pressure p(i)
            //  x(i+1)=D^(-1) [(D-A)x(i)+b]
            p_i_old= V__GetCmp(x,i);
           
            // loop over the term in eq p(i)
	    Len=0;
            for (lj=0;lj<count;lj++) {
	      k1=Q__GetPos(A,i,lj);
	      if(k1<= dof_u){ pos1[Len]=k1; loc_idx[Len]=lj;  indx[k1-1]=Len;  Len++;  }
	    }
            printf(" line %d Len %d \n", i,Len );
         
            Q_Constr(&Mat, "Mat",Len,_LPFalse, Rowws, Normal,_LPTrue);
            V_Constr(&xb, "xb", Len, Normal,_LPTrue);
            V_Constr(&bb, "bb", Len, Normal,_LPTrue);
            V_Constr(&xf, "xf", Len, Normal,_LPTrue);
            V_Constr(&bf, "bf", Len, Normal,_LPTrue);
	    
         for (j2=0;j2<Len;j2++) {
	    val=Q__GetVal(A,i,loc_idx[j2]); 
	      V__SetCmp(&bb,j2+1,val);
          }
            for (j=0;j<Len;j++) {
                il= pos1[j];
// 		   Q_SetEntry(&Mat,j+1,indx[k-1],indx[k-1]+1,1.);
                Lenl=Q__GetLen(A,il);
                sum=V__GetCmp(b,il);
		 Q_SetLen(&Mat,j+1,Len);  for (j1=0;j1<Len;j1++) Q_SetEntry(&Mat,j+1,j1,j1+1,0.);
                for (jl=0;jl<Lenl;jl++) {
                    k=Q__GetPos(A,il,jl);
//                     if(indx[k-1] <0) {sum -= V__GetCmp(x,k)*Q__GetVal(A,il,jl);}
                       if( pos1[j] != k ) {sum -= V__GetCmp(x,k)*Q__GetVal(A,il,jl);}
                    else {
		      Q_SetEntry(&Mat,j+1,indx[k-1],indx[k-1]+1,Q__GetVal(A,il,jl));
		     
		    }
                }
                V__SetCmp(&bf,j+1,sum);  V__SetCmp(&xf,j+1,V__GetCmp(x,il)); V__SetCmp(&xb,j+1,0.);
		printf(" i =  %d bb = %g bbf=%g b=%g \n",j2, V__GetCmp(&bb,j+1),  V__GetCmp(&bf,j+1), il, V__GetCmp(b,il));
            }
             GMRESIter(&Mat,&xf,&bf,10,Dummy,Omega);
            GMRESIter(&Mat,&xb,&bb,10,Dummy,Omega);
            printf("\n Mat \n");
	     ftot=0.;            btot=0.;
            for (lj=0;lj<Len;lj++) {
	      
	        for (k1=0;k1<Len;k1++)  printf(" %g ",Q__GetVal(&Mat,lj+1,k1));
		  printf("  \n");
                ftot +=V__GetCmp(&bb,lj+1)*(V__GetCmp(&xf,lj+1) );
//              btot += V__GetCmp(&bb,lj+1)*V__GetCmp(&xb,lj+1);
                val=Q__GetVal(A,i,loc_idx[lj]);
		val2=V_GetCmp(&bb,lj+1);
		btot += val*val2/Q__GetVal(&Mat,lj+1,lj)    ;//V__GetCmp(&xb,lj+1);  
            }

            // for (lj=0;lj<kk;lj++)
            // pressure   ppp=Omega*p_i_old+(1-Omega)*ftot/btot;
            ppp=-ftot/(btot+1.e-10);
              V__SetCmp(x,i,ppp);
            // velocity
	     
             for (j2=0;j2<Len;j2++){
	        printf(" bb= %g bf= %g xb= %g xf=%g ppp= %g ft=%g btot=%g \n",V__GetCmp(&bb,j2+1),V__GetCmp(&bf,j2+1),
		       V__GetCmp(&xb,j2+1),
		       V__GetCmp(&xf,j2+1),ppp,ftot,btot);
	         val=V__GetCmp(&xf,j2+1)-V__GetCmp(&xb,j2+1)*ppp;
                 V_SetCmp(x,pos1[j2],val);  
	     }

            V_Destr(&xb);
            V_Destr(&xf);
            V_Destr(&bb);
            V_Destr(&bf);
            Q_Destr(&Mat);
//             free(pos1);
        }
        if (Iter < MaxIter)
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
//     }
    }
    free(indx);
    V_Destr(&r);
    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);
// printf("I am VankaNS \n");
// exit(0);
    return (x);
}


/*
// =================================================================
QVector *GaussIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                   PrecondProcType Dummy, _LPDouble Omega) {



  int i,j,k,l,n;
  double max,pivot,coef,s;
  n = Q_GetDim(A);
  // row loop
  for(k=1;k<=n;k++){

  // pivot
    max=fabs(Q__GetVal(A,k,k)); l=k;
    for (i=k+1;i<=n;i++) if(max<fabs(Q__GetVal(A,i,k))) { max=fabs(Q__GetVal(A,i,k));l=i;}
    if (fabs(max) < 1.e-20) printf( " \n Error pivot  =0 \n ");
    // swap line pivot -> k
      if (l!=k){
	pivot=Q__GetCmp(b,k);Q__SetCmp(b,k, Q__GetCmp(b,l)); Q__SetCmp(b,l,pivot);
        for (j=k;j<=n;j++) {
          pivot=Q__GetVal(A,k,j); Q__SetEntry(A,k,j-1,j,Q__GetVal(A,l,j)); Q__SetEntry(A,l,j-1,j,pivot);
        }
      }
        // gauss elimination
      pivot=Q__GetVal(A,k,k);
      for (i=k+1;i<=n;i++) {
        coef=Q__GetVal(A,i,k)/pivot;
	 V_AddCmp(x,i,-coef*Q__GetCmp(b,k));  for (j=k+1;j<n;j++) Q__AddEntry(A,i,j-1,j,-coef*Q__GetVal(A,k,j));
      }
  }
  // back sol
  if (Q__GetVal(A,n,n) < 1.e-20) printf( "\n Error det =0 \n");
     V_SetCmp(x,n,Q__GetCmp(b,n)/Q__GetVal(A,n,n));
    for (i=n-1;i>=1;i--) {
      s=Q__GetCmp(b,i); for (j=i+1;j<=n;j++)  s -=Q__GetVal(A,i,j)*Q__GetCmp(b,j);
     V_SetCmp(x,i,s/Q__GetVal(A,i,i));

  }
  return(0);





}*/
// =================================================================
// =================================================================
int Gauss(double c[], double x[], int n,int ind[]) {
    int i,j,k,p,q,m,itemp;
    double temp,sum,max;
//  for (j=0;j<n-1;j++) ind[j]=j;

    for (j=0;j<n-1;j++) {
        max = fabs(c[j*(n+1)+j]);
        p=j;
        for (m=j+1;m<n;m++) {
            if (fabs(c[m*(n+1)+j]) > max) {
                max = c[m*(n+1)+j];
                p = m;
            }
        }
        if (p != j) {
            itemp=ind[j];
            ind[j]=ind[p];
            ind[p]=itemp;
            for (q=j;q<n+1;q++) {
                temp = c[j*(n+1)+q];
                c[j*(n+1)+q] = c[p*(n+1)+q];
                c[p*(n+1)+q] = temp;
            }
        }
        for (i=j+1;i<n;i++) {
            temp = c[i*(n+1)+j] / c[j*(n+1)+j];
            for (k=j;k<=n;k++) {
                c[i*(n+1)+k] = c[i*(n+1)+k] - (temp * c[j*(n+1)+k]);
            }
        }
    }

    x[n-1] = c[(n-1)*(n+1)+n] / c[(n-1)*(n+1)+(n-1)];

    for (i=n-2;i>=0;i--) {
        sum = 0;
        for (j=i+1;j<n;j++) {
            sum = sum + (c[i*(n+1)+j] * x[j]);
        }
        x[i] = (c[i*(n+1)+n] - sum)/c[i*(n+1)+i];
    }
//     printf("\nSOLUTION OF GIVEN SYSTEM\n");
//     fprintf(fp,"\nSOLUTION OF GIVEN SYSTEM\n");
    for (i=0;i<n;i++)
    {
        printf("\nx%d = %.3f ",i,x[i]);
//         fprintf(fp,"\nx%d = %.3f ",i,x[i]);
    }
//
    printf("\n");
//     fprintf(fp,"\n");
}
// =================================================================
// int Gauss(double t[], double x[], double b[], int n,int indx[]) {
//
// // double tt[10000];
//   int i,j,k,l,krowxn,jrowxn,irowxn;
//   double max,pivot,coef,s,sum1;
//
// //           printf( " \n Gauss  call after %d \n ",n);
// //        for(i=0;i<n;i++){
// //    for(j=0;j<n;j++){
// //      printf( " %g  ", t[i*n+j]);
// //    }
// //    printf( " %g  \n  ",b[i]);
// //  }
//
//
//
//   for (k=0;k<n-1;k++) {
//     // pivot
//     max=fabs(t[indx[k]*n+k]);
//     l=k;
//     for (i=k+1;i<n;i++) {
//       irowxn=indx[i]*n;
//       if (max<fabs(t[irowxn+k])) {
//         max=fabs(t[irowxn+k]);
//         l=i;
//       }
//     }
//     // swap line pivot -> k
//     pivot=indx[l];
//     indx[l]=indx[k];
//     indx[k]=pivot;
//
//     // gauss elimination
//     krowxn=indx[k]*n;
//     pivot=t[krowxn+k];
//     for (i=k+1;i<n;i++) {
//       irowxn=indx[i]*n;
//       coef=t[irowxn+k]/pivot;
//       b[indx[i]] -= coef*b[indx[k]];
//       for (j=k;j<n;j++)  t[irowxn+j] -= coef*t[krowxn+j];
//     }
//   }
//
//   // back sol
//   krowxn=indx[n-1]*n;
//   if (fabs(t[krowxn+ (n-1)]) < 1.e-20) {
//     printf("\n Error det =0 \n");
//   }
//   x[n-1]=b[indx[n-1]]/t[krowxn+ (n-1)];
//   for (i=n-2;i>=0;i--) {
//     s=b[indx[i]];
//     for (j=i+1;j<n;j++)  s -=t[indx[i]*n+j]*x[j];
//     x[i]=s/t[indx[i]*n+i];
//   }
//
//   return (0);
//
//
//
//
//
// };



// =================================================================
QVector *VankaIterF(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                    PrecondProcType Dummy, _LPDouble Omega) {


    int Iter,i,Len,Lenl,j,j1,k1,k,il,jl;
    int *indx,*pos1;
    _LPReal bNorm;
    double sum,aadiag,val;
    size_t Dim;
//   QVector ax;  QVector dd;
    QVector r;
    QVector xx;
    QVector bb;
    QMatrix Mat;
    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);
    Dim = Q_GetDim(A);    // dimension
    // temporary x
// V_Constr(&ax, "ax", Dim, Normal,_LPTrue);  // V_SetAllCmp(&ax,1.);
    // ax=Diag_Q(A)
    // V_Constr(&dd, "dd", Dim, Normal,_LPTrue); // Asgn_VV(&dd, Mul_QV(Diag_Q(A),&ax));
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_SetAllCmp(&r,0.);
    /* r = b - A * x(i) */

    Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
    Iter = 0;
    MaxIter=1;
    bNorm = l2Norm_V(b);
//     double t[20000],xt[300],bt[300];

    indx = (int *) malloc(Dim * sizeof(int));
//     int indd[300];
    for (i=0;i<Dim;i++) indx[i]=-1;

    while (!RTCResult(Iter, l2Norm_V(&r), bNorm, VankaIterId) && Iter < MaxIter) {
        Iter++;




//   while (Iter < MaxIter) {
//     Iter++;

        // Vertex loop
        for (i=1;i<=Dim;i++) {
//         printf(" \n Solving system  %d \n ",i);
//         for (j1=0;j1<Dim;j1++)  indx[j1]=-1;
            //vertex system ++++++++++++++++++++++++++++++++++
            Len=Q__GetLen(A,i);    // dim local system
            pos1 = (int *) malloc(Len * sizeof(int));
            Q_Constr(&Mat, "Mat",Len,_LPFalse, Rowws, Normal,_LPTrue);
            V_Constr(&xx, "xx", Len, Normal,_LPTrue);
            V_Constr(&bb, "bb", Len, Normal,_LPTrue);
            // solution local system x(i+1)=D^(-1) [(D-A)x(i)+b]
            for (j=0;j<Len;j++) {
                il=Q__GetPos(A,i,j)-1;
                pos1[j]=il;
                indx[il]=j;
//           indd[j]=j;
//           xt[j]=V__GetCmp(x,il+1);
                V__SetCmp(&xx,j+1,V__GetCmp(x,il+1));
// for (j1=0;j1<Len;j1++)  t[j*(Len+1)+j1]=0.;
            }

            for (j=0;j<Len;j++) {
                il= pos1[j]+1;
                Lenl=Q__GetLen(A,il);
// 	t[j*Len+j]=1.;
                // sum=(D-A)x(i)+b
                sum=V__GetCmp(b,il);
                Q_SetLen(&Mat,j+1,Len);
                for (jl=0;jl<Lenl;jl++) {
                    k=Q__GetPos(A,il,jl);

                    if ( indx[k-1] <0 ) sum -= V__GetCmp(x,k)*Q__GetVal(A,il,jl);
                    else {
//              t[j*(Len+1)+indx[k-1]]=Q__GetVal(A,il,jl);
                        Q_SetEntry(&Mat,j+1,indx[k-1],indx[k-1]+1,Q__GetVal(A,il,jl));
                        //             t[j*Len+j]=1.;
// 	     printf( " ( %d, %d ) ", j,indx[k-1]);
                    }
                }
//           t[j*(Len+1)+Len]=sum;
                V__SetCmp(&bb,j+1,sum);

// 	 printf( " \n ");
            }
//         for (j=0;j<Len;j++) {
//           indd[j]=j;
//         }
//       t[0*4+0]=1.;   t[1*4+1]=1.;   t[2*4+2]=1.;   t[3*4+3]=1.;
//       bt[0]=1.;bt[1]=1.;bt[2]=1.;bt[3]=1.;
//         Gauss(t,xt,Len,indd);
            JacobiIter(&Mat,&xx,&bb,1,Dummy,Omega);
            for (j=0;j<Len;j++) {
                indx[pos1[j]]=-1;
//           V__SetCmp(x,pos1[j]+1,xt[indd[j]]);
                V__SetCmp(x,pos1[j]+1,V__GetCmp(&xx,j+1));
            }
            // new values
//       for (j=0;j<Len;j++) V__SetCmp(x,Q__GetPos(A,i,j),V__GetCmp(&ax,j+1));
            V_Destr(&xx);
            V_Destr(&bb);
            Q_Destr(&Mat);
            free(pos1);
        }

        if (Iter < MaxIter)   Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
//  double bnorm2=l2Norm_V(&r);
//            printf(" bNorm %g   bnorm2 %g \n ",bNorm ,bnorm2 );

    }
//   V_Destr(&ax);  V_Destr(&dd);
    V_Destr(&r);
    free(indx);
    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}



// =================================================================
QVector *VankaIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                   PrecondProcType Dummy, _LPDouble Omega) {


    int Iter,i,Len,Lenl,j,k,il,jl;
    _LPReal bNorm;
    double sum,aadiag,val;
    size_t Dim;
    QVector ax;
    QVector dd;
    QVector r;
    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);
    Dim = Q_GetDim(A);    // dimension
    // temporary x
    V_Constr(&ax, "ax", Dim, Normal,_LPTrue);
    V_SetAllCmp(&ax,1.);
    // ax=Diag_Q(A)
    V_Constr(&dd, "dd", Dim, Normal,_LPTrue);
    Asgn_VV(&dd, Mul_QV(Diag_Q(A),&ax));
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_SetAllCmp(&r,0.);
    /* r = b - A * x(i) */

    Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
    Iter = 0;
    MaxIter=8;
    bNorm = l2Norm_V(b);
    while (!RTCResult(Iter, l2Norm_V(&r), bNorm, VankaIterId) && Iter < MaxIter) {
        Iter++;





//   while (Iter < MaxIter) {
//     Iter++;

        // Vertex loop
        for (i=1;i<=Dim;i++) {
            //vertex system ++++++++++++++++++++++++++++++++++
            Len=Q__GetLen(A,i);    // dim local system
            // solution local system x(i+1)=D^(-1) [(D-A)x(i)+b]
            for (j=0;j<Len;j++) {
                il=Q__GetPos(A,i,j);
                Lenl=Q__GetLen(A,il);
                // sum=(D-A)x(i)+b
                sum=0.;
                for (jl=0;jl<Lenl;jl++) {
                    k=Q__GetPos(A,il,jl);
                    sum -= V__GetCmp(x,k) *Q__GetVal(A,il,jl);
                    aadiag=V__GetCmp(&dd,il);
                }
                sum += (V__GetCmp(b,il) + V__GetCmp(x,il) *aadiag);
                V__SetCmp(&ax,j+1,sum/aadiag);
            }
            // new values
            for (j=0;j<Len;j++) V__SetCmp(x,Q__GetPos(A,i,j),V__GetCmp(&ax,j+1));
        }

        if (Iter < MaxIter)   Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
//  double bnorm2=l2Norm_V(&r);
//            printf(" bNorm %g   bnorm2 %g \n ",bNorm ,bnorm2 );

    }
    V_Destr(&ax);
    V_Destr(&dd);
    V_Destr(&r);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}


QVector *JacobiIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                    PrecondProcType Dummy, _LPDouble Omega) {
    int Iter;
    _LPReal bNorm;
    size_t Dim;
    QVector r;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        while (!RTCResult(Iter, l2Norm_V(&r), bNorm, JacobiIterId)
                && Iter < MaxIter) {
            Iter++;
            /* x(i+1) = x(i) + Omega * D^(-1) * r */
            AddAsgn_VV(x, Mul_SV(Omega, MulInv_QV(Diag_Q(A), &r)));
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            /* r = b - A * x(i) */
            if (Iter < MaxIter)
                Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        }
    }

    printf("\n Iter %d \n",Iter);
    V_Destr(&r);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *SORForwIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                     PrecondProcType Dummy, _LPDouble Omega) {
    int Iter;
    _LPReal bNorm;
    size_t Dim;
    QVector r;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))    OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        while (!RTCResult(Iter, l2Norm_V(&r), bNorm, SORForwIterId)
                && Iter < MaxIter) {
//       std::cout << bNorm << " = bNorm \n ";
            Iter++;
            /* x(i+1) = x(i) + (D / Omega + L)^(-1) r */
            AddAsgn_VV(x, MulInv_QV(Add_QQ(Mul_SQ(1.0 / Omega, Diag_Q(A)), Lower_Q(A)), &r));
            if (Q_KerDefined(A))        OrthoRightKer_VQ(x, A);
            /* r = b - A * x(i) */
            if (Iter < MaxIter)    Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        }
    }

    V_Destr(&r);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *SORBackwIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                      PrecondProcType Dummy, _LPDouble Omega) {
    int Iter;
    _LPReal bNorm;
    size_t Dim;
    QVector r;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        while (!RTCResult(Iter, l2Norm_V(&r), bNorm, SORBackwIterId)
                && Iter < MaxIter) {
            Iter++;
            /* x(i+1) = x(i) + (D / Omega + U)^(-1) r */
            AddAsgn_VV(x, MulInv_QV(Add_QQ(Mul_SQ(1.0 / Omega, Diag_Q(A)), Upper_Q(A)), &r));
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            /* r = b - A * x(i) */
            if (Iter < MaxIter)
                Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        }
    }

    V_Destr(&r);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *SSORIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                  PrecondProcType Dummy, _LPDouble Omega) {
    int Iter;
    _LPReal bNorm;
    size_t Dim;
    QVector r;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        while (!RTCResult(Iter, l2Norm_V(&r), bNorm, SSORIterId)
                && Iter < MaxIter) {
            Iter++;
            /* x(i+1) = x(i) + (2 - Omega) * (Diag(A) / Omega + Upper(A))^(-1)
                   * Diag(A) / Omega * (Diag(A) / Omega + Lower(A))^(-1) r */
            AddAsgn_VV(x, Mul_SV((2.0 - Omega) / Omega,
                                 MulInv_QV(Add_QQ(Mul_SQ(1.0 / Omega, Diag_Q(A)), Upper_Q(A)),
                                           Mul_QV(Diag_Q(A),
                                                  MulInv_QV(Add_QQ(Mul_SQ(1.0 / Omega, Diag_Q(A)), Lower_Q(A)), &r)))));
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            /* r = b - A * x(i) */
            if (Iter < MaxIter)
                Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        }
    }

    V_Destr(&r);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}


#ifndef _LP_USE_COMPLEX_NUMBERS

QVector *ChebyshevIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                       PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     *  The original algorithm doesn't seem to work correctly.
     *  The three term recursion was therefore in the curent version
     *  adopted after
     *
     *  W. Hackbusch:
     *  Iterative Solution of Large Sparse Systems of Equations,
     *  Springer-Verlag, Berlin, 1994
     *
     */

    int Iter;
    _LPDouble MaxEigenval, MinEigenval;
    _LPDouble c, d, Alpha, Beta;
    _LPDouble T = 0.0, TOld = 0.0, TNew; /* values of Chebyshev polynomials */
    _LPReal bNorm;
    size_t Dim;
    QVector r, p, z;


    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Constr(&z, "z", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        MaxEigenval = GetMaxEigenval(A, PrecondProc, OmegaPrecond);
        MinEigenval = GetMinEigenval(A, PrecondProc, OmegaPrecond);

        c = (MaxEigenval - MinEigenval) / 2.0;
        d = (MaxEigenval + MinEigenval) / 2.0;

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned Chebyshev method */
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, ChebyshevIterId)
                    && Iter < MaxIter) {
                Iter++;
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &z, &r, OmegaPrecond);
                else
                    Asgn_VV(&z, &r);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&z, A);
                if (Iter == 1) {
                    TOld = 1.0;
                    T = d / c;
                    Alpha = 1.0 / d;
                    Asgn_VV(&p, Mul_SV(Alpha, &z));
                } else {
                    TNew = 2.0 * d / c * T - TOld;
                    TOld = T;
                    T = TNew;
                    Alpha = 2.0 / c * TOld / T;
                    Beta = 2.0 * d / c * TOld / T - 1.0;
                    Asgn_VV(&p, Add_VV(Mul_SV(Alpha, &z), Mul_SV(Beta, &p)));
                }
                AddAsgn_VV(x, &p);
                if (Iter < MaxIter)
                    Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
            }
        } else {
            /* plain Chebyshev method (z = r) */
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, ChebyshevIterId)
                    && Iter < MaxIter) {
                Iter++;
                if (Iter == 1) {
                    TOld = 1.0;
                    T = d / c;
                    Alpha = 1.0 / d;
                    Asgn_VV(&p, Mul_SV(Alpha, &r));
                } else {
                    TNew = 2.0 * d / c * T - TOld;
                    TOld = T;
                    T = TNew;
                    Alpha = 2.0 / c * TOld / T;
                    Beta = 2.0 * d / c * TOld / T - 1.0;
                    Asgn_VV(&p, Add_VV(Mul_SV(Alpha, &r), Mul_SV(Beta, &p)));
                }
                AddAsgn_VV(x, &p);
                if (Iter < MaxIter)
                    Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&p);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Destr(&z);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

#else

QVector *ChebyshevIter(QMatrix * /* A */, QVector * /* x */, QVector * /* b */, int /* MaxIter */,
                       PrecondProcType /* PrecondProc */, _LPDouble /* OmegaPrecond */) {
    printf("ERROR: Chebyshev iteration not implemented for complex arithmetic.\n");
    abort();
}

#endif /* _LP_USE_COMPLEX_NUMBERS */



QVector *CGIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */

    int Iter;
    _LPDouble Alpha, Beta, Rho, RhoOld = 0.0;
    _LPReal bNorm;
    size_t Dim;
    QVector r, p, q, z;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    V_Constr(&q, "q", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Constr(&z, "z", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned CG */
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, CGIterId)
                    && Iter < MaxIter) {
                Iter++;
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &z, &r, OmegaPrecond);
                else
                    Asgn_VV(&z, &r);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&z, A);
                Rho = Mul_VV(&r, &z);
                if (Iter == 1) {
                    Asgn_VV(&p, &z);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&p, Add_VV(&z, Mul_SV(Beta, &p)));
                }
                Asgn_VV(&q, Mul_QV(A, &p));
                Alpha = Rho / Mul_VV(&p, &q);
                AddAsgn_VV(x, Mul_SV(Alpha, &p));
                SubAsgn_VV(&r, Mul_SV(Alpha, &q));
                RhoOld = Rho;
            }
        } else {
            /* plain CG (z = r) */
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, CGIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = pow(l2Norm_V(&r), 2.0);
                if (Iter == 1) {
                    Asgn_VV(&p, &r);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&p, Add_VV(&r, Mul_SV(Beta, &p)));
                }
                Asgn_VV(&q, Mul_QV(A, &p));
                Alpha = Rho / Mul_VV(&p, &q);
                AddAsgn_VV(x, Mul_SV(Alpha, &p));
                SubAsgn_VV(&r, Mul_SV(Alpha, &q));
                RhoOld = Rho;
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&p);
    V_Destr(&q);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Destr(&z);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *CGNIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                 PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    int Iter;
    _LPDouble Alpha, Beta, Rho, RhoOld = 0.0;
    _LPReal bNorm;
    size_t Dim;
    QVector r, p, q, s, z;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    V_Constr(&q, "q", Dim, Normal, _LPTrue);
    V_Constr(&z, "z", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Constr(&s, "s", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned CGN */
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, CGNIterId)
                    && Iter < MaxIter) {
                Iter++;
                if (PrecondProc != NULL) {
                    (*PrecondProc)(A, &z, &r, OmegaPrecond);
                    (*PrecondProc)(Transp_Q(A), &q, &z, OmegaPrecond);
                    Asgn_VV(&z, Mul_QV(Transp_Q(A), &q));
                } else {
                    Asgn_VV(&z, &r);
                }
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&z, A);
                Rho = pow(l2Norm_V(&z), 2.0);
                if (Iter == 1) {
                    Asgn_VV(&p, &z);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&p, Add_VV(&z, Mul_SV(Beta, &p)));
                }
                Asgn_VV(&q, Mul_QV(A, &p));
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &s, &q, OmegaPrecond);
                else
                    Asgn_VV(&s, &q);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&s, A);
                Alpha = Rho / pow(l2Norm_V(&s), 2.0);
                AddAsgn_VV(x, Mul_SV(Alpha, &p));
                SubAsgn_VV(&r, Mul_SV(Alpha, &q));
                RhoOld = Rho;
            }
        } else {
            /* plain CGN */
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, CGNIterId)
                    && Iter < MaxIter) {
                Iter++;
                Asgn_VV(&z, Mul_QV(Transp_Q(A), &r));
                Rho = pow(l2Norm_V(&z), 2.0);
                if (Iter == 1) {
                    Asgn_VV(&p, &z);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&p, Add_VV(&z, Mul_SV(Beta, &p)));
                }
                Asgn_VV(&q, Mul_QV(A, &p));
                Alpha = Rho / pow(l2Norm_V(&q), 2.0);
                AddAsgn_VV(x, Mul_SV(Alpha, &p));
                SubAsgn_VV(&r, Mul_SV(Alpha, &q));
                RhoOld = Rho;
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&p);
    V_Destr(&q);
    V_Destr(&z);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Destr(&s);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *GMRESIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                   PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */

    char vName[10];
    int Iter, i, j, k;
    _LPDouble h1, h2, r;
    _LPReal bNorm;
    _LPDouble **h, *y, *s, *c1, *c2;
    size_t Dim;
    _LPBoolean AllocOK;
    QVector *v;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    /* allocation of matrix H and vectors y, s, c1 and c2 */
    AllocOK = _LPTrue;
    h = (_LPDouble **) malloc((GMRESSteps + 1) * sizeof(_LPDouble *));
    if (h == NULL) {
        AllocOK = _LPFalse;
    } else {
        for (i = 1; i <= GMRESSteps; i++) {
            h[i] = (_LPDouble *) malloc((GMRESSteps + 2) * sizeof(_LPDouble));
            if (h[i] == NULL)
                AllocOK = _LPFalse;
        }
    }
    y = (_LPDouble *) malloc((GMRESSteps + 1) * sizeof(_LPDouble));
    if (y == NULL)
        AllocOK = _LPFalse;
    s = (_LPDouble *) malloc((GMRESSteps + 2) * sizeof(_LPDouble));
    if (s == NULL)
        AllocOK = _LPFalse;
    c1 = (_LPDouble *) malloc((GMRESSteps + 1) * sizeof(_LPDouble));
    if (c1 == NULL)
        AllocOK = _LPFalse;
    c2 = (_LPDouble *) malloc((GMRESSteps + 1) * sizeof(_LPDouble));
    if (c2 == NULL)
        AllocOK = _LPFalse;

    /* ... and vectors u */
    Dim = Q_GetDim(A);
    v = (QVector *) malloc((GMRESSteps + 2) * sizeof(QVector));
    if (v == NULL)
        AllocOK = _LPFalse;
    else
        for (i = 1; i <= GMRESSteps + 1; i++) {
            sprintf(vName, "v[%d]", i % 1000);
            V_Constr(&v[i], vName, Dim, Normal, _LPTrue);
        }

    if (!AllocOK)
        LASError(LASMemAllocErr, "GMRESIter", NULL, NULL, NULL);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        /* loop for 'MaxIter' GMRES cycles */
        Iter = 0;
        /* v[1] = r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&v[1], Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&v[1], b);
        }
        while (!RTCResult(Iter, l2Norm_V(&v[1]), bNorm, GMRESIterId)
                && Iter < MaxIter) {
            if (PrecondProc != NULL)
                (*PrecondProc)(A, &v[1], &v[1], OmegaPrecond);
            s[1] = l2Norm_V(&v[1]);
            MulAsgn_VS(&v[1], 1.0 / s[1]);

            /* GMRES iteration */
            i = 0;
            while ((PrecondProc != NULL ? _LPTrue : !RTCResult(Iter, _LPfabs(s[i+1]),
                    bNorm, GMRESIterId)) && i < GMRESSteps && Iter < MaxIter) {
                i++;
                Iter++;
                /* w = v[i+1] */
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &v[i+1], Mul_QV(A, &v[i]), OmegaPrecond);
                else
                    Asgn_VV(&v[i+1], Mul_QV(A, &v[i]));

                /* modified Gram-Schmidt orthogonalization */
                for (k = 1; k <= i; k++) {
                    h[i][k] = InnerProd_VV(&v[i+1], &v[k]);
                    SubAsgn_VV(&v[i+1], Mul_SV(h[i][k], &v[k]));
                }

                h[i][i+1] = l2Norm_V(&v[i+1]);
                MulAsgn_VS(&v[i+1], 1.0 / h[i][i+1]);

                /* Q-R algorithm */
                for (k = 1; k < i; k++) {
                    h1 = c1[k] * h[i][k] + c2[k] * h[i][k+1];
                    h2 = - c2[k] * h[i][k] + c1[k] * h[i][k+1];
                    h[i][k] = h1;
                    h[i][k+1] = h2;
                }
                /* r could be complex-valued!! */

                r = sqrt(h[i][i] * h[i][i] + h[i][i+1] * h[i][i+1]);
                c1[i] = h[i][i] / r;
                c2[i] = h[i][i+1] / r;
                h[i][i] = r;
                h[i][i+1] = 0.0;
                s[i+1] = - c2[i] * s[i];
                s[i] = c1[i] * s[i];
            }

            /* Solving of the system of equations : H y = s */
            for (j = i; j > 0; j--) {
                y[j] = s[j] / h[j][j];
                for (k = j - 1; k > 0; k--)
                    s[k] -= h[j][k] * y[j];
            }

            /* updating solution */
            for (j = i; j > 0; j--)
                AddAsgn_VV(x, Mul_SV(y[j], &v[j]));
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);

            /* computing new residual */
            Asgn_VV(&v[1], Sub_VV(b, Mul_QV(A, x)));
        }
    }

    /* release of vectors u, matrix H and vectors y, s, c1 and c2 */
    if (v != NULL) {
        for (i = 1; i <= GMRESSteps + 1; i++)
            V_Destr(&v[i]);
        free(v);
    }

    if (h != NULL) {
        for (i = 1; i <= GMRESSteps; i++)
            if (h[i] != NULL)
                free(h[i]);
        free(h);
    }
    if (y != NULL)
        free(y);
    if (s != NULL)
        free(s);
    if (c1 != NULL)
        free(c1);
    if (c2 != NULL)
        free(c2);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *BiCGIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                  PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */

    int Iter;
    _LPDouble Alpha, Beta, Rho, RhoOld = 0.0;
    _LPReal bNorm;
    size_t Dim;
    QVector r, r_, p, p_, q, z, z_;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&r_, "r_", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    V_Constr(&p_, "p_", Dim, Normal, _LPTrue);
    V_Constr(&q, "q", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A)) {
        V_Constr(&z, "z", Dim, Normal, _LPTrue);
        V_Constr(&z_, "z_", Dim, Normal, _LPTrue);
    }

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned BiCG */
            Asgn_VV(&r_, &r);
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, BiCGIterId)
                    && Iter < MaxIter) {
                Iter++;
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &z, &r, OmegaPrecond);
                else
                    Asgn_VV(&z, &r);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&z, A);
                if (PrecondProc != NULL)
                    (*PrecondProc)(Transp_Q(A), &z_, &r_, OmegaPrecond);
                else
                    Asgn_VV(&z_, &r_);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&z_, Transp_Q(A));
                Rho = Mul_VV(&z, &r_);
                if (_LPIsZeroNumber(Rho)) {
                    LASError(LASBreakdownErr, "BiCGIter", "Rho", NULL, NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&p, &z);
                    Asgn_VV(&p_, &z_);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&p, Add_VV(&z, Mul_SV(Beta, &p)));
                    Asgn_VV(&p_, Add_VV(&z_, Mul_SV(Beta, &p_)));
                }
                Asgn_VV(&q, Mul_QV(A, &p));
                Alpha = Rho / Mul_VV(&p_, &q);
                AddAsgn_VV(x, Mul_SV(Alpha, &p));
                SubAsgn_VV(&r, Mul_SV(Alpha, &q));
                SubAsgn_VV(&r_, Mul_SV(Alpha, Mul_QV(Transp_Q(A), &p_)));
                RhoOld = Rho;
            }
        } else {
            /* plain BiCG (z = r, z_ = r_) */
            Asgn_VV(&r_, &r);
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, BiCGIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = Mul_VV(&r, &r_);
                if (_LPIsZeroNumber(Rho)) {
                    LASError(LASBreakdownErr, "BiCGIter", "Rho", NULL, NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&p, &r);
                    Asgn_VV(&p_, &r_);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&p, Add_VV(&r, Mul_SV(Beta, &p)));
                    Asgn_VV(&p_, Add_VV(&r_, Mul_SV(Beta, &p_)));
                }
                Asgn_VV(&q, Mul_QV(A, &p));
                Alpha = Rho / Mul_VV(&p_, &q);
                AddAsgn_VV(x, Mul_SV(Alpha, &p));
                SubAsgn_VV(&r, Mul_SV(Alpha, &q));
                SubAsgn_VV(&r_, Mul_SV(Alpha, Mul_QV(Transp_Q(A), &p_)));
                RhoOld = Rho;
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&r_);
    V_Destr(&p);
    V_Destr(&p_);
    V_Destr(&q);
    if (PrecondProc != NULL || Q_KerDefined(A)) {
        V_Destr(&z);
        V_Destr(&z_);
    }

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *QMRIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                 PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */

    int Iter;
    _LPDouble Rho, RhoNew, Xi, Gamma, GammaOld, Eta, Delta, Beta, Epsilon = 0.0,
            Theta, ThetaOld = 0.0;
    _LPReal bNorm;
    size_t Dim;
    QVector r, p, p_, q, y, y_, v, v_, w, w_, z, z_, s, d;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    V_Constr(&p_, "p_", Dim, Normal, _LPTrue);
    V_Constr(&q, "q", Dim, Normal, _LPTrue);
    V_Constr(&v, "v", Dim, Normal, _LPTrue);
    V_Constr(&v_, "v_", Dim, Normal, _LPTrue);
    V_Constr(&w, "w", Dim, Normal, _LPTrue);
    V_Constr(&w_, "w_", Dim, Normal, _LPTrue);
    V_Constr(&s, "s", Dim, Normal, _LPTrue);
    V_Constr(&d, "d", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A)) {
        V_Constr(&y, "y", Dim, Normal, _LPTrue);
        V_Constr(&y_, "y_", Dim, Normal, _LPTrue);
        V_Constr(&z, "z", Dim, Normal, _LPTrue);
        V_Constr(&z_, "z_", Dim, Normal, _LPTrue);
    }

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned QMR (M1 ~ A, M2 = I) */
            Asgn_VV(&v_, &r);
            if (PrecondProc != NULL)
                (*PrecondProc)(A, &y, &v_, OmegaPrecond);
            else
                Asgn_VV(&y, &v_);
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(&y, A);
            RhoNew = l2Norm_V(&y);
            Asgn_VV(&w_, &r);
            Asgn_VV(&z, &w_);
            Xi = l2Norm_V(&z);
            Gamma = 1.0;
            Eta = - 1.0;
            GammaOld = Gamma;
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, QMRIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = RhoNew;
                if (_LPIsZeroNumber(Rho) || _LPIsZeroNumber(Xi)) {
                    LASError(LASBreakdownErr, "QMRIter", "Rho", "Xi", NULL);
                    break;
                }
                Asgn_VV(&v, Mul_SV(1.0 / Rho, &v_));
                MulAsgn_VS(&y, 1.0 / Rho);
                Asgn_VV(&w, Mul_SV(1.0 / Xi, &w_));
                MulAsgn_VS(&z, 1.0 / Xi);
                Delta = Mul_VV(&z, &y);
                if (_LPIsZeroNumber(Delta)) {
                    LASError(LASBreakdownErr, "QMRIter", "Delta", NULL, NULL);
                    break;
                }
                Asgn_VV(&y_, &y);
                if (PrecondProc != NULL)
                    (*PrecondProc)(Transp_Q(A), &z_, &z, OmegaPrecond);
                else
                    Asgn_VV(&z_, &z);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&z_, Transp_Q(A));
                if (Iter == 1) {
                    Asgn_VV(&p, &y_);
                    Asgn_VV(&q, &z_);
                } else {
                    Asgn_VV(&p, Sub_VV(&y_, Mul_SV(Xi * Delta / Epsilon, &p)));
                    Asgn_VV(&q, Sub_VV(&z_, Mul_SV(Rho * Delta / Epsilon, &q)));
                }
                Asgn_VV(&p_, Mul_QV(A, &p));
                Epsilon = Mul_VV(&q, &p_);
                if (_LPIsZeroNumber(Epsilon)) {
                    LASError(LASBreakdownErr, "QMRIter", "Epsilon", NULL, NULL);
                    break;
                }
                Beta = Epsilon / Delta;
                if (_LPIsZeroNumber(Beta)) {
                    LASError(LASBreakdownErr, "QMRIter", "Beta", NULL, NULL);
                    break;
                }
                Asgn_VV(&v_, Sub_VV(&p_, Mul_SV(Beta, &v)));
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &y, &v_, OmegaPrecond);
                else
                    Asgn_VV(&y, &v_);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&y, A);
                RhoNew = l2Norm_V(&y);
                Asgn_VV(&w_, Sub_VV(Mul_QV(Transp_Q(A), &q), Mul_SV(Beta, &w)));
                Asgn_VV(&z, &w_);
                Xi = l2Norm_V(&z);
                Theta = RhoNew / (GammaOld * _LPfabs(Beta));
                Gamma = 1.0 / sqrt(1.0 + Theta * Theta);
                if (_LPIsZeroNumber(Gamma)) {
                    LASError(LASBreakdownErr, "QMRIter", "Gamma", NULL, NULL);
                    break;
                }
                Eta = - Eta * Rho * Gamma * Gamma / (Beta * GammaOld * GammaOld);
                if (Iter == 1) {
                    Asgn_VV(&d, Mul_SV(Eta, &p));
                    Asgn_VV(&s, Mul_SV(Eta, &p_));
                } else {
                    Asgn_VV(&d, Add_VV(Mul_SV(Eta, &p),
                                       Mul_SV(ThetaOld * ThetaOld * Gamma * Gamma, &d)));
                    Asgn_VV(&s, Add_VV(Mul_SV(Eta, &p_),
                                       Mul_SV(ThetaOld * ThetaOld * Gamma * Gamma, &s)));
                }
                AddAsgn_VV(x, &d);
                SubAsgn_VV(&r, &s);
                GammaOld = Gamma;
                ThetaOld = Theta;
            }
        } else {
            /* plain QMR (M1 ~ A, M2 = I, y = y_ = v_, z = z_ = w_) */
            Asgn_VV(&v_, &r);
            RhoNew = l2Norm_V(&v_);
            Asgn_VV(&w_, &r);
            Xi = l2Norm_V(&w_);
            Gamma = 1.0;
            Eta = - 1.0;
            GammaOld = Gamma;
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, QMRIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = RhoNew;
                if (_LPIsZeroNumber(Rho) || _LPIsZeroNumber(Xi)) {
                    LASError(LASBreakdownErr, "QMRIter", "Rho", "Xi", NULL);
                    break;
                }
                Asgn_VV(&v, Mul_SV(1.0 / Rho, &v_));
                MulAsgn_VS(&v_, 1.0 / Rho);
                Asgn_VV(&w, Mul_SV(1.0 / Xi, &w_));
                MulAsgn_VS(&w_, 1.0 / Xi);
                Delta = Mul_VV(&w_, &v_);
                if (_LPIsZeroNumber(Delta)) {
                    LASError(LASBreakdownErr, "QMRIter", "Rho", "Xi", NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&p, &v_);
                    Asgn_VV(&q, &w_);
                } else {
                    Asgn_VV(&p, Sub_VV(&v_, Mul_SV(Xi * Delta / Epsilon, &p)));
                    Asgn_VV(&q, Sub_VV(&w_, Mul_SV(Rho * Delta / Epsilon, &q)));
                }
                Asgn_VV(&p_, Mul_QV(A, &p));
                Epsilon = Mul_VV(&q, &p_);
                if (_LPIsZeroNumber(Epsilon)) {
                    LASError(LASBreakdownErr, "QMRIter", "Epsilon", NULL, NULL);
                    break;
                }
                Beta = Epsilon / Delta;
                if (_LPIsZeroNumber(Beta)) {
                    LASError(LASBreakdownErr, "QMRIter", "Beta", NULL, NULL);
                    break;
                }
                Asgn_VV(&v_, Sub_VV(&p_, Mul_SV(Beta, &v)));
                RhoNew = l2Norm_V(&v);
                Asgn_VV(&w_, Sub_VV(Mul_QV(Transp_Q(A), &q), Mul_SV(Beta, &w)));
                Xi = l2Norm_V(&w_);
                Theta = RhoNew / (GammaOld * _LPfabs(Beta));
                Gamma = 1.0 / sqrt(1.0 + Theta * Theta);
                if (_LPIsZeroNumber(Gamma)) {
                    LASError(LASBreakdownErr, "QMRIter", "Gamma", NULL, NULL);
                    break;
                }
                Eta = - Eta * Rho * Gamma * Gamma / (Beta * GammaOld * GammaOld);
                if (Iter == 1) {
                    Asgn_VV(&d, Mul_SV(Eta, &p));
                    Asgn_VV(&s, Mul_SV(Eta, &p_));
                } else {
                    Asgn_VV(&d, Add_VV(Mul_SV(Eta, &p),
                                       Mul_SV(ThetaOld * ThetaOld * Gamma * Gamma, &d)));
                    Asgn_VV(&s, Add_VV(Mul_SV(Eta, &p_),
                                       Mul_SV(ThetaOld * ThetaOld * Gamma * Gamma, &s)));
                }
                AddAsgn_VV(x, &d);
                SubAsgn_VV(&r, &s);
                GammaOld = Gamma;
                ThetaOld = Theta;
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&p);
    V_Destr(&p_);
    V_Destr(&q);
    V_Destr(&v);
    V_Destr(&v_);
    V_Destr(&w);
    V_Destr(&w_);
    V_Destr(&s);
    V_Destr(&d);
    if (PrecondProc != NULL || Q_KerDefined(A)) {
        V_Destr(&y);
        V_Destr(&y_);
        V_Destr(&z);
        V_Destr(&z_);
    }

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *CGSIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                 PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */

    int Iter;
    _LPDouble Alpha, Beta, Rho, RhoOld = 0.0;
    _LPReal bNorm;
    size_t Dim;
    QVector r, r_, p, p_, q, u, u_, v_;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&r_, "r_", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    V_Constr(&q, "q", Dim, Normal, _LPTrue);
    V_Constr(&u, "u", Dim, Normal, _LPTrue);
    V_Constr(&u_, "u_", Dim, Normal, _LPTrue);
    V_Constr(&v_, "v_", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Constr(&p_, "p_", Dim, Normal, _LPTrue);

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned CGS */
            Asgn_VV(&r_, &r);
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, CGSIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = Mul_VV(&r_, &r);
                if (_LPIsZeroNumber(Rho)) {
                    LASError(LASBreakdownErr, "CGSIter", "Rho", NULL, NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&u, &r);
                    Asgn_VV(&p, &u);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&u, Add_VV(&r, Mul_SV(Beta, &q)));
                    Asgn_VV(&p, Add_VV(&u, Mul_SV(Beta, Add_VV(&q, Mul_SV(Beta, &p)))));
                }
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &p_, &p, OmegaPrecond);
                else
                    Asgn_VV(&p_, &p);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&p_, A);
                Asgn_VV(&v_, Mul_QV(A, &p_));
                Alpha = Rho / Mul_VV(&r_, &v_);
                Asgn_VV(&q, Sub_VV(&u, Mul_SV(Alpha, &v_)));
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &u_, Add_VV(&u, &q), OmegaPrecond);
                else
                    Asgn_VV(&u_, Add_VV(&u, &q));
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&u_, A);
                AddAsgn_VV(x, Mul_SV(Alpha, &u_));
                SubAsgn_VV(&r, Mul_SV(Alpha, Mul_QV(A, &u_)));
                RhoOld = Rho;
            }
        } else {
            /* plain CGS (p_ = p) */
            Asgn_VV(&r_, &r);
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, CGSIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = Mul_VV(&r_, &r);
                if (_LPIsZeroNumber(Rho)) {
                    LASError(LASBreakdownErr, "CGSIter", "Rho", NULL, NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&u, &r);
                    Asgn_VV(&p, &u);
                } else {
                    Beta = Rho / RhoOld;
                    Asgn_VV(&u, Add_VV(&r, Mul_SV(Beta, &q)));
                    Asgn_VV(&p, Add_VV(&u, Mul_SV(Beta, Add_VV(&q, Mul_SV(Beta, &p)))));
                }
                Asgn_VV(&v_, Mul_QV(A, &p));
                Alpha = Rho / Mul_VV(&r_, &v_);
                Asgn_VV(&q, Sub_VV(&u, Mul_SV(Alpha, &v_)));
                Asgn_VV(&u_, Add_VV(&u, &q));
                AddAsgn_VV(x, Mul_SV(Alpha, &u_));
                SubAsgn_VV(&r, Mul_SV(Alpha, Mul_QV(A, &u_)));
                RhoOld = Rho;
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&r_);
    V_Destr(&p);
    V_Destr(&q);
    V_Destr(&u);
    V_Destr(&u_);
    V_Destr(&v_);
    if (PrecondProc != NULL || Q_KerDefined(A))
        V_Destr(&p_);

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

QVector *BiCGSTABIter(QMatrix *A, QVector *x, QVector *b, int MaxIter,
                      PrecondProcType PrecondProc, _LPDouble OmegaPrecond) {
    /*
     *  for details to the algorithm see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */

    int Iter;
    _LPDouble Alpha = 0.0, Beta, Rho, RhoOld = 0.0, Omega = 0.0;
    _LPReal bNorm;
    size_t Dim;
    QVector r, r_, p, p_, v, s, s_, t;

    Q_Lock(A);
    V_Lock(x);
    V_Lock(b);

    Dim = Q_GetDim(A);
    V_Constr(&r, "r", Dim, Normal, _LPTrue);
    V_Constr(&r_, "r_", Dim, Normal, _LPTrue);
    V_Constr(&p, "p", Dim, Normal, _LPTrue);
    V_Constr(&v, "v", Dim, Normal, _LPTrue);
    V_Constr(&s, "s", Dim, Normal, _LPTrue);
    V_Constr(&t, "t", Dim, Normal, _LPTrue);
    if (PrecondProc != NULL || Q_KerDefined(A)) {
        V_Constr(&p_, "p_", Dim, Normal, _LPTrue);
        V_Constr(&s_, "s_", Dim, Normal, _LPTrue);
    }

    if (LASResult() == LASOK) {
        bNorm = l2Norm_V(b);

        Iter = 0;
        /* r = b - A * x(i) */
        if (!_LPIsZeroReal(l1Norm_V(x) / Dim)) {
            if (Q_KerDefined(A))
                OrthoRightKer_VQ(x, A);
            Asgn_VV(&r, Sub_VV(b, Mul_QV(A, x)));
        } else {
            Asgn_VV(&r, b);
        }
        if (PrecondProc != NULL || Q_KerDefined(A)) {
            /* preconditioned BiCGStab */
            Asgn_VV(&r_, &r);
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, BiCGSTABIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = Mul_VV(&r_, &r);
                if (_LPIsZeroNumber(Rho)) {
                    LASError(LASBreakdownErr, "BiCGSTABIter", "Rho", NULL, NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&p, &r);
                } else {
                    Beta = (Rho / RhoOld) * (Alpha / Omega);
                    Asgn_VV(&p, Add_VV(&r, Mul_SV(Beta, Sub_VV(&p, Mul_SV(Omega, &v)))));
                }
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &p_, &p, OmegaPrecond);
                else
                    Asgn_VV(&p_, &p);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&p_, A);
                Asgn_VV(&v, Mul_QV(A, &p_));
                Alpha = Rho / Mul_VV(&r_, &v);
                Asgn_VV(&s, Sub_VV(&r, Mul_SV(Alpha, &v)));
                if (PrecondProc != NULL)
                    (*PrecondProc)(A, &s_, &s, OmegaPrecond);
                else
                    Asgn_VV(&s_, &s);
                if (Q_KerDefined(A))
                    OrthoRightKer_VQ(&s_, A);
                Asgn_VV(&t, Mul_QV(A, &s_));
                Omega = Mul_VV(&t, &s) / pow(l2Norm_V(&t), 2.0);
                AddAsgn_VV(x, Add_VV(Mul_SV(Alpha, &p_), Mul_SV(Omega, &s_)));
                Asgn_VV(&r, Sub_VV(&s, Mul_SV(Omega, &t)));
                RhoOld = Rho;
                if (_LPIsZeroNumber(Omega))
                    break;
            }
        } else {
            /* plain BiCGStab (p_ = p) */
            Asgn_VV(&r_, &r);
            while (!RTCResult(Iter, l2Norm_V(&r), bNorm, BiCGSTABIterId)
                    && Iter < MaxIter) {
                Iter++;
                Rho = Mul_VV(&r_, &r);
                if (_LPIsZeroNumber(Rho)) {
                    LASError(LASBreakdownErr, "BiCGSTABIter", "Rho", NULL, NULL);
                    break;
                }
                if (Iter == 1) {
                    Asgn_VV(&p, &r);
                } else {
                    Beta = (Rho / RhoOld) * (Alpha / Omega);
                    Asgn_VV(&p, Add_VV(&r, Mul_SV(Beta, Sub_VV(&p, Mul_SV(Omega, &v)))));
                }
                Asgn_VV(&v, Mul_QV(A, &p));
                Alpha = Rho / Mul_VV(&r_, &v);
                Asgn_VV(&s, Sub_VV(&r, Mul_SV(Alpha, &v)));
                Asgn_VV(&t, Mul_QV(A, &s));
                Omega = Mul_VV(&t, &s) / pow(l2Norm_V(&t), 2.0);
                AddAsgn_VV(x, Add_VV(Mul_SV(Alpha, &p), Mul_SV(Omega, &s)));
                Asgn_VV(&r, Sub_VV(&s, Mul_SV(Omega, &t)));
                RhoOld = Rho;
                if (_LPIsZeroNumber(Omega))
                    break;
            }
        }
        if (Q_KerDefined(A))
            OrthoRightKer_VQ(x, A);
    }

    V_Destr(&r);
    V_Destr(&r_);
    V_Destr(&p);
    V_Destr(&v);
    V_Destr(&s);
    V_Destr(&t);
    if (PrecondProc != NULL || Q_KerDefined(A)) {
        V_Destr(&p_);
        V_Destr(&s_);
    }

    Q_Unlock(A);
    V_Unlock(x);
    V_Unlock(b);

    return (x);
}

void SetGMRESRestart(int Steps) {
    GMRESSteps = Steps;
}
