// NACA airfoil
// D2Q9
//

#ifndef NACA00XX_A_H
#define NACA00XX_A_H
#endif

#include<math.h>
#include<assert.h>
#include<stdio.h>
#include<stdbool.h>
#define M_PI 3.14159265358979323846

typedef struct AIRFOIL{
// Rigid Body 2D
// 3 DOFs
    double xc; //global coordinates
    double yc; 
    double theta;
    double ux; //velocity
    double uy;
    double wz; //angular velocity
    
    double xx; // NACA00xx
    double c;  // chord
    double xm; // center of mass local
    double ym; // 0 for symmetric
}Airfoil;


double NACA00xx(double x, double c, double xx){
    double yt;
// c: chord length
// x: 0-c
// yt: half thickness
// t: maximum thickness
// 100t=xx
// 0012
    double t;
    t = xx;
    yt=5*t*c*(0.2969*sqrt(x/c)+(-0.1260)*(x/c)+(-0.3516)*pow(x/c,2)+0.2843*pow(x/c,3)+(-0.1015)*pow(x/c,4));

    return yt;
}
double NACA00xxInt(double x, double c, double xx){
    double inty;
// int 2*y(x)*dx
    double t;
    t = xx;
    inty=2*5*t*c*(0.2969/sqrt(c)*1./1.5*pow(x,1.5)+(-0.1260)/c*1./2*pow(x,2)+(-0.3516)/pow(c,2)*1./3*pow(x,3)+0.2843/pow(c,3)*1./4*pow(x,4)+(-0.1015)/pow(c,4)*1./5*pow(x,5));
    
    return inty;
}
double NACA00xxMom(double x, double c, double xx){
    double momx;
// int x*2*y(x)*dx
    double t;
    t = xx;
    momx=2*5*t*c*(0.2969/sqrt(c)*1./2.5*pow(x,2.5)+(-0.1260)/c*1./3*pow(x,3)+(-0.3516)/pow(c,2)*1./4*pow(x,4)+0.2843/pow(c,3)*1./5*pow(x,5)+(-0.1015)/pow(c,4)*1./6*pow(x,6));
    
    return momx;
}
double NACA00xxCent(double c,double xx){
    double xcenter;
// xcenter = int x*2*y(x)*dx / int 2*y(x)*dx
    double intxy;
    double inty;

    intxy = NACA00xxMom(c,c,xx)-NACA00xxMom(0.,c,xx);
    inty = NACA00xxInt(c,c,xx)-NACA00xxInt(0.,c,xx);
    xcenter = intxy/inty;

    return xcenter;
}

double NACA00xxJac(double x, double c,double xx){
    double dydx;
    double t;
    double dx=1e-6;
    t = xx;
    if (x>dx){
    dydx = 5*t*c*(0.2969*1./2/sqrt(c)*pow(x,-1./2)+(-0.1260)/c+(-0.3516)*2/pow(c,2)*x+0.2843*3/pow(c,3)*pow(x,2)+(-0.1015)*4/pow(c,4)*pow(x,3));
    }
    else{
        dydx = (NACA00xx(x+dx,c,xx)-NACA00xx(x,c,xx))/dx;
    }
    return dydx;
}
double OBJ(double t, double x[2], double theta, double c, int n, double xx){
    double obj;
    double v_local[2][9];
    int i;
    int v[2][9] = {{0, 1, 0, -1,  0, 1, -1, -1, 1},
                   {0, 0, 1,  0, -1, 1,  1, -1, -1},};


    for(i=0;i<9;i++){
        v_local[0][i] = v[0][i]*cos(theta) + v[1][i]*sin(theta);
        v_local[1][i] = v[0][i]*-sin(theta) + v[1][i]*cos(theta);
    }
    // NACA00xx: fabs(y) = NACA00xx(x)
    obj = (NACA00xx(x[0]+t*v_local[0][n],c,xx) - fabs(x[1] + t*v_local[1][n]));
    return obj;    
}
double JAC(double t, double x[2], double theta, double c, int n, double xx){
    double jac;
    double v_local[2][9];
    int i;
    int v[2][9] = {{0, 1, 0, -1,  0, 1, -1, -1, 1},
                   {0, 0, 1,  0, -1, 1,  1, -1, -1},};


    for(i=0;i<9;i++){
        v_local[0][i] = v[0][i]*cos(theta) + v[1][i]*sin(theta);
        v_local[1][i] = v[0][i]*-sin(theta) + v[1][i]*cos(theta);
    }
    // obj*obj
    if((x[1]+t*v_local[1][n])>0){
        jac = 2*(NACA00xx(x[0]+t*v_local[0][n],c,xx) - fabs(x[1] + t*v_local[1][n]))*(NACA00xxJac(x[0]+t*v_local[0][n],c,xx)*v_local[0][n] - (v_local[1][n]));
    }else{
        jac = 2*(NACA00xx(x[0]+t*v_local[0][n],c,xx) - fabs(x[1] + t*v_local[1][n]))*(NACA00xxJac(x[0]+t*v_local[0][n],c,xx)*v_local[0][n] + (v_local[1][n]));
    }
    // fabs(obj);
    /*if((NACA00xx(x[0]+t*v_local[0][n],c,xx) - fabs(x[1] + t*v_local[1][n]))>0){
        if((x[1]+t*v_local[1][n])>0){
            jac = (NACA00xxJac(x[0]+t*v_local[0][n],c,xx)*v_local[0][n] - (v_local[1][n]));
        }
        else{
            jac = (NACA00xxJac(x[0]+t*v_local[0][n],c,xx)*v_local[0][n] + (v_local[1][n]));
        }

    }
    else{
        if((x[1]+t*v_local[1][n])>0){
            jac = -(NACA00xxJac(x[0]+t*v_local[0][n],c,xx)*v_local[0][n] - (v_local[1][n]));
        }
        else{
            jac = -(NACA00xxJac(x[0]+t*v_local[0][n],c,xx)*v_local[0][n] + (v_local[1][n]));
        }

    }*/

    return jac;
}
double NR(double x[2], double theta, double c, int n, double xx){
    double t0 = 0.7;
    double t=-1;
    double tol = 1e-12;
    int maxiters = 64;
    int i=0;
    int j;    
    double dt=1e-6;
    double v_local[2][9];
    int v[2][9] = {{0, 1, 0, -1,  0, 1, -1, -1, 1},
                   {0, 0, 1,  0, -1, 1,  1, -1, -1},};

    for(j=0;j<9;j++){
        v_local[0][j] = v[0][j]*cos(theta) + v[1][j]*sin(theta);
        v_local[1][j] = v[0][j]*-sin(theta) + v[1][j]*cos(theta);
    }

    while((fabs(OBJ(t,x,theta,c,n,xx))>tol)&&(OBJ(t,x,theta,c,0,xx)<=0)){
        // obj<=0 to exclude the grid inside airfoil         
        //
        t = t0 - pow(OBJ(t0,x,theta,c,n,xx),2)/JAC(t0,x,theta,c,n,xx);
        t0 = t;
        i++;
        if(i>maxiters){                       
            //printf("iters:%d,t:%f\n",i,t);
            break;
        }
    }
    if((x[0]<=0)&&(x[0]+1*v_local[0][n]>=0)&&(fabs(x[1]+1*v_local[1][n])<=NACA00xx(x[0]+1*v_local[0][n],c,xx))){
        //printf("%f %f\n",x[0],x[1]);

        for(t0=1;(fabs(x[1]+t0*v_local[1][n])<=NACA00xx(x[0]+t0*v_local[0][n],c,xx));t0-=dt){
            if(fabs(OBJ(t0,x,theta,c,n,xx))<1e-3){
                t = t0;
                //printf("%f %f %f\n",x[0],x[1],t);
                break;
            }
        }
    } 
    
    return t;
}

// calculate the wall velocity when node move from solid to fluid
void FlappingWing(double A0, double B0, double freq, double phi, double x0, double y0, double alpha0,double t, Airfoil* pNACA00xx){
    // x(t) = x0 + A0/2*cos(2*pi*f*t)
    // y(t) = y0
    // alpha(t) = alpha0 + B0*sin(2*pi*f*t+phi)
    // r: (x,y,alpha)
    // rdot: (xdot,ydot,alphadot)
    double r[3];
    double rdot[3];

    r[0] = x0 + A0/2*cos(2*M_PI*freq*t);
    r[1] = y0;
    r[2] = alpha0 + B0*sin(2*M_PI*freq*t+phi);

    rdot[0] = A0*M_PI*freq*-sin(2*M_PI*freq*t);
    rdot[1] = 0;
    rdot[2] = B0*2*M_PI*freq*cos(2*M_PI*freq*t+phi);

    pNACA00xx->xc=r[0];
    pNACA00xx->yc=r[1];
    pNACA00xx->wz=r[2];
    pNACA00xx->ux=rdot[0];
    pNACA00xx->uy=rdot[1];
    pNACA00xx->wz=rdot[2];

}

void G2L(int grid[2], double point[2], double xc, double yc, double theta, double xm, double ym){
    // Global grid point converted to local airfoil coordinate
    double x,y;
    int i,j;

    i = grid[0];
    j = grid[1];
 
    x = (i-xc)*cos(theta) + (j-yc)*sin(theta);
    y = (i-xc)*-sin(theta) + (j-yc)*cos(theta);

    point[0] = x+xm;
    point[1] = y+ym;

}

void G2L_F(double t,int grid[2], double point[2], double xm, double ym, double A0, double B0, double freq, double phi, double x0, double y0, double alpha0){
    // Global grid point converted to local airfoil coordinate
    double x,y;
    int i,j;

    double xc;
    double yc;
    double theta;

    i = grid[0];
    j = grid[1];

    xc = x0 + A0/2*cos(2*M_PI*freq*t);
    yc = y0;
    theta = alpha0 + B0*sin(2*M_PI*freq*t+phi);

    x = (i-xc)*cos(theta) + (j-yc)*sin(theta);
    y = (i-xc)*-sin(theta) + (j-yc)*cos(theta);

    point[0] = x+xm;
    point[1] = y+ym;

}


void Wall_F(double t, double t0, double t1, int grid[2], double uw[2], double c, double xx, double xm, double ym, double A0, double B0, double freq, double phi, double x0, double y0, double alpha0){
    double point[2];
    double point0[2];
    double point1[2];
    double pointm[2];

    double tl = t0;
    double tr = t1;
    double tm = 0;

    int iteration=0;
    int maxiter=1e5;
    double h = (t1-t0)/maxiter;
    bool maxiter_flag = 0;
    bool success_flag = 0;

    double xc;
    double yc;
    double theta;
    double ux;
    double uy;
    double wz;
    double r;

    G2L_F(t0,grid,point0,xm,ym,A0,B0,freq,phi,x0,y0,alpha0);
    G2L_F(t1,grid,point1,xm,ym,A0,B0,freq,phi,x0,y0,alpha0);
    // t0: node inside solid
    // t1: node outside solid
    printf("i: %d j: %d\n",grid[0],grid[1]);
    printf("0: %f %f 1: %f %f\n",point0[0],point0[1],point1[0],point1[1]);
    printf("0: %f 1: %f \n",NACA00xx(point0[0],c,xx)-fabs(point0[1]),NACA00xx(point1[0],c,xx)-fabs(point1[1]));
    //printf("OBJ 0: %f 1: %f\n", OBJ(0,point0,0,c,0,xx), OBJ(0,point1,0,c,0,xx)); 
    
    assert(NACA00xx(point0[0],c,xx)-fabs(point0[1])>0);
    if((point0[0]>=0)&&(point1[0]>=0)){
    assert(NACA00xx(point1[0],c,xx)-fabs(point1[1])<0);
    
    // binary search
    while(tl<tr){
        tm = (tl+tr)*0.5; 
        G2L_F(tm,grid,pointm,xm,ym,A0,B0,freq,phi,x0,y0,alpha0);
        if(NACA00xx(pointm[0],c,xx)-fabs(pointm[1])>0){
            tl=tm;
        }else{
            tr=tm;
        }
        iteration ++;
	if(fabs(NACA00xx(pointm[0],c,xx)- pointm[1])<h){
		success_flag = 1;
		break;
	}
        if(iteration>maxiter){
            printf("maximum iteration %d reached\n",maxiter);
	    printf("%f %f %f\n",pointm[0],pointm[1],NACA00xx(pointm[0],c,xx)- pointm[1]);
	    maxiter_flag = 1;
            break;
        }
    }
    }
    if((point0[0]>=0&&point1[0]<0)||(1==maxiter_flag)){ // increasing search
	for(tm=t0;tm<t1;tm+=h*1e-2){
        	G2L_F(tm,grid,pointm,xm,ym,A0,B0,freq,phi,x0,y0,alpha0);
		if((NACA00xx(pointm[0],c,xx)- pointm[1])<h){
			success_flag = 1;
			break;
		}
	}
    }
    			
    //assert(fabs(NACA00xx(pointm[0],c,xx)- pointm[1])<h);
    if(1==success_flag){
    t = tm;
    point[0] = pointm[0];
    point[1] = pointm[1];
    }else{
	printf("%f %f %f\n",pointm[0],pointm[1],NACA00xx(pointm[0],c,xx)- pointm[1]);
	t = (t0+t1)*0.5;	
        G2L_F(t,grid,point,xm,ym,A0,B0,freq,phi,x0,y0,alpha0);
	printf("searching failed, use middle time instead\n");
    }
    
    xc = x0 + A0/2*cos(2*M_PI*freq*t);
    yc = y0;
    theta = alpha0 + B0*sin(2*M_PI*freq*t+phi);

    ux = A0*M_PI*freq*-sin(2*M_PI*freq*t);
    uy = 0;
    wz = B0*2*M_PI*freq*cos(2*M_PI*freq*t+phi);

    r = sqrt(pow(point[0]-xc,2)+pow(point[1]-yc,2));
    uw[0] = ux + wz*r*-sin(theta);
    uw[1] = uy + wz*r*cos(theta);

    
}





        


                

    

    
    



