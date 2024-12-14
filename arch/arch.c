#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include<time.h>
#include"../linear-algebra/matrix.c"
#include"../linear-algebra/lu.c"
#include"../newton-raphson/nr.c"
#include"../elasticity/element.c"
#include"../elasticity/structure.c"
#include"catenary.c"

double E;

double h;	
double l;

int n;
	
double width;
double arch_thick;
double upper_thick;

double small_section;
double big_section;

double rho;

void set_parameters()
{
	E=20e6;
	
	h=60;
	n=8;
	
	l=30;
	width=8;
	arch_thick=l/11;
	upper_thick=l/10;
	
	big_section=M_PI*pow(0.30,2);
	small_section=M_PI*pow(0.20,2);
	
	rho=7500;
	g=0;
	
	N=14;
}

void phi(struct matrix *x,struct matrix *v)
{
	double th=90;
	double ph=90;
	double f=8000*pow(getm(x,0,2)/h,4);
	zeros(v,3,1);
	setm(v,0,0,cos(ph)*sin(th)*f);
	setm(v,1,0,sin(ph)*sin(th)*f);
	setm(v,2,0,cos(th)*f);
}

double tau(double epsilon)
{
	return E*epsilon;
}

double dtau(double epsilon)
{
	return E;
}

void xsquare(struct matrix *A,struct matrix *B,struct matrix *C,struct matrix *D)
{
	inizialize_substructure();
	_append_column(&subX,A);
	_append_column(&subX,B);
	_append_column(&subX,C);
	_append_column(&subX,D);
	
	substr.n=4;
	substr.nelements=6;
	substr.elements=(struct element*)malloc(6*sizeof(struct element));
	
	struct matrix AB,BC,CD,DA,AC,BD;
	difference(A,B,&AB);
	difference(B,C,&BC);
	difference(C,D,&CD);
	difference(D,A,&DA);
	difference(A,C,&AC);
	difference(B,D,&BD);
	
	struct element e,f;
	
	e.A=big_section;
	e.tau=tau;
	e.dtau=dtau;
	
	f.A=small_section;
	f.tau=tau;
	f.dtau=dtau;
	
	e.e1=0;
	e.e2=1;
	e.L=norm(&AB);
	e.m=e.L*e.A*rho;
	copye(&e,substr.elements);
	
	e.e1=1;
	e.e2=2;
	e.L=norm(&BC);
	e.m=e.L*e.A*rho;
	copye(&e,substr.elements+1);
	
	e.e1=2;
	e.e2=3;
	e.L=norm(&CD);
	e.m=e.L*e.A*rho;
	copye(&e,substr.elements+2);
	
	e.e1=3;
	e.e2=0;
	e.L=norm(&DA);
	e.m=e.L*e.A*rho;
	copye(&e,substr.elements+3);
	
	f.e1=0;
	f.e2=2;
	f.L=norm(&AC);
	f.m=f.L*f.A*rho;
	copye(&f,substr.elements+4);
	
	f.e1=1;
	f.e2=3;
	f.L=norm(&BD);
	f.m=f.L*f.A*rho;
	copye(&f,substr.elements+5);
	
	add_substructure();
	clear_substructure();
	
	clearm(&AB);
	clearm(&BC);
	clearm(&CD);
	clearm(&DA);
	clearm(&AC);
	clearm(&BD);
}

void xcube(struct matrix *A,struct matrix *B,struct matrix *C,struct matrix *D,struct matrix *E,struct matrix *F,struct matrix *G,struct matrix *H)
{
	xsquare(A,B,C,D);
	xsquare(A,B,F,E);
	xsquare(B,C,G,F);
	xsquare(A,D,H,E);
	xsquare(C,D,G,H);
	xsquare(E,F,G,H);
}

void arch_circle()
{
	int n=1;
	double theta=0.5*M_PI;
	
	struct matrix A,B,C,D,E,F,G,H;
	
	double width=3;
	double r,R;
	
	r=10;
	R=12;
	
	for(int i=-n;i<n;i++)
	{
		double phi=(double)i/(double)n*theta;
		
		loadm(&A,3,1);
		setm(&A,0,0,R*sin(phi));
		setm(&A,1,0,0);
		setm(&A,2,0,R*cos(phi));
	
		loadm(&B,3,1);
		setm(&B,0,0,r*sin(phi));
		setm(&B,1,0,0);
		setm(&B,2,0,r*cos(phi));
	
		loadm(&C,3,1);
		setm(&C,0,0,r*sin(phi));
		setm(&C,1,0,width);
		setm(&C,2,0,r*cos(phi));
	
		loadm(&D,3,1);
		setm(&D,0,0,R*sin(phi));
		setm(&D,1,0,width);
		setm(&D,2,0,R*cos(phi));
		
		phi=(double)(i+1)/(double)n*theta;
	
		loadm(&E,3,1);
		setm(&E,0,0,R*sin(phi));
		setm(&E,1,0,0);
		setm(&E,2,0,R*cos(phi));
	
		loadm(&F,3,1);
		setm(&F,0,0,r*sin(phi));
		setm(&F,1,0,0);
		setm(&F,2,0,r*cos(phi));
	
		loadm(&G,3,1);
		setm(&G,0,0,r*sin(phi));
		setm(&G,1,0,width);
		setm(&G,2,0,r*cos(phi));
	
		loadm(&H,3,1);
		setm(&H,0,0,R*sin(phi));
		setm(&H,1,0,width);
		setm(&H,2,0,R*cos(phi));
	
		xcube(&A,&B,&C,&D,&E,&F,&G,&H);
	}
}

void upper_bridge(double shift)
{	
	double a=find_a(h,l);
	
	double time=-1+(n-1)*2.0/(2.0*(double)n-1.0);
	
	double x=a*asinh(time*sinh(l/a));
	double z=h+a*(1-cosh(x/a));
	double length=x;
	time+=2.0/(2.0*(double)n-1.0);
	x=a*asinh(time*sinh(l/a));
	length=x-length;
	
	int m=ceil((l-x)/length);
		
	struct matrix A,B,C,D,E,F,G,H;
	
	for(int i=m;i>=-m;i--)
	{
		double X=x+i*length;
		
		loadm(&A,3,1);
		setm(&A,0,0,X);
		setm(&A,1,0,shift);
		setm(&A,2,0,z+upper_thick);
	
		loadm(&B,3,1);
		setm(&B,0,0,X);
		setm(&B,1,0,shift);
		setm(&B,2,0,z);
	
		loadm(&C,3,1);
		setm(&C,0,0,X);
		setm(&C,1,0,shift+width);
		setm(&C,2,0,z);
	
		loadm(&D,3,1);
		setm(&D,0,0,X);
		setm(&D,1,0,shift+width);
		setm(&D,2,0,z+upper_thick);
		
		X=x+(i-1)*length;
		
		loadm(&E,3,1);
		setm(&E,0,0,X);
		setm(&E,1,0,shift);
		setm(&E,2,0,z+upper_thick);
	
		loadm(&F,3,1);
		setm(&F,0,0,X);
		setm(&F,1,0,shift);
		setm(&F,2,0,z);
	
		loadm(&G,3,1);
		setm(&G,0,0,X);
		setm(&G,1,0,shift+width);
		setm(&G,2,0,z);
	
		loadm(&H,3,1);
		setm(&H,0,0,X);
		setm(&H,1,0,shift+width);
		setm(&H,2,0,z+upper_thick);
		
		xcube(&A,&B,&C,&D,&E,&F,&G,&H);
	
		if(i==m)
		{
			fix_node(&A);
			fix_node(&B);
			fix_node(&C);
			fix_node(&D);
		}
		if(i==-m)
		{
			fix_node(&E);
			fix_node(&F);
			fix_node(&G);
			fix_node(&H);
		}
	}
}


void arch_catenary(double shift)
{
	double a=find_a(h,l);

	double time=-1;
	
	for(int i=-n;i<n;i++)
	{
		double x=a*asinh(time*sinh(l/a));
		time+=2.0/(2.0*(double)n-1.0);
		double z=h+a*(1-cosh(x/a));
		
		double t=sinh(x/a);
		double s=t/sqrt(1+t*t);
		double c=sqrt(1-s*s);
		double m=l-s*arch_thick;
		
		double u=x-arch_thick*s;
		double w=z-arch_thick*c;
		
		struct matrix A,B,C,D,E,F,G,H;
		
		loadm(&A,3,1);
		setm(&A,0,0,x);
		setm(&A,1,0,shift);
		setm(&A,2,0,z);
	
		loadm(&B,3,1);
		setm(&B,0,0,u);
		setm(&B,1,0,shift);
		setm(&B,2,0,w);
	
		loadm(&C,3,1);
		setm(&C,0,0,u);
		setm(&C,1,0,shift+width);
		setm(&C,2,0,w);
	
		loadm(&D,3,1);
		setm(&D,0,0,x);
		setm(&D,1,0,shift+width);
		setm(&D,2,0,z);
		
		if(i==-1) i++;
		
		x=a*asinh(time*sinh(l/a));
		z=h+a*(1-cosh(x/a));
		
		t=sinh(x/a);
		s=t/sqrt(1+t*t);
		c=sqrt(1-s*s);
		m=l-s*arch_thick;
		
		u=x-arch_thick*s;
		w=z-arch_thick*c;
	
		loadm(&E,3,1);
		setm(&E,0,0,x);
		setm(&E,1,0,shift);
		setm(&E,2,0,z);
	
		loadm(&F,3,1);
		setm(&F,0,0,u);
		setm(&F,1,0,shift);
		setm(&F,2,0,w);
	
		loadm(&G,3,1);
		setm(&G,0,0,u);
		setm(&G,1,0,shift+width);
		setm(&G,2,0,w);
	
		loadm(&H,3,1);
		setm(&H,0,0,x);
		setm(&H,1,0,shift+width);
		setm(&H,2,0,z);
	
		xcube(&A,&B,&C,&D,&E,&F,&G,&H);

		if(i==-n)
		{
			fix_node(&A);
			fix_node(&B);
			fix_node(&C);
			fix_node(&D);
		}
		if(i==n-1)
		{
			fix_node(&E);
			fix_node(&F);
			fix_node(&G);
			fix_node(&H);
		}
	}
}

int main()
{
	omp_set_num_threads(4);

	set_parameters();

	inizialize_structure();
	
	arch_catenary(0);
	//upper_bridge(0);
	//arch_catenary(width);
	//upper_bridge(width);
	
	inizialize_force();
	add_weight();
	add_referential_force(phi);
	
	clock_t tic;
	tic=omp_get_wtime();
	equilibrium(1e-6);
	printf("Time elapsed: %f\n",(double)(omp_get_wtime()-tic));
	
	export_elements("e.dat");
	fprintm("X.dat",&X,"w");
	fprintm("x.dat",&x,"w");
	export_story("xs.dat");
}